const std = @import("std");

pub const TabulationHashFunction = struct {
    mask: u64 = 0xFF,
    T1: [256]u64,
    T2: [256]u64,
    T3: [256]u64,
    T4: [256]u64,
    T5: [256]u64,
    T6: [256]u64,
    T7: [256]u64,
    T8: [256]u64,

    pub inline fn init(rand: *std.rand.DefaultPrng) TabulationHashFunction {
        var hashd = TabulationHashFunction{ .T1 = undefined, .T2 = undefined, .T3 = undefined, .T4 = undefined, .T5 = undefined, .T6 = undefined, .T7 = undefined, .T8 = undefined };

        for (0..256) |i| {
            hashd.T1[i] = rand.next();
            hashd.T2[i] = rand.next();
            hashd.T3[i] = rand.next();
            hashd.T4[i] = rand.next();
            hashd.T5[i] = rand.next();
            hashd.T6[i] = rand.next();
            hashd.T7[i] = rand.next();
            hashd.T8[i] = rand.next();
        }
        return hashd;
    }

    pub fn hash(self: *TabulationHashFunction, x: u64) u64 {
        var hashvalue = self.T1[x & self.mask] ^ self.T2[(x >> 8) & self.mask] ^ self.T3[(x >> 16) & self.mask] ^ self.T4[(x >> 24) & self.mask] ^
            self.T5[(x >> 32) & self.mask] ^ self.T6[(x >> 40) & self.mask] ^ self.T7[(x >> 48) & self.mask] ^ self.T8[(x >> 56) & self.mask];
        return hashvalue;
    }
};

pub const TabulationHashFunction32 = struct {
    mask: u64 = 0xFF,
    T1: [256]u64,
    T2: [256]u64,
    T3: [256]u64,
    T4: [256]u64,

    pub inline fn init(rand: *std.rand.DefaultPrng) TabulationHashFunction {
        var hashd = TabulationHashFunction{
            .T1 = undefined,
            .T2 = undefined,
            .T3 = undefined,
            .T4 = undefined,
        };

        for (0..256) |i| {
            hashd.T1[i] = rand.next();
            hashd.T2[i] = rand.next();
            hashd.T3[i] = rand.next();
            hashd.T4[i] = rand.next();
        }
        return hashd;
    }

    pub fn hash(self: *TabulationHashFunction, x: u64) u64 {
        return self.T1[x & self.mask] ^ self.T2[(x >> 8) & self.mask] ^ self.T3[(x >> 16) & self.mask] ^ self.T4[(x >> 24) & self.mask];
    }
};


pub const DartHash = struct {
	pub const Dart = struct {
		index: u64,
		rank: f64,
	};

	t: u64,
	T_nu: TabulationHashFunction32,
	T_rho: TabulationHashFunction32,
	T_w: TabulationHashFunction32,
	T_r: TabulationHashFunction32,

	T_i: TabulationHashFunction,
	T_p: TabulationHashFunction,
	T_q: TabulationHashFunction,
	F: TabulationHashFunction,
	M: TabulationHashFunction,

	powers_of_two: std.ArrayList(f64),
	negative_powers_of_two: std.ArrayList(f64),
	poisson_cdf: std.ArrayList(f64),

	darts: std.ArrayList(Dart),

	pub fn init(allocator: std.mem.Allocator, t: u64, rng: *std.rand.DefaultPrng) !DartHash {
		var powers_of_two = std.ArrayList(f64).initCapacity(allocator,1000);
		var negative_powers_of_two = std.ArrayList(f64).initCapacity(allocator,1000);
    var p: f64 = 1.0;
    var q: f64 = 1.0;
		
		for(0..1000) |_| {
			try powers_of_two.append(p);
			p = 2.0*p;
			try negative_powers_of_two.append(q);
			q = 0.5*q;
		}

		var poisson_cdf = std.ArrayList(f64).initCapacity(allocator,100);

    var pdf:f64 = std.math.exp(-1.0);
		var cdf:f64 = pdf;
		for(0..100) |i| {
			poisson_cdf.push_back(cdf);
			pdf = pdf/(i + 1);
			cdf += pdf;
		}

		var darts= std.ArrayList(Dart).initCapacity(allocator,3*t);

		return DartHash {
			.t = t,

			.T_nu = TabulationHashFunction32.init(rng),
			.T_rho = TabulationHashFunction32.init(rng),
			.T_w = TabulationHashFunction32.init(rng),
			.T_r = TabulationHashFunction32.init(rng),

			.T_i = TabulationHashFunction.init(rng),
			.T_p = TabulationHashFunction.init(rng),
			.T_q = TabulationHashFunction.init(rng),
			.F = TabulationHashFunction.init(rng),
			.M = TabulationHashFunction.init(rng),

			.powers_of_two = powers_of_two,
			.negative_powers_of_two = negative_powers_of_two,
			.poisson_cdf = poisson_cdf,

			.darts = darts,
		};
	}

		 pub fn toUnit(x: u64) f64 {
			 return @intToFloat(f64,x/0xFFFFFFFFFFFFFFFF);
		 }

		pub fn toUnits(x: u64) struct{f64,f64} {
			 return .{@intToFloat(f64,(x>>32)/0xFFFFFFFF), @intToFloat(f64,(x&0xFFFFFFFF)/0xFFFFFFFF)};
		 }



	pub fn hash(self: *DartHash, comptime IteratorType: type, iter: IteratorType, total_weight: f64, theta: f64) !*std.ArrayList(Dart) {
		self.darts.clearRetainingCapacity();

		var max_rank: f64 = theta/total_weight;
		var t_inv:f64 = 1.0/self.t;
		var RHO:u32 = @floor(std.math.log2(1.0 + max_rank));
		

		while(iter.next()) |it| {
			const index: u64 = it.index;
			const weight: f64 = it.weight;

			const index_hash: u64 = self.T_i.hash(index);
			
			const NU:u32 = @floor(std.math.log2(1.0+self.t*weight));
			

			for(0..NU+1) |nu| {
				const nu_hash = self.T_nu.hash(nu);

				for(0..RHO+1) |rho| {
					const region_hash:u64 = nu_hash ^ self.T_rho.hash(rho);
					
					const two_nu:f64 = self.powers_of_two.items[nu];
					const two_rho:f64 = self.powers_of_two.items[rho];
					
					const W:f64 = (two_nu-1.0)*t_inv;
					const R:f64 = two_rho-1.0;
					
					const delta_nu:f64 = two_nu * t_inv * self.negative_powers_of_two[rho];
					
					const delta_rho:f64 = two_rho * self.negative_powers_of_two[nu];
					
					var w0:f64 = W;
					
					const w_max:u32 = if(rho < 32) (1<<rho) else (1<<31);
					

					for(0..w_max) |w| {
						if(weight < w0) break;
						const w_hash:u64 = self.T_w.hash(w);
						var r0:f64 = R;
						
						const r_max: u32 = if(nu < 32) (1<<nu) else (1<<31);
						for(0..r_max) |r| {
							if(max_rank < r0) break;
							const area_hash:u64 = w_hash^self.T_r.hash(r);
							const z:u64 = index_hash ^ region_hash ^ area_hash;

							const p_z:f64 = DartHash.toUnit(self.T_p.hash(z));
							var p:u8 = 0;
							while(p_z > self.poisson_cdf.items[p]) {
								p+=1;
							}

							var q:u64 = 0;
							while(q < p) {
								const z_q:u64 = z ^ (q<<56) ^ (q<<48) ^ (q<<40) ^ (q<<32) ^ (q<<24) ^ (q<<16) ^ (q<<8) ^ q;

								const uniform_weight_rank:f64 = DartHash.toUnits(self.T_q.hash(z_q));
								const dart_weight:f64 = w0 + delta_nu*uniform_weight_rank.@"0";
								const rank:f64 = r0 + delta_rho*uniform_weight_rank.@"1";
								if(dart_weight < weight and rank < max_rank) {
									try self.darts.append(Dart {
										.index = self.F.hash(z_q),
										.rank = dart_weight,
									});
								}
								q+=1;
							}

							r0 += delta_rho;
						}
						w0 += delta_nu;
					}
				}
			}
		}
		return &self.darts;
	}

};


pub const DartMinHash = struct {
	k: u64,
	pseudo_hash: TabulationHashFunction,
	dart: DartHash,

	min_hash: []u8,
	ranks: []f64,

	pub fn init(allocator: std.mem.Allocato, rng: *std.rand.DefaultPrng, k: u64) DartMinHash {
		var dart = DartHash.init(allocator,k*std.math.log2_int(u64,k)+2*k, rng);
		var pseudo_hash = TabulationHashFunction.init(rng);
		
		var ranks = try allocator.alloc(f64,k);
		var min_hash = try allocator.alloc(u8,k);

		return DartMinHash {
			.k = k,
			.pseudo_hash = pseudo_hash,
			.dart = dart,
			.ranks = ranks,
			.min_hash = min_hash,
		};
	}


	pub fn hash(self: *DartMinHash, comptime IteratorType: type, iter: IteratorType, total_weight: f64) ![]u8 {
		@memset(self.ranks, std.math.f64_max);

		var theta: f64 = 1.0;
		var hashed:u32 = 0;
		while(hashed < self.k) {
			const darts = try self.dart.hash(IteratorType, iter, total_weight, theta);
			for(darts) |*dart| {
				const j = self.hash.hash(dart.index)%self.k;

				if(dart.rank < self.ranks[j]) {
					if(self.ranks[j] == std.math.f64_max) {
						hashed+=1;
					}
					self.ranks[j] = dart.rank;
					self.min_hash[j] = @truncate(u8,dart.index);
				}

				theta += 0.5;
				iter.reset();
			}
			return self.min_hash;
		}
	}
	
};
