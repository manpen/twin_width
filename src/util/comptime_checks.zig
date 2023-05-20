pub fn checkIfIsCompatibleInteger(comptime T: type) bool {
    comptime if (T != u8 and T != u16 and T != u32) {
        return false;
    };
    return true;
}

pub fn log2Type(comptime T: type) type {
    return switch (T) {
        u8 => u3,
        u16 => u4,
        u32 => u5,
        u64 => u6,
        else => unreachable,
    };
}
