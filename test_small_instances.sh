#/usr/bin/env bash
set -e
echo "Compile ..."
zig build main_exact -Doptimize=ReleaseSafe -Dcpu=haswell

echo "Compile list ..."
#rm -rf list_test_instances
for d in instances/small-random/n*; do
    echo $d;
    ls -1 $d/*cir*.gr >> list_test_instances

done
wc -l list_test_instances

echo "Run ..."
cat list_test_instances | \
    parallel --shuf zig-out/bin/main_exact ::: 2>&1 | \
    tee test.log | grep -E "FAILED|panic" 
    
    
