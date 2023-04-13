#!/usr/bin/env bash
for sol in instances/heuristic-public/*.solution; do
    graph=$(echo "$sol" | perl -pe "s/.gr.*.solution\$/.gr/")
    
    if ./pace_verifier.py $graph $sol ; then
        # verified -> let keep that solution
        hash=$(grep -v 'c' $sol | sha1sum -)
        archived="$graph.vsolution-${hash:0:8}"
        
        if [ ! -f $archived ]; then
            echo "Found new solutiuon $archived"
            cp $sol $archived
        fi
    fi
done
