#!/bin/bash
for f in ../lore_count_spacers/*_R1_001_lore_library_count.csv
    do
            SAMPLE=$(basename $f _R1_001_lore_library_count.csv)
            echo "GUIDE,${SAMPLE}" | cat - $f | sed 's///' > "$SAMPLE".txt
    done            
