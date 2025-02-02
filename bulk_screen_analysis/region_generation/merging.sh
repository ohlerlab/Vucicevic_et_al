# This script uses bedtools to merge overlapping or nearby intervals in a BED file.
# The intervals are merged if they are within 500 base pairs of each other.
# The fourth column of the merged intervals will contain a comma-separated list of values from the original intervals.
# The input file is 'sig_window_bins.bed' and the output file is 'merged_windowed_500.bed'.
#!/bin/bash
bedtools merge -d 500 -c 4 -o collapse -i sig_window_bins.bed > merged_windowed_500.bed 