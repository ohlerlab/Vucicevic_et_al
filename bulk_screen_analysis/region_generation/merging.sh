#!/bin/bash
bedtools merge -d 500 -c 4 -o collapse -i sig_window_bins.bed > merged_windowed_500.bed 