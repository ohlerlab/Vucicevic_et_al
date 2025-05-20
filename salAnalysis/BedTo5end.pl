#!/usr/bin/perl

use strict;
use warnings;
use feature ':5.10';
use Cwd;
my $dir = getcwd;

my @commands = @ARGV;
open BED, "<$commands[0]" or die "Cannot open input file.\n";
open OUT, ">$dir/$commands[0]_true5end.bed" or die "Cannot open output file.\n";

while(<BED>)
{
	chomp;
	my @d = split;
	if ($d[5] eq "+")
	{
		my $end = $d[1] + 1;
		print OUT "$d[0]\t$d[1]\t$end\t$d[3]\t$d[4]\t$d[5]\n";
	}
	if ($d[5] eq "-")
	{
		my $start = $d[2] - 1;
		print OUT "$d[0]\t$start\t$d[2]\t$d[3]\t$d[4]\t$d[5]\n";
	}
}
