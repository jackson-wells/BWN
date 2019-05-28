#!/usr/bin/perl -w
use strict;
use POSIX qw(ceil floor);
use Storable qw(dclone);
use Getopt::Long;
use Data::Dumper;


my @aminoAcids = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y");
my $searchString = $aminoAcids[int rand(20)];

for(my $i = 0; $i < 25; $i++)
{
        $searchString = $aminoAcids[int rand(20)].$searchString;
}


my $fh;

open($fh, '>', "longDistanceAnalysis$$.txt") or die $!;

for(my $i = 0; $i < 25; $i++)
{
	my $avgTime = 0;
	for(my $j = 0; $j < 20; $j++)
	{
		my $outputTime = `(time -f '%e' bwp-search -s $searchString -a 1) 2>&1`;
		$outputTime =~ s/\n//g;
		$avgTime += $outputTime;
	}
	$avgTime = $avgTime/20;
	my $stringLength = $i + 26;
	print $fh "$stringLength\t$avgTime\n";
	$searchString = $aminoAcids[int rand(20)].$searchString;
}

print $fh "$searchString\n";

close($fh);
