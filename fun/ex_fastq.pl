#!/usr/bin/perl
#
## usage:  extractSeqByID.pl SEQ123 < huge.fsa > my.fsa
#
use warnings;
use strict;

my $lookup = shift @ARGV;  # ID to extract
my $lookupf;
open($lookupf, "<", $lookup);
my @lookup;
while(<$lookupf>){
        chomp;
        push @lookup, $_;
}

local $/ = "\n>";  # read by FASTA record

while (my $seq = <>) {
        chomp $seq;
        foreach my $i (@lookup){
                if($seq =~ /$i/){
                        #$seq =~ s/^>*.+\n//;  # remove FASTA header
                        #$seq =~ s/\n//g;  # remove endlines
			print ">";
                        print "$seq";
			print "\n";
                }
        }
}

