#!/usr/bin/env perl

#  remove_duplicate_fasta.pl
#  perl_diplicate_remover
#
#  Created by Richard Gunning on 8/6/14.
#  Copyright (c) 2014 Richard Gunning. All rights reserved.

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;              # avoid redundant &Usage()

my %sequences;

my $outfile = './output';
my $sequencefile = '';
my $linenumber = 0;
my $opt_help;
my $fastaid;

GetOptions ('file|f=s' => \$sequencefile,
'outfile|o=s' => \$outfile,
'help|h!' => \$opt_help,
) or pod2usage(-verbose => 1) && exit;;
pod2usage(-verbose => 1) && exit if defined $opt_help;


die "file doesn't exist\n" if (!-e $sequencefile);
open (my $outfh, '>', $outfile) or die "can't write output to file";
open (my $seqfh, '<', $sequencefile) or die "Can't Open File\n";
while (my $line = <$seqfh>)
{
    $linenumber++;
    chomp $line;
    die "not sure how you managed this" if $linenumber == 0;
    if ($linenumber == 1) {
        if (substr($line,0,1) eq '>') {
                $fastaid = $line;
        } else { die "multi-line fasta at line $." }
    } elsif ($linenumber == 2) {
        $sequences{$line} = $fastaid;
        $linenumber = 0;
    } else {
        die "Are You some kind of magician"
    }
}

while (my $id = each %sequences) {
    print $outfh "$sequences{$id}\n$id\n";
}


close $seqfh;
close $outfh;

=head1 NAME
 
 perl_duplicate_remover.pl
 
=head1 SYNOPSIS
 
 perl perl_duplicate_remover.pl [Options]
 
 remove Duplicates sequences from Fasta file
 Keep second occurence of sequence
 
=head1 OPTIONS
 
 -f,--file       input file
 -o,--outfile    output file
 -h,--help
 
=cut