#!/usr/bin/env perl

#  transcript_to_gene.pl
#  transcript_to_gene.
#
#  Created by Richard Gunning on 8/6/14.
#  Copyright (c) 2014 Richard Gunning. All rights reserved.

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;              # avoid redundant &Usage()
#use Bio::Tools::GFF;


my %tgmhash;
my %genereads;
my %genekmers;
my %genetpm;
my %generpkm;
my %genekpkm;

my $outfile = './output';
my $sailfish = '/Users/rg12/Documents/Perl/transcript_to_gene/job6quant.sf';
my $tgm = '/Users/rg12/Documents/Perl/transcript_to_gene/gencode.vM3.annotation.transcript.ERCC.gff3';

my $opt_help;
my $fastaid;

GetOptions ('tgm|t=s' => \$tgm,
'sailfish|s=s' => \$sailfish,
'outfile|o=s' => \$outfile,
'help|h!' => \$opt_help,
) or pod2usage(-verbose => 1) && exit;;
pod2usage(-verbose => 1) && exit if defined $opt_help;


die "file doesn't exist\n" if (!-e $sailfish);
die "file doesn't exist\n" if (!-e $tgm);

open (my $outfh, '>', $outfile) or die "can't write output to file";
open (my $sailfh, '<', $sailfish) or die "Can't Open File\n";
open (my $tgmfh, '<', $tgm) or die "Can't Open File\n";

#read tgm to hash
while (my $line = <$tgmfh>){
    chomp $line;
    # read attributes from column 9
    my @columns = split(/\t/,$line);
    my @attributes = split(";",$columns[8]);
    # find transcript_id and gene_id
    my ($rest1, $gene_id) = split(/=/,$attributes[2]);
    my ($gene,$junk1) = split(/\./,$gene_id);
    #print "$gene \n";
    my ($rest, $transcript_id) = split(/=/,$attributes[3]);
    my ($transcript,$junk) = split(/\./,$transcript_id);
    #print "$transcript \n";
    $tgmhash{$transcript}=$gene;
}

#foreach my $transcript (keys %tgmhash) {
#    print "$transcript: $tgmhash{$transcript}\n";
#}
my $failurecount = 0;
while (my $line = <$sailfh>)
{
    next if ($line =~ /^#/);
    chomp $line;
    my @columns = split(/\t/,$line);
#    $columns[0] =~ /(^\w+.\d)/;
    my ($transcript_id,$junk) = split(/\|/,$columns[0],2);
    my ($transcript,$morejunk) = split(/\./,$transcript_id,2);
    if (exists $tgmhash{$transcript}) {
        my $gene_id = $tgmhash{$transcript};
        if (exists $genetpm{$gene_id}) {
            $genetpm{$gene_id} += $columns[2];
            $generpkm{$gene_id} += $columns[3];
            $genekpkm{$gene_id} += $columns[4];
            $genekmers{$gene_id} += $columns[5];
            $genereads{$gene_id} += $columns[6];
        }else{
            $genetpm{$gene_id} = $columns[2];
            $generpkm{$gene_id} = $columns[3];
            $genekpkm{$gene_id} = $columns[4];
            $genekmers{$gene_id} = $columns[5];
            $genereads{$gene_id} = $columns[6];
        }
    }else{
        if (exists $genetpm{"NOT_in_tgm"}) {
            $genetpm{"NOT_in_tgm"} += $columns[2];
            $generpkm{"NOT_in_tgm"} += $columns[3];
            $genekpkm{"NOT_in_tgm"} += $columns[4];
            $genekmers{"NOT_in_tgm"} += $columns[5];
            $genereads{"NOT_in_tgm"} += $columns[6];
        }else{
            $genetpm{"NOT_in_tgm"} = $columns[2];
            $generpkm{"NOT_in_tgm"} = $columns[3];
            $genekpkm{"NOT_in_tgm"} = $columns[4];
            $genekmers{"NOT_in_tgm"} = $columns[5];
            $genereads{"NOT_in_tgm"} = $columns[6];
        }
        $failurecount ++;
    }
}

print "$failurecount:\tTRANSCRIPTS FOUND THAT ARE NOT IN TGM\n";

print $outfh "Gene\tTPM\tRPKM\tKPKM\tEstimatedNumKmers\tEstimatedNumReads\n";
foreach my $gene_id (keys %genetpm) {
    print $outfh "$gene_id\t$genetpm{$gene_id}\t$generpkm{$gene_id}\t$genekpkm{$gene_id}\t$genekmers{$gene_id}\t$genereads{$gene_id}\n"
}

close $tgmfh;
close $sailfh;
close $outfh;

=head1 NAME
 
 transcript_to_gene.pl
 
=head1 SYNOPSIS
 
 perl transcript_to_gene.pl [Options]
 
 remove Duplicates sequences from Fasta file
 Keep second occurence of sequence
 
=head1 OPTIONS
 
 -t,-tgm        tgm input file
 -s -sailfish   sailfish Quant file
 -o,-outfile    output file
 -h,-help
 
=cut