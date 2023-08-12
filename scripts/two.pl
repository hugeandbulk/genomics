#!/usr/bin/perl
use warnings;
use strict;


print 'ENTER THE FILENAME OF THE DNA SEQUENCE:= ';
chomp( my $DNAfilename = <STDIN> );

open my $DNAFILE, $DNAfilename or die qq[Cannot open file "$DNAfilename" because: $!];

local $/;
( my $DNA = uc <$DNAFILE> ) =~ tr/ACGT//cd;

print "\nThe original DNA file is:\n$DNA\n";

my %codon2aa = qw(
    TCA  S  TCC  S  TCG  S  TCT  S  TTC  F  TTT  F  TTA  L  TTG  L
    TAC  Y  TAT  Y  TAA  _  TAG  _  TGC  C  TGT  C  TGA  _  TGG  W
    CTA  L  CTC  L  CTG  L  CTT  L  CCA  P  CCC  P  CCG  P  CCT  P
    CAC  H  CAT  H  CAA  Q  CAG  Q  CGA  R  CGC  R  CGG  R  CGT  R
    ATA  I  ATC  I  ATT  I  ATG  M  ACA  T  ACC  T  ACG  T  ACT  T
    AAC  N  AAT  N  AAA  K  AAG  K  AGC  S  AGT  S  AGA  R  AGG  R
    GTA  V  GTC  V  GTG  V  GTT  V  GCA  A  GCC  A  GCG  A  GCT  A
    GAC  D  GAT  D  GAA  E  GAG  E  GGA  G  GGC  G  GGG  G  GGT  G
    );

my $protein = '';
while ( $DNA =~ /(...)/g ) {
    exists $codon2aa{ $1 } or die qq[Bad codon "$1"!!\n];
    $protein .= $codon2aa{ $1 };
    }

print "The translated protein is :\n$protein\n";
<STDIN>;
