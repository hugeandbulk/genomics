#!/usr/bin/env perl

=head1 NAME

    find_internal_stops.pl

=head1 SYNOPSIS
 
    find_internal_stops.pl input_fasta          
        where input_fasta is the input fasta file of protein translations.

=head1 DESCRIPTION

    This script takes an input fasta file of protein translations, and prints out
    a list of all proteins that have internal stop codons (not including the stop
    codon at the end of the translation).

=head1 VERSION
  
    Perl script last edited 3-Jul-2013.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut

# 
# Perl script find_internal_stops.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 3-Jul-13.
# Last edited 3-Jul-2013.
# SCRIPT SYNOPSIS: find_internal_stops.pl: given an input fasta file of protein translations, prints a list of proteins that have internal stop codons (assumes * is STOP).
#
#------------------------------------------------------------------#

# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:

use strict;
use warnings;

# xxx
# BEGIN { 
#     unshift (@INC, '/nfs/users/nfs_a/alc/Documents/git/helminth_scripts/modules'); 
# }

use HelminthGenomeAnalysis::AvrilGffUtils;
use HelminthGenomeAnalysis::AvrilFileUtils;
use HelminthGenomeAnalysis::AvrilFastaUtils;

my $num_args               = $#ARGV + 1;
if ($num_args != 1)
{
    print "Usage of find_internal_stops.pl\n\n";
    print "perl find_internal_stops.pl <input_fasta>\n";
    print "where <input_fasta> is the input fasta file of protein translations\n";
    print "For example, >perl find_internal_stops.pl final.translations.fa\n";
    exit;
}

# FIND THE PATH TO THE INPUT FASTA FILE:                     

my $input_fasta            = $ARGV[0];

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

&run_main_program($input_fasta);

print STDERR "FINISHED.\n";

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

sub run_main_program
{
   my $input_fasta         = $_[0]; # THE INPUT FASTA FILE
   my $input_fasta_obj;             # OBJECT FOR THE INPUT FASTA FILE 
   my $protname2internalstops;      # HASH TABLE OF THE NUMBER INTERNAL STOP CODONS IN SEQUENCES. 
   my $name;                        # PROTEIN NAME
   my $internal_stops;              # NUMBER OF INTERNAL STOP CODONS IN A PROTEIN
 
   # READ IN THE SEQUENCES IN THE INPUT FASTA FILE:
   $input_fasta_obj        = HelminthGenomeAnalysis::AvrilFastaUtils->new(fasta_file => $input_fasta); 

   # CHECK IF ANY OF THE SEQUENCES HAVE INTERNAL STOP CODONS:  
   $protname2internalstops = HelminthGenomeAnalysis::AvrilFastaUtils::count_internal_stops($input_fasta_obj);
   foreach $name (keys %{$protname2internalstops})
   {
      $internal_stops      = $protname2internalstops->{$name};
      print STDERR "WARNING: $name has $internal_stops internal stop codons\n";
   }
}

#------------------------------------------------------------------#
