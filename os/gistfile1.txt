#!/usr/bin/env perl

=head1 NAME

    split_up_fasta.pl

=head1 SYNOPSIS
 
    split_up_fasta.pl input_fasta num_seqs_per_file prefix outputdir
        where input_fasta is the input fasta file,
              num_seqs_per_fasta is the number of sequences to put in each fasta file,
              prefix is the prefix to use for the output files,
              outputdir is the output directory for writing output files. 

=head1 DESCRIPTION

    This script takes an input fasta file, and splits it up into new files, with 
    num_seqs_per_file sequences in each new file. The new files are written in the
    directory outputdir, and have prefix 'prefix' (are called prefix1, prefix2, etc.)
    
=head1 VERSION
  
    Perl script last edited 5-Sept-2012.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut

# 
# Perl script split_up_fasta.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 5-Sept-12.
# Last edited 5-Sept-2012.
# SCRIPT SYNOPSIS: split_up_fasta.pl: splits up an input fasta file into smaller files with n sequences each
#
#------------------------------------------------------------------#

# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:

use strict;
use warnings;
use POSIX; # HAS THE ceil() FUNCTION

# 
# BEGIN { 
#    unshift (@INC, '/nfs/users/nfs_a/alc/Documents/git/helminth_scripts/modules'); 
#}

use HelminthGenomeAnalysis::AvrilFastaUtils;

my $num_args               = $#ARGV + 1;
if ($num_args != 4)
{
    print "Usage of split_up_fasta.pl\n\n";
    print "perl split_up_fasta.pl <input_fasta> <num_seqs_per_file> <prefix> <outputdir>\n";
    print "where <input_fasta> is the input fasta file,\n";
    print "      <num_seqs_per_file> is the number of sequences to put in each fasta file,\n";
    print "      <prefix> is the prefix to use for the output files,\n";
    print "      <outputdir> is the output directory for writing output files\n";
    print "For example, >perl split_up_fasta.pl velvet.k55.04.sspace_with_454_8kb.assembly.fa 100\n";
    print "velvet /lustre/scratch108/parasites/alc/ReaprRepeats/Plasmodium\n";
    exit;
}

# FIND THE PATH TO THE INPUT FASTA FILE:                     

my $input_fasta            = $ARGV[0];

# FIND THE NUMBER OF SEQUENCES TO PUT IN EACH FASTA FILE:

my $num_seqs_per_file      = $ARGV[1];

# FIND THE PREFIX TO USE FOR THE OUTPUT FILES:

my $prefix                 = $ARGV[2];

# FIND THE DIRECTORY TO USE FOR OUTPUT FILES:      

my $outputdir              = $ARGV[3];

#------------------------------------------------------------------#

# TEST SUBROUTINES: 

my $PRINT_TEST_DATA        = 0;   # SAYS WHETHER TO PRINT DATA USED DURING TESTING.
&test_print_error;
&test_run_main_program($outputdir);
print STDERR "Finished running tests, running main code now...\n";

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

&run_main_program($outputdir,$input_fasta,$num_seqs_per_file,$prefix,0);

print STDERR "FINISHED.\n";

#------------------------------------------------------------------#

# TEST &run_main_program

sub test_run_main_program
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN.
   my $random_number;               # RANDOM NUMBER TO USE IN TEMPORARY FILE NAMES.
   my $fasta;                       # INPUT FASTA FILES.
   my $prefix;                      # PREFIX TO USE FOR THE OUTPUT FILES   
   my $i;                           # 
   my $output;                      # NAME OF THE OUTPUT FILE 
   my $seqcnt;                      # NUMBER OF SEQUENCES IN A FILE. 
   my $line;                        # 
   my $errorcode;                   # RETURNED AS 0 IF THERE IS NO ERROR.
   my $errormsg;                    # RETURNED AS 'none' IF THERE IS NO ERROR.

   $random_number            = rand();
   $fasta                    = $outputdir."/tmp".$random_number;
   open(FASTA,">$fasta") || die "ERROR: test_split_up_input_file: cannot open $fasta\n";
   print FASTA ">seq1\n";
   print FASTA "AAAAA\n";
   print FASTA ">seq2\n";
   print FASTA "TTTTT\n";
   print FASTA ">seq3\n";
   print FASTA "GGGGG\n";
   print FASTA ">seq4\n";
   print FASTA "CCCCC\n";
   print FASTA ">seq5\n";
   print FASTA "ACTGA\n";
   print FASTA ">seq6\n";
   print FASTA "ATCAG\n";
   print FASTA ">seq7\n";
   print FASTA "AAAAA\n";
   close(FASTA); 
   $random_number            = rand();
   $prefix                   = "tmp".$random_number."file";
   ($errorcode,$errormsg)    = &run_main_program($outputdir,$fasta,3,$prefix,1);
   # CHECK THAT THERE ARE THE EXPECTED NUMBER OF OUTPUT FILES, AND THEY HAVE THE EXPECTED NUMBER OF SEQUENCES EACH:
   for ($i = 1; $i <= 3; $i++)
   {
      $output                = $prefix.$i; 
      $output                = $outputdir."/".$output;
      if (!(-e "$output"))
      {
         print STDERR "ERROR: test_split_up_input_file: failed test1a (output $output does not exist)\n";
         exit;
      }
      # CHECK THAT THE FILE HAS THE EXPECTED NUMBER OF SEQUENCES: 
      $seqcnt                = 0;
      open(OUTPUT,"$output") || die "ERROR:  test_split_up_input_file: cannot open $output\n";
      while(<OUTPUT>)
      {
         $line               = $_;
         chomp $line;
         if (substr($line,0,1) eq ">") { $seqcnt++;}
      }
      close(OUTPUT);
      if ($i == 1)
      {
         if ($seqcnt != 2) { print STDERR "ERROR: test_split_up_input_file: failed test1b (seqcnt = $seqcnt)\n"; exit;}
      }
      elsif ($i == 2)
      {
         if ($seqcnt != 3) { print STDERR "ERROR: test_split_up_input_file: failed test1c (seqcnt = $seqcnt)\n"; exit;}
      }
      elsif ($i == 3)
      {
         if ($seqcnt != 2) { print STDERR "ERROR: test_split_up_input_file: failed test1d (seqcnt = $seqcnt)\n"; exit;}
      }
      system "rm -f $output";
   }
   system "rm -f $fasta\n";

}

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

sub run_main_program
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN.
   my $input_fasta         = $_[1]; # THE INPUT FASTA FILE
   my $num_seqs_per_file   = $_[2]; # THE NUMBER OF SEQUENCES TO PUT IN EACH OUTPUT FILE
   my $prefix              = $_[3]; # THE PREFIX TO USE FOR THE OUTPUT FILES. 
   my $testing             = $_[4]; # SAYS WHETHER THIS WAS CALLED BY A TESTING SUBROUTINE.
   my $errorcode;                   # RETURNED AS 0 IF THERE IS NO ERROR.
   my $errormsg;                    # RETURNED AS 'none' IF THERE IS NO ERROR. 
   my $input_fasta_obj;             # OBJECT FOR THE INPUT FASTA FILE 
   my $contigs2len;                 # HASH TABLE OF LENGTHS OF SCAFFOLDS IN THE ASSEMBLY
   my $input_fasta_contigs2seq;     # HASH TABLE OF SEQUENCES IN $input_fasta 
   my @scaffolds;                   # ARRAY OF SCAFFOLD NAMES
   my $filenum;                     # FILE NUMBER TO PUT A CERTAIN SEQUENCE INTO 
   my $num_files;                   # NUMBER OF OUTPUT FILES THAT WE WILL NEED
   my $output;                      # OUTPUT FILE 
   my $i;                           # 
   my $scaffold;                    # A SCAFFOLD 
   my $seq;                         # SCAFFOLD SEQUENCE
   my $returnvalue;                  # RETURNVALUE FROM A FUNCTION
   
   # READ IN THE SEQUENCES IN THE INPUT FASTA FILE FOR THE ASSEMBLY:
   $input_fasta_obj        = HelminthGenomeAnalysis::AvrilFastaUtils->new(fasta_file => $input_fasta); 
   $input_fasta_contigs2seq= $input_fasta_obj->_contigs2seq;

   # GET THE LENGTHS OF THE SEQUENCES IN THE INPUT FASTA FILE FOR THE ASSEMBLY:
   $contigs2len            = HelminthGenomeAnalysis::AvrilFastaUtils::build_contigs2len($input_fasta_obj),

   # GET AN ARRAY OF ALL THE SCAFFOLD NAMES:
   @scaffolds              = keys %{$contigs2len};
   # SORT BY SCAFFOLD LENGTH:
   @scaffolds              = sort { $contigs2len->{$a} <=> $contigs2len->{$b} } @scaffolds;
  
   # FIND OUT HOW MANY OUTPUT FILES WE WILL NEED:
   $num_files              = ($#scaffolds + 1)/$num_seqs_per_file;
   $num_files              = ceil($num_files); # ROUND UP
   
   # FIGURE OUT WHICH OUTPUT FILE EACH SCAFFOLD WILL GO INTO:
   for ($i = 1; $i <= $#scaffolds + 1; $i++)
   {
      $scaffold            = $scaffolds[($i-1)];
      $filenum             = ($i % $num_files) + 1;
      $output              = $prefix.$filenum;
      $output              = $outputdir."/".$output; 
      if (-e $output)
      {
         $seq              = $input_fasta_contigs2seq->{$scaffold};
         $returnvalue      = HelminthGenomeAnalysis::AvrilFastaUtils::print_seq_to_fasta($output,$seq,$scaffold,'yes');
      }
      else
      {
         $seq              = $input_fasta_contigs2seq->{$scaffold};
         if ($testing == 0) { print STDERR "Opening $output...\n";}
         $returnvalue      = HelminthGenomeAnalysis::AvrilFastaUtils::print_seq_to_fasta($output,$seq,$scaffold,'no');
      }
   } 

}

#------------------------------------------------------------------#

# TEST &print_error

sub test_print_error
{
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE WAS NO ERROR
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE WAS NO ERROR

   ($errormsg,$errorcode)  = &print_error(45,45,1);
   if ($errorcode != 12) { print STDERR "ERROR: test_print_error: failed test1\n"; exit;}

   ($errormsg,$errorcode)  = &print_error('My error message','My error message',1);
   if ($errorcode != 11) { print STDERR "ERROR: test_print_error: failed test2\n"; exit;}

   ($errormsg,$errorcode)  = &print_error('none',45,1);
   if ($errorcode != 13) { print STDERR "ERROR: test_print_error: failed test3\n"; exit;} 

   ($errormsg,$errorcode)  = &print_error('My error message', 0, 1);
   if ($errorcode != 13) { print STDERR "ERROR: test_print_error: failed test4\n"; exit;}
}

#------------------------------------------------------------------#

# PRINT OUT AN ERROR MESSAGE AND EXIT.
 
sub print_error
{
   my $errormsg            = $_[0]; # THIS SHOULD BE NOT 'none' IF AN ERROR OCCURRED.
   my $errorcode           = $_[1]; # THIS SHOULD NOT BE 0 IF AN ERROR OCCURRED.
   my $called_from_test    = $_[2]; # SAYS WHETHER THIS WAS CALLED FROM test_print_error OR NOT

   if ($errorcode =~ /[A-Z]/ || $errorcode =~ /[a-z]/) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 11; $errormsg = "ERROR: print_error: the errorcode is $errorcode, should be a number.\n"; # ERRORCODE=11
         return($errormsg,$errorcode);
      }
      else 
      { 
         print STDERR "ERROR: print_error: the errorcode is $errorcode, should be a number.\n"; 
         exit;
      }
   }

   if (!($errormsg =~ /[A-Z]/ || $errormsg =~ /[a-z]/)) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 12; $errormsg = "ERROR: print_error: the errormessage $errormsg does not seem to contain text.\n"; # ERRORCODE=12
         return($errormsg,$errorcode);
      }
      else
      {
         print STDERR "ERROR: print_error: the errormessage $errormsg does not seem to contain text.\n"; 
         exit;
      }
   }

   if    ($errormsg eq 'none' || $errorcode == 0) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 13; $errormsg = "ERROR: print_error: errormsg $errormsg, errorcode $errorcode.\n"; # ERRORCODE=13
         return($errormsg,$errorcode);
      }
      else 
      {
         print STDERR "ERROR: print_error: errormsg $errormsg, errorcode $errorcode.\n"; 
         exit;
      }
   }
   else                                           
   { 
      print STDERR "$errormsg"; 
      exit;                                                      
   } 

   return($errormsg,$errorcode);
}

#------------------------------------------------------------------#


