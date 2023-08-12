#!/usr/local/bin/perl

=head1 NAME

    translate_treefam_cigars_to_hmms.pl

=head1 SYNOPSIS
 
    translate_treefam_cigars_to_hmm.pl treefam_version hmms_output outputdir cigars alntype family hmmer_bin alns_output map_output

=head1 DESCRIPTION

    For a TreeFam family of interest (<family>), this script reads in the cigar-format
    alignment from a file with cigar-format TreeFam alignments (<cigars>), for 
    a particular version of the TreeFam database, and translates the cigar alignment
    for the family of interest into a fasta format alignment, and then builds a HMM for 
    the alignment using HMMER. The argument <alntype> says whether "full" or "seed"  
    alignments should be used. The HMM is written in the output file <hmms_output>. 
    The multiple alignment is written in the output file <alns_output>. The hmmbuild
    map output is written to file <map_output>. <hmmer_bin> is the directory containing 
    the hmmer executables.

=head1 VERSION
  
    Perl script last edited 23-Oct-2012.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut

# 
# Perl script translate_treefam_cigars_to_hmms.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 23-Oct-12.
# Last edited 23-Oct-2012.
# SCRIPT SYNOPSIS: translate_treefam_cigars_to_hmms.pl: reads in a cigar-format alignment for a TreeFam familiy, and makes a HMM for the family
#
#------------------------------------------------------------------#

# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:

use strict;
use warnings;
use DBI;
use Devel::Size qw(size total_size);

my $num_args               = $#ARGV + 1;
if ($num_args != 9)
{
    print "Usage of translate_treefam_cigars_to_hmms.pl\n\n";
    print "perl translate_treefam_cigars_to_hmms.pl <treefam_version> <hmms_output> <outputdir> <cigars> <alntype> <family> <hmmer_bin> <alns_output> <map_output>\n";
    print "where <treefam_version> is the version of TreeFam to use (eg. 7),\n";
    print "      <hmms_output> is the output file of HMMs,\n";
    print "      <outputdir> is the directory for writing output files,\n";
    print "      <cigars> is the file with cigar-format alignments,\n";
    print "      <family> is the family of interest,\n";
    print "      <hmmer_bin> is the directory containing the hmmer executable\n";
    print "For example, >perl translate_treefam_cigars_to_hmms.pl 7 myhmms\n";
    print "/nfs/users/nfs_a/alc/Documents/StrongyloidesTreeFam alignments full TF101000\n";
    print "/software/hmmer/bin/ myalns mymap\n";
    exit;
}

# FIND THE TREEFAM VERSION TO USE:                           

my $treefam_version        = $ARGV[0];

# FIND THE OUTPUT FILE NAME OF HMMs:

my $hmms_output            = $ARGV[1];

# FIND THE DIRECTORY TO WRITE OUTPUT FILES TO:

my $outputdir              = $ARGV[2];

# FIND THE FILE WITH CIGAR-FORMAT ALIGNMENTS:

my $cigars                 = $ARGV[3];

# FIND THE TYPE OF ALIGNMENT TO USE:

my $alntype                = $ARGV[4];

# FIND THE FAMILY OF INTEREST:            

my $family                 = $ARGV[5];

# FIND THE DIRECTORY WITH THE HMMER EXECUTABLES:

my $hmmer_bin              = $ARGV[6];

# FIND THE OUTPUT MULTIPLE ALIGNMENT FILE:

my $alns_output            = $ARGV[7];

# FIND THE hmmbuild MAP OUTPUT FILE:

my $map_output             = $ARGV[8];

# TEST SUBROUTINES: 

my $PRINT_TEST_DATA        = 0;   # SAYS WHETHER TO PRINT DATA USED DURING TESTING.
my $PRINT_WARNINGS         = 0;   # SAYS WHETHER TO PRINT WARNINGS ABOUT MISSING DATA. 
&test_print_error;
&test_get_fasta_aln_for_seq;
&test_make_hmm($outputdir);
&test_add_sequences_to_cigarfile($outputdir);
&test_remove_nonunique_seqs_from_cigarfile($outputdir);

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

&run_main_program($treefam_version,$hmms_output,$outputdir,$cigars,$alntype,$family,$hmmer_bin,$alns_output,$map_output);

print STDERR "FINISHED.\n";

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

sub run_main_program
{
   my $version             = $_[0]; ## THE TREEFAM VERSION TO USE. 
   my $output              = $_[1]; ## THE OUTPUT FILE.     
   my $outputdir           = $_[2]; ## DIRECTORY TO WRITE OUTPUT FILES TO.
   my $cigarfile           = $_[3]; ## THE FILE WITH TREEFAM CIGAR-FORMAT ALIGNMENT.
   my $alntype             = $_[4]; ## FIND THE TYPE OF ALIGNMENT TO USE
   my $family              = $_[5]; ## THE FAMILY OF INTEREST       
   my $hmmer_build         = $_[6]; ## DIRECTORY CONTAINING THE HMMER EXECUTABLES
   my $alns_output         = $_[7]; ## OUTPUT MULTIPLE ALIGNMENT FILE 
   my $map_output          = $_[8]; ## hmmbuild MAP OUTPUT FILE
   my $errorcode;                   ## RETURNED AS 0 IF THERE IS NO ERROR.
   my $errormsg;                    ## RETURNED AS 'none' IF THERE IS NO ERROR.
   my $new_cigarfile;               ## CIGAR FILE WITH SEQUENCES WITH NON-UNIQUE DISPLAY_IDs REMOVED
   my $new_cigarfile2;              ## CIGAR FILE WITH DISPLAY_IDs REPLACED BY SEQUENCE IDs 
   my $new_cigarfile3;              ## CIGAR FILE WITH THE SEQUENCES ADDED 

   # REPLACE THE DISPLAY_IDs IN THE CIGAR-FORMAT ALIGNMENT FILE WITH THE SEQUENCE IDENTIFIERS:
   ($new_cigarfile,$errorcode,$errormsg) = &replace_display_ids_with_seq_ids($cigarfile,$version,$alntype,$family,$outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }

   # WRITE THE SEQUENCES IN THE CIGAR-FORMAT ALIGNMENT FILE: 
   ($new_cigarfile2,$errorcode,$errormsg) = &add_sequences_to_cigarfile($new_cigarfile,$version,$outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   system "rm -f $new_cigarfile";

   # JUST TAKE SEQUENCES WITH UNIQUE DISPLAY_IDs IN THE CIGAR-FORMAT ALIGNMENT FILE:
   ($new_cigarfile3,$errorcode,$errormsg) = &remove_nonunique_seqs_from_cigarfile($new_cigarfile2); 
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   system "rm -f $new_cigarfile2";

   # CONVERT THE CIGAR-FORMAT ALIGNMENT FILE TO A HMM FILE:
   ($errorcode,$errormsg)  = &make_hmm_from_cigar_alignment($new_cigarfile3,$hmmer_bin,$output,$family,$alns_output,$map_output);  
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   system "rm -f $new_cigarfile3";
}

#------------------------------------------------------------------#

# TEST &add_sequences_to_cigarfile

sub test_add_sequences_to_cigarfile
{
   my $outputdir            = $_[0]; # DIRECTORY FOR WRITING OUTPUT FILES TO
   my $random_number;                # NUMBER TO USE IN TEMPORARY FILE NAMES
   my $cigarfile;                    # THE CIGAR FILE
   my $new_cigarfile;                # NAME OF THE NEW CIGAR FILE
   my $expected_cigarfile;           # FILE CONTAINING EXPECTED CONTENTS OF $new_cigarfile 
   my $differences;                  # DIFFERENCES BETWEEN $new_cigarfile AND $expected_cigarfile
   my $length_differences;           # LENGTH OF $differences 
   my $errorcode;                    # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg;                     # RETURNED AS 'none' IF THERE IS NO ERROR
   my $line;                         #  

   # OPEN A NEW FILE FOR WRITING THE NEW CIGAR-FORMAT ALIGNMENT:
   $random_number           = rand();
   $cigarfile               = $outputdir."/tmp".$random_number;
   open(CIGARFILE,">$cigarfile") || die "ERROR: test_add_sequences_to_cigarfile: cannot open $cigarfile\n";
   print CIGARFILE "BGIOSIFCE007851.1_species BGIOSIFCE007851.1 50D21M1D1M1D8M1D7M3D2M31D4M3D5M38D6M68D11M15D16M2D27M5D16M4D1M1D21M6D4M6D4M10D46M20D75M1D35M3D24M1D27M23D29M158D\n";   
   print CIGARFILE "Bm1_12035A_species Bm1_12035A 12D78M3D89M66D16M11D17M2D27M4D17M4D1M1D21M7D3M4D6M10D17M2D11M1D2M1D12M20D28M2D2M1D23M3D51M4D3M1D20M1D25M28D26M76D28M3D11M40D\n";
   print CIGARFILE "ZK507.6.1_species ZK507.6.1 26D47M1D21M31D12M24D20M67D12M14D17M2D24M8D16M4D23M4D16M10D17M2D11M1D2M1D12M20D28M2D2M1D23M3D51M4D3M1D20M1D25M26D28M76D28M3D27M24D\n";
   close(CIGARFILE);
   # OPEN A FILE CONTAINING THE EXPECTED CONTENTS OF $new_cigarfile:
   $random_number           = rand();
   $expected_cigarfile      = $outputdir."/tmp".$random_number;
   open(EXPECTED,">$expected_cigarfile") || die "ERROR: test_add_sequences_to_cigarfile: cannot open $expected_cigarfile\n";
   print EXPECTED "BGIOSIFCE007851.1_species BGIOSIFCE007851.1 50D21M1D1M1D8M1D7M3D2M31D4M3D5M38D6M68D11M15D16M2D27M5D16M4D1M1D21M6D4M6D4M10D46M20D75M1D35M3D24M1D27M23D29M158D MDSIMEPYVADLLADDITASMVELLSGDGGAAQMDVGVLDAYLRAIGALPAHPAAPGADLAAAAEVESMASNDDTNGNWDTKVDAKVPSAFLPPPPGFPPLPVPALANEPVYAAPVDEGDAIRAFMQQLEWSEQYNGDDDAPAPDDSMASRPQLCAPYDDDIDANLRAMEKDAAERPSPDYLDTVHNGQISAASRASLVAWMGRLTHRYELAAGTLHRAVSYFDRFLSARALPSYTEHQLSLVGATAVYTAAKYEDQGTVFKLDAREIASYGEFASAQEVLAMEREMMAALGYRLGGPNAETFVEHFTRYSKGKEELRVQRLARHIADRSLESYGCLGYLPSVVAAAVISIARWTLNPPGALPWSSELHELTGYSSQHISSCVLTVLNTQ\n";   
   print EXPECTED "Bm1_12035A_species Bm1_12035A 12D78M3D89M66D16M11D17M2D27M4D17M4D1M1D21M7D3M4D6M10D17M2D11M1D2M1D12M20D28M2D2M1D23M3D51M4D3M1D20M1D25M28D26M76D28M3D11M40D MNGASDENAENRFGEENEENEQFLRRKGQRKSRPHALAYRNTNISGSVSKKYEEKGVYISIIEKSLLLDKKRRRESFQVLFIFLCQRNQEPPTRSRWLPTGISAFTVFCDEDNGQEVGSSKDPEKNDQQQGEGSKRSGEMSTQYLQIPRKIRRPLATLSPNRDDSDSESNRPLARNLSPNHEESDLESSTPSDANLQTDEGEDSYHSRNSTTSIDSSDISANELRPDDEVVFVVDSVSNNVRRAFKSVPDISSQAATTSSSIIYTTEEHEDRVRTDPVYDADIYLYMRYRELKLCPSGDFLEFQEEICGEVRYLLVEWICDSAREFNLSTESLHMAVSIVDRVLNTLQCPREKLQLLGAAALLLAAKFEEIYPPEVKEFARITAYTFHETEIIRMERIVFARVEYQLLAPTSWWFASRLVRMAHVPRIIYCMMRYLLELALLDHTFLDFRPSVIGAASFFLSNVIFKSDYLLLTTETDIGVNELRSPARKMLELFHEVPNKDYCSVFEKYATEKYEQVSCLLLPDNFSELDVTR\n";
   print EXPECTED "ZK507.6.1_species ZK507.6.1 26D47M1D21M31D12M24D20M67D12M14D17M2D24M8D16M4D23M4D16M10D17M2D11M1D2M1D12M20D28M2D2M1D23M3D51M4D3M1D20M1D25M26D28M76D28M3D27M24D MRSALSLKPSNGNAAKSQAVNNKNVIKNAPLGGKLTRQIGTSNLLQQALPSKKIDESPIIKIDAKDSFKVFEDQEPEKENSSENVDATEKDSNVIPAEDNNMIHELERKMEEKSRAEKLKFKFMQTRDNSDITSRFSEPPSEFSVLCDDDDCDKVSVASSTFTTSVRATFSSFHFDENQRKKEFGKEEAVKKIQKKAAKEARDDSMFSSEEFFPDIIKYMLHRQTKNRASHECFDIQSQVNEEMRTILIDWFSDVVKEYNFQKETFHLAVSLVDRALSMFNIDKMRFQLVGTTSMMIAVKYEEIFPPEIEDFALITDNTYRVPDILLMERFLLGKFDFVVAMPTSSWFGTCFAKRMNFTKKMRNTVHYLLELSLIDVHFLRYRPSDIAAAACCFANLQADVESWPQKMVDDTGISTEDFVDVLRDLHRMYLNASTADFKSIFYNYSETAQMEVALLPAPTDKLRSMFPSIFVTAPKSSNDSSSPQ\n";
   close(EXPECTED);
   # RUN &add_sequences_to_cigarfile:
   ($new_cigarfile,$errorcode,$errormsg) = &add_sequences_to_cigarfile($cigarfile,"7",$outputdir);
   $differences              = "";
   if ($new_cigarfile eq '') { print STDERR "ERROR: test_add_sequences_to_cigarfile: failed test1 (new_cigarfile $new_cigarfile)\n"; exit;}
   open(TEMP,"diff $new_cigarfile $expected_cigarfile |");
   while(<TEMP>)
   {
      $line                  = $_;
      $differences           = $differences.$line;
   }
   close(TEMP);  
   $length_differences       = length($differences);
   if ($length_differences != 0 || $errorcode != 0) { print STDERR "ERROR: test_add_sequences_to_cigarfile: failed test1\n"; exit;}
   # DELETE TEMPORARY FILES:
   system "rm -f $new_cigarfile";
   system "rm -f $expected_cigarfile";
   system "rm -f $cigarfile";
}

#------------------------------------------------------------------#

# WRITE THE SEQUENCES IN THE CIGAR-FORMAT ALIGNMENT FILE: 

sub add_sequences_to_cigarfile
{
   my $cigarfile            = $_[0]; ## THE CIGAR FILE
   my $version              = $_[1]; ## VERSION OF THE TREEFAM DATABASE TO USE
   my $outputdir            = $_[2]; ## THE DIRECTORY TO WRITE OUTPUT FILES TO
   my $new_cigarfile;                ## NEW CIGAR FILE WITH THE SEQUENCES OF INTEREST
   my $id;                           ## SEQUENCE IDENTIFIER FOR THE DISPLAY ID
   my $errorcode            = 0;     ## RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg             = "none";## RETURNED AS 'none' IF THERE IS NO ERROR
   my @temp;                         ##  
   my $line;                         ##  
   my $st;                           ##  
   my $sth;                          ##  
   my $rv;                           ##  
   my @array;                        ##  
   my $database;                     ## THE TREEFAM DATABASE TO USE
   my $dbh;                          ##  
   my $sequence;                     ## THE SEQUENCE FOR THE IDENTIFIER
   my $cigar;                        ## THE CIGAR FOR THE IDENTIFIER 
   my $table_w;                      ## THE DATABASE TABLE TO USE
   my $random_number;                ## RANDOM NUMBER TO USE IN TEMPORARY FILE NAME
   my $display_id;                   ## THE DISPLAY_ID 

   # OPEN A NEW FILE FOR WRITING THE NEW CIGAR-FORMAT ALIGNMENT:
   $random_number           = rand();
   $new_cigarfile           = $outputdir."/tmp".$random_number;
   open(NEW_CIGARFILE,">$new_cigarfile") || die "ERROR: add_sequences_to_cigarfile: cannot open $new_cigarfile\n";

   # CONNECT TO THE DATABASE AND GET ALL FAMILIES OF THIS TYPE:
   $database                = "treefam_".$version;
   $dbh                     = DBI->connect("dbi:mysql:$database:db.treefam.org:3308", 'anonymous', '') || return;

   # GET THE ALIGNMENT FOR EACH FAMILY FROM THE DATABASE:
   $table_w                 = "aa_seq";
   open(CIGARFILE,"$cigarfile") || die "ERROR: cannot open $cigarfile\n";
   while(<CIGARFILE>)
   {
      $line                 = $_;
      chomp $line;
      @temp                 = split(/\s+/,$line);
      $display_id           = $temp[0];
      $id                   = $temp[1];
      $cigar                = $temp[2];
      $st                   = "SELECT SEQ from $table_w WHERE ID=\'$id\'"; 
      $sth                  = $dbh->prepare($st) or die "Cannot prepare $st: $dbh->errstr\n";
      $rv                   = $sth->execute or die "Cannot execute the query: $sth->errstr";
      $sequence             = "none";
      if ($rv >= 1)
      {
         while ((@array) = $sth->fetchrow_array)
   {
	    $sequence       = $array[0];
         }
      }
      print NEW_CIGARFILE "$display_id $id $cigar $sequence\n";
   }
   close(CIGARFILE);
   close(NEW_CIGARFILE);
   $dbh->disconnect();

   return($new_cigarfile,$errorcode,$errormsg);

}

#------------------------------------------------------------------#

# REPLACE THE DISPLAY_IDs IN THE CIGAR-FORMAT ALIGNMENT FILE WITH THE SEQUENCE IDENTIFIERS:
# NOTE: THIS IS HARD TO TEST AS IT EXPECTS THAT $new_idfile AND $new_cigarfile WILL HAVE THE SAME EXACT IDENTIFIERS,
# AND I WOULD NEED A VERY SMALL FAMILY TO USE AS A TEST CASE.
# SUBROUTINE SYNOPSIS: replace_display_ids_with_seq_ids(): replaces display IDs in a cigar-format alignment for a TreeFam family with sequence IDs

sub replace_display_ids_with_seq_ids
{
   my $cigarfile            = $_[0]; ## THE CIGAR FILE
   my $version              = $_[1]; ## VERSION OF THE TREEFAM DATABASE TO USE
   my $aln_type             = $_[2]; ## TYPE OF ALIGNMENT TO USE
   my $family               = $_[3]; ## TREEFAM FAMILY 
   my $outputdir            = $_[4]; ## DIRECTORY TO WRITE OUTPUT FILES TO
   my $errorcode            = 0;     ## RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg             = "none";## RETURNED AS 'none' IF THERE IS NO ERROR
   my $new_cigarfile        = "none";## THE NEW CIGAR FILE, SORTED BY DISPLAY_ID
   my $new_cigarfile2       = "none";## THE NEW CIGAR FILE, WITH ID INSTEAD OF DISPLAY_ID
   my $random_number;                ## RANDOM NUMBER TO USE IN TEMPORARY FILE NAMES 
   my $idfile;                       ## TEMPORARY FILE FOR PUTTING THE IDENTIFIERS FOR THE FAMILY
   my $new_idfile;                   ## A SORTED VERSION OF $idfile 
   my $new_idfile_lineno;            ## NUMBER OF LINES IN $new_idfile
   my $new_cigarfile_lineno;         ## NUMBER OF LINES IN $new_cigarfile
   my $cmd;                          ## COMMAND TO RUN 

   # SORT THE CIGAR FILE BY THE DISPLAY_ID:
   $random_number           = rand();
   $new_cigarfile           = $outputdir."/tmp".$random_number;
   system "sort -k1,1 $cigarfile > $new_cigarfile";

   # GET THE IDENTIFIERS FOR THE FAMILY IN A TEMPORARY FILE:
   ($idfile,$errorcode,$errormsg) = &get_ids_for_family($family,$version,$aln_type,$outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }

   # SORT THE IDENTIFIERS FILE BY DISPLAY_ID:
   $random_number           = rand();
   $new_idfile              = $outputdir."/tmp".$random_number;
   system "sort -k1,1 $idfile > $new_idfile";
   # WE CAN NOW DELETE $idfile:
   system "rm -f $idfile";
   
   # CHECK THAT $new_cigarfile AND $new_idfile HAVE THE SAME NUMBER OF LINES:
   open(TMP,"wc -l $new_idfile | cut -d\" \" -f1 |");
   while(<TMP>)
   {
      $new_idfile_lineno    = $_;
      chomp $new_idfile_lineno;
   }
   close(TMP); 
   open(TMP,"wc -l $new_cigarfile | cut -d\" \" -f1 |");
   while(<TMP>)
   {
      $new_cigarfile_lineno = $_;
      chomp $new_cigarfile_lineno;
   }
   close(TMP);
   if ($new_idfile_lineno != $new_cigarfile_lineno)
   {
      $errormsg             = "ERROR: replace_display_ids_with_seq_ids: new_idfile_lineno $new_idfile_lineno ($new_idfile) new_cigarfile_lineno $new_cigarfile_lineno ($new_cigarfile)\n";
      $errorcode            = 2; # ERRORCODE=2
      return($new_cigarfile2,$errorcode,$errormsg);
   }

   # MAKE A NEW CIGAR FILE WITH THE DISPLAY_ID, ID COLUMN FROM $new_idfile_lineno AND THE CIGAR COLUMN FROM $new_cigarfile_lineno:
   $random_number           = rand();
   $new_cigarfile2          = $outputdir."/tmp".$random_number;
   $cmd                     = "paste $new_idfile $new_cigarfile | awk '{print \$1,\$2,\$4}' > $new_cigarfile2";
   system "$cmd";
   # WE CAN NOW DELETE $new_cigarfile:
   system "rm -f $new_cigarfile";
   system "rm -f $new_idfile";

   return($new_cigarfile2,$errorcode,$errormsg);
}

#------------------------------------------------------------------#

# TEST &remove_nonunique_seqs_from_cigarfile

sub test_remove_nonunique_seqs_from_cigarfile
{
   my $outputdir            = $_[0]; # DIRECTORY FOR WRITING OUTPUT FILES TO
   my $random_number;                # NUMBER TO USE IN TEMPORARY FILE NAMES
   my $cigarfile;                    # THE CIGAR FILE
   my $new_cigarfile;                # NAME OF THE NEW CIGAR FILE
   my $expected_cigarfile;           # FILE CONTAINING EXPECTED CONTENTS OF $new_cigarfile 
   my $differences;                  # DIFFERENCES BETWEEN $new_cigarfile AND $expected_cigarfile
   my $length_differences;           # LENGTH OF $differences 
   my $errorcode;                    # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg;                     # RETURNED AS 'none' IF THERE IS NO ERROR
   my $line;                         #  

   # OPEN A NEW FILE FOR WRITING THE NEW CIGAR-FORMAT ALIGNMENT:
   $random_number           = rand();
   $cigarfile               = $outputdir."/tmp".$random_number;
   open(CIGARFILE,">$cigarfile") || die "ERROR: test_add_sequences_to_cigarfile: cannot open $cigarfile\n";
   print CIGARFILE "BGIOSIFCE007851.1_species BGIOSIFCE007851.1 50D21M1D1M1D8M1D7M3D2M31D4M3D5M38D6M68D11M15D16M2D27M5D16M4D1M1D21M6D4M6D4M10D46M20D75M1D35M3D24M1D27M23D29M158D AAA\n";   
   print CIGARFILE "Bm1_12035A_species Bm1_12035A 12D78M3D89M66D16M11D17M2D27M4D17M4D1M1D21M7D3M4D6M10D17M2D11M1D2M1D12M20D28M2D2M1D23M3D51M4D3M1D20M1D25M28D26M76D28M3D11M40D GGG\n";
   print CIGARFILE "ZK507.6.1_species ZK507.6.1 26D47M1D21M31D12M24D20M67D12M14D17M2D24M8D16M4D23M4D16M10D17M2D11M1D2M1D12M20D28M2D2M1D23M3D51M4D3M1D20M1D25M26D28M76D28M3D27M24D CCC\n";
   print CIGARFILE "Bm1_12035A_species Bm1_12035A 12D78M3D89M66D16M11D17M2D27M4D17M4D1M1D21M7D3M4D6M10D17M2D11M1D2M1D12M20D28M2D2M1D23M3D51M4D3M1D20M1D25M28D26M76D28M3D11M40D TTT\n";
   close(CIGARFILE);
   # OPEN A FILE CONTAINING THE EXPECTED CONTENTS OF $new_cigarfile:
   $random_number           = rand();
   $expected_cigarfile      = $outputdir."/tmp".$random_number;
   open(EXPECTED,">$expected_cigarfile") || die "ERROR: test_add_sequences_to_cigarfile: cannot open $expected_cigarfile\n";
   print EXPECTED "BGIOSIFCE007851.1_species BGIOSIFCE007851.1 50D21M1D1M1D8M1D7M3D2M31D4M3D5M38D6M68D11M15D16M2D27M5D16M4D1M1D21M6D4M6D4M10D46M20D75M1D35M3D24M1D27M23D29M158D AAA\n";   
   print EXPECTED "ZK507.6.1_species ZK507.6.1 26D47M1D21M31D12M24D20M67D12M14D17M2D24M8D16M4D23M4D16M10D17M2D11M1D2M1D12M20D28M2D2M1D23M3D51M4D3M1D20M1D25M26D28M76D28M3D27M24D CCC\n";
   close(EXPECTED);
   # RUN &remove_nonunique_seqs_from_cigarfile:
   ($new_cigarfile,$errorcode,$errormsg) = &remove_nonunique_seqs_from_cigarfile($cigarfile);
   $differences              = "";
   open(TEMP,"diff $new_cigarfile $expected_cigarfile |");
   while(<TEMP>)
   {
      $line                  = $_;
      $differences           = $differences.$line;
   }
   close(TEMP);  
   $length_differences       = length($differences);
   if ($length_differences != 0 || $errorcode != 0) { print STDERR "ERROR: test_remove_nonunique_seqs_from_cigarfile: failed test1\n"; exit;}
   # DELETE TEMPORARY FILES:
   system "rm -f $new_cigarfile";
   system "rm -f $expected_cigarfile";
   system "rm -f $cigarfile";

   # RUN &remove_nonunique_seqs_from_cigarfile:
   ($new_cigarfile,$errorcode,$errormsg) = &remove_nonunique_seqs_from_cigarfile("");
   if ($errorcode != 3) { print STDERR "ERROR: test_remove_nonunique_seqs_from_cigarfile: failed test2\n"; exit;}
   
}

#------------------------------------------------------------------#

# JUST TAKE SEQUENCES WITH UNIQUE DISPLAY_IDs IN THE CIGAR-FORMAT ALIGNMENT FILE:

sub remove_nonunique_seqs_from_cigarfile
{
   my $cigarfile            = $_[0]; ## THE FILE WITH THE CIGAR-FORMAT ALIGNMENT
   my %CNT                  = ();    ## COUNT OF TIMES A SEQUENCE APPEARS IN THE CIGAR-FORMAT FILE 
   my $display_id;                   ## DISPLAY_ID FOR A SEQUENCE IN THE CIGAR-FORMAT FILE
   my $line;                         ##  
   my @temp;                         ##  
   my $new_cigarfile        = "none";## THE NEW CIGAR FILE
   my $random_number;                ## RANDOM NUMBER TO USE IN TEMPORARY FILE NAMES 
   my $errorcode            = 0;     ## RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg             = "none";## RETURNED AS 'none' IF THERE IS NO ERROR 
   my $id;                           ## IDENTIFIER FOR A SEQUENCE 

   # FIND OUT HOW MANY TIMES EACH DISPLAY_ID APPEARS IN $cigarfile:
   if ($cigarfile eq '') 
   { 
      $errormsg             = "ERROR: remove_nonunique_seqs_from_cigarfile: cigarfile $cigarfile\n";
      $errorcode            = 3; # ERRORCODE=3
      return($new_cigarfile,$errorcode,$errormsg);
   }
   open(CIGARFILE,"$cigarfile") || die "ERROR: cannot open $cigarfile\n";
   while(<CIGARFILE>)
   {
      $line                 = $_;
      @temp                 = split(/\s+/,$line);
      if ($#temp == 3)
      {
         $display_id        = $temp[0];
         if (!($CNT{$display_id})) { $CNT{$display_id} = 1;}
         else                      { $CNT{$display_id}++;  }
      }
   }
   close(CIGARFILE);

   # WRITE A NEW FILE WITH JUST THOSE DISPLAY_IDs THAT APPEARED ONCE IN $cigarfile:
   $random_number           = rand();
   $new_cigarfile           = $outputdir."/tmp".$random_number;
   open(NEW_CIGARFILE,">$new_cigarfile") || die "ERROR: remove_nonunique_seqs_from_cigarfile: cannot open $new_cigarfile\n";
   open(CIGARFILE,"$cigarfile") || die "ERROR: cannot open $cigarfile\n";
   while(<CIGARFILE>)
   {
      $line                 = $_;
      @temp                 = split(/\s+/,$line);
      if ($#temp == 3)
      {
         $display_id        = $temp[0]; 
         if (!($CNT{$display_id}))
         {
            $errormsg       = "ERROR: remove_nonunique_seqs_from_cigarfile: do not now count for $display_id\n";
            $errorcode      = 1; # ERRORCODE=1
            return($new_cigarfile,$errorcode,$errormsg);
         }
         if ($CNT{$display_id} == 1) 
         { 
            print NEW_CIGARFILE "$line";
         } 
      }
      # IGNORE THE FIRST LINE OF $cigarfile (FAMILY NAME) AND LAST LINE (JUST SAYS '#END')
   }
   close(CIGARFILE);
   close(NEW_CIGARFILE);

   return($new_cigarfile,$errorcode,$errormsg);
}

#------------------------------------------------------------------#

# READ IN THE CIGAR-FORMAT ALIGNMENTS FOR THE FAMILIES:

sub make_hmm_from_cigar_alignment
{
   my $cigarfile            = $_[0]; ## THE FILE WITH THE CIGAR-FORMAT ALIGNMENTS
   my $hmmer_bin            = $_[1]; ## DIRECTORY CONTAINING THE HMMER EXECUTABLES
   my $output               = $_[2]; ## THE OUTPUT FILE  
   my $family               = $_[3]; ## THE FAMILY
   my $mfafile              = $_[4]; ## MULTIPLE SEQUENCE ALIGNMENT FILE 
   my $map_output           = $_[5]; ## hmmbuild MAP OUTPUT FILE
   my $errorcode            = 0;     ## RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg             = "none";## RETURNED AS 'none' IF THERE IS NO ERROR
   my $line;                         ##  
   my @temp;                         ##    
   my $cigar;                        ## CIGAR STRING FOR A SEQUENCE
   my $mfa_line_length      = "none";## LENGTH OF THE ENTRIES IN THE MULTIPLE ALIGNMENT 
   my $sequence;                     ## SEQUENCE FOR A IDENTIFIER
   my $id;                           ## IDENTIFIER
   my $fasta;                        ## FASTA FORMAT ALIGNMENT FOR SEQUENCE 
   my $random_number;                ## RANDOM NUMBER FOR TEMPORARY FILE NAMES 
   my $display_id;                   ## DISPLAY_ID FOR A SEQUENCE 
   my $core_species_genes   = 0;     ## SAYS WHETHER GENES ARE TAKEN FROM THE CORE SPECIES 

   # OPEN A FILE FOR STORING THE FASTA-FORMAT MULTIPLE ALIGNMENT: 
   open(MFAFILE,">$mfafile") || die "ERROR: make_hmm_from_cigar_alignment: cannot open $mfafile\n";

   # READ IN THE CIGAR-FORMAT ALIGNMENTS:
   open(CIGARS,"$cigarfile") || die "ERROR: read_cigar_alignments: cannot open $cigarfile\n";
   while(<CIGARS>)
   {
      $line                 = $_;
      chomp $line;
      @temp                 = split(/\s+/,$line);
      if ($#temp != 3) 
      {
         $errormsg          = "ERROR: make_hmm_from_cigar_alignment: temp $#temp - line $line in $cigarfile\n";
         $errorcode         = 4;
         return($errorcode,$errormsg);
      }
      $display_id           = $temp[0];
      $id                   = $temp[1]; 
      $cigar                = $temp[2];
      $sequence             = $temp[3];
      # CHECK IF THIS IS ONE OF THE TREEFAM 'CLEAN' SPECIES:
      if ($display_id =~ /_CAEEL/ || $display_id =~ /_CAEBR/ || $display_id =~ /_BRUMA/ || 
          $display_id =~ /_DROME/ || $display_id =~ /_DROPE/ || $display_id =~ /_ANOGA/ || $display_id =~ /_AEDAE/ ||
          $display_id =~ /_SCHMA/ || $display_id =~ /_NEMVE/ || $display_id =~ /_CIOIN/ ||
          $display_id =~ /_HUMAN/ || $display_id =~ /_PANTR/ || $display_id =~ /_MACMU/ ||
          $display_id =~ /_RAT/   || $display_id =~ /_MOUSE/ || $display_id =~ /_BOVIN/ || $display_id =~ /_CANFA/ ||
          $display_id =~ /_MONDO/ || $display_id =~ /_CHICK/ || $display_id =~ /_BRARE/ || $display_id =~ /_XENTR/ ||
          $display_id =~ /_TETNG/ || $display_id =~ /_GASAC/ || $display_id =~ /_YEAST/ || $display_id =~ /_SCHPO/ ||
          $display_id =~ /_ARATH/ || $display_id =~ /_STRPU/ || $display_id =~ /_DICDI/)
      {
         $core_species_genes = 1;
         # CONVERT THE CIGAR FOR $seq TO A FASTA FORMAT ALIGNMENT:
         ($fasta,$errorcode,$errormsg) = &get_fasta_aln_for_seq($sequence,$cigar);
         if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
         if ($fasta eq 'none') 
         {
            print STDERR "WARNING: $family failed (cigar does not match sequence)\n";
            return($errorcode,$errormsg);
         }
         if ($mfa_line_length eq 'none') { $mfa_line_length = length($fasta);}
         else
         {
            if ($mfa_line_length != length($fasta))
            {
               $errormsg    = "ERROR: read_cigar_alignments: length of $fasta is not expected length for family $family\n";
               $errorcode   = 15; # ERRORCODE=15
               return($errorcode,$errormsg);
            }
         }
         print MFAFILE ">$id\n";
         print MFAFILE "$fasta\n";
      }
   }
   close(CIGARS); 
   close(MFAFILE);
   if ($core_species_genes == 0)
   {
      print STDERR "WARNING: $family failed: no core species genes for multiple alignment\n";
   }
   else 
   {
      # MAKE A HMM FOR THIS FAMILY:
      ($errorcode,$errormsg) = &make_hmm($mfafile,$hmmer_bin,$output,$map_output);
      if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   }

   return($errorcode,$errormsg);
}

#------------------------------------------------------------------#

# TEST &make_hmm

sub test_make_hmm
{
   my $outputdir            = $_[0]; ## DIRECTORY FOR WRITING OUTPUT FILES TO 
   my $errorcode;                    ## RETURNED AS 0 BY A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                     ## RETURNED AS 'none' BY A FUNCTION IF THERE IS NO ERROR
   my $random_number;                ## RANDOM NUMBER FOR TEMPORARY FILE NAME
   my $mfafile;                      ## FILE FOR WRITING MULTIPLE ALIGNMENT TO 
   my $hmmfile;                      ## FILE FOR WRITING THE HMM TO 

   $random_number           = rand();
   $mfafile                 = $outputdir."/tmp".$random_number;
   $random_number           = rand();
   $hmmfile                 = $outputdir."/tmp".$random_number;
   open(MFAFILE,">$mfafile") || die "ERROR: test_make_hmm: cannot open $mfafile\n";
   print MFAFILE ">seq1\nACGTACC\n";
   print MFAFILE ">seq2\n-CCTACC\n";
   print MFAFILE ">seq3\nAC-TACC\n";
   close(MFAFILE);
   ($errorcode,$errormsg) = &make_hmm($mfafile,"/software/hmmer/bin/",$hmmfile);
   if ($errorcode != 0) { print STDERR "ERROR: test_make_hmm: failed test1\n"; exit;}
   system "rm -f $mfafile";
   system "rm -f $hmmfile";
}

#------------------------------------------------------------------#

# MAKE A HMM FOR THIS FAMILY:

sub make_hmm
{
   my $mfafile              = $_[0]; ## FILE CONTAINING THE MULTIPLE ALIGNMENT FOR THE FAMILY
   my $hmmer_bin            = $_[1]; ## DIRECTORY CONTAINING THE HMMER EXECUTABLES
   my $output               = $_[2]; ## NAME OF THE OUTPUT HMM FILE 
   my $map_output           = $_[3]; ## NAME OF THE hmmbuild MAP OUTPUT FILE
   my $errorcode            = 0;     ## RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg             = "none";## RETURNED AS 'none' IF THERE IS NO ERROR
   my $cmd;                          ## COMMAND TO RUN 

   # RUN HMMER2 TO MAKE THE HMM:
   
   $cmd                      = $hmmer_bin."/hmmbuild -o $map_output $output $mfafile >& /dev/null";
   system "$cmd";
   $cmd                      = $hmmer_bin."/hmmcalibrate $output >& /dev/null";
   system "$cmd";

   return($errorcode,$errormsg); 
}

#------------------------------------------------------------------#

# TEST &get_fasta_aln_for_seq 

sub test_get_fasta_aln_for_seq
{
   my $errorcode;                    # RETURNED AS 0 BY A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                     # RETURNED AS 'none' BY A FUNCTION IF THERE IS NO ERROR
   my $fasta;                        # FASTA FORMAT ALIGNMENT FOR A SEQUENCE

   ($fasta,$errorcode,$errormsg) = &get_fasta_aln_for_seq("ABCDEF","3D6M");
   if ($fasta ne '---ABCDEF' || $errorcode != 0) { print STDERR "ERROR: test_get_fasta_aln_for_seq: failed test1\n";}

   ($fasta,$errorcode,$errormsg) = &get_fasta_aln_for_seq("ABCDEF","1D6M1D");
   if ($fasta ne '-ABCDEF-' || $errorcode != 0) { print STDERR "ERROR: test_get_fasta_aln_for_seq: failed test2\n";}

   ($fasta,$errorcode,$errormsg) = &get_fasta_aln_for_seq("ABCDEF","3M3D3M");
   if ($fasta ne 'ABC---DEF' || $errorcode != 0) { print STDERR "ERROR: test_get_fasta_aln_for_seq: failed test3\n";}

   ($fasta,$errorcode,$errormsg) = &get_fasta_aln_for_seq("ABCDEF","3M3I3M");
   if ($errorcode != 14) { print STDERR "ERROR: test_get_fasta_aln_for_seq: failed test4\n";}

   ($fasta,$errorcode,$errormsg) = &get_fasta_aln_for_seq("ABCDEF","3M");
   if ($errorcode != 0 || $fasta ne 'none') { print STDERR "ERROR: test_get_fasta_aln_for_seq: failed test5\n";}

   ($fasta,$errorcode,$errormsg) = &get_fasta_aln_for_seq("ABCDEF","20M");
   if ($errorcode != 0 || $fasta ne 'none') { print STDERR "ERROR: test_get_fasta_aln_for_seq: failed test6\n";}

}

#------------------------------------------------------------------#

# CONVERT THE CIGAR FOR A SEQUENCE TO A FASTA FORMAT ALIGNMENT:
# SUBROUTINE SYNOPSIS: get_fasta_aln_for_seq(): converts a cigar-format alignment for a TreeFam family into fasta-format

sub get_fasta_aln_for_seq
{
   my $sequence             = $_[0]; # PROTEIN SEQUENCE FOR A TREEFAM GENE
   my $cigar                = $_[1]; # CIGAR FORMAT ALIGNMENT FOR THE SEQUENCE
   my $errorcode            = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg             = "none";# RETURNED AS "none" IF THERE IS NO ERROR
   my $fasta                = "";    # FASTA FORMAT ALIGNMENT FOR THE SEQUENCE
   my @tmp;                          #   
   my @tmp2;                         #  
   my $i;                            # 
   my $cigar_type;                   # TYPE OF ELEMENT IN THE CIGAR STRING
   my $cigar_len;                    # LENGTH OF ELEMENT IN THE CIGAR STRING
   my $j;                            # 
   my $sequence_index;               # INDEX IN THE SEQUENCE $sequence
   my $aa;                           # AN AMINO ACID IN SEQUENCE $sequence 
   my $sequence_length;              # LENGTH OF $sequence 
   my $last_letter;                  # LAST LETTER OF $sequence 
   
   # REMOVE * FROM THE END OF $sequence:
   if ($sequence eq '')
   {
      $errormsg             = "ERROR: get_fasta_aln_for_seq: sequence $sequence\n";
      $errorcode            = 3;
      return($fasta,$errorcode,$errormsg);
   }
   $last_letter             = substr($sequence,length($sequence)-1,1);
   if ($last_letter eq '*') { chop($sequence);}
   # REMOVE . FROM THE END OF $sequence (OCCURS FOR SOME SPECIES):
   $last_letter             = substr($sequence,length($sequence)-1,1);
   if ($last_letter eq '.') { chop($sequence);}
   @tmp		            = split(/\d+/,$cigar);
   @tmp2		    = split(/\D+/,$cigar);
   # LOOK AT EACH ELEMENT IN THE CIGAR STRING
   $sequence_index          = -1;
   $sequence_length         = length($sequence);
   for ($i = 0; $i <= $#tmp2; $i++)
   {
      $cigar_type           = $tmp[($i+1)];
      $cigar_len            = $tmp2[$i];
      if ($cigar_type eq 'D')
      {
         for ($j = 1; $j <= $cigar_len; $j++) 
         { 
	    $fasta          = $fasta."-";
         }
      }
      elsif ($cigar_type eq 'M')
      {
         for ($j = 1; $j <= $cigar_len; $j++)
         {
            $aa             = "*";
            while($aa eq '*' || $aa eq '.') # IGNORE *s OR .s IN THE SEQUENCE (SOME SEQUENCES HAVE INTERNAL '*' OR '.'):
            { 
               $sequence_index++;
               if ($sequence_index > $sequence_length - 1)
               {
                  $fasta    = "none";
                  if ($PRINT_WARNINGS == 1) { print STDERR "WARNING: get_fasta_aln_for_seq: sequence_length $sequence_length but sequence_index $sequence_index\n";}
                  return($fasta,$errorcode,$errormsg);
               } 
               $aa          = substr($sequence,$sequence_index,1);
            }
            $fasta          = $fasta.$aa;
         } 
      }
      else
      {
         $errormsg          = "ERROR: get_fasta_aln_for_seq: cigar_type $cigar_type\n";
         $errorcode         = 14; # ERRORCODE=14
         return($fasta,$errorcode,$errormsg);
      }
   }
   if ($sequence_index != $sequence_length - 1) # THE CIGAR ALIGNMENT AND SEQUENCE ARE NOT COMPATIBLE:
   {
      $fasta                = 'none';
      if ($PRINT_WARNINGS == 1) { print STDERR "WARNING: get_fasta_aln_for_seq: sequence_length $sequence_length but sequence_index $sequence_index\n"; }
      return($fasta,$errorcode,$errormsg);
   }
 
   return($fasta,$errorcode,$errormsg);
}

#------------------------------------------------------------------#

# GET THE IDENTIFIERS FOR THE SEQUENCES IN THE FAMILY:
# SUBROUTINE SYNOPSIS: get_ids_for_family(): get the identifiers for the sequences in a TreeFam family

sub get_ids_for_family
{
   my $family               = $_[0]; ## TREEFAM FAMILY IDENTIFIER
   my $version              = $_[1]; ## VERSION OF THE TREEFAM DATABASE TO USE
   my $alntype              = $_[2]; ## TYPE OF ALIGNMENTS TO USE (SEED/FULL)
   my $outputdir            = $_[3]; ## DIRECTORY FOR WRITING OUTPUT FILES TO
   my $errorcode            = 0;     ## RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg             = "none";## RETURNED AS 'none' IF THERE IS NO ERROR
   my $database;                     ## DATABASE TO CONNECT TO
   my $dbh;                          ## 
   my $table_w;                      ## DATABASE TABLE THAT WE WANT TO LOOK AT
   my $st;                           ##  
   my $sth;                          ##  
   my $rv;                           ##  
   my @array;                        ##  
   my $id;                           ## TAXONOMY ID.
   my $disp_id;                      ## DISPLAY ID. 
   my $idfile;                       ## TEMPORARY FILE WITH THE IDENTIFIERS FOR A FAMILY 
   my $random_number;                ## RANDOM NUMBER TO USE IN TEMPORARY OUTPUT FILE 

   # SORT THE CIGAR FILE BY THE DISPLAY_ID:
   $random_number           = rand();
   $idfile                  = $outputdir."/tmp".$random_number;
   open(IDFILE,">$idfile") || die "ERROR: get_ids_for_family: cannot open $idfile\n";

   # CONNECT TO THE DATABASE:
   $database                = "treefam_".$version;
   $dbh                     = DBI->connect("dbi:mysql:$database:db.treefam.org:3308", 'anonymous', '') || return;

   # CHECK THAT WE ARE USING EITHER SEED OR FULL ALIGNMENTS:
   if ($alntype ne 'seed' && $alntype ne 'full')
   {
      $errormsg             = "ERROR: get_ids_for_family: alntype $alntype\n";
      $errorcode            = 16; # ERRORCODE=16
      $dbh->disconnect();
      return($idfile,$errorcode,$errormsg); 
   }
 
   # GET THE SEQUENCES FOR THEFAMILY FROM THE DATABASE:
   if ($alntype eq 'full')
   {
      $table_w              = "aa_full_align";
   }
   else
   {
      $table_w              = "aa_seed_align";
   }
   $st                      = "SELECT ID, DISP_ID from $table_w WHERE AC=\'$family\'"; 
   $sth                     = $dbh->prepare($st) or die "Cannot prepare $st: $dbh->errstr\n";
   $rv                      = $sth->execute or die "Cannot execute the query: $sth->errstr";
   if ($rv >= 1)
   {
      while ((@array) = $sth->fetchrow_array)
      {
         $id                = $array[0];
         $disp_id           = $array[1];
         print IDFILE "$disp_id $id\n";
       }
   }
   else
   {
      print STDERR "WARNING: $family failed (no sequences found)\n";
      if ($PRINT_WARNINGS == 1) { print STDERR "WARNING: get_ids_for_family: no sequences found for family $family\n";}
      $dbh->disconnect(); 
      return($idfile,$errorcode,$errormsg);
   }
   $dbh->disconnect();

   return($idfile,$errorcode,$errormsg); 

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
