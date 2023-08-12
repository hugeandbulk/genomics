#!/usr/local/bin/perl

$fasta = $ARGV[0];
$gff = $ARGV[1];

# read in the gff file to find which sequences to take:
%TAKE = ();
$num_to_take = 0;
open(GFF,"$gff");
while(<GFF>)
{  
   $line = $_;
   chomp $line;
   if (substr($line,0,1) ne '#')
   {  
      # seq0    LTRharvest      LTR_retrotransposon     227145  239120  .       -       .       ID=LTR_retrotransposon1;Parent=repeat_region
1;ltr_similarity=98.06;seq_number=0
      # seq0    LTRharvest      LTR_retrotransposon     1847908 1862843 .       -       .       ID=LTR_retrotransposon2;Parent=repeat_region
2;ltr_similarity=87.74;seq_number=0
      # seq0    LTRharvest      LTR_retrotransposon     3847365 3852453 .       +       .       ID=LTR_retrotransposon3;Parent=repeat_region
3;ltr_similarity=96.35;seq_number=0
      @temp = split(/\t+/,$line);
      $seq = $temp[0]; # eg. seq0
      $feature = $temp[2];
      if ($feature eq 'LTR_retrotransposon')
      {  
         $start = $temp[3];
         $end = $temp[4];
         $seq = $seq."_".$start."_".$end;
         $TAKE{$seq} = 1;
         $num_to_take++;
      }
   }
}
close(GFF);

# read in the fasta file of sequences, and print out those we want to take:
$take = 0;
$num_found = 0;
open(FASTA,"$fasta");
while(<FASTA>)
{  
   $line = $_;
   chomp $line;
   if (substr($line,0,1) eq '>') # >chromosome:WBcel235:I:1:15072434:1 chromosome I (dbseq-nr 0) [65157,67865]
   {
      @temp = split(/\s+/,$line);
      $seqno = $temp[$#temp-1]; # eg. 0)
      chop($seqno); # eg. 0
      $pos = $temp[$#temp]; # eg. [65157,67865]
      $pos = substr($pos,1,length($pos)-2);
      @temp = split(/\,/,$pos);
      $start = $temp[0];
      $end = $temp[1];
      $seq = "seq".$seqno."_".$start."_".$end;
      if ($TAKE{$seq}) { $take = 1; $num_found++;} else { $take = 0;}
   }    
   if ($take == 1) { print "$line\n";}
}
close(FASTA);
if ($num_found != $num_to_take) { print STDERR "ERROR: num_found $num_found num_to_take $num_to_take\n"; exit;}
print STDERR "FINISHED\n";
