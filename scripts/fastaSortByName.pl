#!/usr/bin/perl -w

# Copyright 2016, Naoki Takebayashi <ntakebayashi@alaska.edu>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# Version: 20190211

my $usage="\nUsage: $0 [-hrg] [fastaFileName1 ...]\n".
    "  -h: help\n".
    "  -q: input is FASTQ instead of FASTA\n".
    "  -r: reverse\n" .
    "  -g: remove gaps '-' from the sequence\n".
    "Sort FASTA sequences alphabetically by names.  If multiple files are \n".
    "given, sequences in all files are marged before sorting.  If no \n".
    "argument is given, it will take STDIN as the input.\n" .
    "Note that the entire sequence label including spaces is used as\n".
    "the name.\n";

# Changelog
# Version 20190211
# - STDIN as input was broken, so it is fixed
# - renamed from fastaSortByName.pl fastxSortByName.pl
# Version 20160717
# - New feature: -q option to work with fastq files

our($opt_h, $opt_g, $opt_r, $opt_q);

use Bio::SeqIO;

use Getopt::Std;
getopts('hgrq') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

die "$usage\nERROR: -g with fastq file is not implemented\n" 
	if (defined ($opt_g) && defined ($opt_q));
my $format = (defined $opt_q) ? "fastq" : "fasta";
my @seqArr = ();

@ARGV = ('-') unless @ARGV;
while (my $file = shift) {
    my $seqio_obj = ($file eq '-') ?
	Bio::SeqIO->new(-fh => \*STDIN, -format => $format) :
	Bio::SeqIO->new(-file => $file, -format => $format);
    while (my $seq = $seqio_obj->next_seq()) {
       # need to deal with spaces
	$seq->desc( $seq->id . " ". $seq->desc);

	push(@seqArr, $seq);
    }
}

if (defined($opt_r)) {
    @seqArr = sort { - ($a->desc() cmp $b->desc()) } @seqArr;
} else {
    @seqArr = sort { $a->desc() cmp $b->desc() } @seqArr;
}

my $seqOut = Bio::SeqIO->new(-fs => \*STDOUT, -format => $format);
foreach my $s (@seqArr) {
    # prints "id desc", and desc was modified, returning it to original
    my $thisDesc = $s->desc;
    $thisDesc =~ s/^\S+ //; # remove the first word.
    $s->desc($thisDesc);

    if(defined($opt_g)) {
	my $tmp = $s->seq();
	$tmp =~ s/-//g;
	$s->seq($tmp);
    }
    $seqOut->write_seq($s);
}

exit;
