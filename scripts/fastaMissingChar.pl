#!/usr/bin/perl -w

# Copyright 2013, Naoki Takebayashi <ntakebayashi@alaska.edu>
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

# Version: 20130612

my $usage="\nUsage: $0 [-h] [-m char] [fastaFileName1 ...]\n".
    "  -h: help\n".
    "  -m: missing character\n".
    "Print out the name of sequences with characters other than ATGC-.\n".
    "If -m is specified, the ambiguous characters are repleced with the\n".
    "specified character.  e.g. -m '?' will place ? to the ambigous characters.\n" .
    "If multiple files are given, sequences in all files are marged.  If no \n".
    "argument is given, it will take STDIN as the input\n";

our($opt_h, $opt_m);

use Bio::SeqIO;

use Getopt::Std;
getopts('hm:') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

my $format = "fasta";
my @seqArr = ();

@ARGV = ('-') unless @ARGV;
while (my $file = shift) {
    my $seqio_obj = Bio::SeqIO->new(-file => $file, -format => $format);
    while (my $seq = $seqio_obj->next_seq()) {
	push(@seqArr, $seq);
    }
}

#@seqArr = sort { $a->id() cmp $b->id() } @seqArr;

foreach my $s (@seqArr) {
    my $thisSeq = $s->seq();
    my $ambig = AmbiguousChar($thisSeq);
    if ($ambig ne "") {
	print STDERR $s->id(), "\t$ambig\n";
	if (defined($opt_m)) {
	    $thisSeq = ReplaceAmbiguousChar($thisSeq, $opt_m);
	    $s->seq($thisSeq);
	}
    }
}

if (defined($opt_m)) {
    my $seqOut = Bio::SeqIO->new(-fs => \*STDOUT, -format => $format);
    foreach my $s (@seqArr) {
	$seqOut->write_seq($s);
    }
}
exit;

sub AmbiguousChar {
    my $string = shift;
    $string =~ s/[ATGC-]//g;

    $string =~ s/\s+//g;
    return $string;
}

sub ReplaceAmbiguousChar {
    my ($string, $char) = @_;
    $string =~ s/[^ATGC-]/$char/g;
    return $string;
}
