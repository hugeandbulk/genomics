#! perl
use strict;
use warnings;

my (%seqs, $id, $dna);

while (my $line = <>)
{
    chomp $line;

    if ($line =~ / ^ > (.+) /x)
    {
        $seqs{$id} = $dna if defined $id;
        $id        = $1;
        $dna       = '';
    }
    else
    {
        $dna      .= $line;
    }
}

$seqs{$id} = $dna if defined $id;

for my $key (sort { length $seqs{$a} <=>
                    length $seqs{$b} } keys %seqs)
{
    printf "%s:%d\n", $key, length $seqs{$key};
}
