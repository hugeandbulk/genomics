#!/usr/bin/perl

# get the accession numbers from FASTA file

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

while (<>){
    chomp;
    next unless (/^>/);
    s/^>\s*//;

    my @line = split /\s+/;

    my $first = shift (@line);
    my @numbers = split /\|/, $first;
    
    $accNum = $numbers[3];
    $accNum =~ s/\.\d+$//;  # remove version numbers

    print "$accNum\t\t# " . $_ . "\n";
}
