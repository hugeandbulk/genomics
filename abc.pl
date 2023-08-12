#!/usr/bin/perl
if (! -d abc) {
    mkdir "abc", "0777" or die "Can't mkdir abc: $!";
}
