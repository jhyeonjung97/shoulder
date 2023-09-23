#!/usr/bin/env perl
#;-*- Perl -*-

use Cwd;

$nbands = `grep NBANDS OUTCAR|cut -c 105-108`;
$lines = $nbands+3;
$lines2 = $nbands+4;
$spin1 = `grep -A $lines 'spin component 1' OUTCAR | tail -n $lines2`;
$spin2 = `grep -A $lines 'spin component 2' OUTCAR | tail -n $lines2`;

open FH, ">", "upspin" or die "$!\n";
print FH $spin1;
close FH;

open FH, ">", "downspin" or die "$!\n";
print FH $spin2;
close FH;
