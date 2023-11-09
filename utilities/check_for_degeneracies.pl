#!/usr/bin/perl

use strict;
use warnings;


# Check for degeneracies between two states in batch job directory
#
# Arguments: [system name] [state # 1] [state # 2]
my $usage = "[system name] [state # 1] [state # 2]";

if (scalar @ARGV != 3) {
    die "Usage: $usage\n";
}
my $sysname = $ARGV[0];
my $state1  = $ARGV[1];
my $state2  = $ARGV[2];

# Read in the states file and check the degeneracy of the states of interest
open(FH, "<", "${sysname}/${sysname}_states.data") or die $!;
my @energies = grep(/\n/i, <FH>);
close(FH);
my $ediff = $energies[$state2 - 1] - $energies[$state1 - 1];
if ($ediff > 0.00001) {
    print "States not degenerate! ediff = $ediff\n";
    print @energies;
    print "\n";
}

