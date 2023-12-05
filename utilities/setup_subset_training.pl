#!/usr/bin/perl

use strict;
use warnings;

use Cwd qw(getcwd);
use List::Util qw(shuffle);


# Get command line arguments
my ($fname, $nbatches) = @ARGV;
if ( ! defined $fname || ! defined $nbatches) {
    die "Usage: setup_subset_training.pl [training file] [num batches]\n";
}

# Read in and shuffle training data file
open(FILE, "<", $fname);
my @data = shuffle <FILE>;
close(FILE);

# Print new training data files
my $npts_subset = int((scalar @data) / $nbatches);
for (my $i = 1; $i <= $nbatches; $i++) {
    my $start = ($i - 1) * $npts_subset;
    my $final = $i * $npts_subset - 1;

    if ($i == $nbatches) {$final = $#data;}

    open(FILE, ">", "training.dat_${i}");
    print FILE @data[$start..$final];
    close(FILE);
}



