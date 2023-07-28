#!/usr/bin/perl

use strict;
use warnings;

use Cwd qw(getcwd);

#
# Program to compute penalty function
#

# F = w_1 * density + w_2 * excited states + w_3 * ZPE
#

# Process command line arguments
my ($w_den, $w_xst, $w_gst) = @ARGV;
if ( ! defined $w_den || ! defined $w_xst || ! defined $w_gst ) {
    die "Please provide weights for function.\n";
}

# Read in each RMSE value
open(FILE, "<", "density.rmse.dat") or die "Could not open density.rmse.dat!\n";
chomp(my @dens_error = <FILE>);
close(FILE);
open(FILE, "<", "excited_states.rmse.dat") or die "Could not open excited_states.rmse.dat!\n";
chomp(my @xst_error = <FILE>);
close(FILE);
open(FILE, "<", "zero_point.rmse.dat") or die "Could not open zero_point.rmse.dat!\n";
chomp(my @gst_error = <FILE>);
close(FILE);

my $rmse_val = $w_den * $dens_error[0] + $w_xst * $xst_error[0] + $w_gst * $gst_error[0];


# Read parameters
open(FILE, "<", "var.dat") or die "Could not open file var.dat!\n";
chomp(my @var = <FILE>);
close(FILE);

# Print final output
open(FILE, ">", "density.rmse.dat") or die "Could not open file to write density.rmse.dat!\n";
for (my $i = 0; $i <= $#var; $i++) {
    printf FILE " %10.5f", $var[$i];
}
printf FILE " %15.8f\n", $rmse_val;
close(FILE);
