#!/usr/bin/perl

use strict;
use warnings;

use Cwd qw(getcwd);

#
# Program to compute RMSE for a given set of parameters
#

#
# Molecular systems
#
my @molec_sys = (["hcn", 96, 3, 4], ["hehhe", 96, 3, 4], ["fhf", 96, 3, 4]);


# Ensure we have directories for each system
for (my $i = 0; $i <= $#molec_sys; $i++) {
    if (! -d $molec_sys[$i][0]) { print "No directory: $molec_sys[$i][0]!\n"};
}

my $rmse_val = 0.0;
# Compute RMSE for each system
for (my $i = 0; $i <= $#molec_sys; $i++) {

    # Enter directory
    if (-e "$molec_sys[$i][0]/$molec_sys[$i][0].input" ) {
        chdir "$molec_sys[$i][0]";
        my $mname = $molec_sys[$i][0];
        my $jtype = $molec_sys[$i][2];
        my $npts  = $molec_sys[$i][1];
        my $nstate= $molec_sys[$i][3];
        
        # Read in ground states. 
        open(FILE, "<", "${mname}_states.data");
        chomp(my @qce = <FILE>);
        close(FILE);
        open(FILE, "<", "${mname}_states.fgh.data");
        chomp(my @ref = <FILE>);
        close(FILE);
        
        # Get excited states in cm-1
        for (my $i = 1; $i < $nstate; $i++) {
            $qce[$i] = ($qce[$i] - $qce[0]) * 219474.63;
            $ref[$i] = ($ref[$i] - $ref[0]) * 219474.63;
        }
        
        my $rmse = 0.0;
        # Compute RMSE
        for (my $i = 1; $i < $nstate; $i++) {
            my $diff = $qce[$i] - $ref[$i];
            $rmse = $rmse + ($diff * $diff);
        }
        $rmse = $rmse / $nstate;
        $rmse = sqrt($rmse);
        
        # Add RMSE to RMSE_VAL
        $rmse_val = $rmse_val + $rmse;
        
        # Leave directory
        chdir "../";
    }
    #else
    #{
    #    print "WARNING: Directory not found! $molec_sys[$i][0]\n";
    #}
}

# Read parameters
open(FILE, "<", "var.dat") or die "Could not open file var.dat!\n";
chomp(my @var = <FILE>);
close(FILE);

# Print final output
open(FILE, ">", "excited_states.rmse.dat") or die "Could not open file to write excited_states.rmse.dat!\n";
printf FILE " %15.8f\n", $rmse_val;
close(FILE);

 
