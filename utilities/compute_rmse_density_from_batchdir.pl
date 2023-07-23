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
my @molec_sys = (["hcn", 64, 2], ["hehhe", 64, 2]);


# Ensure we have directories for each system
for (my $i = 0; $i <= $#molec_sys; $i++) {
    if (! -d $molec_sys[0][0]) { die "No directory: $molec_sys[$i][0]!\n"};
}

my $rmse_val = 0.0;
# Compute RMSE for each system
for (my $i = 0; $i < $#molec_sys; $i++) {

    # Enter directory
    chdir "$molec_sys[$i][0]";
    my $mname = $molec_sys[$i][0];
    my $jtype = $molec_sys[$i][2];
    my $npts  = $molec_sys[$i][1];
    
    system("\$MLBASOPT/bin/cubedata_into_training_data.x den_p_0.cube ${jtype} > ${mname}_dens.data");
    system("\$MLBASOPT/bin/compute_rmse_density.x ${mname}_dens.data ${mname}_dens.fgh.data ${npts} > ${mname}.rmse");
    
    open(FILE, "<", "${mname}.rmse") or die "Could not open file ${mname}.rmse!\n";
    my $rmse = <FILE>;
    chomp($rmse);
    close(FILE);
    
    # Add RMSE to RMSE_VAL
    $rmse_val = $rmse_val + $rmse;

    # Leave directory
    chdir "../";
}

# Read parameters
open(FILE, "<", "var.dat") or die "Could not open file var.dat!\n";
chomp(my @var = <FILE>);
close(FILE);

# Print final output
open(FILE, ">", "rmse.dat") or die "Could not open file to write rmse.dat!\n";
for (my $i = 0; $i <= $#var; $i++) {
    printf FILE " %10.5f", $var[$i];
}
printf FILE " %15.8f\n", $rmse_val;
close(FILE);

 
