#! /usr/bin/perl -w

# Runs a bunch of PAML input files through PAML codeml program by writing 
# a new control file (codeml.ctl) for the codeml program for each file. 
# Renames output of each run so that results aren't overwritten.
#
# Note: Input files should be named *.nuc and treefiles should be named *.tree
# and present in directory called NUCs in the directory where codeml is located. 
# Output files are named *.out

##  Updated October 2007 ##


$masterfilename = $ARGV[0];

open (MASTER, $masterfilename);
while (<MASTER>) {

    # grab filename from line
    $input_file = $_;

    # modify information for outputting etc.
    chomp($input_file);
    $trunc_filename = $input_file;
    chop($trunc_filename);
    chop($trunc_filename);
    chop($trunc_filename);

    $output_filename = $trunc_filename . "out";
   # $treefile = $trunc_filename . "tree";
    
    print "$input_file\n";

    # initialize a new PAML control file
    open(NF, ">yn00.ctl");

    print NF "seqfile = $input_file\n";
    print NF "outfile = $output_filename\n";
  #  print NF "noisy = 9\n";
    print NF "verbose = 1\n";
  #  print NF "runmode = 0\n";
  #  print NF "seqtype = 1\n"; 
  #  print NF "CodonFreq = 2\n"; 
  #  print NF "aaDist = 0\n"; 
  #  print NF "aaRatefile = wag.dat\n";
  #  print NF "model = 0\n";   ############# 
  #  print NF "NSsites = 7\n"; ############# Change this part to specify different models (M7, M8, etc.)
    print NF "icode = 0\n";
    print NF "weighting = 0\n";
    print NF "commonf3x4 = 0\n";
  #  print NF "Mgene = 0\n"; 
  #  print NF "fix_kappa = 0\n";
  #  print NF "kappa = 2\n";
  #  print NF "fix_omega = 0\n"; 
  #  print NF "omega = .4\n";
  #  print NF "fix_alpha = 1\n"; 
  #  print NF "alpha = 0.\n";
  #  print NF "Malpha = 0\n"; 
  #  print NF "ncatG = 10\n";  ############ 
  #  print NF "clock = 0\n"; 
  #  print NF "getSE = 1\n"; 
  #  print NF "RateAncestor = 1\n";
  #  print NF "Small_Diff = .5e-6\n";
  #  print NF "method = 1\n";
  #  print NF "cleandata = 1\n";
  #  print NF "treefile = NUCs/$treefile\n";

    close(NF);

    # execute codeml program
    system("yn00");
    print "$input_file Done\n";
}    
