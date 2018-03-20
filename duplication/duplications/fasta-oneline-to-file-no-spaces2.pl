#!/usr/bin/perl -w

# converts
# >GeneA
# ATGGCCG
# GTCCGTT
# to
# >GeneA
# ATGGCCGGTCCGTT

### AND PRINTS this to a file named after the original file ##########

$file = $ARGV[0];
chomp($file);
$newfilename = $file . ".oneline";

open(NF,">$newfilename");

open(FILE, $file);
while (<FILE>) {
    		
    if (/^>/) {
	$count++;
	if($count==1) {
	    print NF "$_";
	} else {
	    print NF "\n$_";
	}
    }
    if (/^[a-zA-Z\*]/i) {
	s/\s//g;
	chomp($_);
	print NF $_;
    }

 }    


  
