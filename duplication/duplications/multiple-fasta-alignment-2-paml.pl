#!/usr/bin/perl -w

# Converts a multiple FASTA alignment to PAML format
#
# A "multiple FASTA alignment"
# >gene1 NNNNNNNNNNNNNNNNNNNNNN
# >gene2 NNNNNNNNNNN---NNNNNNNN
# >gene3 -----NNNNNNNNNNNNNNNNN

# Note: To process many multiple alignments use in conjunction w/ run-a-script.pl

# _________________________________________________________________________________
#
# Copyright (C) 2004 Cristian I. Castillo-Davis
# ccastillo-davis@stat.harvard.edu
#
# If you use this script in published work, please cite the reference below...
# i.e. "using scripts obtained/derived from Castillo-Davis et al. (2004)"
#
# Castillo-Davis, C. I., R. J. Kulathinal, F. A. Kondrashov, D. L. Hartl. 2004.
# The functional genomic distribution of protein divergence in two animal phyla:
# co-evolution, genomic conflict, and constraint. Genome Research. 4(5):802-11.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# _________________________________________________________________________________


# a file with FASTA alignments
$file = $ARGV[0];

#intitialize
$too_short_files = "";
$non_matching_length = "";


$count = 0;

open(FILE, "$file");
while(<FILE>) {

    if (/^\>(.*)\n/) {
	$count++; 
	$genename = $1;
	
	# put name into an array for use later
	$namearray[$count] = $genename;
    }	
    if ((/^(A|T|G|C|N|\-)(A|T|G|C|N|\-)(A|T|G|C|N|\-)(.*)/i)) {
	chomp($_);
	$seq = $_;
	
        # put sequence into a hash with name as key
	$nthash{$genename} = $seq;

    }
}

for($i=1;$i<=$count;$i++) {

    $seqlength = 0;
    
    #grab genename for each seq
    $genename = $namearray[$i];

    #grab seq for each gene
    $seq = $nthash{$genename};

    #chop off last codon because PAML can't deal with stop codons
    #even at end of sequence
    chop($seq);
    chop($seq);
    chop($seq);

    $_ = $seq;
	
    #replace every instance of N with a dash
    tr/(N)/-/;
    
    #count how long the sequence is using transliteration trick
    $seqlength = tr/A|C|T|G|\-//;
    
    print "\nseq length is $seqlength\t";
    
    # Add dashes to end of sequence if not ending in multiple of codons
    $modulus_of_count = ($seqlength % 3);   #divide by three and check mod
    #print "modulus 3 of count is: $modulus_of_count\n";	
    
    if ($modulus_of_count == 1) {
	print "add \-\-\n";
	$addition ="\-\-";
	$seqlength = $seqlength + 2;
    } 
    if ($modulus_of_count == 2) {
	print "add \-\n";
	$addition = "\-";
	$seqlength = $seqlength + 1;
    }
    if ($modulus_of_count == 0) {
	print "add nothing\n";
	    $addition = "";
    }
    
    #save clean sequences
    $clean_nthash{$genename} = $seq;

    $seqlengthhash{$genename} = $seqlength;

}

for($i=1;$i<=$count;$i++) {
    
    #grab genename for each seq
    $genename = $namearray[$i];
    $value = $clean_nthash{$genename};
    print "$genename\n$value\n";
}



################## Print to PAML file #################################

## check that seqs are all of same length
$sum = 0;
for($i=1;$i<=$count;$i++) {
    
    #grab genename for each seq
    $genename = $namearray[$i];
    #grab seq length for each gene
    $sum = $sum + $seqlengthhash{$genename};
}
# divide sum by number of genes (count)... should be the same as one sequence 
$oneseq = $seqlengthhash{$genename};
if (($sum/$count) == $oneseq) {         
    
    $newfilename = $file . ".nuc";
    print "$newfilename\n";
    
    # PAML requires at least 333bp
    if ($oneseq > 333) {
	open (NF, ">$newfilename");
	print NF "   $count   $oneseq   I\n\n";

	for($i=1;$i<=$count;$i++) {
    	    #grab genename for each gene
	    $genename = $namearray[$i];
	    
	    print NF "$genename\n";           
	}
	print NF "\n";
	for($i=1;$i<=$count;$i++) {
    	    #grab seq for each gene
	    $genename = $namearray[$i];
	    $seq = $clean_nthash{$genename};
    	    
	    print NF "$seq\n";           
	}
	close(NF);

    } else {
	$too_short_files = $too_short_files .  "$file\n";
    }
} else {
    $non_matching_length = $file;
    print "Warning $non_matching_length has sequences of NON-MATCHING LENGTH\n";
}
open(SHORT, ">>too_short_for_paml");
print SHORT $too_short_files;
close(SHORT);

