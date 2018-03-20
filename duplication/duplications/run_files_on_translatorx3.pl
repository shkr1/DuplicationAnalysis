#!/usr/bin/perl -w

my @aln_file = <*.txt>;



foreach my $file (@aln_file){
	
	$outfilename =$file;
	$outfilename =~ s/\.new/out/;

	open (HERVEINPUT, ">herve.input");
	print HERVEINPUT "$file\n";
	print HERVEINPUT "1\n";
	print HERVEINPUT "$outfilename\n";
	print HERVEINPUT "y\n";
	print HERVEINPUT "1\n";
	print HERVEINPUT "m\n";
	close HERVEINPUT;
	system "perl translatorx3.pl < herve.input";
}