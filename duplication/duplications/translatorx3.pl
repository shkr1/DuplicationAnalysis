#!/usr/bin/perl -w

use strict;

# Variables.
my @nt_seq;
my $nt_seq;
my $nt_name;
my %nt_nameseq;
my @nt_spnames;
my $code;
my $codehash;
my %namecode;
my $namecode_name;
my $namecode_code;
my @tempseqarray;
my $temptriplet;
my @tempaminoacidseq;
my $pauser;
my $clustaloutfile;
my $aa_name;
my @aa_seq;
my $aa_seq;
my %aa_nameseq;
my @aa_spnames;
my $shifted_nts;
my $fastaoutfile;
my $file1; 
my $clustalinfile;
my $pirfirstline;
my $pir_seq;
my @pir_seq;
my $pir_name;
my %pir_nameseq;
my @pir_spnames;
my $muscleinfile;
my $muscleoutfile;
my $samecodeYN; 
my $spcode; 
my $namecodefile; 
my $xcounter; 
my $discardline;
my $aa_posn; 
my @aa_array; 
my @nt_array; 
my $aa_length;
my @finalnt_array; 
my $finalaoutfile;
my $firstline;
my $secondline;
my $t_coffeeoutname;
my $filetype;
my $function;
my $tempnt =();
my $tempaa =();
my $tempntlength = ();
my $tempaalength = ();

# These are the codon table hashes
my %Univcodons = (
"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
"TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "X",
"TGT" => "C", "TGC" => "C", "TGA" => "X", "TGG" => "W",
"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
"ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
"AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");

my %VertMt = (
"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
"TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "X",
"TGT" => "C", "TGC" => "C", "TGA" => "W", "TGG" => "W",
"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
"ATT" => "I", "ATC" => "I", "ATA" => "M", "ATG" => "M",
"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
"AGT" => "S", "AGC" => "S", "AGA" => "X", "AGG" => "X",
"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");

my %BranchLanMt = (
"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
"TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "X",
"TGT" => "C", "TGC" => "C", "TGA" => "W", "TGG" => "W",
"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
"ATT" => "I", "ATC" => "I", "ATA" => "M", "ATG" => "M",
"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
"AGT" => "S", "AGC" => "S", "AGA" => "G", "AGG" => "X",
"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");

my %BranchFlMt = (
"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
"TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "X",
"TGT" => "C", "TGC" => "C", "TGA" => "W", "TGG" => "W",
"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
"ATT" => "I", "ATC" => "I", "ATA" => "M", "ATG" => "M",
"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
"AGT" => "S", "AGC" => "S", "AGA" => "S", "AGG" => "X",
"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");

my %YeastMt= (
"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
"TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "X",
"TGT" => "C", "TGC" => "C", "TGA" => "W", "TGG" => "W",
"CTT" => "T", "CTC" => "T", "CTA" => "T", "CTG" => "T",
"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
"ATT" => "I", "ATC" => "I", "ATA" => "M", "ATG" => "M",
"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
"AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");

my %MoldDipoloMt = (
"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
"TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "X",
"TGT" => "C", "TGC" => "C", "TGA" => "W", "TGG" => "W",
"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
"ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
"AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");

my %InvertMt = (
"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
"TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "X",
"TGT" => "C", "TGC" => "C", "TGA" => "W", "TGG" => "W",
"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
"ATT" => "I", "ATC" => "I", "ATA" => "M", "ATG" => "M",
"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
"AGT" => "S", "AGC" => "S", "AGA" => "S", "AGG" => "S",
"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");

my %CiliateNuc = (
"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
"TAT" => "Y", "TAC" => "Y", "TAA" => "Q", "TAG" => "Q",
"TGT" => "C", "TGC" => "C", "TGA" => "X", "TGG" => "W",
"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
"ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
"AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");

my %EchPlatyMt = (
"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
"TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "X",
"TGT" => "C", "TGC" => "C", "TGA" => "W", "TGG" => "W",
"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
"ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
"AAT" => "N", "AAC" => "N", "AAA" => "N", "AAG" => "K",
"AGT" => "S", "AGC" => "S", "AGA" => "S", "AGG" => "S",
"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");

my %EuplotidNuc = (
"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
"TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "X",
"TGT" => "C", "TGC" => "C", "TGA" => "C", "TGG" => "W",
"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
"ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
"AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");

my %AltYeastNuc = (
"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
"TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "X",
"TGT" => "C", "TGC" => "C", "TGA" => "X", "TGG" => "W",
"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "S",
"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
"ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
"AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");

my %AscidianMt = (
"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
"TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "X",
"TGT" => "C", "TGC" => "C", "TGA" => "W", "TGG" => "W",
"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
"ATT" => "I", "ATC" => "I", "ATA" => "M", "ATG" => "M",
"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
"AGT" => "S", "AGC" => "S", "AGA" => "G", "AGG" => "G",
"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");

my %BelpharismaNuc = (
"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
"TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "Q",
"TGT" => "C", "TGC" => "C", "TGA" => "X", "TGG" => "W",
"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
"ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
"AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");

my %MetazMtUniv2 = (
"TTT" => "-", "TTC" => "-", "TTA" => "-", "TTG" => "-",
"TCT" => "1", "TCC" => "1", "TCA" => "1", "TCG" => "-1",
"TAT" => "-", "TAC" => "-", "TAA" => "X", "TAG" => "X",
"TGT" => "-", "TGC" => "-", "TGA" => "-", "TGG" => "-",
"CTT" => "-", "CTC" => "-", "CTA" => "-", "CTG" => "-",
"CCT" => "-", "CCC" => "-", "CCA" => "-", "CCG" => "-",
"CAT" => "-", "CAC" => "-", "CAA" => "-", "CAG" => "-",
"CGT" => "-", "CGC" => "-", "CGA" => "-", "CGG" => "-",
"ATT" => "-", "ATC" => "-", "ATA" => "-", "ATG" => "-",
"ACT" => "-", "ACC" => "-", "ACA" => "-", "ACG" => "-",
"AAT" => "-", "AAC" => "-", "AAA" => "-", "AAG" => "-",
"AGT" => "2", "AGC" => "2", "AGA" => "2", "AGG" => "2",
"GTT" => "-", "GTC" => "-", "GTA" => "-", "GTG" => "-",
"GCT" => "-", "GCC" => "-", "GCA" => "-", "GCG" => "-",
"GAT" => "-", "GAC" => "-", "GAA" => "-", "GAG" => "-",
"GGT" => "-", "GGC" => "-", "GGA" => "-", "GGG" => "-", "---" => "-");

my %ArthropAGGMt = (
"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
"TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "X",
"TGT" => "C", "TGC" => "C", "TGA" => "W", "TGG" => "W",
"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
"ATT" => "I", "ATC" => "I", "ATA" => "M", "ATG" => "M",
"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
"AGT" => "S", "AGC" => "S", "AGA" => "S", "AGG" => "K",
"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");

my %HemichordMt = (
"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
"TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "X",
"TGT" => "C", "TGC" => "C", "TGA" => "W", "TGG" => "W",
"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
"ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
"AGT" => "S", "AGC" => "S", "AGA" => "S", "AGG" => "S",
"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");




	
# #########################################################################
# 
# Load up unaligned nucleotide file
# 
# 
# #########################################################################

print "Enter name of file of nucleotides to be aligned (FASTA or NBRF/PIR formats ONLY!!)\n";
$file1 = <STDIN>;
chomp $file1;


print "\n\nIf your file has Mac line endings TranslatorX is now converting them to UNIX line endings

To convert a Unix file back into a Mac file type the following at the command prompt

perl -pi -e 's/\\n\\r?/\\r/g' YOURFILENAME\n\n";

system "perl -pi -e 's/\\r\\n?/\\n/g' $file1";

if (! $file1){
	print "unable to open file\n"; 
	exit 1;
	} 
		

# Read in nucleotide file and put it into a hash keyed by the species name	

open (FILE1, $file1)  || die; # put dcse file into hash

$/ = '>'; # alter record separator for fasta/nbrf file

while (<FILE1>) {
	chomp;
	if ($_ =~ m/\w/i) {
		if ($_ =~ m/\w\w;/ and $_ =~ /\*/i) { # This means it is an nbrf/pir file
			($firstline, $nt_name, @nt_seq) = split "\n"; # name into $name, several seq lines into @Sequence
		}
		else { # its a FASTA file
			($nt_name, @nt_seq) = split "\n"; # name into $name, several seq lines into @Sequence
		}
		$nt_seq = join ('',@nt_seq); #now join up elements in @n_seq to create $n_seq.
		$nt_name =~ s/\s/_/g;
		$nt_name =~ s/[\(\)]//g;# getting rid of brackets
		$nt_seq =~ s/\s//g;
		$nt_name =~ s/[\[\]]//g;# not getting rid of square brackets
		$nt_nameseq{$nt_name} = $nt_seq;
		push @nt_spnames, $nt_name;
	}
	else {
		next;
	}
}
	
	

# #########################################################################
# 
# 				Decide how to use the software
# 
# 
# #########################################################################


$/ = "\n";

print "1	Start with unaligned nucleotides and no amino acid file
2	Start with aligned amino acids and unaligned nucleotides\n\n";


$function = <STDIN>;
chomp $function;

$function =~ m/[1-2]/ || die print "\n\n!!Must specify choice from 1 or 2!!\n\n";


# #########################################################################
# 
# 			Translate nucleotides, align the amino acids and 
# 			then go to the end bit where they are recombined
# 
# #########################################################################

if ($function == 1) {

	print "\n\nAligned nucleotides will have suffix .fasta.\nClustal output files will have suffixes .aln .dnd and .pir (aligned aa seqs)\nT_coffee output files will have suffixes .dnd and .pir (aligned aa seqs)\n\nOutput file prefix....\n";
	$clustalinfile = <STDIN>;
	chomp $clustalinfile;

	
	
# Find out what genetic code is to be used
	
	$/ = "\n";
	print "Do your taxa share the same genetic code? (Y/N)\n";
	$samecodeYN = <STDIN>;
	chomp $samecodeYN;
	
	if ($samecodeYN =~ m/y/i) {
		print "Which code?\n\n1\tUniversal\n2\tVertebrate mitochondrial\n3\tYeast mitochondrial\n4\tMold-Diploblast\n5\tInvertebrate mitochondrial\n6\tCiliate nuclear\n7\tEchinoderem-Rhabditophoran platyhelminth\n8\tEuplotid nuclear\n9\tAlternative yeast nuclear\n10\tAscidian mitochondrial\n11\tBelpharisma nuclear\n12\tBranchiostoma floridae\n13\tBranchiostoma lancelolatum\n14\tArthropodAGG=K\n15\tHemichordate\n";
		$code = <STDIN>;
		if ($code == 1) {
			$codehash = \%Univcodons;
			}
		elsif ($code == 2) {
			$codehash = \%VertMt;
			}
		elsif ($code == 3) {
			$codehash = \%YeastMt;
			}
		elsif ($code == 4) {
			$codehash = \%MoldDipoloMt;
			}
		elsif ($code == 5) {
			$codehash = \%InvertMt;
			}
		elsif ($code == 6) {
			$codehash = \%CiliateNuc;
			}
		elsif ($code == 7) {
			$codehash = \%EchPlatyMt;
			}
		elsif ($code == 8) {
			$codehash = \%EuplotidNuc;
			}
		elsif ($code == 9) {
			$codehash = \%AltYeastNuc;
			}
		elsif ($code == 10) {
			$codehash = \%AscidianMt;
			}
		elsif ($code == 11) {
			$codehash = \%BelpharismaNuc;
			}
		elsif ($code == 12) {
		$codehash = \%BranchFlMt;
		}
		elsif ($code == 13) {
		$codehash = \%BranchLanMt;
		}
		elsif ($code == 14) {
		$codehash = \%ArthropAGGMt;
		}
		elsif ($code == 15) {
		$codehash = \%HemichordMt;
		}
		else {
			print "not a valid choice\n";
			exit 1;
		}
	}
	
# If codes vary between taxa they can be downloaded from a file or input one by one.
	
	$spcode = 0;
	$namecodefile = '';
	$/ = "\n";
	
	if ($samecodeYN =~/[^y]/i) {
		print "Enter name of file: format \"Species_names TAB genetic_code\". Or press return\n";
		$namecodefile = <STDIN>;
		chomp $namecodefile;
	
		if ($namecodefile eq '') {
			print "Which code?\n\n1\tUniversal\n2\tVertebrate mitochondrial\n3\tYeast mitochondrial\n4\tMold-Diploblast\n5\tInvertebrate mitochondrial\n6\tCiliate nuclear\n7\tEchinoderem-Rhabditophoran platyhelminth\n8\tEuplotid nuclear\n9\tAlternative yeast nuclear\n10\tAscidian mitochondrial\n11\tBelpharisma nuclear\n12\tBranchiostoma floridae\n13\tBranchiostoma lancelolatum\n14\tArthropodAGG=K\n15\tHemichordate\n";
			foreach (@nt_spnames) {
				$spcode = 0;
				until ($spcode < 12 and $spcode > 0) {
					print "$_\t";
					$spcode = <STDIN>;
					chomp $spcode;
					$namecode{$_} = $spcode;
				}
			}
		}
		else  {
			open (NAMECODEFILE, $namecodefile) || die "unable to open $namecodefile\n";
			while (<NAMECODEFILE>) {
				chomp;
				($namecode_name, $namecode_code) = split "\t";
				$namecode_name =~ s/\s/_/g;
				$namecode_name =~ s/[\(\)]//g;
				$namecode_code =~ s/\s//g;
				$namecode{$namecode_name} = $namecode_code;
				print "name is $namecode_name code is $namecode_code\n";
			}
			close NAMECODEFILE;
		}
	}
	
	
# Open up the file to be sent to clustal.  
# This has to have each nucleotide translated into amino acids according to its appropriate code.
	$muscleinfile=$clustalinfile . '.min';
	$xcounter =0;
	open (CLUSTALOUT, ">$clustalinfile");
	open (muscleOUT, ">$muscleinfile");
	foreach (@nt_spnames) {
			if ($samecodeYN =~/n/i) {
				if (! $namecode{$_}){
					print "Crashed because no genetic code set for $_\n"; 
					exit 1;
				}
				$code = $namecode{$_};
				if ($code == 1) {
					$codehash = \%Univcodons;
				}
				elsif ($code == 2) {
					$codehash = \%VertMt;
				}
				elsif ($code == 3) {
					$codehash = \%YeastMt;
				}
				elsif ($code == 4) {
					$codehash = \%MoldDipoloMt;
				}
				elsif ($code == 5) {
					$codehash = \%InvertMt;
				}
				elsif ($code == 6) {
					$codehash = \%CiliateNuc;
				}
				elsif ($code == 7) {
					$codehash = \%EchPlatyMt;
				}
				elsif ($code == 8) {
					$codehash = \%EuplotidNuc;
				}
				elsif ($code == 9) {
					$codehash = \%AltYeastNuc;
				}
				elsif ($code == 10) {
					$codehash = \%AscidianMt;
				}
				elsif ($code == 11) {
					$codehash = \%BelpharismaNuc;
				}
				elsif ($code == 12) {
					$codehash = \%BranchFlMt;
				}
				elsif ($code == 13) {
					$codehash = \%BranchLanMt;
				}
				elsif ($code == 14) {
				$codehash = \%ArthropAGGMt;
				}
				elsif ($code == 15) {
				$codehash = \%HemichordMt;
				}
				else {
					print "not a valid choice\n";
					exit 1;
				}
			}
			else {
			}
		splice @tempseqarray;
		$temptriplet = '';
		splice @tempaminoacidseq;
		@tempseqarray = $nt_nameseq{$_} =~ /(.{3})/g;
		foreach  $temptriplet (@tempseqarray) {
			$temptriplet = (uc $temptriplet);
			if ($$codehash{$temptriplet}){
				if ($$codehash{$temptriplet} =~/x/i) {
					$xcounter++;
				}
			  push  (@tempaminoacidseq, $$codehash{$temptriplet})
			}
			else {
			  push (@tempaminoacidseq, "X");
			}
		}
		if ($xcounter >1) {
			open (ERRORLOG, ">>error.log");
			print ERRORLOG "in file $file1 Species $_ has $xcounter termination codons.  Suggest check correct code used and that in frame\n"
		}
		print CLUSTALOUT">P1;$_\n$_\n",@tempaminoacidseq,"*\n";
		$xcounter =0;
		print muscleOUT">$_\n",@tempaminoacidseq, "\n";
	}			
	
	close ERRORLOG;
	close CLUSTALOUT;
	close muscleOUT;
	close FILE1;
	
	
	
	$clustaloutfile = $clustalinfile . '.pir';
	$muscleoutfile=$clustalinfile . '.mout';
	
	
	
	print "For clustal alignment (quicker) press C\nFor T-Coffee alignment (slower but better) press T\nFor muscle alignment (fast AND furious) press M\n";
	$pauser = <STDIN>;
	chomp $pauser;
	
# Make a file that will send instructions to Clustal with appropriate name of input file
	
	if ($pauser =~/[C]/i) {
	
		open (CLUSTALINPUT, ">clustal.input");
		print CLUSTALINPUT "1\n$clustalinfile\n2\n9\n1\n2\n\n1\n\n\nx\n\nx\n\nx\n\n";
		close CLUSTALINPUT;
		
		# Start up CLustalW and run it according to the instruction file just made.
		
		system "clustalW <clustal.input";
	}
	
	elsif ($pauser =~/[T]/i) {
		$t_coffeeoutname=$clustalinfile . '.pir';
		system "t_coffee -infile=$clustalinfile -output=pir_aln -outfile=$t_coffeeoutname";
	}
	
	elsif ($pauser =~/[M]/i) {
		$muscleoutfile=$clustalinfile . '.mout';	
		system "./muscle -in $muscleinfile -out $muscleoutfile";
	}
	
	
	else {
	print "not a valid choice\n";
			exit 1;
	}


# Open up the clustal alignment output file and put the ALIGNED amino acids into a hash keyed by the species name
	if ($pauser =~/[TC]/i) {
	$clustaloutfile = $clustalinfile . '.pir';
	}
	elsif ($pauser =~/[M]/i) {
	$clustaloutfile = $clustalinfile . '.mout';
}
	
}


# #########################################################################
# 
# 					Get the right file names for options 1 and 2 
# 								
# 
# #########################################################################


if ($function == 1){
	$fastaoutfile = $clustalinfile . '.fasta';
}

if ($function == 2) {
	print "\n\n Name of your pir/nbrf or fasta format aligned amino acids file\n\n";
	$clustaloutfile = <STDIN>;
	chomp $clustaloutfile;
	
	system "perl -pi -e 's/\\r\\n?/\\n/g' $clustaloutfile";

	
	print "\n\nOutput file prefix....\n";
	$fastaoutfile = <STDIN>;
	chomp $fastaoutfile;
	
	$fastaoutfile = $fastaoutfile . '.output';
}


# #########################################################################
# 
# Insert dashes into the unaligned nucleotide file according to the aligned 
# 								amino acids
# 
# #########################################################################

open (CLUSTALOUTFILE, $clustaloutfile) || die "unable to open $clustaloutfile\n"; # put nbrf file into hash

$/ = '>'; # alter record separator for fasta/nbrf file

while (<CLUSTALOUTFILE>) {
	chomp;
	if ($_ =~ m/\w/i) {
		if ($_ =~ m/\w\w;/ and $_ =~ /\*/i) { # This means it is an nbrf/pir file
			($aa_name, $discardline, @aa_seq) = split "\n"; # name into $name, several seq lines into @Sequence
		}
		else { # its a FASTA file
			($aa_name, @aa_seq) = split "\n"; # name into $name, several seq lines into @Sequence
		}
		$aa_seq = join ('',@aa_seq); #now join up elements in @n_seq to create $n_seq.
		$aa_name =~ s/\s/_/g;
		$aa_name =~ s/P1;//;
		$aa_name =~ s/[\(\)]//g; # not getting rid of underscores
		$aa_seq =~ s/\s//g;
		$aa_seq =~ s/\*//g; # get rid of * at end of NBRF file
		$aa_nameseq{$aa_name} = $aa_seq;
		push @aa_spnames, $aa_name;
		$discardline='';
	}
	else {
		next;
	}
}

# #########################################################################
# 
# Check each nt seq is 3x length of equivalent aa seq
# 							
# 
# #########################################################################

foreach (@aa_spnames) {
	$tempnt = $nt_nameseq{$_};
	$tempaa = $aa_nameseq{$_};
	$tempnt =~ s/\-//g;
	$tempaa =~ s/\-//g;
	$tempntlength = length $tempnt;
	$tempntlength = int ($tempntlength/3);
	$tempaalength = length $tempaa;
	if ($tempntlength ne $tempaalength) {
		print "length ($tempaalength) of amino acid seq from $_ does not correspond to 1/3 of nucelotide seq ($tempntlength)\n";
		exit;
	}
	else {
		next;
	}
	$tempnt =();
	$tempaa =();
	$tempntlength = ();
	$tempaalength = ();
}




# Now compare original NUC file with newly aligned AMINO ACID file and put dashes in as appropriate.

foreach (@nt_spnames) {
	if (! $aa_nameseq{$_}) {
		print "Warning:	Sequence $_ not found in amino acid file\n"
	}
}


foreach (@aa_spnames) {
	
	if ($nt_nameseq{$_}) {
		$aa_posn = 0;
		@aa_array = split(//, $aa_nameseq{$_});
		@nt_array = $nt_nameseq{$_} =~ /(.{3})/g;
		$aa_length = scalar(@aa_array);
		@finalnt_array = ();
			while  ($aa_posn < $aa_length) {
				if ($aa_array[$aa_posn] =~ /\-/) {
				push @finalnt_array, "---";
				$aa_posn++;
				}
				elsif ($aa_array[$aa_posn] =~ /[A-Za-z]/) {
				$shifted_nts = shift @nt_array;
				push @finalnt_array, $shifted_nts;
				$aa_posn++;
				}
				else {
				$aa_posn++;
				}
			}
		delete $nt_nameseq{$_};
	}

	else {
		print "Warning:	Sequence $_ not found in nucleotide file\n"
	}
	
open (FASTAOUT, ">>$fastaoutfile");
print FASTAOUT ">$_\n",@finalnt_array,"\n";

}
close FASTAOUT;
