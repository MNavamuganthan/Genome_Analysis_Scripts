#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

my $PWD = getcwd();
my $PATH = "/Users/MCN/Documents/UConn/Genomes/Yersinia/Submission_Preparation/ENA_Validate/";
my $outfile = $PATH . "test_out.txt";

# Pre-execution safety checking
if ( ! @ARGV) {
    die "Script aborted: Must supply an input fasta as a script variable!";
}
open (INFILE, "$ARGV[0]") or die "Cannot open input file\n";
open (OUTFILE, ">$outfile") or die "Cannot create output file\n";

my $peptide;
my $thisline;
my $line_number;
my $start;
while ($thisline = <INFILE>) {             # Begin reading in the input file, line by line
    $line_number++;
    $thisline =~ s/\s+$//;                 # Remove any trailing whitespaces and end of line characters
    # print OUTFILE "$line_number: ";                # Print the current line number to STDOUT
    # If the current line begins with FT (Feature Table), and feature is in 5'->3' orientation
    if ($thisline =~ /^FT {19}\/translation=("[A-Z]+")$/) {    # If the line contains a translation, then
        $peptide = $1;                     # Place the amino acid sequence into the string peptide
        print OUTFILE "FT                   /translation=", substr($peptide, 0, 46), "\n";
        my $pep_length = length($peptide);
        $start = 46;
        while ($pep_length > 59) {
            print OUTFILE "FT                   ", substr($peptide, $start, 59), "\n";
            $pep_length = $pep_length - 59;
            $start = $start + 59;
        }
        if ($pep_length > 0) {
            print OUTFILE "FT                   ", substr($peptide, $start), "\n";            
        }
    }
    else{
        print OUTFILE "$thisline\n";
    }
}