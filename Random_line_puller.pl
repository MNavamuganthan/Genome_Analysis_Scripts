#!/usr/bin/perl
use strict; use warnings;

my $PATH = "/Users/MCN/Documents/UConn/Genomes/Yersinia/Submission_Preparation/ENA_Validate/";
my $log = $PATH . "Log.txt";

# Pre-execution safety checking
if ( ! @ARGV) {
    die "Script aborted: Must supply an input fasta as a script variable!";
}
open (INFILE, "$ARGV[0]") or die "Cannot open input file\n";

my $thisline;                      # The current line that we're at in the file
my $line_number = 0;               # Lets us know what line we're at
my $badlines = 0;

while ($thisline = <INFILE>) {             # Begin reading in the input file, line by line
    $thisline =~ s/\s+$//;                 # Remove any trailing whitespaces and end of line characters
    $line_number++;                        # Increment line number by 1
    my $length = length($thisline);
    if ($thisline =~ /[A-Z]+[0-9]+/) {
        if ($thisline =~ /locus_tag=/) {
            next;
        }
        elsif ($thisline =~ /note=/) {
            next;
        }
        $badlines++;
        print "$line_number: $thisline\n";
    }
}

print "\nThere are $badlines line left still to fix :(\n"