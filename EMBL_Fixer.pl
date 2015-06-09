#!/usr/bin/perl
use strict; 
use warnings;
use Cwd;

my $PWD = getcwd();
my $PATH = "/Users/MCN/Documents/UConn/Genomes/Yersinia/Submission_Preparation/ENA_Validate/";
my $outfile = $PATH . $ARGV[1];
my $log = $PATH . "Log.txt";

# Pre-execution safety checking
if ( ! @ARGV) {
    die "Script aborted: Must supply an input fasta as a script variable!";
}
open (INFILE, "$ARGV[0]") or die "Cannot open input file\n";
open (OUTFILE, ">$outfile") or die "Cannot create output file\n";
open (LOG, ">$log") or die "Cannont create log file\n";

print LOG "Your input file was $ARGV[0]\n";
print LOG "The output file will be named $outfile\n";
print LOG "The log file is named Log.txt\n";

# Initialize global use variables
my $line_number = 0;               # Lets us know what line we're at
my $locus_tag_prefix = "$ARGV[2]";   # The locus tag prefix
my $protein_id_prefix = "";        # The protein ID prefix; Rarely added
my $locus_tag_count = 0;           # Current count of added locus tags

# Initialize state variables for tracking current position in table file
my $thisline;                      # The current line that we're at in the file
my $is_key;                        # Define if the current line as either a feature key (CDS, RNA, etc.)
my $is_qual;                       # Or as a qualifier for the feature (locus_tag, product, etc.)
my $current_key;                   # Define what the current key type is
my $current_loc;                   # Define the current location
my $current_qual;                  # Define the current qualifier
my $current_val;                   # Define the value
my @master_loc;                    # Track the mater location of all features
my @tokens;                        # 

clearflags ();                # Clear all variables before running through file, shouldn't be necessary but safety first

# Things that need fixing in current file examples:
# -------------------------------------------------
# 1. gene features need to be added for every CDS, rRNA, tRNA feature  
# 2. Locus tags need to be added to genes and features, incrementing by 5 each time
# 3. Amino acid translations need to be wrapped at 80 characters per line.
# 4. 
#
#
#

while ($thisline = <INFILE>) {             # Begin reading in the input file, line by line
    $thisline =~ s/\s+$//;                 # Remove any trailing whitespaces and end of line characters
    $line_number++;                        # Increment line number by 1
    #print "$line_number: ";                # Print the current line number to STDOUT

# If the current line begins with FT (Feature Table), and feature is in 5'->3' orientation
    if ($thisline =~ /^FT   (\S+) +(\d+\.\.\d+)$/) {     # \1 = the key term, \2 = feature locus range
        $is_key = 1;                       # This line denotes a feature key so we need to capture that
        $is_qual = 0;                      # It is not a qualifier, those lines comes next so reset 
        $current_key = $1;                 # Set current key to be the name of the feature term
        $current_loc = "$2";               # Set current location to the locus range
        @master_loc = ();                  # Keep track of the overall location
        push (@master_loc, $current_loc);  # Append the location into the master tracking array
        flushline ();                      # Run flushine to print feature
    }
# If the current line begins with FT (Feature Table), and feature is in 3'->5' orientation
    elsif ($thisline =~ /^FT   (\S+) +(complement\(\d+\.\.\d+\))/) {
        #print STDOUT "ORIG: $thisline\n";
        # We have a feature in 3'->5' direction
        $is_key = 1;
        $is_qual = 1;
        $current_key = $1;
        $current_loc = "$2";
        @master_loc = ();
        push (@master_loc, $current_loc);
        flushline ();
    } 
    # If the current line begins with FT (Feature Table), but does not give a location
    elsif ($thisline =~ /^FT.*$/) {
        # We have a qualifier with data, now what?
        print OUTFILE "$thisline\n";  # Print it I guess
    }
    elsif ($thisline =~ '^$') {   # If line is empty then skip; often occurs when file has blank line at the end
        next;
    }
    else {                    # If line not capture by either of the above conditions then just print it
          print OUTFILE "$thisline\n";
    }
}

# Subroutine for clearing the state variables back to null/0 values
sub clearflags {                   
  $thisline = "";
  $is_key = 0;
  $is_qual = 0;
  $current_key = "";
  $current_loc = "";
  $current_qual = "";
  $current_val = "";
  @master_loc = ();
}

# subroutine for printing feature key and location
sub printkey {                 # @_ is the array of variables passed to the subroutine for 
  my $thiskey = shift (@_);    # Initialize and set current key value
  my $is_first = 1;            # Set that this is the first instance
  foreach my $thisloc (@master_loc) {  # For each location in the master list, do...
    if ($is_first == 1) {              # If the value of is_first is 1 then
        if ($thiskey eq "CDS") {       # If the key value is CDS then print the feature key line to file
            print OUTFILE "FT   $thiskey             $thisloc\n";
        } elsif ($thiskey eq "rRNA" || $thiskey eq "tRNA" || $thiskey eq "gene") {
            print OUTFILE "FT   $thiskey            $thisloc\n";
        } elsif ($thiskey eq "sig_peptide") {
            print OUTFILE "FT   $thiskey     $thisloc\n";
        } elsif ($thiskey eq "repeat_region") {
            print OUTFILE "FT   $thiskey   $thisloc\n";
        } elsif ($thiskey eq "misc_feature") {
            print OUTFILE "FT   $thiskey    $thisloc\n";
        } elsif ($thiskey eq "source") {
            print OUTFILE "FT   $thiskey          $thisloc\n";
        }
    } 
    $is_first = 0;
  }
}

# subroutine for adding locus tags by feature key
sub printfeat {
  my $thiskey = shift (@_);

  if ($thiskey eq "CDS" || $thiskey eq "rRNA" || $thiskey eq "tRNA" || $thiskey eq "ncRNA" || $thiskey eq "tmRNA") {
    printkey ("gene");
    my $locus_tag_string = "0000";
    $locus_tag_count = $locus_tag_count + 5;
    $locus_tag_string = "0000";
    if ($locus_tag_count > 999) {
      $locus_tag_string = $locus_tag_prefix . "_" . $locus_tag_count;
    } elsif ($locus_tag_count > 99) {
      $locus_tag_string = $locus_tag_prefix . "_0" . $locus_tag_count;
    } elsif ($locus_tag_count > 9) {
      $locus_tag_string = $locus_tag_prefix . "_00" . $locus_tag_count;
    } else {
      $locus_tag_string = $locus_tag_prefix . "_000" . $locus_tag_count;
    }     
    print OUTFILE "FT                   /locus_tag=\"$locus_tag_string\"\n";
    printkey ($thiskey);
    print OUTFILE "FT                   /locus_tag=\"$locus_tag_string\"\n";
    return;
    }

  printkey ($thiskey);
}

# subroutine for filtering by qualifier and value
sub printqual {
  my $thisqual = shift (@_);
  my $thisval = shift (@_);
  my $thiskey = shift (@_);
  print "\n We're running the printqual subroutine\n";
  print OUTFILE "\t\t\t$thisqual\t$thisval\n";
}

# Subroutine to print feature key / location / qualifier lines
sub flushline {
    if ($is_key == 1) {  # If all filter subroutine for feature keys
        printfeat ($current_key);
    } elsif ($is_qual == 1) {
        if ($current_val eq "") {
            print OUTFILE "\t\t\t$current_qual\n";
        } else {
            # call filter subroutine for qualifiers
            printqual ($current_qual, $current_val, $current_key);
        }
      }
}

# Close input and output files for safety purposes
close (INFILE);
close (OUTFILE);
close (LOG);