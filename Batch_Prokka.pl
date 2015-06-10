#!/usr/bin/env perl
#
# Script to take a list of genome sequences and run prokka as a batch
# The input list file is a tab-delimited, one entry per line format.
# Genome files should be in fasta/multi-fasta format.
# Output for each genome is put into individual directories named according to the ID given by the user in the input list file
# User specifies the kmer size during execution, which will be appended to the sorted output filename (_Nmer.txt)
#
# NB: This script serves as a batch wrapper for prokka which does all of the actual work.
#     The prokka executable must be present somewhere in the users PATH in order for this script to work.
#
# Script skeleton and some tasty bits taken from prokka (https://github.com/tseemann/prokka)
#
# Created by: Michael C. Nelson
# Version: 1.1
# Created on: 2015-06-03
# Revised on: 2015-06-10
# License: GPL3
######################################################################################################

use strict;
#use warnings;
use FindBin;
use File::Copy;
use Time::Piece;
use Time::Seconds;
use Getopt::Long;
use Scalar::Util qw(openhandle);

###### Initial variable instancing and options setting ######
my $input;
my $output;
my $outdir;
my $lazy;
my $logfile;
my $genus;
my $center = 'UCONNMISEQ';
my $cpus = 4;
my $starttime = localtime;
my $fp = find_exe("prokka");
err("FATAL ERROR: Can't find Prokka in your \$PATH!") if !$fp;
my $procmd = "prokka --addgenes --rawproduct ";

my @Options;
setOptions();

###### Begin initial checks and building of prokka command ######
if ($input) {
    open(IN, $input) or err("ERROR: Could not open input file.\n");    
}
else {
    err("ERROR: Input file not provided")
}

$procmd .= "--cpus $cpus ";

if ($lazy) {
    $procmd .= "--debug ";
}

if ($genus) {
    $procmd .= "--usegenus ";
}

if ($outdir) {
    if (-d $outdir){
        err("Output directory already exists, choose a new name for --outdir that is not $outdir");
    }
    else {
        runcmd("mkdir -p \Q$outdir\E")
    }
    $logfile = "$outdir/Batch_Prokka.log";    
}
else {
    $logfile = "Batch_Prokka.log";
}

###### ACTUAL WORK GETS DONE HERE ######

open LOG, '>', $logfile or err("Can't open logfile");
msg("Began running Batch_Prokka.pl at $starttime");
if ($outdir) {
        msg("Output directory for each genome will be put into $outdir\/.");
    }
    else {
        msg("Output directory for each genome will be put into current directory.");
    }
msg("Will use maximum of $cpus cores.");
msg("Writing log to: $logfile");
msg("Using $input as the input file.");

while (<IN>) {
    chomp;
    my @line = split(/\t/);
    my $gFP = $line[0];
    my $gGenus = $line[1];
    my $gSpecies = $line[2];
    my $gID = $line[3];
    my $gLocus = $line[4];
    my $cmd = $procmd;
    ## Finalize the prokka command using table data
    if ($outdir) {
        $output = $outdir."/".$gID;
    }
    else {
        $output = $gID;
    }
    $cmd .= "--outdir $output --genus $gGenus --species $gSpecies --strain $gID --locustag $gLocus --prefix $gID --centre $center $gFP";

    ## Now run prokka
    msg("Annotating genome $gID using the file $gFP:");
    print "The prokka command that would be run is:\n\n$cmd\n\n";
    #runcmd($cmd);
}


my $endtime = localtime;
my $walltime = $endtime - $starttime;
my $pretty = sprintf "%.2f minutes", $walltime->minutes;
msg("Finished processing. Total time taken was: $pretty");

###### Sub-routines ######
sub find_exe {
    my($bin) = shift;
    for my $dir (File::Spec->path) {
        my $exe = File::Spec->catfile($dir, $bin);
        return $exe if -x $exe;
    }
    return;
}

sub runcmd {
    system(@_)==0 or err("Could not run command:", @_);
}

sub msg {
    my $t = localtime;
    my $line = "[".$t->hms."] @_\n";
    print LOG $line if openhandle(\*LOG);
    print STDERR $line;
}

sub err {
  msg(@_);
  exit(2);
}

sub setOptions {
    @Options = (
    'Mandatory:',
    {OPT=>"input=s", VAR=>\$input, DESC=>"The input table of genomes to annotate."},
    'Options:',
    {OPT=>"outdir=s", VAR=>\$outdir, DESC=>"Directory where results will be put. [DEFAULT=$outdir]\n\t\t    NOTE: Each annotated genome will be put in its own directory named according to the Strain ID provided in the input file."},
    {OPT=>"usegenus!", VAR=>\$genus, DESC=>"Use custom database for initial annotation.\n\t\t    NOTE: Genus name provided in input file must match a valid database name or error will occur."},
    {OPT=>"center=s", VAR=>\$center, DEFAULT=>$center, DESC=>"Sequencing center name."},
    {OPT=>"cpus=i", VAR=>\$cpus, DESC=>"Numer of CPUs to use. [DEFAULT=$cpus]"},
    {OPT=>"lazy!", VAR=>\$lazy, DESC=>"Don't delete intermediate files."},    
    'Help:',
    {OPT=>"help", VAR=>\&usage, DESC=>"Print this help message."},
    {OPT=>"example", VAR=>\&example, DESC=>"Print an example input table."},
    );
    
    (!@ARGV) && (usage());
    
    &GetOptions(map {$_->{OPT}, $_->{VAR}} grep { ref } @Options) || usage();
    
    # Now setup default values.
    foreach (@Options) {
        if (ref $_ && defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
            ${$_->{VAR}} = $_->{DEFAULT};
        }
    }
}

sub usage {
    print STDERR
    "\nBatch_prokka.pl: A prokka wrapper for batch annotation of multiple genome sequences.\n",
    "Usage:\tBatch_prokka.pl [options] --input genomes_list.txt\n\n";
    foreach (@Options) {
        if (ref) {
            my $def = defined($_->{DEFAULT}) ? " [DEFAULT='$_->{DEFAULT}']" : "";
            my $opt = $_->{OPT};
            $opt =~ s/!$//;
            $opt =~ s/=s$/ <xxx>/;
            $opt =~ s/=i$/ <\#>/;
            printf STDERR "  --%-15s %s%s\n", $opt, $_->{DESC}, $def;
        }
        else {
            print STDERR "$_\n";
        }
    }
    print "\n";
    exit(1);
}

sub example {
    print STDERR
    "\nThe input file format is tab-delimited as follows, without the header line.
Note that the genus name must match that of one of the custom databases available if using the --genus flag.
Installed databases can be checked by running prokka --listdb.
If actual values are not known they should be filled in with placeholders (e.g. unknown for Species)\n
    Contig files    Genus      Species     StrainID/Output directory  Locus tag  
    contigs1.fasta  Aeromonas  hydrophila  Ah1                        ALO05        
    contigs2.fasta  Aeromoans  veronii     Hm21                       M001
    ...\n
";
    exit(0);
}