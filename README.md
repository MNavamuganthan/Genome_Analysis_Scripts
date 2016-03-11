# Genome Analysis Scripts
A collection of scripts for annotating and analyzing genomes. Also includes scripts used for modifying/fixing files for use in annotated genome submissions to an INSDC repository.

## Batch_Prokka.pl
Wrapper script to perform batch annotation on a set of genomes using [Prokka](https://github.com/tseemann/prokka). It's highly recommended that you understand how Prokka works and it's options before using this batch script.

Takes a tabular input table listing the location of the contigs and other information and runs prokka, using a custom annotation database if one is available and requested.

__Options:__
```shell
Mandatory:
  --input <xxx>     The input table of genomes to annotate.
Options:
  --outdir <xxx>    Directory where results will be put. [DEFAULT=]
		    NOTE: Each annotated genome will be put in its own directory named according to the Strain ID provided in the input file.
  --usegenus        Use custom database for initial annotation.
		    NOTE: Genus name provided in input file must match a valid database name or error will occur.
  --center <xxx>    Sequencing center name. [DEFAULT='UCONNMISEQ']
  --cpus <#>        Numer of CPUs to use. [DEFAULT=4]
  --quiet           Don't display prokka output.
  --lazy            Don't delete intermediate files.
Help:
  --help            Print this help message.
  --example         Print an example input table.
```

__Useage:__
```shell
Batch_Prokka.pl --input bath_input.txt --usegenus --center GrafLab --cpus 8
```
##






