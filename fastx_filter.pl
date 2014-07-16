#!/usr/bin/perl

#===============================================================================
#   Author: Alejandro Manzano Marin
#
#   File: fastx_filter.pl
#   Date: 30-11-2012
#   Version: 1.0
#
#   Usage:
#      perl fastx_filter.pl --infile|-i inFile.fasta --list|-l listFile.lst [--informat|-if 'fasta'|'fastq'] [options]
#
#      Check out 'perl fastx_filter.pl -h' for short usage manual and info on the software.
#
#    Description: This program is designed to be a general purpose filtering tool for both FASTA and
#                 FASTQ formats. It takes as input a FAST(Q/A) formatted file and can perform both
#                 filtering and reverse filtering (selecting all that are not on your list) and perform
#                 some reformatting mainly useful for creating files to use with NGS assemly, mapping,
#                 etc..
#                 
#
#    Contact: Contact the author at alejandro.manzano@uv.es using 'fastx_filter: ' as
#             as begining for subject for bug reporting, feature request or whatever
#             (of course realted to the software).
#
#    COPYRIGHT: Copyright (C) 2012  Alejandro Manzano Marin.
#
#    LICENCE: This program is free software: you can redistribute it and/or modify it under the terms
#             of the GNU General Public License as published by the Free Software Foundation, either
#             version 3 of the License, or (at your option) any later version.
#             This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#             without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#             See the GNU General Public License for more details.
#             You should have received a copy of the GNU General Public License along with this program.
#             If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


# Load modules
use Getopt::Long;
use Pod::Usage;


# Define subroutines
sub printVersion {
	my $software= $_[0];
	my $version= $_[1];
	my $fileHandle= $_[2];
	print $fileHandle "$software v$version\n";
	exit (0);
}


# Variable definition

## Define other variables
my $file= "";
my $line= "";
my $id= "";
my $flag= 0;
my %reads= ();
my $inHead= ">";

## General variables
my $PROGRAMNAME= "fastx_filter";
my $VERSION= "1.0";

## Define options default values
my $opt_inFiles= "";
my $opt_listFile= "";

my $opt_inFormat="fasta";
my $opt_outSingle= 0;
my $opt_noHdQual= 0;

my $opt_exactMatch= 0;
my $opt_revFilter= 0;

my $opt_verbose= 0;
my $opt_man= 0;
my $opt_help= 0;
my $opt_printVersion= 0;

## Define options hash
GetOptions(\%opts, 
	'infile|i=s' => \$opt_inFile, 
	'list|l=s' => \$opt_listFile, 
	'informat|if:s'=> \$opt_inFormat, 
	'outSingle!'=> \$opt_outSingle, 
	'noHdQual!' => \$opt_noHdQual, 
	'exact!' => \$opt_exactMatch, 
	'revFilter!' => \$opt_revFilter, 
	'verbose|v!' => \$opt_verbose, 
	'help|h!' => \$opt_help, 
	'man!'  => \$opt_man, 
	'version!' => \$opt_printVersion) || pod2usage(-exitval => 1,  -verbose => 2);

if ($opt_help){
	pod2usage(-exitval => 1,  -verbose => 1);
}
if ($opt_man){
	pod2usage(-exitval => 0, -verbose => 2);
}
if ($opt_printVersion){
	&printVersion($PROGRAMNAME, $VERSION, \*STDERR);
}


# Script documetation

=pod

=head1 NAME

fastx_filter

=head1 VERSION

fastx_filter v1.0

=head1 SYNOPSIS

perl fastx_filter.pl --infile|-i inFile.fasta --list|-l listFile.lst [--informat|-if 'fasta'|'fastq']
[--outSingle] [--noHdQual] [--exact] [--revFilter] [--verbose|-v] [--help|-h] [--man] [--version]

=head1 DESCRIPTION

This program is designed to be a general purpose filtering tool for both FASTA and FASTQ formats. It
takes as input a FAST(Q/A) formatted file and can perform both filtering and reverse filtering (selecting
all that are not on your list) and perform some reformatting mainly useful for creating files to use with
NGS assemly, mapping, etc..

=head1 OPTIONS

=head2 INPUT

=over 8

=item B<-i> | B<--infile> <string> (mandatory)

File containing sequences to be filtered in FAST(A/Q) format. Repeated as many times as input files.

=item B<-l> | B<--list> <string> (mandatory)

File containing in the first column sequence identifiers (Remember they cannot contain spaces, usually
the case with NGS data).

=back

=head2 FORMAT

=over 8

=item B<-if> | B<--informat> <string> (default: "fasta")

Input file format, either "fasta" or "fastq".

=item B<--outSingle> <boolean> (default: 0)

Output sequence IDs in single-end format (especially useful when filtering 454 or Illumina reads from
paired- or mate-pair runs, takes out the .f,.r,/1,/2 endning of a read).

=item B<--noHdQual> <boolean> (default: 0)

Output fastq files without the header in the quality header (really good for reducing file size).

=back

=head2 FILTERING

=over 8

=item B<--exact> <boolean> (default: 0)

Only do exact matching for sequence IDs (especially usefull when you only want a mate of a pair selected,
leave off if interested in collecting both mates even if only one is on the list).

=item B<--revFilter> <boolean> (default: 0)

If on, perform reverse filtering (as in select all that are not in the list).

=back

=head2 INFO AND HELP

=over 8

=item B<-v> | B<--verbose> <boolean> (default: 0)

Prints status and info messages while processing.

=item B<-h> | B<--help> <boolean>

Print useful help on using this script.

=item B<--man> <boolean>

Print the full documentation.

=item B<--version> <boolean>

Print program version.

=back

=head1 AUTHOR

Alejandro Manzano Marin, C<< <alejandro_dot_manzano_at_uv_dot_es> >>

=head1 BUGS

If you find a bug please email me at C<< <alejandro_dot_manzano_at_uv_dot_es> >> so that I can keep making
fastx_filter better.

=head1 COPYRIGHT

Copyright (C) 2012  Alejandro Manzano Marin.

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see
<http://www.gnu.org/licenses/>.

=cut


# Assign values to variables dependant on options

## Format
if ($opt_inFormat=~ m/fastq/){
	$inHead="@";
}


# Check options
if (!$opt_inFile || !$opt_listFile || ($opt_inFormat && $opt_inFormat!~ m/(fasta|fastq)/)){
	print STDERR "ERROR:\n";
	if (!$opt_inFile){
		print STDERR "FAST(A/Q) infile missing\n";
	}
	if (!$opt_listFile){
		print STDERR "List file missing please check usage manual\n";
	}
	if ($opt_inFormat && $opt_inFormat!~ m/(fasta|fastq)/){
		print STDERR "Input file format neither 'fasta' nor 'fastq'\n";
	}
	print STDERR "Please check usage manual\n";
	pod2usage(-exitval => 1,  -verbose => 1);
}


# Main program

## Read list of accepted reads
if ($opt_verbose){
	print STDERR "Reading list for filtering " . $opt_listFile . "..";
}

open (LISTIN, "$opt_listFile") || die "Unable to open file $opt_listFile for reading\n$!\n";;
while ($line=<LISTIN>){
	if ($line =~ m/^(\S+)/){
		$id= $1;
		$id=~ s/^[@\+>]//;
		if (!$opt_exactMatch){
			$id=~ s/[\.\/][fr12]$//;
		}
		$reads{$id}=1;
	}
}
close (LISTIN);

if ($opt_verbose){
	print STDERR "DONE\n";
}

## Read fastx files and write output to STDOUT
if ($opt_verbose){
	print STDERR "Filtering FASTX file " . $opt_inFile . "..";
}

open (FASTXIN, "$opt_inFile") || die "Unable to open file $opt_inFile for reading\n$!\n";
while ($line= <FASTXIN>){
	chomp $line;
	if ($flag && $opt_inFormat=~ m/fasta/ && $line !~ m/^$inHead/){
		print STDOUT $line . "\n";
		next;
	}
	if ($line =~ m/^$inHead(\S+)/){
		$id= $1;
		if (!$opt_exactMatch){
			$id=~ s/[\.\/][12fr]$//;
		}
		$id=~ s/^$inHead//;
		if ($opt_outSingle){
			$line=~ s/\s+.*$//;
			$line=~ s/[\.\/][fr12]$//; #for single-end illumina files
			$line=~ s/[\.]fn$//;
		}
		if (!$opt_revFilter){ # If filter IS NOT set to reverse
			if ($reads{$id}){
				print STDOUT $line . "\n";
				if ($opt_inFormat=~ m/fasta/){
					$flag= 1;
				}
				else {
					$line= <FASTXIN>;
					print STDOUT $line;
					$line= <FASTXIN>;
					if ($opt_noHdQual){
						$line = "+\n";
					}
					print STDOUT $line;
					$line= <FASTXIN>;
					print STDOUT $line;
				}
			}
			else{
				if ($opt_inFormat=~ m/fasta/){
					$flag= 0;
				}
				else {
					$line= <FASTXIN>;
					$line= <FASTXIN>;
					$line= <FASTXIN>;
				}
			}
		}
		else { # If filter IS set to reverse
			if (!$reads{$id}){
				print STDOUT $line . "\n";
				if ($opt_inFormat=~ m/fasta/){
					$flag= 1;
				}
				else {
					$line= <FASTXIN>;
					print STDOUT $line;
					$line= <FASTXIN>;
					if ($opt_noHdQual){
						$line = "+\n";
					}
					print STDOUT $line;
					$line= <FASTXIN>;
					print STDOUT $line;
				}
			}
			else{
				if ($opt_inFormat=~ m/fasta/){
					$flag=0;
				}
				else {
					$line= <FASTXIN>;
					$line= <FASTXIN>;
					$line= <FASTXIN>;
				}
			}
		}
	}
	else {
		next;
	}
}
close (FASTXIN);

if ($opt_verbose){
	print STDERR "DONE\n";
}


exit (0);
