#!/usr/bin/perl

#===============================================================================
#   Author: Alejandro Manzano Marin
#
#   File: gbkTools.pl
#   Date: 19-08-2013
#   Version: 0.1
#
#   Usage:
#      perl gbkTools.pl --infile|-i inFile.gbk --program|-p gbkToolsProgram [options]
#
#      Check out 'perl gbkTools.pl -h' for short usage manual and info on the software.
#
#    Description: This program will consist of a series of subprograms which will perform specific
#                 tasks on Genbank (*.gbk or *.gb) files to rapidly extract information in popular
#                 formats for use as input in other programs.
#                 
#
#    Contact: Contact the author at alejandro.manzano@uv.es using 'gbkTools: ' as
#             as begining for subject for bug reporting, feature request or whatever
#             (of course realted to the software).
#
#    COPYRIGHT: Copyright (C) 2013  Alejandro Manzano-Marin.
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
use Data::Dumper;


# Define subroutines
sub printVersion {
	# Usage: printVersion(<programName>, <version>, <typeglobRef>);
	# Example: printVersion('program_name', '3.2', \*STDERR);
	# Capture arguments
	my $software= $_[0];
	my $version= $_[1];
	my $fileHandle= $_[2];
	
	print $fileHandle "$software v$version\n";
	exit (0);
}

sub fullGbk2Hash {
	# Capture Arguments
	my $gbkFile= $_[0];
	my $program= $_[1];
	my $captFeatures= $_[2];
	
	# Define variables
	my $fileHandle;
	my $locusID= '';
	my $featCount= 0;
	my $featType= '';
	my $posString= '';
	my $featStrand= 1;
	my $featPartialStart= 0;
	my $featStart= 1;
	my $featPartialEnd= 0;
	my $featEnd= 0;
	my %featHash= ();
	my %flag=('printFeature'=> 0,
		'printQualifier'=>0
	);
	my $featQual= '';
	my %locusSeq= ();
	my $i= 0;
	my @tempArray= ();
	
	open ($fileHandle, '<', $gbkFile) || die "Unable to open file $gbkFile for reading\nERROR: $!\n";
		while (my $line= <$fileHandle>){
			chomp $line;
			if ($line=~ m/^LOCUS\s+(\S+)\s+\d+ bp\s+/){
				$locusID= $1;
				$featCount= 0;
				next;
			}
			if ($line=~ m/^\s+(\S+)\s+((join|complement)\()?(((complement|join)\()?<?(\d+)\.\.>?(\d+)\)?,?)+\)?/){
				$featStrand= 1;
				$featType= $1;
				$locationMod= $3;
				$posString= $4;
				if ($locationMod eq 'complement'){
					$featStrand= -1;
				}
				if ($featCount == 0){
					@{$featHash{$locusID}}= ();
				}
				if ($captFeatures eq 'all' || $captFeatures=~ m/$featType/ || ($captFeatures=~ m/pseudo/ && $featType eq 'gene')){
					$flag{'printFeature'}= 1;
					$featCount++;
					%{$featHash{$locusID}[$featCount-1]}= ();
					while ($posString=~ m/((complement|join)\(|)?(<?)(\d+)\.\.(>?)(\d+)\)?/g){
						if ($1 eq 'complement('){
							$featStrand= -1;
						}
#						elsif (){
#							$featStrand= 1;
#						}
						if ($3 eq '<'){
							$featPartialStart= 1;
						}
						else {
							$featPartialStart= 0;
						}
						$featStart= $4;
						if ($5 eq '<'){
							$featPartialEnd= 1;
						}
						else {
							$featPartialEnd= 0;
						}
						$featEnd= $6;
						$featHash{$locusID}[$featCount-1]{'featType'}= $featType;
						push(@{$featHash{$locusID}[$featCount-1]{'pos'}}, [$featStart, $featEnd, $featStrand, $featPartialStart, $featPartialEnd]);
					}
				}
				else {
					$flag{'printFeature'}= 0;
				}
				next;
			}
			if ($flag{'printFeature'} && $line!~ m/^ORIGIN/){
				if ($line=~ m/^\s+\/([^=]+)=?(.+)?/){
					$featQual= $1;
					$featQualDescr= $2;
					if ($featQual eq 'pseudo'){
						$featQual= 'pseudogene';
						$featQualDescr= 'unknown';
					}
					elsif ($featQual && !$featQualDescr){
						$featQualDescr= 1;
					}
					$featQualDescr=~ s/^\"//;
					$featQualDescr=~ s/\"$//;
					if (defined($featHash{$locusID}[$featCount-1]{$featQual})){
						$featHash{$locusID}[$featCount-1]{$featQual}.= ';' . $featQualDescr;
					}
					else {
#					if ($opt_qualifiers eq 'all' || $opt_qualifiers=~ m/$featQual/){
						$featHash{$locusID}[$featCount-1]{$featQual}= $featQualDescr;
#						$featHash{$locusID}[$featCount-1]{$featQual}=~ s/^\"//;
						$featHash{$locusID}[$featCount-1]{$featQual}=~ s/\"$//;
						$flag{'printQualifier'}= 1;
#					}
#					else {
#						$flag{'printQualifier'}= 0;
#					}
					}
				}
				elsif ($flag{'printQualifier'} && $line=~ m/^\s+(\S+.*)/){
#					if ($featQual eq 'translation'){
#						$featHash{$locusID}[$featCount-1]{$featQual}.= $1;
#					}
#					else {
						$featHash{$locusID}[$featCount-1]{$featQual}.= ' ' . $1;
#					}
					$featHash{$locusID}[$featCount-1]{$featQual}=~ s/\"$//;
				}
				next;
			}
			if ($line=~ m/^ORIGIN/){
				$flag{'printFeature'}= 0;
				$flag{'printQualifier'}= 0;
				if ($program eq 'gbk2fasta' || $program eq 'gbk2circos'){
					$locusSeq{$locusID}= '';
					while ($line!~ m/^\/\//){
						$line= <$fileHandle>;
						chomp $line;
						if ($line=~ m/\s+\d+(( [a-z]+)+)/){
							$locusSeq{$locusID}.= $1;
							$locusSeq{$locusID}=~ s/\s+//g;
						}
					}
				}
			}
		}
	close ($fileHandle);
	
#	if ($captFeatures=~ m/pseudo/ && $captFeatures!~ m/gene/){
		foreach $locusID (sort {lc($a) cmp lc($b)} keys %featHash){
			for ($i=0; $i<scalar(@{$featHash{$locusID}}); $i++){
				if ($captFeatures=~ m/pseudo/ && $captFeatures!~ m/gene/ && $featHash{$locusID}[$i]{'featType'} eq 'gene' && !defined($featHash{$locusID}[$i]{'pseudogene'})){
					splice (@{$featHash{$locusID}}, $i, 1);
					$i--;
					next;
				}
				### Flatten multiposition feature and remove partial
				for ($j=0; $j<scalar(@{$featHash{$locusID}[$i]{'pos'}}); $j++){
					splice (@{$featHash{$locusID}}, $i+$j+1, 0, $featHash{$locusID}[$i]);
					$featHash{$locusID}[$i+$j+1]{'start'}= $featHash{$locusID}[$i]{'pos'}[$j][0];
					$featHash{$locusID}[$i+$j+1]{'end'} = $featHash{$locusID}[$i]{'pos'}[$j][1];
					$featHash{$locusID}[$i+$j+1]{'strand'}= $featHash{$locusID}[$i]{'pos'}[$j][2];
					$featHash{$locusID}[$i+$j+1]{'partialFive'}= $featHash{$locusID}[$i]{'pos'}[$j][3];
					$featHash{$locusID}[$i+$j+1]{'partialThree'}= $featHash{$locusID}[$i]{'pos'}[$j][4];
					delete ($featHash{$locusID}[$i+$j+1]{'pos'});
				}
				splice (@{$featHash{$locusID}}, $i, 1);
				$i+=$j-1;
#				if (defined($featHash{$locusID}[$i+$j+1]{'translation'})){
#					$featHash{$locusID}[$i+$j+1]{'translation'}=~ s/\s+//g;
#				}
			}
			if (scalar(@{$featHash{$locusID}}) == 0){
				delete ($featHash{$locusID});
			}
		}
#	}
	
#	print Dumper %featHash; $pene=<STDIN>;
	if ($program eq 'gbk2tab'){
		return (\%featHash);
	}
	if ($program eq 'gbk2fasta'){
		return (\%featHash, \%locusSeq);
	}
}

sub printFasta {
	my $header= $_[0];
	my $sequence= $_[1];
	my $charPerLine= $_[2];
	my $fileHandle= $_[3];
	
	my $i= 0;
	
	print $fileHandle ">" . $header . "\n";
	if ($charPerLine < 0){
		print STDERR "printFasta ERROR: Illegal number of characters per line, number must be >=0\n";
		return (1);
	}
	if ($charPerLine == 0){
		print $fileHandle $sequence . "\n";
	}
	else {
		for ($i=0; $i<length($sequence); $i+=70){
			if (length(substr($sequence, $i, length($sequence)-$i)) < $charPerLine){
				print $fileHandle substr($sequence, $i, length($sequence)-$i)."\n";
			}
			else {
				print $fileHandle substr($sequence, $i, $charPerLine)."\n";
			}
		}
	}
	return (0);
}

sub revCompRef {
	my $seqRef= $_[0];
	$$seqRef= scalar(reverse($$seqRef));
	$$seqRef=~ tr/atgcATGC/tacgTACG/;
	return (0);
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
my $PROGRAMNAME= "gbkTools";
my $VERSION= "1.0";

## Define options default values
my @opt_inFiles= ();
my $opt_program= '';

my $opt_features='all';
my $opt_qualifiers= 'locusID;locus_tag;gene,locus_tag;description;location;strandString';
my $opt_noHdQual= 0;

my $opt_exactMatch= 0;
my $opt_revFilter= 0;

my $opt_verbose= 0;
my $opt_man= 0;
my $opt_help= 0;
my $opt_printVersion= 0;

## Define options hash
GetOptions(\%opts, 
	'inFile|i=s' => \$opt_inFile, 
	'program|p=s' => \$opt_program, 
	'features:s'=> \$opt_features, 
	'qualifiers:s'=> \$opt_qualifiers, 
	'writeSeq:s'=> \$opt_writeSeq, 
	'outPrefix|o'=> \$opt_outPrefix, 
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

gbkTools

=head1 VERSION

gbkTools v1.0

=head1 SYNOPSIS

gbkTools.pl --infile|-i inFile.gbk --program|-p gbkToolsProgram [options] [--verbose|-v] [--help|-h]
[--man] [--version]

=head1 DESCRIPTION

This program will consist of a series of subprograms which will perform specific tasks on Genbank (*.gbk or *.gb)
files to rapidly extract information in popular formats for use as input in other programs.

=head1 OPTIONS

=head2 INPUT

=over 8

=item B<-i> | B<--infile> <string> (mandatory)

File containing genbank (or multigenbak) to processed.

=item B<-p> | B<--program> <string> (mandatory)

Subprogram to be used to process Genbank file (gbk2tab, gbk2fasta, gbk2ctg).

=back

=head2 gbk2tab

=over 8

=item B<--features> <string> (default: 'all')

String of aliases* or names (http://www.insdc.org/files/feature_table.html) of features to be printed.
List must be separated by ';' and delimited by 'list' (single-quotation marks).

* Aliases

	all: Print ALL features present in the genbank, including 'source', 'assembly_gap', etc...
	locusID: Print info for LOCUS of the genbank (You can think of them as genomes, contigs, scaffolds, etc..).
	pseudo: Print ALL gene qualifiers with /pseudo or /pseudogene qualifier.

=item B<--qualifiers> <string> (default: 'locus_tag,gene;locus_tag;description;location;strandString')

String of the aliases* or names (http://www.insdc.org/files/feature_table.html) of qualifiers to print
for all features selected to print out. List of names of qualifiers must be separated by ';' and
delimited as 'list' (single-quotation marks). A ',' may be used to indicate to print a certain qualifier if
present, if not print the next one ('qual1;qual2,altQual2,..;qual3;..').

* Aliases

	description: If 'product' qualifier available print it, if not print note. For pseudogenes print note only.
	location: Print qualifier of location in 'artemis' format 'start:end'.
	locusID: Print info for the LOCUS that a feature belongs to. You can think of them as genomes, contigs, scaffolds, etc...
	size: print size of feature in bp.
	strandString: Prinf strand of feature as 'forward'/'reverse'.
	strandNumber: Prinf strand of feature as '1'/'-1'.

=back

=head2 gbk2fasta

=over 8

=item B<--features> <string> (default: 'all')

String of aliases* or names (http://www.insdc.org/files/feature_table.html) of features to be printed.
List must be separated by ';' and delimited by 'list' (single-quotation marks).

* Aliases

	all: Print ALL features present in the genbank, including 'source', 'assembly_gap', etc...
	pseudo: Print ALL gene qualifiers with /pseudo or /pseudogene qualifier.

=item B<--qualifiers> <string> (default: 'locus_tag,gene;locus_tag;description;location;strandString')

String of the aliases* or names (http://www.insdc.org/files/feature_table.html) of qualifiers to print
for all features selected to print out. List of names of qualifiers must be separated by ';' and
delimited as 'list' (single-quotation marks). A ',' may be used to indicate to print a certain qualifier
if present, if not print the next one ('qual1,qual2|otherQual2|..,qual3,..').

* Aliases

	description: If 'product' qualifier available print it, if not print note. For pseudogenes print note only.
	location: Print feature location in 'artemis' format 'start:end'.
	strand: Prinf strand of feature as '1'/'-1'.
	strandString: Prinf strand of feature as 'forward'/'reverse'.
	

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

Alejandro Manzano-Marin, C<< <alejandro_dot_manzano_at_uv_dot_es> >>

=head1 BUGS

If you find a bug please email me at C<< <alejandro_dot_manzano_at_uv_dot_es> >> so that I can keep making
gbkTools better.

=head1 COPYRIGHT

Copyright (C) 2013  Alejandro Manzano Marin.

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut


# Check options
if (!$opt_inFile || !$opt_program){
	print STDERR "ERROR:\n";
	if (!$opt_inFile){
		print STDERR "Genbank infile missing\n";
	}
	if (!$opt_program){
		print STDERR "Program to process genbank not specified\n";
	}
	print STDERR "Please check usage manual\n";
	pod2usage(-exitval => 1,  -verbose => 1);
}


# Assign values to variables dependant on options


# Main program

## Read list of accepted reads
if ($opt_verbose){
	print STDERR "Reading Genbank and storing information" . $opt_listFile . "..";
}

if ($opt_program eq 'gbk2tab' || $opt_program eq 'gbk2fasta'){
	### Define variables for subscript
#	@featArray= split (',', $opt_features);
	my $featHashRef= 0;
	my $locusSeqRef= 0;
	my @qualArray= split (';', $opt_qualifiers);
	my $locusID= '';
	my $i= 0;
	my $qualString= '';
	my $j= 0;
	my $qualName= '';
	my @atlQual= ();
	my $k= 0;
	my $featSeq= '';
	
	### Read and process Genbank file
	if ($opt_program eq 'gbk2tab'){
		($featHashRef)= fullGbk2Hash($opt_inFile, $opt_program, $opt_features);
	}
	if ($opt_program eq 'gbk2faa'){
		($featHashRef)= fullGbk2Hash($opt_inFile, $opt_program, $opt_features);
	}
	elsif ($opt_program eq 'gbk2fasta'){
		($featHashRef, $locusSeqRef)= fullGbk2Hash($opt_inFile, $opt_program, $opt_features);
	}

	### Substitute some aliases in @qualArray
	for ($i=0; $i<scalar(@qualArray); $i++){
		if ($qualArray[$i] eq 'description'){
			$qualArray[$i]= 'product,note';
		}
	}
	
	### Print tabular to STDOUT
	foreach $locusID (sort {lc($a) cmp lc($b)} keys %{$featHashRef}){
		for ($i=0; $i<scalar(@{$featHashRef->{$locusID}}); $i++){
			$qualString= '';
#			print "Feature Number $i\n";
#			print Dumper %{$featHashRef->{$locusID}[$i]};
			for ($j=0; $j<scalar(@qualArray); $j++){
#				print "Original Qualifier: $qualArray[$j]\n";
				if ($j){
					$qualString.= "\t";
				}
				$qualName= '';
				@atlQual= ();
				if ($qualArray[$j]=~ m/,/){
					@altQual= split(/,/, $qualArray[$j]);
					for ($k=0; $k<scalar(@altQual); $k++){
						$qualName= $altQual[$k];
#						print "\t-Testing $qualName\n";
						if (defined($featHashRef->{$locusID}[$i]{$qualName})){
							last;
						}
					}
				}
				else {
					$qualName= $qualArray[$j];
				}
#				print "Chosen Qualifier: $qualName\n";
				if ($qualName eq 'location'){
					$qualString.= $featHashRef->{$locusID}[$i]{'start'} . ':' . $featHashRef->{$locusID}[$i]{'end'};
#					print "Qualifier value: " . $featHashRef->{$locusID}[$i]{'start'} . ':' . $featHashRef->{$locusID}[$i]{'end'} . "\n";# $pene=<STDIN>;
				}
				elsif ($qualName eq 'strandString'){
					if ($featHashRef->{$locusID}[$i]{'strand'}==1){
						$qualString.= 'forward';
#						print "Qualifier value: forward\n"; #$pene=<STDIN>;
					}
					else {
						$qualString.= 'reverse';
#						print "Qualifier value: reverse\n"; #$pene=<STDIN>;
					}
				}
				elsif ($qualName eq 'sizeNuc'){
					$qualString.= ($featHashRef->{$locusID}[$i]{'end'}-$featHashRef->{$locusID}[$i]{'start'}+1);
				}
				elsif ($qualName eq 'locusID'){
					$qualString.= $locusID;
				}
				else {
					$qualString.= $featHashRef->{$locusID}[$i]{$qualName};
#					print "Qualifier value: " . $featHashRef->{$locusID}[$i]{$qualName} . "\n"; #$pene=<STDIN>;
				}
			}
			if (defined($featHashRef->{$locusID}[$i]{'pseudo'})){
				$qualString.= "\t" . 'pseudogene';
			}
			elsif (defined($featHashRef->{$locusID}[$i]{'pseudogene'})){
				$qualString.= "\t" . 'pseudogene=' . $featHashRef->{$locusID}[$i]{'pseudogene'};
			}
			
			if ($opt_program eq 'gbk2tab'){
				print STDOUT $qualString . "\n";
			}
			elsif ($opt_program eq 'gbk2fasta'){
				$featSeq= substr ($locusSeqRef->{$locusID}, ($featHashRef->{$locusID}[$i]{'start'}-1), ($featHashRef->{$locusID}[$i]{'end'}-$featHashRef->{$locusID}[$i]{'start'}+1));
				if ($featHashRef->{$locusID}[$i]{'strand'}==-1){
					revCompRef(\$featSeq);
				}
				printFasta($qualString, $featSeq, 70, \*STDOUT);
			}
			elsif ($opt_program eq 'gbk2faa'){
				printFasta($qualString, $featHashRef->{$locusID}[$i]{'translation'}, 70, \*STDOUT);
			}
		}
	}
}

#if ($opt_program eq 'gbk2circos'){
#	
#	### Define variables for subscript
#	my $featHashRef= 0;
#	my $locusSeqRef= 0;
#	my @qualArray= split (';', $opt_qualifiers);
#	my $locusID= '';
#	my $i= 0;
#	my $qualString= '';
#	my $j= 0;
#	my $qualName= '';
#	my @atlQual= ();
#	my $k= 0;
#	my $featSeq= '';
#	
#	### Read and process Genbank file
#	($featHashRef, $locusSeqRef)= fullGbk2Hash($opt_inFile, $opt_program, $opt_features);
#	
#	## Create files and folders
#	if (-d "$opt_outPrefix"){
#		unlink (glob("$opt_outPrefix/tracks/*"));
#		rmdir ("$opt_outPrefix/tracks"); 
#		rmdir ("$opt_outPrefix");
#	}
#	
#	mkdir ("$opt_outPrefix");
#	mkdir ("$opt_outPrefix/tracks");

#	open (CHROUT, ">$opt_outPrefix/tracks/chrs_karyotype.txt") || die "Unable to open file $opt_outPrefix/tracks/chrs_karyotype.txt for writing\n$!\n";

#	open (CDSFWDOUT, ">$opt_outPrefix/tracks/chrs_CDS_fwd.txt") || die "Unable to open file $opt_outPrefix/tracks/chrs_CDS_fwd.txt for writing\n$!\n";
#	open (CDSREVOUT, ">$opt_outPrefix/tracks/chrs_CDS_rev.txt") || die "Unable to open file $opt_outPrefix/tracks/chrs_CDS_rev.txt for writing\n$!\n";
#	open (RNAFWDOUT, ">$opt_outPrefix/tracks/chrs_RNA_fwd.txt") || die "Unable to open file $opt_outPrefix/tracks/chrs_RNA_fwd.txt for writing\n$!\n";
#	open (RNAREVOUT, ">$opt_outPrefix/tracks/chrs_RNA_rev.txt") || die "Unable to open file $opt_outPrefix/tracks/chrs_RNA_rev.txt for writing\n$!\n";
#	open (RPTOUT, ">$opt_outPrefix/tracks/chrs_RPT.txt") || die "Unable to open file $opt_outPrefix/tracks/chrs_RPT.txt for writing\n$!\n";

#	if ($opt_extraFeat){
#		open (EXTFWDOUT, ">$opt_outPrefix/tracks/chrs_EXT_fwd.txt") || die "Unable to open file $opt_outPrefix/tracks/chrs_EXT_fwd.txt for writing\n$!\n";
#		open (EXTREVOUT, ">$opt_outPrefix/tracks/chrs_EXT_rev.txt") || die "Unable to open file $opt_outPrefix/tracks/chrs_EXT_rev.txt for writing\n$!\n";
#	}
#	if ($opt_gcSkew){
#		open (GCSOUT, ">$opt_outPrefix/tracks/chrs_GCS.txt") || die "Unable to open file $opt_outPrefix/tracks/chrs_pos.txt for writing\n$!\n";
#	}
#	
#	### Print circos track files
#	foreach $locusID (sort {lc($a) cmp lc($b)} keys %{$featHashRef}){
#		for ($i=0; $i<scalar(@{$featHashRef->{$locusID}}); $i++){
#			$qualString= '';
##			print "Feature Number $i\n";
##			print Dumper %{$featHashRef->{$locusID}[$i]};
##			for ($j=0; $j<scalar(@qualArray); $j++){
###				print "Original Qualifier: $qualArray[$j]\n";
##				if ($j){
##					$qualString.= "\t";
##				}
##				$qualName= '';
##				@atlQual= ();
##				if ($qualArray[$j]=~ m/,/){
##					@altQual= split(/,/, $qualArray[$j]);
##					for ($k=0; $k<scalar(@altQual); $k++){
##						$qualName= $altQual[$k];
###						print "\t-Testing $qualName\n";
##						if (defined($featHashRef->{$locusID}[$i]{$qualName})){
##							last;
##						}
##					}
##				}
##				else {
##					$qualName= $qualArray[$j];
##				}
###				print "Chosen Qualifier: $qualName\n";
##				if ($qualName eq 'location'){
##					$qualString.= $featHashRef->{$locusID}[$i]{'start'} . ':' . $featHashRef->{$locusID}[$i]{'end'};
###					print "Qualifier value: " . $featHashRef->{$locusID}[$i]{'start'} . ':' . $featHashRef->{$locusID}[$i]{'end'} . "\n";# $pene=<STDIN>;
##				}
##				elsif ($qualName eq 'strandString'){
##					if ($featHashRef->{$locusID}[$i]{'strand'}==1){
##						$qualString.= 'forward';
###						print "Qualifier value: forward\n"; #$pene=<STDIN>;
##					}
##					else {
##						$qualString.= 'reverse';
###						print "Qualifier value: reverse\n"; #$pene=<STDIN>;
##					}
##				}
##				elsif ($qualName eq 'sizeNuc'){
##					$qualString.= ($featHashRef->{$locusID}[$i]{'end'}-$featHashRef->{$locusID}[$i]{'start'}+1);
##				}
##				elsif ($qualName eq 'locusID'){
##					$qualString.= $locusID;
##				}
##				else {
##					$qualString.= $featHashRef->{$locusID}[$i]{$qualName};
###					print "Qualifier value: " . $featHashRef->{$locusID}[$i]{$qualName} . "\n"; #$pene=<STDIN>;
##				}
##			}
#			if ($qual){
#			
#			}
#			if (defined($featHashRef->{$locusID}[$i]{'pseudo'})){
#				$qualString.= "\t" . 'pseudogene';
#			}
#			elsif (defined($featHashRef->{$locusID}[$i]{'pseudogene'})){
#				$qualString.= "\t" . 'pseudogene=' . $featHashRef->{$locusID}[$i]{'pseudogene'};
#			}
#			
#			if ($opt_program eq 'gbk2tab'){
#				print STDOUT $qualString . "\n";
#			}
#			elsif ($opt_program eq 'gbk2fasta'){
#				$featSeq= substr ($locusSeqRef->{$locusID}, ($featHashRef->{$locusID}[$i]{'start'}-1), ($featHashRef->{$locusID}[$i]{'end'}-$featHashRef->{$locusID}[$i]{'start'}+1));
#				if ($featHashRef->{$locusID}[$i]{'strand'}==-1){
#					revCompRef(\$featSeq);
#				}
#				printFasta($qualString, $featSeq, 70, \*STDOUT);
#			}
#		}
#	}
#}

if ($opt_verbose){
	print STDERR "DONE\n";
}


exit (0);
