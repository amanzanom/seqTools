# Script that takes an aminoacid seuqence FASTA file and creates all versions with a single point mutation
﻿%aa20=("A" => 1,
	"R" => 1,
	"N" => 1,
	"D" => 1,
	"C" => 1,
	"E" => 1,
	"Q" => 1,
	"G" => 1,
	"H" => 1,
	"I" => 1,
	"L" => 1,
	"K" => 1,
	"M" => 1,
	"F" => 1,
	"P" => 1,
	"S" => 1,
	"T" => 1,
	"W" => 1,
	"Y" => 1,
	"V" => 1);
%seq=();
$id='';
open (IN, "$ARGV[0]");
while ($line=<IN>){
	chomp $line;
	if ($line =~ m/^>(\S+)/ || eof){
		if (length($seq{$id})>0){
			@seqArray= split(//, $seq{$id});
			$seqTemp= $seq{$id};
			for ($i=0; $i<scalar(@seqArray); $i++){
				if ($i == 0){ # Comment if wish to change first amino acid
					next;
				}
				%aaTemp=%aa20;
				delete $aaTemp{$seqArray[$i]};
				foreach $aminoAcid (sort {lc($a) cmp lc($b)} keys %aaTemp){
					print STDOUT ">$id|$seqArray[$i]". ($i+1) . "$aminoAcid\n";
					$seqTemp= substr($seqTemp, $i, 1, $aminoAcid);
					for ($j=0; $j<length($seqTemp); $j+=70){
						if (length(substr($seqTemp, $j, length($seqTemp)-$j))<70){
							print STDOUT substr($seqTemp, $j, length($seqTemp)-$j)."\n";
						}
						else {
							print STDOUT substr($seqTemp, $j, 70)."\n";
						}
					}
				}
			}
		}
		$id= $1;
		$seq{$id}='';
		next;
	}
	$seq{$id}.= $line;
}
close (IN);
