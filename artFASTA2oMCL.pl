#!/usr/local/bin/perl
# Change artemis' FASTA file and make a FASTA compliant with orthoMCL

use Getopt::Long;

GetOptions(\%opts, "infile|i=s@", "org:s", "db=s", "MCL!", "list|l!", "hyp!", "help|h!");
        sub usage(){
                die "USAGE :: head_trans.pl -infile string -outfile string -org string\n\n
                infile|i\t\tFile with artemis gene prediction [String]\n
                org\t\tName of organism (only use if db=NCBI) [String]\n
                db\t\tString identifying database fotmat, either NCBI or Artemis output[String]\n
                MCL\t\tMake an orthoMCL fasta format to file\n
                list|l\t\tMake a list format file\n
                help\t\tThis information\n\n
                ";
        }

if ($opts{'help'} || !$opts{'infile'} || !$opts{'db'}){
	&usage;
}

#####Capture all gene ID's with gene start and gene end positions including strand#####
$flag=0;

if ($opts{'MCL'}){
	open (MCL, ">$opts{'org'}.fasta") or die ("Unable to open file for writting: $opts{'org'}.fasta\n$!\n");
	if ($opts{'hyp'}){
		open (HYPMCL, ">$opts{'org'}.fasta.hyp") or die ("Unable to open file for writting: $opts{'org'}.fasta.hyp\n$!\n");
	}
}
if ($opts{'list'}){
	open (LIST, ">$opts{'org'}.lst") or die ("Unable to open file for writting: $opts{'list'}\n$!\n");
	if ($opts{'hyp'}){
		open (HYPLIST, ">$opts{'org'}.lst.hyp") or die ("Unable to open file for writting: $opts{'list'}.hyp\n$!\n");
	}
}

foreach $file (@{$opts{'infile'}}){
	$flag{'hyp'}=0;
	open (IN, "$file") or die ("Unable to open file for reading: $file\n$!\n");
	while($line=<IN>){
		chomp $line;
		if ($line=~ /^#/){	#Skip commented lines
			next;
		}
		if ($line=~/^>/){
			if ($opts{'hyp'} && ($line=~ /hypothetical/i || $line=~ /conserved protein/i)){	#Skip commented lines
				$flag{'hyp'}=1;
			}
			elsif ($opts{'hyp'} && $line!~ /hypothetical/){
				$flag{'hyp'}=0;
			}
	#		$line=~s/>//;
			#####Use for NCBI RefSeq format#####
			if ($opts{'db'} eq "NCBI"){
				$line=~s/\s\[.*\]//;
				@temp=split(/\s+/, $line);
				$gi=shift(@temp);
				$gi=~ s/^>gi\|//;
				$gi=~ s/\|.*//;
				if ($opts{'MCL'} && !$flag{'hyp'}){
					print MCL ">$opts{'org'}|$gi\n";
				}
				elsif ($opts{'MCL'} && $flag{'hyp'}){
					print HYPMCL ">$opts{'org'}|$gi\n";
				}
				if ($opts{'list'} && !$flag{'hyp'}){
					print LIST "$opts{'org'}|$gi|".join("_", @temp)."\n";
				}
				elsif ($opts{'list'} && $flag{'hyp'}){
					print HYPLIST "$opts{'org'}|$gi|".join("_", @temp)."\n";
				}
				undef @temp;
			}
			#####Use for artemis output of fasta sequences#####
			else {
				$line=~ s/\s+\d+:\d+.*$//;
				@temp=split(/\s+/, $line);
				$id=shift(@temp);
				$id=~s/^>//;
				$gene_name=shift(@temp);
				if ($opts{'MCL'} && !$flag{'hyp'}){	
					print MCL ">$opts{'org'}|$id\n";
				}
				elsif ($opts{'MCL'} && $flag{'hyp'}){	
					print HYPMCL ">$opts{'org'}|$id\n";
				}
				if ($opts{'list'} && !$flag{'hyp'}){
					print LIST "$opts{'org'}|$id|$gene_name|".join("_", @temp)."\n";
				}
				elsif ($opts{'list'} && $flag{'hyp'}){
					print HYPLIST "$opts{'org'}|$id|$gene_name|".join("_", @temp)."\n";
				}
				undef @temp;
				undef $id;
				undef $gene_name;
			}
		}
		if ($opts{'MCL'} && $line!~/^>/ && !$flag{'hyp'}){
			print MCL $line."\n";
		}
		elsif ($opts{'MCL'} && $line!~/^>/ && $flag{'hyp'}){
			print HYPMCL $line."\n";
		}
	}
	close(IN);
}


if ($opts{'MCL'}){
	close (MCL);
	if ($opts{'hyp'}){
		close (HYPMCL);
	}
}
if ($opts{'list'}){
	close (LIST);
	if ($opts{'hyp'}){
		close (HYPLIST);
	}
}

