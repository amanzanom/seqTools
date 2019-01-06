#!/usr/local/bin/perl -w

use Getopt::Long;

GetOptions(\%opts, "infile|i=s@", "format|f=s", "help|h!");
        sub usage(){
                die "USAGE :: catAln.pl -infile|i aln.file -format format_ext [-help]\n\n
                -infile|i\t\tAlignment file(s) (IMPORTANT: files must be input in concatenation order) [String]\n
                -format|f\t\tFormat of the alignment (fasta, clustal, etc...) [String]\n
                -help|h\t\tPrint usage manual\n\n
                ";
        }

sub printFasta {
	my $header= $_[0];
	my $seq= $_[1];
	my $i= 0;
	print ">$header\n";
	for ($i=0; $i<length($seq); $i+=70){
		if (length(substr($seq, $i, length($seq)-$i))<70){
			print substr($seq, $i, length($seq)-$i)."\n";
		}
		else {
			print substr($seq, $i, 70)."\n";
		}
	}
	return (0);
}

sub captureFastaHash {
	local (*FILE)= $_[0];
	my %seq= %{$_[1]};
	my $header= '';
	my $line='';
	while ($line=<FILE>){
		chomp $line;
		if ($line=~ /^>(.+)/){
			$header=$1;
			if (!$seq{$header}){
				$seq{$header}='';
			}
			next;
		}
		$seq{$header}.=$line;
	}
	return (\%seq);
}

if ($opts{'help'} || !$opts{'infile'} || !$opts{'format'}){
	if (!$opts{'infile'}){
		print "Alignment file missing, please check usage manual\n";
	}
	if (!$opts{'format'}){
		print "Format of alignment not specified, please check usage manual\n";
	}
	&usage;
}

$header='';
$seqRef='';
%seq=();
foreach $file (@{$opts{'infile'}}){
	open (IN, "$file") || die ("Unable to open: Group_genes.txt\n$!\n");
	if ($opts{'format'} eq "fasta"){
		$seqRef= &captureFastaHash (*IN, \%seq);
	}
	%seq=%$seqRef;
	close (IN);
}

foreach $header (sort {lc($a) cmp lc($b)} keys %seq){
	&printFasta($header, $seq{$header});
}
