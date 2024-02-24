## filter the fastq files based on the matched reads in a sam file
## for filtering out the mouse reads in a pdx samples

use strict;
use Bio::SeqIO;

# get the file name, somehow 
# my $file = shift; 
# my $file_reads = shift;
my $file_ids  = shift;

# for test
$file_reads 	= '../exonseq/Project_07813_X/JAX_0424/Sample_DS_bla_221_X_IGO_07813_X_1/DS_bla_221_X_IGO_07813_X_1_S37_R1_001.fastq.gz';
$file_ids   	= './r_fang/initFiles/UCC15X/UCC15X_ucc15x1s6_hsa_id.tsv';
$file_out	= $file_reads =~ s/(_R1_.*.fastq.gz)/_FIL_\1/;

my $fin; 
my %hashid;
open($fin, '<', $file_ids) or die('ID file: $file_ids can not opened');
while(<>){
	chomp;
	$hashid{$_} = 1
}

my $seqin = Bio::SeqIO->new(-fh => *STDIN, -format => 'fastq'); 
my $seqout = Bio::SeqIO->new(-fh => *STDOUT, -format => 'fastq'); 

while(my $seqi = $seqio -> next()){
	if($hashid{$seqi -> display_id} > 0){
		$seqout -> write_seq($seqi);
	}
}
