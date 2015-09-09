#!/usr/bin/perl
#
# Usage: make_pscores.pl SEQFILE SCOREFILE
#
# Runs BLAST and writes output to a file "SEQFILE.out".
# This file is parsed into lines of the form "seqname1 \t seqname2 \t bitscore"
# (suitable for input to POA), which are written to SCOREFILE.
#
# NB: The bitscore increases with increasing sequence similarity.
#

$seq_file = $ARGV[0];
$pscore_file = $ARGV[1];

open(PSCORE_OUT, ">$pscore_file");
system("./formatdb -i $seq_file -p T");
system("./blastall -p blastp -d $seq_file -i $seq_file -M BLOSUM80 -o $seq_file.out");

open(BLAST_OUT, "<$seq_file.out");
while(<BLAST_OUT>){
     @my_parse = split(/\s/, $_);
     if ($my_parse[0] =~ /^Query=/){
	 $seq_name1 = $my_parse[1];
	 while(<BLAST_OUT>){         
	    if ($_ =~ />/){
		last;
            }
            if ($_ =~ /bits/){
		while(<BLAST_OUT>){
		    if ($_ =~ />/){
			last;
		    }
                    @other_parse = split(/\s+/, $_);
                    $seq_name2 = $other_parse[0]; 
	            $bit_score = $other_parse[1];
                    if ($seq_name2 ne ""){
                       printf PSCORE_OUT "$seq_name1\t$seq_name2\t$bit_score.0\n";
		   }
                }
               if ($_ =~ />/){
		last;
               } 
            }
	 }
     }
}
close(BLAST_OUT);
close(PSCORE_OUT);
