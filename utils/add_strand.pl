#!/usr/bin/perl
use strict;
use Getopt::Long;
use Bio::DB::Fasta;

my ($LIBRARY_TYPE) = ('fr-secondstrand');

#Reads options from command line
GetOptions( 
	    '-library-type=s' => \$LIBRARY_TYPE,
	    
	    
	   );
die "cat [SAM] | $0 -library-type [fr-firststrand|fr-secondstrand]\n" if($LIBRARY_TYPE ne 'fr-firststrand' and $LIBRARY_TYPE ne 'fr-secondstrand');


while(<STDIN>) {
  chomp;
  if($_=~/^\@/) {
     print $_,"\n";
  } else {
     my $line1 = $_;
     my ($id1,$flag1,$ref1,$pos1,$mapq1,$cigar1,$mrn1,$mpos1,$tlen1,$seq1,$qual1,@opt1) = split /\t/,$line1;
     my $line2 = <STDIN>;
     chomp $line2;
     my ($id2,$flag2,$ref2,$pos2,$mapq2,$cigar2,$mrn2,$mpos2,$tlen2,$seq2,$qual2,@opt2) = split /\t/,$line2;     
     
     next if($ref1 ne $ref2 or $pos1 == $pos2);
     
     my $strand;
     if($pos1 < $mpos1 and $pos2 > $mpos2) {
     	$strand = ($LIBRARY_TYPE ne 'fr-secondstrand') ? '+' : '-';
     } elsif($pos1 > $mpos1 and $pos2 < $mpos2) {
     	$strand = ($LIBRARY_TYPE ne 'fr-secondstrand') ? '-' : '+';
     } else {
     	die "\n",$line1,"\n",$line2;
     }
     unshift @opt1,"XS:A:$strand";
     unshift @opt2,"XS:A:$strand";
     
     print join("\t",($id1,$flag1,$ref1,$pos1,$mapq1,$cigar1,$mrn1,$mpos1,$tlen1,$seq1,$qual1,@opt1)),"\n";
     print join("\t",($id2,$flag2,$ref2,$pos2,$mapq2,$cigar2,$mrn2,$mpos2,$tlen2,$seq2,$qual2,@opt2)),"\n";

  }
} 
