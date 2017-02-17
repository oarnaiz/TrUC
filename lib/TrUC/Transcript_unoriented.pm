package TrUC::Transcript_unoriented;
use strict;
use base 'TrUC::Transcript';

use TrUC::Config;
use File::Basename;
use FindBin qw($Bin);

use Bio::GFF3::LowLevel qw/ gff3_parse_feature /;

=head1 NAME

 TrUC TSS module - Predict Transcription Start Sites (Need CapSeq mapping)

=head1 AUTHORS

2017, I2BC - CNRS , GNU GPL3

=cut

# ALL PARAMETERS FOR THIS MODULE
my %PARAMETERS = ();


=head2 new

 Title   : new
 Usage   : my $factory = PARTIE::MIRET->new({....})
 Function: Create a factory object
 Returns : the factory object

=cut

sub new {
  my ($class,$ref) = @_;
  
  # Check the class
  $class = ref($class) || $class;

  # Link between the object and the class
  my $self = bless {},$class;
  
  # init parent class
  $self = $self->SUPER::new($ref);
  
  # load parameters
  foreach my $param (sort keys %PARAMETERS) {
     $self->{$param} = ((defined $ref->{"-".lc($param)} and $ref->{"-".lc($param)} ne '') and (ref($ref->{"-".lc($param)}) ne 'ARRAY' or scalar @{$ref->{"-".lc($param)}} !=0)) ? $ref->{"-".lc($param)} : $PARAMETERS{$param}->{DEFAULT}; 
	die "$param should be TRUE or FALSE not : $self->{$param}\n" if($PARAMETERS{$param}->{TYPE} eq 'BOOLEAN' and $self->{$param} ne 'TRUE' and $self->{$param} ne 'FALSE');
  }   
  return $self;
}



=head2 calculate

 Title   : Calculate function
 Usage   : $factory->calculate($seq);
 Function: 
 Returns : Nothing
 Args    : Bio::Seq object

=cut


sub calculate {
  my ($self,$seq) = @_;  
  
   my $seq_id = $seq->id;  
   $self->stderr("Calculation $seq_id\n");
   
   my $seq_length = $seq->length;
   my $sequence = $seq->seq;


   my %INTRON_POSITIONS;
   my %INTRON_COVERAGE;
   my %INSERT_COVERAGE;
   
   
   my @SAM;
   foreach my $bam_file (@{$self->{BAM}}) {
      my $fname = basename($bam_file);
      $self->stderr("Read bam file $fname for $seq_id ... \n");      
      push @SAM,  Bio::DB::Sam->new(-bam  =>$bam_file, -fasta=> $self->{GENOME} ); 
   }
    
   foreach my $sam (@SAM) {   
      foreach my $pair  ($sam->features(-type   => 'read_pair', -seq_id => $seq_id
      #,-start=>197000,-end=>210000
      )) {
       
         my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
	 next if(!$first_mate or !$second_mate);
	 next if(!$self->_is_unique($first_mate,$second_mate)) ;
	 
	 my @pos = sort {$a<=>$b} ($first_mate->start,$first_mate->end,$second_mate->start,$second_mate->end);
	 my ($insert_start,$insert_end) = ($pos[0],$pos[$#pos]);
	 next if(abs($insert_start-$insert_end) > $self->{MAX_LENGTH_BTW_PAIRS});
	 
	 my $strand = '.';
	 for(my $i=$insert_start;$i<=$insert_end;$i++) {
            $INSERT_COVERAGE{$strand}->{$i}++;
         }
	 
	 foreach my $read ($first_mate,$second_mate) {
	    my $cigar_str = $read->cigar_str;
	    # find introns positions
	    ########################
            if($cigar_str=~/^(\d+)M(\d+)N(\d+)M$/) {
	    
               my ($left_match,$gap,$right_match)= ($1,$2,$3);
               if($gap >= $self->{MIN_INTRON_LENGTH} and $gap <= $self->{MAX_INTRON_LENGTH})  {
	    	  my ($intron_start,$intron_end) = ($read->start+$left_match,$read->start+$left_match+$gap-1);
		  
		  my $intron_consensus = $INTRON_POSITIONS{$strand}->{"$intron_start-$intron_end"}->{INTRON_CONSENSUS};
		  
		  $INTRON_POSITIONS{$strand}->{"$intron_start-$intron_end"}->{COVERAGE}++;
		  $INTRON_POSITIONS{$strand}->{"$intron_start-$intron_end"}->{INTRON_CONSENSUS} = $intron_consensus;
		  
		  for(my $i=$intron_start;$i<=$intron_end;$i++) {
                     $INTRON_COVERAGE{$strand}->{$i}++;
                  }				     
	       }
            }
         }
      } 			     
   }

   # FIND INTRONS
   my @INTRONS;
   foreach my $strand (keys %INTRON_COVERAGE) {
      my $start_covtig;      
      my @intron_positions = keys %{$INTRON_POSITIONS{$strand}};      
      for(my $i=1; $i<=($seq_length+1); $i++) {
         my $cov = $INTRON_COVERAGE{$strand}->{$i} ? $INTRON_COVERAGE{$strand}->{$i} : 0; 
	 if($cov < 1 and $start_covtig) {
	    my $end_covtig = $i-1;
	    my %intron;
	    foreach my $intron_position (@intron_positions) {
	       my ($intron_start,$intron_end) = split /-/,$intron_position;
	       next if($intron_end < $start_covtig or $intron_start > $end_covtig);
	       if(!%intron 
	       	    or $intron{COVERAGE} < $INTRON_POSITIONS{$strand}->{$intron_position}->{COVERAGE}
		    or ($intron{INTRON_CONSENSUS} == 0 and $INTRON_POSITIONS{$strand}->{$intron_position}->{INTRON_CONSENSUS} == 1)) {		    
		   
		   my $intron_seq = uc(substr($sequence,$intron_start-1,$intron_end-$intron_start+1));
		   
		   my $deduced_strand = $strand;
		   if($self->{INTRON_CONSENSUS} eq 'TRUE') {
		      if($intron_seq=~/^GT.+AG$/) {
		         $deduced_strand="+";
		      } elsif($intron_seq=~/^CT.+AC$/) {
		         $deduced_strand="-";
		      }
		   }
		   
		   
		   %intron = ( POSITIONS => $intron_position, STRAND => $deduced_strand, START => $intron_start , END =>$intron_end,
		   		COVERAGE =>$INTRON_POSITIONS{$strand}->{$intron_position}->{COVERAGE},
		   		INTRON_CONSENSUS =>$INTRON_POSITIONS{$strand}->{$intron_position}->{INTRON_CONSENSUS},
				SPLICING_RATE=>'NA');
	       }
	       
	    }
	    
	    
	    if($self->{MIN_SPLICING_RATE} ne '') {
	       my ($nb_reads_unspliced_intron) = (0,0);
	       foreach my $sam (@SAM) { 
	          foreach my $read ($sam->get_features_by_location(-seq_id => $seq_id, -start  => $intron{START}, -end    => $intron{END})) {
		     next if($read->end < $intron{START} or $read->start > $intron{END});
		     
		     if($read->cigar_str =~/^\d+M$/) {
		        $nb_reads_unspliced_intron++;
		     }
		  }
	       }
	       $intron{SPLICING_RATE} = ($intron{COVERAGE}/($intron{COVERAGE}+$nb_reads_unspliced_intron));
	       
	    }
	    #push @INTRONS,\%intron if($intron{COVERAGE} >= $self->{MIN_INTRON_COVERAGE} and ($intron{SPLICING_RATE} eq 'NA' or $intron{SPLICING_RATE}>=$self->{MIN_SPLICING_RATE}));
	    push @INTRONS,\%intron if($intron{COVERAGE} >= $self->{MIN_INTRON_COVERAGE} and ($intron{SPLICING_RATE} eq 'NA' or $intron{SPLICING_RATE}>=$self->{MIN_SPLICING_RATE}));
	    
	    $start_covtig='';
	 } elsif($cov >= 1 and !$start_covtig) {
   	    $start_covtig = $i;
   	 }
      }
   } 
   
   my %RESULTS;
   $RESULTS{$seq_id}->{INTRONS} = \@INTRONS;
  
   foreach my $type (qw(transcription_start_site transcription_end_site)) {
      next if(!$self->{FEATURES}->{$type}->{$seq_id});
      foreach my $strand (keys %{$self->{FEATURES}->{$type}->{$seq_id}}) {
         foreach my $pos (keys %{$self->{FEATURES}->{$type}->{$seq_id}->{$strand}}) {
	    my $in_intron = 0;
	    foreach my $intron (sort {$a->{START}<=>$b->{START}} @INTRONS) {
	       next if($in_intron or $intron->{STRAND} ne $strand);
	       next if($intron->{COVERAGE} < $self->{MIN_INTRON_COVERAGE});
	       if($pos >= $intron->{START} and $pos <= $intron->{END}) {
	          $in_intron = 1;
	       }
	    }  
	    $INSERT_COVERAGE{$strand}->{$pos} = 0 if(!$in_intron);
	 }
      }
   }  
   
  
   my $MAX_EXTEND = 10;
  
   my @TRANSCRIPTS;
   my @COVTIGS;
   foreach my $strand (keys %INSERT_COVERAGE) {
      my $start_covtig;      
      my $sum_coverage=0;  
      for(my $i=1; $i<=($seq_length+1); $i++) {
         my $cov = $INSERT_COVERAGE{$strand}->{$i} ? $INSERT_COVERAGE{$strand}->{$i} : 0; 
	 if($cov < $self->{MIN_COVERAGE} and $start_covtig) {
	    my $end_covtig = $i-1;
     	    my $length = $end_covtig-$start_covtig;
	    if($length > $self->{MIN_LENGTH}) {
	    
	       
	       my $coverage = sprintf('%.2f',$sum_coverage/$length);
	       if($coverage >= $self->{MIN_SCORE}) {
	       

	       
	          my %covtig = (STRAND => $strand, START => $start_covtig , END =>$end_covtig, COVERAGE => $coverage);
	          push @COVTIGS,\%covtig;
	       
	       
	       
	          my @exons_positions = ($start_covtig);
	          foreach my $intron (sort {$a->{START}<=>$b->{START}} @INTRONS) {
	             #next if($intron->{STRAND} ne $strand);
	             next if($intron->{END} < $start_covtig or $intron->{START} > $end_covtig);
		     next if($intron->{COVERAGE} < $self->{MIN_INTRON_COVERAGE});
		     next if($intron->{SPLICING_RATE} ne 'NA' and $intron->{SPLICING_RATE} < $self->{MIN_SPLICING_RATE});
     	             push @exons_positions,$intron->{START}-1,$intron->{END}+1;
	          }
   	          push @exons_positions,$end_covtig;
                  my $nb_exons = (scalar @exons_positions) / 2;
     	          die "Pb with $seq_id:$start_covtig .. $end_covtig : nb exons" if((scalar @exons_positions) % 2 != 0);          
	          my %transcript = (STRAND => $strand, START => $start_covtig , END =>$end_covtig, COVERAGE => $coverage, EXONS => \@exons_positions);
	          push @TRANSCRIPTS,\%transcript;
	       } 
	    }
	    $start_covtig='';
	 } elsif($cov >= $self->{MIN_COVERAGE} and !$start_covtig) {
   	    $start_covtig = $i;
	    $sum_coverage=0;
   	 }
	 $sum_coverage+=$cov;
      }
   } 
   $RESULTS{$seq_id}->{COVTIGS} = \@COVTIGS;
   $RESULTS{$seq_id}->{TRANSCRIPTS} = \@TRANSCRIPTS;
   

   return %RESULTS;
}



1;
