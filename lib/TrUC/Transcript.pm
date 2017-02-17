package TrUC::Transcript;
use strict;
use base 'TrUC::Root';

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
my %PARAMETERS = (
			INTRON_CONSENSUS => {
				MANDATORY=>0, DEFAULT=>'TRUE', TYPE=>'BOOLEAN', RANK => 3,
				DESCRIPTION=>"Restrict intron detection to the GT..AG intron consensus"
				},
			MIN_INTRON_LENGTH => {
				MANDATORY=>1, DEFAULT=>15, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Minimum intron length"
				},				
			MAX_INTRON_LENGTH => {
				MANDATORY=>1, DEFAULT=>10000, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Maximum intron length"
				},				
			MIN_INTRON_COVERAGE => {
				MANDATORY=>1, DEFAULT=>1, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Minimum intron coverage"
				},
			NO_OVERLAP => {
				MANDATORY=>0, DEFAULT=>'TRUE', TYPE=>'BOOLEAN', RANK => 3,
				DESCRIPTION=>"No overlap detection"
				},
			MIN_SPLICING_RATE => {
				MANDATORY=>0, DEFAULT=>0.5, TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Minimum splicing to detect as intron"
				},		
			MIN_LENGTH => {
				MANDATORY=>1, DEFAULT=>100, TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"Minimum length of the predicted unit"
				},
			TSS => {
				MANDATORY=>0, DEFAULT=>'', TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Transcription Start Site GFF3 file"
				},	
			TTS => {
				MANDATORY=>0, DEFAULT=>'', TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Transcription End Site GFF3 file"
				},
			UNORIENTED => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 4,
				DESCRIPTION=>"Unoriented sequecing data"
				},	
				
		);


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

#############################################
# PRIVATE FUNCTIONS
#############################################

sub get_parameters {
  my ($self) = @_;  
  return \%PARAMETERS;
}

=head2 _check_mandatory_parameters

 Title   : "Private" function _check_mandatory_parameters
 Usage   : $factory->_check_mandatory_parameters
 Function: Check the mandatory parameters or deduce them with -auto
 Returns : 1 (success) or 0 (error)

=cut

sub _check_mandatory_parameters {
  my ($self) = @_; 
  
  return 0 if(!$self->SUPER::_check_mandatory_parameters); 
  
  if($self->{MIN_SPLICING_RATE} and ($self->{MIN_SPLICING_RATE} > 1 or $self->{MIN_SPLICING_RATE} < 0)) {
     print STDERR "Splicing rate should be a float value between 0 and 1";
     return 0;
  }  
  
  return 1;
}

#############################################
# PUBLIC FUNCTIONS
#############################################



=head2 init

 Usage   : $factory->init
 Function: Initiate processes before multi thread calculation
 Returns : Nothing
 Args    : Nothing

=cut

sub init {
  my ($self) = @_;  
  $self->SUPER::_init;
  
  if($self->{TSS}) {
     $self->stderr("Read ".basename($self->{TSS})." \n");
     $self->{FEATURES}->{transcription_start_site} = TrUC::Utils->read_gff_file_by_position($self->{TSS});
  }
  if($self->{TTS}) {
     $self->stderr("Read ".basename($self->{TTS})." \n");
     $self->{FEATURES}->{transcription_end_site} = TrUC::Utils->read_gff_file_by_position($self->{TTS});
  }
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
	 
	 my $strand = $first_mate->get_tag_values('XS');
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
		  if($INTRON_POSITIONS{$strand}->{"$intron_start-$intron_end"}->{INTRON_CONSENSUS} eq '') {
		     my $intron_seq = uc(substr($sequence,$intron_start-1,$intron_end-$intron_start+1));
		     $intron_seq = TrUC::Utils->revcomp($intron_seq) if($strand eq '-');
		     $intron_consensus = ($intron_seq=~/^GT.+AG$/) ? 1 : 0;
		  }
		  next if($self->{INTRON_CONSENSUS} eq 'TRUE' and !$intron_consensus);
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
		   %intron = ( POSITIONS => $intron_position, STRAND => $strand, START => $intron_start , END =>$intron_end,
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
		     
		     next if($read->get_tag_values('XS') ne $strand);
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
	       
	          if($strand eq '+') {
		     my $five_prime_start = $start_covtig;
		     for(my $p=($start_covtig-$MAX_EXTEND);$p<$start_covtig;$p++) {
		        $five_prime_start = $p if($self->{FEATURES}->{transcription_start_site}->{$seq_id}->{$strand}->{$p} and $p < $five_prime_start);
		     }
		     $start_covtig = $five_prime_start;
		     
		     my $three_prime_end = $end_covtig;
		     for(my $p=$end_covtig+1;$p<($end_covtig+$MAX_EXTEND);$p++) {
		        $three_prime_end = $p if($self->{FEATURES}->{transcription_end_site}->{$seq_id}->{$strand}->{$p} and $p > $three_prime_end);
		     }
		     $end_covtig = $three_prime_end;
		     
		     
		     
		  } else {
		     my $five_prime_start = $end_covtig;
		     for(my $p=$end_covtig+1;$p<($end_covtig+$MAX_EXTEND);$p++) {
		        $five_prime_start = $p if($self->{FEATURES}->{transcription_start_site}->{$seq_id}->{$strand}->{$p} and $p > $five_prime_start);
		     }
		     $end_covtig = $five_prime_start;
		     
		     my $three_prime_end = $start_covtig;
		     for(my $p=($start_covtig-$MAX_EXTEND);$p<$start_covtig;$p++) {
		        $three_prime_end = $p if($self->{FEATURES}->{transcription_end_site}->{$seq_id}->{$strand}->{$p} and $p < $three_prime_end);
		     }
		     $start_covtig = $three_prime_end;
		  }
		  
		  
	       
	          my %covtig = (STRAND => $strand, START => $start_covtig , END =>$end_covtig, COVERAGE => $coverage);
	          push @COVTIGS,\%covtig;
	       
	       
	          my @exons_positions = ($start_covtig);
	          foreach my $intron (sort {$a->{START}<=>$b->{START}} @INTRONS) {
	             next if($intron->{STRAND} ne $strand);
		     next if($intron->{COVERAGE} < $self->{MIN_INTRON_COVERAGE});
		     next if($intron->{SPLICING_RATE} ne 'NA' and $intron->{SPLICING_RATE} < $self->{MIN_SPLICING_RATE});
		     if($intron->{END} < $end_covtig and $intron->{START} > $start_covtig) {
     	               push @exons_positions,$intron->{START}-1,$intron->{END}+1;
		     }
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

=head2 save

 Usage   : $factory->save
 Function: 

 Returns : Nothing
 Args    : Nothing

=cut

sub save {
   my ($self,$data) = @_;    
   foreach my $seq_id (keys %$data) {
      next if(!$data->{$seq_id});
      foreach my $feat_type (qw(INTRONS COVTIGS TRANSCRIPTS)) {
         next if(!$data->{$seq_id}->{$feat_type});
	 foreach my $feat (sort {$a->{START} <=> $b->{START}} @{$data->{$seq_id}->{$feat_type}}) {
	    my %saved_feat = (START=>$feat->{START}, END=>$feat->{END}, STRAND=>$feat->{STRAND}, COVERAGE=>$feat->{COVERAGE});
	    
	    if($feat->{EXONS}) {
	       foreach my $p (@{$feat->{EXONS}}) {
	          push @{$saved_feat{EXONS}},$p;
	       }
	    } elsif($feat_type eq 'INTRONS') {
	       $saved_feat{SPLICING_RATE} = $feat->{SPLICING_RATE};
	    }
	    push @{$self->{RESULTS}->{$feat_type}->{$seq_id}},\%saved_feat;
	 }
      }
   }  
}

=head2 write_results

 Usage   : $factory->write_results
 Function: 

 Returns : Nothing
 Args    : Nothing

=cut

sub write_results {
   my ($self) = @_;  
   
   my %GFF_TYPES = (INTRONS => 'intron', COVTIGS => 'covtig' , TRANSCRIPTS=> 'transcript');
   
   my %FEATURES;
   foreach my $feat_type (qw(INTRONS COVTIGS TRANSCRIPTS)) {
      my $out_file=$self->{OUT_DIR}.'/'.lc($feat_type).'.gff3';
      $self->stderr("Write $out_file \n");
      open(OUT,">$out_file") or die "Can not open $out_file";
      foreach my $seq_id (keys %{$self->{RESULTS}->{$feat_type}}) {
         foreach my $feat (@{$self->{RESULTS}->{$feat_type}->{$seq_id}}) {
            my ($start,$end,$strand,$coverage) = ($feat->{START}, $feat->{END}, $feat->{STRAND}, $feat->{COVERAGE});
	    my $id = "TRUC:".uc($GFF_TYPES{$feat_type}).":$seq_id:$start..$end";
   	    my @attributes = ("ID=$id");
	    push @attributes,"splicing_rate=".$feat->{SPLICING_RATE} if($feat->{SPLICING_RATE} and $feat->{SPLICING_RATE} ne 'NA');
     	    print OUT join("\t",($seq_id,'TRUC',$GFF_TYPES{$feat_type},$start,$end,$coverage,$strand,'.',join(";",@attributes))),"\n";
	    if($feat->{EXONS}) {
	       my @exons_positions = @{$feat->{EXONS}};
	       for(my $i=0;$i<@exons_positions;$i+=2) {
   	          my $exon_id = "TRUC:EXON:$seq_id:$exons_positions[$i]..$exons_positions[$i+1]";
   	          my @exon_atts = ("ID=$exon_id","Name=$exon_id","Parent=$id");
   	          print OUT join("\t",($seq_id,'TRUC','exon',$exons_positions[$i],$exons_positions[$i+1],'.',$strand,'.',join(";",@exon_atts))),"\n";
   	 
   	       }
	       
	    }
	    $FEATURES{$seq_id}->{$id} = { START=>$start, END=>$end, COVERAGE=>$coverage, LENGTH=>($end-$start) } if($feat_type eq 'TRANSCRIPTS' and $self->{NO_OVERLAP} eq 'TRUE');
	    
         }
      }
      close OUT;
   }
   
   
   if($self->{NO_OVERLAP} eq 'TRUE') {
      
      my %SKIP_TRANSCRIPTS;
      
      foreach my $seq_id (keys %FEATURES) {
         my @ids = keys %{$FEATURES{$seq_id}};
         foreach my $id1 (@ids) {
            my ($start1,$end1,$cov1,$length1) = ($FEATURES{$seq_id}->{$id1}->{START},$FEATURES{$seq_id}->{$id1}->{END},$FEATURES{$seq_id}->{$id1}->{COVERAGE},$FEATURES{$seq_id}->{$id1}->{LENGTH});
            foreach my $id2 (@ids) {
               next if($id1 eq $id2);
               my ($start2,$end2,$cov2,$length2) = ($FEATURES{$seq_id}->{$id2}->{START},$FEATURES{$seq_id}->{$id2}->{END},$FEATURES{$seq_id}->{$id2}->{COVERAGE},$FEATURES{$seq_id}->{$id2}->{LENGTH});
               my $min_length = ($length1 < $length2) ? $length1 : $length2;
	       next if($cov2 > $cov1);
	          
		  # full overlap
	       if( ($start2 > $start1 and $end2 < $end1)
	          #partial overlapping
	          or (TrUC::Utils->overlap($start1,$end1,$start2,$end2)/$min_length) >= 0.7
	       ) {
		  $SKIP_TRANSCRIPTS{$id2} = 1;
	       } 
	    
	    }
	 }
      }

      my $feat_type = 'TRANSCRIPTS';      
      my $out_no_overlapping_file=$self->{OUT_DIR}.'/'.lc($feat_type).'_no_overlap.gff3';
      my $out_overlapping_file=$self->{OUT_DIR}.'/overlapping.gff3';
      my $in_file=$self->{OUT_DIR}.'/'.lc($feat_type).'.gff3';
      $self->stderr("Write $out_no_overlapping_file \n");
      open(OUT,">$out_no_overlapping_file") or die "Can not open $out_no_overlapping_file";
      open(SKIP,">$out_overlapping_file") or die "Can not open $out_overlapping_file";
      open(IN,"$in_file") or die "Can not open $in_file";
      my $id;
      while(<IN>) {
        chomp;
        next if /^#/;
        my $feat = gff3_parse_feature($_);
	if($_!~/Parent/) {
	   ($id) = @{$feat->{attributes}->{ID}};
	}
	print OUT $_,"\n" if(!$SKIP_TRANSCRIPTS{$id});
	print SKIP $_,"\n" if($SKIP_TRANSCRIPTS{$id});
	
      }     
      close IN;
      close OUT;
      close SKIP;

   }
   


}


=head2 finish

 Usage   : $factory->finish
 Function: Finish all the procedure and/or comput all results
           Details :
 	      - compare the current retention scores to a control dataset if provided (CONTROL)
 	      - calls R source to compute statistical tests and returns significance (T/F)
 	      - write the final GFF file

 Returns : Nothing
 Args    : Nothing

=cut

sub finish {
   my ($self, $results)=@_;
}

#############################################
# PRIVATE FUNCTIONS
#############################################


1;
