package TrUC::Intron;
use strict;
use base 'TrUC::Root';

use TrUC::Config;
use File::Basename;
use FindBin qw($Bin);

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
   my %INTRONS;
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
				SPLICING_RATE=>'NA', NB_SPLICED_READS =>0, NB_UNSPLICED_READS => 0);
	       }
	       
	    }
	    
	    
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
	       $intron{NB_SPLICED_READS} = $intron{COVERAGE};
	       $intron{NB_UNSPLICED_READS} =$nb_reads_unspliced_intron;
	       $intron{SPLICING_RATE} = ($intron{COVERAGE}/($intron{COVERAGE}+$nb_reads_unspliced_intron));
	   
	    push @{$INTRONS{$seq_id}},\%intron if($intron{COVERAGE} >= $self->{MIN_INTRON_COVERAGE} and ($intron{SPLICING_RATE} eq 'NA' or $intron{SPLICING_RATE}>=$self->{MIN_SPLICING_RATE}));
	    
	    
	    $start_covtig='';
	 } elsif($cov >= 1 and !$start_covtig) {
   	    $start_covtig = $i;
   	 }
      }
   } 





   return %INTRONS;
   

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
      foreach my $feat (sort {$a->{START} <=> $b->{START}} @{$data->{$seq_id}}) {
         my %saved_feat = (START=>$feat->{START}, END=>$feat->{END}, STRAND=>$feat->{STRAND}, COVERAGE=>$feat->{COVERAGE}, 
	 			NB_SPLICED_READS=> $feat->{NB_SPLICED_READS}, NB_UNSPLICED_READS=> $feat->{NB_UNSPLICED_READS}, SPLICING_RATE => $feat->{SPLICING_RATE});
         
         push @{$self->{RESULTS}->{$seq_id}},\%saved_feat;
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
   
   my $out_file=$self->{OUT_DIR}.'/introns.gff3';
   $self->stderr("Write $out_file\n");
   
   my $feat_type ='intron';
   open(OUT,">$out_file") or die "Can not open $out_file";
      foreach my $seq_id (keys %{$self->{RESULTS}}) {
         foreach my $feat (@{$self->{RESULTS}->{$seq_id}}) {
            my ($start,$end,$strand,$coverage) = ($feat->{START}, $feat->{END}, $feat->{STRAND}, $feat->{COVERAGE});
	    my $id = "TRUC:".uc($feat_type).":$seq_id:$start..$end";
   	    my @attributes = ("ID=$id");
	    push @attributes,"nb_spliced_reads=".$feat->{NB_SPLICED_READS};
	    push @attributes,"nb_unspliced_reads=".$feat->{NB_UNSPLICED_READS};
	    push @attributes,"splicing_rate=".$feat->{SPLICING_RATE};
     	    print OUT join("\t",($seq_id,'TRUC',$feat_type,$start,$end,$coverage,$strand,'.',join(";",@attributes))),"\n";
	    #push @{$FEATURES{$seq_id}->{$coverage}},$feat if($feat_type eq 'TRANSCRIPTS' and $skip_antisens);
	    
         }
      }
   close OUT;

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
