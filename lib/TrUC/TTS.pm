package TrUC::TTS;
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
			NB_MIN_A => {
				MANDATORY=>1, DEFAULT=>5, TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Minimum number of Adenine"
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
  
  return $self->SUPER::_check_mandatory_parameters;
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
   
   my %POSITIONS;
   my %COVERAGE;
   
   foreach my $bam_file (@{$self->{BAM}}) { 
      my $fname = basename($bam_file);
      $self->stderr("# Read bam file $fname for $seq_id ... \n");
      my $sam =  Bio::DB::Sam->new(-bam  =>$bam_file, -fasta=> $self->{GENOME} );
      foreach my $pair  ($sam->features(-type   => 'read_pair', -seq_id => $seq_id)) {
       
         my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
	 next if(!$first_mate or !$second_mate);
	 next if(!$self->_is_unique($first_mate,$second_mate)) ;
	 
	 my @pos = sort {$a<=>$b} ($first_mate->start,$first_mate->end,$second_mate->start,$second_mate->end);
	 my ($insert_start,$insert_end) = ($pos[0],$pos[$#pos]);
	 next if(abs($insert_start-$insert_end) > $self->{MAX_LENGTH_BTW_PAIRS});
	 
	 my $strand = $first_mate->get_tag_values('XS');
	 
	 foreach my $read ($first_mate,$second_mate) {
	    my $cigar_str = $read->cigar_str;
	    my $read_seq = $read->query->dna;;
	    my ($tts,$nb_A) = (0,0);
	    if($cigar_str=~/^(\d+)M(\d+)S$/) {
	       my ($match,$space) = ($1,$2);
	       my $space_seq = substr($read_seq,$match);
	       if($space_seq=~/^A+$/ or $space_seq=~/^T+$/) {
	    	  $nb_A=$space;
	    	  $tts = $read->end;
	       } 
	    } elsif($cigar_str=~/^(\d+)S(\d+)M$/) {
	       my ($space,$match) = ($1,$2);
	       my $space_seq = substr($read_seq,0,$space);
	       if($space_seq=~/^A+$/ or $space_seq=~/^T+$/) {
	    	  $nb_A=$space;
	    	  $tts = $read->start;
	       } 
	    }
	    if($tts and $nb_A >= $self->{NB_MIN_A}) {
	    
	       $POSITIONS{$strand}->{$tts}->{COVERAGE}++;
	       $POSITIONS{$strand}->{$tts}->{REPLICATES}->{$fname}=1;
	       for(my $i=$insert_start;$i<=$insert_end;$i++) {
                  $COVERAGE{$strand}->{$i}++;
               }
	    }
	 }
	 
      } 			     
   }   
   
   
   
   my %TTS;
   my $MIN_COVERAGE = 1 ;
   foreach my $strand (keys %COVERAGE) {
      my $start_covtig;
      for(my $i=1; $i<=($seq_length+1); $i++) {
         my $cov = $COVERAGE{$strand}->{$i} ? $COVERAGE{$strand}->{$i} : 0; 
	 if($cov < $MIN_COVERAGE and $start_covtig) {
	    my $end_covtig = $i-1;
	    my %tts;
	    my $sum_coverage = 0;
	    my $nb_replicates = 0;
	    for(my $p=$start_covtig;$p<=$end_covtig;$p++) {
	       next if(!$POSITIONS{$strand}->{$p});
	       
	       if(!%tts
		    or $tts{COVERAGE} < $POSITIONS{$strand}->{$p}->{COVERAGE} 
		    or ($tts{COVERAGE} == $POSITIONS{$strand}->{$p}->{COVERAGE} and $strand eq '+')) {
		  %tts = (POSITION => $p, STRAND => $strand, COVERAGE =>$POSITIONS{$strand}->{$p}->{COVERAGE} );
	       
	       }
	       $nb_replicates = scalar keys %{$POSITIONS{$strand}->{$p}->{REPLICATES}} if($nb_replicates < scalar keys %{$POSITIONS{$strand}->{$p}->{REPLICATES}});
	       $sum_coverage+=$POSITIONS{$strand}->{$p}->{COVERAGE} ;
	    }
	    $tts{SUM_COVERAGE}=$sum_coverage;
	    $tts{NB_REPLICATES}=$nb_replicates;
	    push @{$TTS{$seq_id}},\%tts if($sum_coverage >= $self->{MIN_COVERAGE});
	    $start_covtig='';
	 } elsif($cov >= $MIN_COVERAGE and !$start_covtig) {
   	    $start_covtig = $i;
   	 }
	 
      }
   }   
   
   return %TTS;
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
      foreach my $tts (@{$data->{$seq_id}}) {     
	 $self->{RESULTS}->{$seq_id}->{$tts->{POSITION}} = { strand =>$tts->{STRAND} , coverage =>$tts->{SUM_COVERAGE}, nb_replicates => $tts->{NB_REPLICATES}}; 
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
   
   my $out_file=$self->{OUT_DIR}.'/transcription_end_site.gff3';
   $self->stderr("Write $out_file \n");
   open(OUT,">$out_file") or die $out_file;
   foreach my $seq_id (keys %{$self->{RESULTS}}) {
      foreach my $pos (sort {$a<=>$b} keys %{$self->{RESULTS}->{$seq_id}}) {
     	 my ($strand,$coverage,$nb_replicates) = ($self->{RESULTS}->{$seq_id}->{$pos}->{strand},$self->{RESULTS}->{$seq_id}->{$pos}->{coverage},$self->{RESULTS}->{$seq_id}->{$pos}->{nb_replicates});
     	 if($nb_replicates >= $self->{NB_REPLICATES} and $coverage >= $self->{MIN_SCORE}) {  
	    my $id = "TRUC:TTS:$seq_id:$pos";
     	    my @attributes = ("ID=$id","nb_replicates=$nb_replicates");
     	    print OUT join("\t",($seq_id,'TRUC','transcription_end_site',$pos,$pos,$coverage,$strand,'.',join(";",@attributes))),"\n";
         }
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
