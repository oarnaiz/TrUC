package TrUC::UTR;
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
			ANNOTATION => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Gene annotation file (GFF3)"
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


  if($self->{ANNOTATION} !~/\.gff3$/i and $self->{ANNOTATION} !~/\.gff$/i)  {
     print STDERR "Annotation file (-annotation) should be a GFF3 file";
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
    
  $self->stderr("Read ".basename($self->{ANNOTATION})." \n");
  $self->{FEATURES} = TrUC::Utils->read_gff_file_by_id($self->{ANNOTATION});

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
   
   my $MAX_INSERT_SIZE=$self->{MAX_LENGTH_BTW_PAIRS};
   
   my %UTRs;
   foreach my $bam_file (@{$self->{BAM}}) { 
      my $fname = basename($bam_file);
      $self->stderr("Read bam file $fname for $seq_id ... \n");
      my $sam =  Bio::DB::Sam->new(-bam  =>$bam_file, -fasta=> $self->{GENOME} );

      foreach my $gene_id (keys %{$self->{FEATURES}->{$seq_id}}) {
         my ($gene,@feats) = @{$self->{FEATURES}->{$seq_id}->{$gene_id}};
	 my ($cds_start,$cds_end) = _get_cds_position(@feats);
	 next if(!$cds_start  or !$cds_end);
	 my ($strand) = ($gene->{strand});
         foreach my $pos ($cds_start,$cds_end) {
	    foreach my $pair  ($sam->features(-type   => 'read_pair', -seq_id => $seq_id ,
	 						-start=> ($pos < $MAX_INSERT_SIZE) ? 1 : ($pos - $MAX_INSERT_SIZE),
							-end=>  ($pos + $MAX_INSERT_SIZE) )) {
       
               my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
	       next if(!$first_mate or !$second_mate);	 
	       my @pos = sort {$a<=>$b} ($first_mate->start,$first_mate->end,$second_mate->start,$second_mate->end);
	       my ($insert_start,$insert_end) = ($pos[0],$pos[$#pos]);
	       next if(abs($insert_start-$insert_end) > $MAX_INSERT_SIZE);
	    
	       #next if(!$self->_is_unique($first_mate,$second_mate)) ;
	    
	       if( ($cds_start <= $first_mate->start and $cds_end >= $first_mate->end) or ($cds_start <= $second_mate->start and $cds_end >= $second_mate->end)) {
	          $UTRs{$seq_id}->{$gene_id}->{MAX_START} = $insert_start if(!$UTRs{$seq_id}->{$gene_id}->{MAX_START} or $UTRs{$seq_id}->{$gene_id}->{MAX_START} > $insert_start);
	          $UTRs{$seq_id}->{$gene_id}->{MAX_END} = $insert_end if(!$UTRs{$seq_id}->{$gene_id}->{MAX_END} or $UTRs{$seq_id}->{$gene_id}->{MAX_END} < $insert_end);
	       
	       }
	    }
	 }
  
      }
   }

   return %UTRs;

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
      foreach my $gene_id (keys %{$data->{$seq_id}}) {
         foreach my $key (keys %{$data->{$seq_id}->{$gene_id}}) {
            $self->{RESULTS}->{$seq_id}->{$gene_id}->{$key} = $data->{$seq_id}->{$gene_id}->{$key};
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
   
   my $annot_file = $self->{ANNOTATION};
   my $out_file=$self->{OUT_DIR}.'/'.basename($annot_file).'_with_UTRs.gff3';
   $self->stderr("Write $out_file \n"); 
   open(OUT,">$out_file") or die $out_file; 
   $self->stderr("Add UTRs to $annot_file \n");
   open(GFF,$annot_file) or die $annot_file;
   my @feats;
   my $id;
   while(<GFF>) {
     chomp;
     print OUT $_,"\n" if(/^#/);
     next if /^#/;
     my $feat = gff3_parse_feature($_);
     if(!$feat->{attributes}->{Parent}) {
        print OUT join("\n",_write_gff_features($self->{RESULTS}->{$feats[0]->{seq_id}}->{$id},@feats)),"\n" if(@feats);
	($id) = @{$feat->{attributes}->{ID}};
        @feats=();
     } 
     push @feats,$feat;
     
     push @{$self->{FEATURES}->{$feat->{seq_id}}->{$id}},$feat;
     

   }	
   print OUT join("\n",_write_gff_features($self->{RESULTS}->{$feats[0]->{seq_id}}->{$id},@feats)),"\n" if(@feats);
   close GFF;
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
sub _get_cds_position {
   my @feats = @_;
   my ($start_cds,$end_cds);
   foreach my $feat (@feats) {
      next if($feat->{type} ne 'CDS');
      $start_cds = $feat->{start} if(!$start_cds or $start_cds > $feat->{start});
      $end_cds = $feat->{end} if(!$end_cds or $end_cds < $feat->{end});
      
   }
   return ($start_cds,$end_cds);
}

sub _write_gff_features {
   my ($utr,@feats) = @_;
   my $gene = $feats[0];
   
   my $MAX_UTR_SIZE = 200;
   my ($gene_start,$gene_end,$strand) = ($gene->{start},$gene->{end},$gene->{strand});
   my $new_gene_start = ($utr->{MAX_START} and $utr->{MAX_START} < $gene_start and (abs($utr->{MAX_START} - $gene_start) < $MAX_UTR_SIZE) ) ? $utr->{MAX_START} : $gene_start;
   my $new_gene_end = ($utr->{MAX_END} and $utr->{MAX_END} > $gene_end and (abs($utr->{MAX_END} - $gene_end) < $MAX_UTR_SIZE)) ? $utr->{MAX_END} : $gene_end;
   
#   print STDERR "\n";
#   foreach my $feat (@feats) {
#      print STDERR TrUC::Utils->gff_line($feat),"\n";
#   
#   }
#   
#   print STDERR "\n";
#   
   
   my @gff_lines;
   
   my $print_utr = 0;
   my $mrna_id;
   foreach my $feat (@feats) {
      $feat->{start} = $new_gene_start if($feat->{start} == $gene_start and $feat->{type} ne 'CDS');
      $feat->{end} = $new_gene_end if($feat->{end} == $gene_end and $feat->{type} ne 'CDS');
      if($feat->{type} eq 'CDS' and $new_gene_start != $gene_start and !$print_utr) {
         
	 my $cds = $feat;
         my ($cds_id) = @{$cds->{attributes}->{ID}};
         my ($utr_id) = split /:/,$cds_id;
         $utr_id=~s/C(\d+)$/U$1/;
         $utr_id.=":".($new_gene_start).'..'.($gene_start-1);
         my $utr = { seq_id => $cds->{seq_id}, source => $cds->{source}, 
      			type => ($strand eq '+') ? 'five_prime_UTR' : 'three_prime_UTR' , 
			start => ($new_gene_start) , end => ($gene_start-1), score => '.', strand=>$strand, phase=>'.',
      			attributes => { ID => [$utr_id], Name => [$utr_id], Parent => [$mrna_id]}};
      
      
         push @gff_lines, TrUC::Utils->gff_line($utr);
      
      
	 
	 $print_utr = 1;
      } elsif($feat->{type} eq 'mRNA' ) {
         ($mrna_id) = @{$feat->{attributes}->{ID}};
      }
      
      push @gff_lines, TrUC::Utils->gff_line($feat);
   }   
   if($new_gene_end != $gene_end) {
      my $cds = $feats[$#feats];
      my ($cds_id) = @{$cds->{attributes}->{ID}};
      my ($utr_id) = split /:/,$cds_id;
      $utr_id=~s/C(\d+)$/U$1/;
      $utr_id.=":".($gene_end+1).'..'.$new_gene_end;
      my $utr = { seq_id => $cds->{seq_id}, source => $cds->{source}, 
      			type => ($strand eq '+') ? 'three_prime_UTR' : 'five_prime_UTR', 
			start => ($gene_end+1) , end => $new_gene_end, score => '.', strand=>$strand, phase=>'.',
      			attributes => { ID => [$utr_id], Name => [$utr_id], Parent => [$mrna_id]}};
      
      
      push @gff_lines, TrUC::Utils->gff_line($utr);
      
      
      
   }
   return @gff_lines;
   

}
1;
