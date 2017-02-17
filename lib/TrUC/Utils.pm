package TrUC::Utils;
use strict;


use Bio::GFF3::LowLevel qw/ gff3_parse_feature /;
use Statistics::R;
use File::Basename;
use Bio::SeqIO;
use List::MoreUtils qw(uniq); 


#############################################
# PUBLIC FUNCTIONS
#############################################


  
=head2 read_gff_file_by_id

 Title   : Read gff file 
 Usage   : Read gff file and store the objects by seq_ids
 Function: read_gff_file($gff_file)
 Returns : Hash of arrays 
 Args    : GFF file

=cut

sub read_gff_file_by_id {
   my ($self,$gff_file) = @_;
   
   my %features;
   my $id;
   open(FILE,$gff_file) or die "No gff file : $gff_file";
   while(<FILE>) {
      chomp;
      next if /^#/;
      my $feat = gff3_parse_feature($_);
      if(!$feat->{attributes}->{Parent}) {
        ($id) = @{$feat->{attributes}->{ID}};
      }
      push @{$features{$feat->{seq_id}}->{$id}},$feat;
   }
   close FILE;
   return \%features;

}
  
  
=head2 read_gff_file_by_position

 Title   : Read gff file 
 Usage   : Read gff file and store the objects by seq_ids
 Function: read_gff_file($gff_file)
 Returns : Hash of arrays 
 Args    : GFF file

=cut

sub read_gff_file_by_position {
   my ($self,$gff_file) = @_;
   
   my %features;
   open(FILE,$gff_file) or die "No gff file : $gff_file";
   while(<FILE>) {
      chomp;
      next if /^#/;
      my $feat = gff3_parse_feature($_);
      my $strand = ($feat->{strand} > 0 or $feat->{strand} eq '+') ? '+' : '-';
      $features{$feat->{seq_id}}->{$strand}->{$feat->{start}} = $feat;

   }
   close FILE;
   return \%features;

}
   
  
  
=head2 gff_line

 Title   : write a gff line
 Usage   : gff_line($obj)
 Function: format the object in GFF3 format
 Returns : String
 Args    : Bio::GFF3::LowLevel object

=cut

sub gff_line {
   my ($self,$feat) = @_;
   my ($id) = @{$feat->{attributes}->{ID}};
   my @attributes = ("ID=$id");
   foreach my $key (sort keys %{$feat->{attributes}}) {
      next if($key eq 'ID');
      foreach my $value (uniq @{$feat->{attributes}->{$key}}) {
         push @attributes,"$key=$value";
      }
   }
   my $score = ($feat->{score} eq '') ? '.' : $feat->{score};
   my $strand = ($feat->{strand} eq '') ? '.' : $feat->{strand};
   my $phase = ($feat->{phase} eq '') ? '.' : $feat->{phase};
   return join("\t",($feat->{seq_id},$feat->{source},$feat->{type},$feat->{start},$feat->{end},$score,$strand,$phase,join(";",@attributes)));
}

=head2 revcomp

 Title   : reverse complement of the sequence
 Usage   : revcomp($seq)
 Function: 
 Returns : DNA String
 Args    : DNA String

=cut

sub revcomp {
   my ($self,$nt_seq) = @_;
   $nt_seq = reverse $nt_seq;
   $nt_seq=~tr/ATGC/TACG/;
   return $nt_seq;
}

=head2 overlap

 Title   : 
 Usage   : overlap($start1,$end1, $start2, $end2)
 Function: 
 Returns : Boolean
 Args    : Positions

=cut

sub overlap {
   my ($start1,$end1, $start2, $end2) = @_;
   
   return 0 if($end1 < $start2 or $end2 < $start1);
   
   if($start2 >= $start1 and $start2 <= $end1) {
      return abs($end1-$start2);
   } elsif($start1 >= $start2 and $start1 <= $end2) {
      return abs($end2-$start1);
   } else { die "$start1,$end1, $start2, $end2"; }
   
}

1;
