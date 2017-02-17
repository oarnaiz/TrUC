package TrUC::Root;
use strict;
use FileHandle;
use File::Basename;

=head1 NAME

 TrUC Root module

=head1 AUTHORS

2017, I2BC - CNRS , GNU GPL3

=cut


# COMMON PARAMETERS
my %ROOT_PARAMETERS = (
			GENOME => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"Reference genome"
				},
			OUT_DIR => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"Out directory"
				},
			BAM => {
				MANDATORY=>1, DEFAULT=>[], TYPE=>'MULTIPLE', RANK => 1,
				DESCRIPTION=>"BAM mapping file"
				},								

			MIN_COVERAGE => {
				MANDATORY=>0, DEFAULT=>'10', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"Minimum coverage"
				},
			MIN_SCORE => {
				MANDATORY=>0, DEFAULT=>'10', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"Minimum score"
				},
			MAX_LENGTH_BTW_PAIRS => {
				MANDATORY=>0, DEFAULT=>'1000', TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Maximum length btw pairs"
				},
			NB_REPLICATES => {
				MANDATORY=>0, DEFAULT=>'1', TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Minimum number of replicate to valid the prediction"
				},
			SEQ_ID => {
				MANDATORY=>0, DEFAULT=>'', TYPE=>'VALUE', RANK => 4,
				DESCRIPTION=>"Calculate only on this seq_id"
				},
			THREADS => {
				MANDATORY=>0, DEFAULT=>'6', TYPE=>'VALUE', RANK => 0,
				DESCRIPTION=>"Number of threads"
				},
			VERBOSE => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 0,
				DESCRIPTION=>"Verbose"
				},			
			QUIET => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 0,
				DESCRIPTION=>"Quiet mode"
				},
		);
		

my $STDLOG_FILE = "stdout.log";
my $STDERR_FILE = "stderr.log";

=head2 new

 Title   : new
 Usage   : my $factory = PARTIE::Root->new({....})
 Function: Create a factory object
 Returns : the factory object

=cut

sub new {
  my ($class,$ref) = @_;
  
  # Check the class
  $class = ref($class) || $class;
  
  # Link between the object and the class
  my $self = bless {},$class;
  
  my %root_parameters = %{$self->get_root_parameters};
  foreach my $param (keys %root_parameters) {
     $self->{$param} = ((defined $ref->{"-".lc($param)} and $ref->{"-".lc($param)} ne '') and (ref($ref->{"-".lc($param)}) ne 'ARRAY' or scalar @{$ref->{"-".lc($param)}} !=0)) ? $ref->{"-".lc($param)} : $root_parameters{$param}->{DEFAULT}; 
	die "$param should be TRUE or FALSE not : $self->{$param}\n" if($root_parameters{$param}->{TYPE} eq 'BOOLEAN' and $self->{$param} ne 'TRUE' and $self->{$param} ne 'FALSE');
     
  }
  $self->{VERSION} = $ref->{"-version"};
  
  
  return $self;
}

=head2 finish

 Usage   : $factory->finish
 Function: Finish all the procedure and/or comput all results
 Returns : Nothing
 Args    : Nothing

=cut

sub finish {
  my ($self) = @_;  
  
  $self->{LOG}->close;

}


###################################
# PRIVATE FUNCTIONS
###################################

sub _init {
  my ($self) = @_;  
  
  
  $self->stderr("Initiation\n");
  # create out directory
  my $path = $self->{OUT_DIR}."/";
  
  $self->stderr("Create directory $path\n");
  system("mkdir -p $path");
  $self->{PATH}=$path;

  $self->{LOG} = FileHandle->new(">$self->{PATH}/".$self->get_mode.".$STDLOG_FILE"); 
  $self->stdlog("# COMMAND LINE (TrUC version ".$self->{VERSION}.") :\n".join(" \\\n",$self->_command_line)."\n\n");

  
}

# Display command line
sub _command_line {
   my ($self) = @_;

   my @command = "$0 ".$self->get_mode;

   my %all_parameters=(%{$self->get_root_parameters}, %{$self->get_parameters});
   
   foreach my $param (sort keys %all_parameters) {
      my $type = $all_parameters{$param}->{TYPE} ;   
      if($type eq 'VALUE') {
        my $param_value = $self->{$param};
	$param_value = "'$param_value'" if($param_value=~/ /);
	next if($param_value eq '');
        push @command," -".lc($param)." ".$param_value; 
      } elsif($type eq 'MULTIPLE') {
         foreach my $param_value (@{$self->{$param}}) {
	    $param_value = "'$param_value'" if($param_value=~/ /);
            push @command, " -".lc($param)." ".$param_value;	
	 }
     } elsif($type eq 'BOOLEAN') {
        push @command," -".lc($param) if($self->{$param} eq 'TRUE');    
     }
   }   
   return @command;      
}



=head2 _check_mandatory_parameters

 Title   : "Private" function _check_mandatory_parameters
 Usage   : $factory->_check_mandatory_parameters
 Function: Check the mandatory parameters or deduce them with -auto
 Returns : 1 (success) or 0 (error)

=cut

sub _check_mandatory_parameters {
  my ($self) = @_;  
  my %all_parameters=(%{$self->get_root_parameters}, %{$self->get_parameters});
  foreach my $param (keys %all_parameters) {
     if($all_parameters{$param}->{MANDATORY} and ($self->{$param} eq '' or (ref($self->{$param}) eq 'ARRAY' and scalar @{$self->{$param}} ==0 ))) {
        print STDERR "\n# ERROR -",lc($param)," IS MANDATORY\n\n";
        return 0;
     }
  }
  
  if($self->{GENOME} !~/\.fa$/i and $self->{GENOME} !~/\.fasta$/i)  {
     print STDERR "Genome file (-genome) should be a FASTA file";
     return 0;
  }
  

    
  return 1;
}


# Print usage
sub _usage {
  my ($self) = @_;
  
  
  my %all_parameters=(%{$self->get_root_parameters}, %{$self->get_parameters});
  print STDERR "USAGE : \n$0 ".$self->get_mode."\n";
  my $last_rank;
  foreach my $param (sort {$all_parameters{$a}->{RANK} <=> $all_parameters{$b}->{RANK} } keys %all_parameters) {
     my $mandatory = '(MANDATORY)' if($all_parameters{$param}->{MANDATORY});
     my $type = $all_parameters{$param}->{TYPE};
     my $desc = $all_parameters{$param}->{DESCRIPTION};
     
     print STDERR "\n" if($last_rank ne $all_parameters{$param}->{RANK});
     if($type eq 'VALUE') {
        my $value = ($self->{$param} ne '') ? $self->{$param} : $all_parameters{$param}->{DEFAULT};
        print STDERR "\t-".lc($param)." [".$value."]\t".$desc." $mandatory\n";
     } elsif($type eq 'MULTIPLE') {
        print STDERR "\t{\n";
	my @values = ($self->{$param}) ? @{$self->{$param}} : @{$all_parameters{$param}->{DEFAULT}};
	print STDERR "\t\t-".lc($param)." []\n" if(scalar @values == 0);
        foreach my $value (@values) {
	   print STDERR "\t\t-".lc($param)." [".$value."]\n";
	}
	print STDERR "\t} (MULTIPLE VALUES) $mandatory\t".$desc."\n";	
     } elsif($type eq 'BOOLEAN') {
        my $value = ($self->{$param}) ? $self->{$param} : $all_parameters{$param}->{DEFAULT};
	print STDERR "\t(-".lc($param)." $mandatory (default $all_parameters{$param}->{DEFAULT}))\t$desc\n";   
     }
     else { die $param,' ',$type; }
     $last_rank = $all_parameters{$param}->{RANK}; 
  }   
  print STDERR "\n";
 
  exit 0;
}


sub _is_unique {
   my ($self,$first_mate,$second_mate) = @_;
   
   die "Need paired-end mapping : ".$first_mate->qname." AND ".$second_mate->qname."\n" if($first_mate->qname ne $second_mate->qname);
   
   my $first_mate_xt = $first_mate->get_tag_values('XT') ? $first_mate->get_tag_values('XT') : 'U';
   my $second_mate_xt = $second_mate->get_tag_values('XT')? $second_mate->get_tag_values('XT') : 'U';
   
   return 0 if($first_mate_xt eq 'R' and $second_mate_xt eq 'R');
   
   return 1 if($self->{NOT_STRANDED});
   
   my $xs1 = $first_mate->get_tag_values('XS');
   my $xs2 = $second_mate->get_tag_values('XS');
   
   if($xs1 and $xs1 eq $xs2 and ($xs1 eq '+' or $xs1 eq '-')) {
      return 1;
   } else {
      $self->stdlog("Pb with orientation $xs1 VS $xs2 for ".$first_mate->qname."\n");
   }
   return 0;
   
   
}

###################################
# PROTOTYPES
###################################


=head2 calculate

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut


sub calculate {
  my ($self,$seq) = @_;  
  die "Not managed for $self";
}


=head2 save

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub save {
  my ($self,$data) = @_;  
  die "Save not manage for $self";
  return $self;
}

=head2 write_results

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub write_results {
  my ($self) = @_;  
  die "write_results not manage for $self";
  return $self;
}


###################################
# GENERIC/UTILS FUNCTIONS
###################################

sub get_root_parameters {
  my ($self) = @_;  
  return \%ROOT_PARAMETERS;
}


sub get_mode {
  my ($self) = @_;  
  my $module_name = ref($self);
  if( $module_name =~/::(\S+)$/) {
     return $1;
  }
  die $module_name;
}


sub verbose {
   my ($self) = @_;
   return ($self->{VERBOSE} eq 'TRUE') ? 1 : 0;
}

sub stdlog {
   my ($self,$message) = @_;
   my $fh = $self->{LOG};
   print $fh "$message";
}

sub stderr {
   my ($self,$message) = @_;
   return if($self->{QUIET} eq 'TRUE');
   my $mode = $self->get_mode;
   $mode=~s/_unoriented//;
   print STDERR "# [".$mode."] : $message" if($self->verbose);
}
1;
