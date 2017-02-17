package TrUC::Config;

use strict;


use File::Which qw(which);
use Bio::GFF3::LowLevel qw/ gff3_parse_feature /;
use Bio::SeqIO;
use Bio::DB::Sam;

use TrUC::Utils;
use TrUC::TSS;
use TrUC::TTS;
use TrUC::Transcript;
use TrUC::Transcript_unoriented;
use TrUC::UTR;
use TrUC::Intron;



=head2 get_modes

 Title   : Get all availables modes
 Usage   : get_modes;
 Function: 
 Returns : Hash with all modes
 Args    : Nothing

=cut

sub get_modes {
   my %modes = ( 
		TSS 		=> { N=> 1, DESC=>'Predict Transcription Start Sites (Need CapSeq mapping)'},
		TTS 		=> { N=> 2, DESC=>'Predict Transcription Terminal Site (polyA mRNA sequencing)'},
		Transcript	=> { N=> 3, DESC=>'Predict Transcription Units'},
		UTR 		=> { N=> 4, DESC=>'Add UTRs to a GFF3 file'},
		Intron 		=> { N=> 5, DESC=>'Predict Introns'},
               );
   return %modes;
} 

1;
