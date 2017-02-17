# TrUC
Transcription Units by Coverage

# Description

# Install
An installation example is provided at the end of the user manual, as well as in the "INSTALL" file. 
Once all the dependencies are installed, use the "check" file to make sure everything is ok.
```bash
 ./check
```

# Usage
To run TrUC use the following command
```bash
truc (version 1.00) [MODE] : Transcription Units by Coverage
 TSS : Predict Transcription Start Sites (Need CapSeq mapping)
 TTS : Predict Transcription Terminal Site (polyA mRNA sequencing)
 Transcript : Predict Transcription Units
 UTR : Add UTRs to a GFF3 file
 Intron : Predict Introns

```

# License
This software is distributed under the GNU GPL v3 license. See the "LICENSE" file for details.

# How to cite
If you use this software please cite the following publication

Olivier Arnaiz, Erwin Van Dijk, Mireille Bétermier, Maoussi Lhuillier-Akakpo, Augustin De Vannsay, Sandra Duharcourt, Erika Sallet, Jérôme Gouzy and Linda Sperling 
Improved methods and resources for Paramecium genomics: Transcription Units, Gene Annotation and Gene Expression. 
Submited

# Example
The example directory contains the following files :

example/reference.fa : A refence sequence
example/CapSeq.bam : 5'CapSeq data
example/mRNAseq.bam : mRNA seq data
example/mRNAseq_partial_mapping.bam : mRNA seq data allowing partial mapping
example/annotation.gff3 : GFF3 annotation file

TSS : Predict Transcription Start Sites (Need CapSeq mapping)
```bash
OUT_DIR=TEST
THREADS=6
MIN_COVERAGE=15
MIN_SCORE=15

./truc TSS -v -genome example/reference.fa -out_dir $OUT_DIR -threads $THREADS -min_coverage $MIN_COVERAGE -min_score $MIN_SCORE -bam example/CapSeq.bam


```



TTS : Predict Transcription Terminal Site (polyA mRNA sequencing)
```bash
OUT_DIR=TEST
THREADS=6
MIN_COVERAGE=5
MIN_SCORE=10
NB_MIN_A=5

./truc TTS -v -genome example/reference.fa -out_dir $OUT_DIR -threads $THREADS -min_coverage $MIN_COVERAGE -min_score $MIN_SCORE -nb_min_A $NB_MIN_A -bam example/mRNAseq_partial_mapping.bam


```

Intron : Predict introns
```bash
OUT_DIR=TEST
THREADS=6
MIN_COVERAGE=5
MIN_SCORE=10
MIN_INTRON=15
MAX_INTRON=100
MIN_INTRON_COVERAGE=15
MIN_SPLICING_RATE=0.7

./truc Intron -v -genome example/reference.fa -out_dir $OUT_DIR -threads $THREADS -min_coverage $MIN_COVERAGE -min_score $MIN_SCORE -intron_consensus -min_intron_length $MIN_INTRON -max_intron_length $MAX_INTRON -min_intron_coverage $MIN_INTRON_COVERAGE -bam example/mRNAseq.bam


```


UTR : Add UTRs to a GFF3 file
```bash
OUT_DIR=TEST
THREADS=6
MIN_COVERAGE=5
MIN_SCORE=10
ANNOTATION=example/annotation.gff3

./truc UTR -v -genome example/reference.fa -out_dir $OUT_DIR -threads $THREADS -min_coverage $MIN_COVERAGE -min_score $MIN_SCORE -annotation $ANNOTATION -bam example/mRNAseq.bam


```




Transcript : Predict Transcription Units
```bash
OUT_DIR=TEST
THREADS=6
MIN_COVERAGE=5
MIN_SCORE=10
MIN_LENGTH=300
MIN_INTRON=15
MAX_INTRON=100
MIN_INTRON_COVERAGE=15
MIN_SPLICING_RATE=0.7
TSS=$OUT_DIR/transcription_start_site.gff3
TTS=$OUT_DIR/transcription_end_site.gff3

./truc Transcript -v -genome example/reference.fa -out_dir $OUT_DIR -threads $THREADS -min_coverage $MIN_COVERAGE -min_score $MIN_SCORE -min_length $MIN_LENGTH \
 -intron_consensus -min_intron_length $MIN_INTRON -max_intron_length $MAX_INTRON -min_intron_coverage $MIN_INTRON_COVERAGE \
 -no_overlap -tss $TSS -tts $TTS \
 -bam example/mRNAseq.bam


```





