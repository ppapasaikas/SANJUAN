SANJUAN -- *S*plicing *AN*alysis & *JU*nction *AN*notation
version 1.0 beta


0. CONTENTS
===========
1. SHORT INTRODUCTION
2. DEPENDENCIES
3. INSTALLATION
4. ADDING GENOMES
5. ADDING MAPPING INDEXES
6. RUNNING SANJUAN
7. OUTPUT OF SANJUAN
8. LICENSE
9. CONTACT / BUG REPORTS


1. SHORT INTRODUCTION
=====================
SANJUAN determines differential inclusion levels of
alternative splicing events between two conditions,
e.g., knock down versus wild type.

SANJUAN is not designed to determine inclusion levels of 
alternative splicing events under only one condition.

SANJUAN relies mostly on the RNAseq data to de-novo detect
splicing junctions and alternative splicing events and
only minimally on available annotation data of transcript
structure for refining results.

SANJUAN comes with the capability of making available 
new genomes for splicing analysis to your
local SANJUAN installation with a few easy steps. 

The following genomes are pre-compiled and ready-to-use.
1.  hg19       (human)
2.  hg38       (human)
3.  mm10       (mouse)
4.  dm6        (fly)
5.  danRer10   (zebrafish)
6.  ce11       (worm)
7.  ath10      (arabidopsis)

SANJUAN takes as input either
1. RNAseq data 
- from paired- or single-end experiments
- raw FASTQ, gzipped -.gz- or bzipped -.bz2-
or directly
2. mapped reads 
- when the use wants to do pre-processing of 
  RNAseq data and mapping by himself
- BAM files

When you supply RNAseq read data, SANJUAN relies
on the STAR aligner using it for adapter removal, 
read trimming and mapping.

SANJUAN is designed to run on Linux-like platforms
including Mac OS.

SANJUAN is currently a beta version. You use it
at your own risk.


2. DEPENDENCIES
===============
SANJUAN is a Perl pipeline and was tested under
Perl v5.10.1 and v5.18.2. It depends on the
following programs:

1. samtools v.>=1.1
2. bedtools v.>=2.26
3. awk                     (+)
4. command line tool sort  (+)
5. STAR aligner v.>=2.4.0  (*)

(+): These programs are normally already
installed on Linux-like platforms.
(*): These programs need only be installed if you
intend to use the pre-processing and mapping
functionality of SANJUAN.


3. INSTALLATION
===============
Checkout SANJUAN from its repository on GitHub.com

> git clone https://github.com/ppapasaikas/SANJUAN.git 

or download and un-compress the zip archive.

The SANJUAN (SANJUAN-master) directory should contain 
the following sub-directories:
1. bin
2. lib
3. perllib
4. db
5. indexes

1. bin: contains the main script sanjuan.pl which should
         be used to run SANJUAN. 
2. lib: contains miscellaneous scripts used internally
         by SANJUAN.
3. perllib: contains two Perl modules (Text and
         Statistics) which are used internally by SANJUAN.
4. db: contains annotation data used 
         internally by SANJUAN.
5. mapping_indexes: contains STAR indexes compiled by 
         the user. They are used only by the pre-processing
         and mapping routine of SANJUAN. Can be left empty
         if you don't use this functionality of SANJUAN.

IMPORTANT: contents of the db folder is necessary for
the  splicing analysis. Adding genomes will populate
this folder with data. See section ADDING GENOMES for
installation details.

IMPORTANT: contents of the indexes folder are used by
the pre-processing routines of SANJUAN (essentially mapping).
If the user wants to use the pre-processing routine of SANJUAN, 
the user needs to create the necessary indexes. See section 
ADDING MAPPING INDEXES for details. These indexes are unnecessary 
if you want to do splicing analysis only. Because these indexes
are bulky and unnecessary for the splicing analysis, they are 
not part of SANJUAN.

If necessary make the file sanjuan.pl in the bin directory 
executable with

> chmod 755 sanjuan.pl

Now you can run sanjuan through its full
path

> /full/path/to/SANJUAN/bin/sanjuan.pl

We recommend you to define a link to sanjuan.pl in any of
your PATH directories to allow easier execution.

> ln -s /full/path/to/SANJUAN/bin/sanjuan.pl /one/PATH/directory/sanjuan

Having set this link, you can call SANJUAN simply by

> sanjuan

To run the splicing analysis, you have to add genomes.
For details see section ADDING GENOMES.


4. ADDING GENOMES
===========================
Before running splicing analysis, the user needs to
add annotation data for genomes.
We offer these data ready-to-use for the following genomes:
1.  hg19       (human)
2.  hg38       (human)
3.  mm10       (mouse)
4.  dm6        (fly)
5.  danRer10   (zebrafish)
6.  ce11       (worm)
7.  ath10      (arabidopsis)

To see an overview of all available annotation files 
of splicing events, go to:
https://s3.amazonaws.com/PAN/SANJUAN/aws_S3_index.html

Installation steps:
1. Go to the main SANJUAN installation folder with
name SANJUAN.
2. Download the annotation data for genomes you are
interested in
> wget https://s3.amazonaws.com/PAN/SANJUAN/SANJUAN_db_hg19.tar.gz
> wget https://s3.amazonaws.com/PAN/SANJUAN/SANJUAN_db_hg38.tar.gz
> wget https://s3.amazonaws.com/PAN/SANJUAN/SANJUAN_db_mm10.tar.gz
> wget https://s3.amazonaws.com/PAN/SANJUAN/SANJUAN_db_danRer10.tar.gz
> wget https://s3.amazonaws.com/PAN/SANJUAN/SANJUAN_db_ce11.tar.gz
> wget https://s3.amazonaws.com/PAN/SANJUAN/SANJUAN_db_dm6.tar.gz
> wget https://s3.amazonaws.com/PAN/SANJUAN/SANJUAN_db_ath10.tar.gz
3. un-compress all downloaded files
> ls SANJUAN_db_*.tar.gz | xargs -i tar xfzv {}
4. Delete or keep the downloaded archives.
5. Run SANJUAN with option -g to see available genomes and
corresponding short names. When you have decided to download
all annotation data, the output should look like
> sanjuan -g
   available genomes: 
   ID          FOR MAPPING     FOR SPLICING ANALYSIS
   ath10                no                       yes
   ce11                 no                       yes
   danRer10             no                       yes
   dm6                  no                       yes
   hg19                 no                       yes
   hg38                 no                       yes
   mm10                 no                       yes

If you want to add more genomes, SANJUAN offers convenient
Perl scripts helping you. Please find more details in file
ADDING_MORE_GENOMES.txt .


5. ADDING MAPPING INDEXES
=========================
We provide more details on the installation
of mapping indexes as part of the manual
for adding new genomes to SANJUAN. Details
are given in the file ADDING_MORE_GENOMES.txt.


6. RUNNING SANJUAN
==================
IMPORTANT: SANJUAN was developed at the CRG, Barcelona,
and initially tailored to using local cluster resources
at the CRG, i.e., sending automatically jobs to the 
CRG cluster. For use without access to the CRG cluster, 
SANJUAN needs to be run with the option -noqsub to prevent
that SANJUAN sends jobs to the CRG cluster. Nevertheless,
it is possible to run SANJUAN with option -noqsub 
as a single job on any cluster you have access to.

Once SANJUAN is installed, the user might call

> sanjuan

to obtain a help message with further explanations.

The parameters of SANJUAN can be specified by
1. arguments of the SANJUAN command line call
2. a text file containing parameter definitions

Both ways of defining parameters are almost identical
with one important exception regarding the way of
defining the input files

If input files are given as arguments to the command
line call, they have to be specified with full path 
and be given in a predefined order. The file names
do not have to follow any specific naming convention.

If input files are given through the SANJUAN parameter
file, they are not specified individually but given
by a link to a directory which contains all input 
files. These files are then assigned to the two 
groups / conditions automatically by SANJUAN and, 
if they are RNAseq FASTQ files, they get paired 
automatically with respect to the paired-RNAseq data.
Therefor, these files need to follow a specific naming 
convention. Details on all parameters, naming conventions,
and the  two different ways to bypass parameters to 
SANJUAN, are given within the help message of SANJUAN.


IMPORTANT: when doing pre-processing and mapping by yourself,
please note the following important points,

1. the XS attribute field in BAM files is necessary
as SANJUAN relies among st other things on the 
XS attribute field for junction reads in the BAM files. 
When using the STAR aligner the XS attribute for (canonical)
junction reads can be generated for all types of libraries
by setting the "--outSAMstrandField" switch to "intronMotif".
The STAR aligner manual gives more details:
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
This attribute field is present by default when mapping 
with Tophat2. 

2. The chromosome ids might have the prefix "chr" or not.
Chromosome ids in your BAM files need to have this prefix 
if the genome annotation files used internally by SANJUAN 
use this prefix. Call sanjuan -g to see if installed genome 
annotations for the splicing analysis use this "chr" 
prefix or not.

3. The BAM files need to be sorted and have a header section.


7. OUTPUT OF SANJUAN
====================
SANJUAN applies a set of constraints on differential 
splicing junctions and reports only those which fulfill
the constraints.
The constraints are grouped according to different
levels of confidence on being differentially spliced.
Confidence levels:
1. VHC (very high confidence)
2.  HC (high confidence, default)
3.  MC (medium confidence)
4.  LC (low confidence)
5.  NC (no confidence)

When running SANJUAN, you may select on of these 
confidence levels. These confidence levels imply
into the following constraints on splicing junctions.

                                          VHC       HC       MC       LC
min. DeltaPSI among conditions:           0.2     0.15      0.1        0              
min. ln.fold change among conditions:    0.15      0.1     0.05    0.005
min. number of junction reads:              9        7        5        3
min. Neighb.Junction/Junction reads:      0.1     0.05     0.01    0.001
min. Junction/Neighb.Junction reads:    0.005    0.004    0.002    0.001
max. length of junction:               100000   100000   100000   200000
min. length of junction:                   50       50       50       50
max. p value of Hypergeometric test:   0.0001    0.001     0.01      0.3

The confidence level NC implies no constraints at all.

Depending on the user's choice of the confidence level,
differentially spliced junctions are reported in file
Annotated_Diff_Junctions.txt.
Differentially retained introns are reported in file
Annotated_Diff_Introns.txt.

Columns of Annotated_Diff_Junctions.txt are as follows.
Here, we assume that the label for the first group is
"cntr" and for the second group "ko".
1.  INCL_COORDs : Genomic coordinates delimiting the inclusive interval of the Alternative Splicing Event.
2.  Gene_Name(s) : Official gene symbol or gene ID if gene symbol is not available
3.  High_Confidence_Junction: ID of Differentially used junction. Junction ID is of the form chrName_starCoord_endCoord_strand 
4.  Competing_Junction: ID of identified Competing Junction.  
5.  minJDist: Minimum distance between the boundaries of the competing junctions. For alt5'/3'ss or Cassette Exons this is 0.
6.  COMPET_TYPE: Type of competition between the two junctions. See section Competition Type below for details.
7.  HCJ_5'ss: Annotation status of the HC Junction 5'ss. Either known or novel.
8.  HCJ_3'ss: Annotation status of the HC Junction 3'ss. Either known or novel.
9.  HCJ_Junc: Annotation status of the HC Junction. Either known or novel.
10. CompJ_5'ss: Annotation status of the Competin Junction 3'ss. Either known or novel.
11.CompJ_3'ss: Annotation status of the Competing Junction 5'ss. Either known or novel.
12. CompJ_Junc: Annotation status of the Competing Junction. Either known or novel.
13. HCJ_LR(cond2/cond1): Log (natural) ratio of HC Junction usage (e.g ratio in ko/cntr).
14. CompJ_LogRatio(cond2/cond1): Log (natural) ratio of Comp Junction usage (e.g ratio in ko/cntr).
15. HCJ_Delta(cond2-cond1): Delta PSI of HC Junction.
16. CompJ_Delta(cond2-cond1): Delta PSI of Comp Junction.
17. HCJ_Pval: p-value for differential usage of the HC Junction.
18. CompJ_Pval: p-value for differential usage of the Comp. Junction.
19. HCJ_Eff_cond1: HC Junction Efficiency* in condition 1 (e.g cntr). 
20. HCJ_Eff_cond2: HC Junction Efficiency* in condition 2 (e.g ko). 
21. HCJ_N_cond1: Number of reads supporting HC Junction in condition 1 (e.g cntr). 
22. HCJ_N_cond2: Number of reads supporting HC Junction in condition 2 (e.g ko). 
23. HCJ_PSI_cond1: HC Junction PSI in condition 1 (e.g cntr). 
24. HCJ_PSI_cond2: HC Junction PSI in condition 2 (e.g ko). 
25. CompJ_Eff_cond1: Comp Junction Efficiency* in condition 1 (e.g cntr). 
26. CompJ_Eff_cond2: Comp Junction Efficiency* in condition 2 (e.g ko). 
27. CompJ_N_cond1: Number of reads supporting Comp Junction in condition 1 (e.g cntr). 
28. CompJ_N_cond2: Number of reads supporting Comp Junction in condition 2 (e.g ko). 
29. CompJ_PSI_cond1: Comp Junction PSI in condition 1 (e.g cntr). 
30. CompJ_PSI_cond2: Comp Junction PSI in condition 2 (e.g ko). 

*Efficiency quantifies relative junction usage in a "locus" by taking into account neighboring junctions. It is a less intuitive, yet more robust measure, of junction usage.

Columns in Annotated_Diff_Introns.txt are as follows.
1.  INCL_COORDs : Genomic coordinates delimiting the inclusive interval of the Alternative Splicing Event.
2.  Gene_Name(s) : Official gene symbol or gene ID if gene symbol is not available
3.  High_Confidence_Junction: ID of Differentially used junction. Junction ID is of the form chrName_starCoord_endCoord_strand 
4.  COMPET_TYPE: Type of competition between the two junctions. See section Competition Type below for details.
5.  HCJ_5'ss: Annotation status of the HC Junction 5'ss. Either known or novel.
6.  HCJ_3'ss: Annotation status of the HC Junction 3'ss. Either known or novel.
7.  HCJ_Junc: Annotation status of the HC Junction. Either known or novel.
8.  IRLR: Intron retention (natural) log ratio.
9.  PvalIR: p-value for differential intron retention.
10. HCJ_Delta: Differential Intron Retention. In testing status, NOT recommended for general use.
11. HCJ_Pval. p-value for differential intron retention.
12. HCJ_N_ko: Number of reads supporing the junction in condition 1 (e.g ko)
13. HCJ_N_cntr: Number of reads supporting the junction in condition 2 (e.g cntr)
14. HCJ_PSI_ko: A PSI-like index for intron retention in condition 1. In testing status, NOT recommended for general use.
15. HCJ_PSI_cntr: A PSI-like index for intron retention in condition 2. In testing status, NOT recommended for general use.

IMPORTANT:
Currently SANJUAN mainly relies on fold changes for calling differentially retained introns.
A more accurate calculation of intron retention PSI values is being implemented and will be available in near-future versions.

In addition, raw information are written into files
1. Diff_Junctions_<VHC,HC,MC,LC,NC>.txt
2. Diff_Junctions_LC.txt
3. Diff_Junctions_NC.txt
where the first file contains information for differential
splicing junctions according to the user's choice of
confidence level. The other two files always get created
for confidence levels LC and NC.


8. LICENSE
==========
SANJUAN is free software: you can redistribute it and 
modify under the terms of the 

GNU General Public License version 3 

or (at your option) any later version as published by the
Free Software Foundation.


9. CONTACT / BUG REPORTS
========================
Panagiotis Papasaikas developed SANJUAN: panagiotis.papasaikas@fmi.ch
Andre Gohr extended and added functionality: andre.gohr@crg.eu (Support)
Claudia Vivori (testing and support): claudia.vivori@crg.eu 
