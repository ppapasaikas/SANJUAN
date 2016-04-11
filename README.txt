Last update of this file: April 11th, 2016

We are currently re-working SANJAUN to make
integration of new genomes possible.
A first complete version of SANJUAN shall be available
within a few days.

SANJUAN -- *S*plicing *AN*alysis & *JU*nction *AN*notation
version 1.0 beta


0. CONTENTS
===========
1. SHORT INTRODUCTION
2. DEPENDENCIES
3. INSTALLATION
4. ADDING GENOMES / SPECIES
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

It is not designed to determine inclusion levels of 
alternative splicing events under only one condition.

SANJUAN detects in a de-novo manner splicing events
but also relies on annotation data of splicing events.
SANJUAN comes with the capability of adding new 
genomes / species by the user. Therefore, SANJUAN
offers convenient scripts. In addition, the following
genomes / species are available ready-to-use.
1. human      (hg19)
2. human      (hg38)
3. mouse      (mm10)
4. fly        (dm6)
5. zebrafish  (danRer10)
6. worm       (ce11)


As input SANJUAN takes either
1. RNAseq data (FASTQ, FASTQGZ files, only paired-end)
2. mapped reads (BAM files)

For adapter removal, trimming and mapping,
SANJUAN relies on Cutadapt, Trim-Galore and Tophat2.
The user might decide to do these pre-processing
setps and mapping with other programs by herself
and apply SANJUAN on the resulting BAM files,
essentially skipping the build-in pre-processing
of SANJUAN.

SANJUAN is designed to run on "Linux" like platforms.
SANJUAN is currently a beta version. You use it
on your own risk.


2. DEPENDENCIES
===============
SANJUAN is a Perl pipeline and was tested under
Perl v5.10.1 and v5.18.2. It relies on the
following programs.

1. trim_galore            (*)
2. cutadapt               (*)
3. tophat2 (& botwie2)    (*)
4. samtools
5. bedtools
6. overlapSelect
7. awk                    (+)
8. command line tool sort (+)


(*): These marked prorgams need to be installed only
if you will use the pre-processing (adapter removal, 
trimming) and mapping routine of SANJUAN.
(+): These programs / tools are normally already
installed on "Linux" like platforms.


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
4. db: contains annotations data for the 
         annotation of alternative splicing events used 
         internally by SANJUAN.
5. mapping_indexes: contains bowtie2 and tophat2 indexes as 
         well as GTF gene annotations. These indexes and 
         files are used only by the pre-processing and 
         mapping routine of SANJUAN. 

IMPORTANT: contents of the db folder is necessary for the
splicing analysis. Adding genomes / species will populate
this folder with data. See section ADDING GENOMES / SPECIES
for installation details.

IMPORTANT: contents of the indexes folder are used by
the pre-processing routines of SANJUAN (essentially mapping).
If the user wants to use the pre-processing routine of SANJUAN, 
the user needs to create the necessary indexes. See section 
ADDING MAPPING INDEXES for details. These indexes are unnecessary 
if you want to do splicing analysis only. Because these indexes
take much space and are unnecessary for the splicing analysis,
they are not part of SANJUAN.

If necessary make the file sanjuan.pl in the bin directory 
executable with

> chmod 770 sanjuan.pl

Now you can run sanjuan through its full
path

> /full/path/to/SANJUAN/bin/sanjuan.pl

You might define a link to sanjuan.pl in any of
your PATH directories to allow easier execution.

> ln -s /full/path/to/SANJUAN/bin/sanjuan.pl /one/of/my/directories/inPATH/sanjuan

Having set this link, you might call SANJUAN simply by

> sanjuan

4. ADDING GENOMES / SPECIES
===========================
Before running splicing analysis, the user needs to
add annotation data for genomes / species.
We offer these data ready-to-use for the species:
1. human      (hg19)
2. human      (hg38)
3. mouse      (mm10)
4. fly        (dm6)
5. zebrafish  (danRer10)
6. worm       (ce11)

To install them, go through the following steps.
1. go to the main SANJUAN installation folder with
name SANJUAN
2. download annotation data by
wget https://s3.amazonaws.com/PAN/SANJUAN/SANJUAN_db.tar.gz
3. uncompress
tar -xvzf SANJUAN_db.tar.gz
4. delete or keep the file SANJUAN_db.tar.gz 
5. run SANJUAN with option -g to see available species and
corresponding short names

If you want to add more genomes, SANJUAN offers convenient
Perl scripts helping you. Please find more details in file
ADDING_MORE_GENOMES.txt .


5. ADDING MAPPING INDEXES
=========================
We provide more details on the installation
of mapping indexes as part of the manual
for adding / installing new genomes / species.
Please find details in ADDING_MORE_GENOMES.txt .


6. RUNNING SANJUAN
==================

IMPORTANT: running without access to CRG cluster 
resources SANJUAN was developed at the CRG, Barcelona, 
and initially tailored for using local cluster resources
at the CRG, i.e., sending automatically jobs to the 
CRG cluster. For use without access to the cluster 
resources of the CRG, SANJUAN needs to be run with
the option -noqsub to prevent that SANJUAN sends 
automatically jobs to the CRG cluster. Nevertheless,
it is possible to run SANJUAN with option -noqsub 
as a single job on any cluster you have access to.


Once SANJUAN is installed, the user might call

> sanjuan

to obtain a help message with further explanations.

Definition of parameters: 
SANJUAN allows the following two ways.
1. by arguments of the command line call
2. by a text file containing parameter definitions

Both ways of defining parameters are almost identical
with one important exception regarding the way of
defining the input files.

If input files are given as arguments to the comman 
line call, they have to be specified with full path 
and be given in a predefined order. The file names
do not have to follow any specific naming convention.

If input files are given through the SANJUAN parameters
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


IMPORTANT: when doing pre-processing and mapping by yourself
1. the XS attribute field in BAM files is necessary
Skipping SANJUAN's pre-processing and mapping routine, the
user can do adapter removal, trimming and mapping of paired 
RNAseq reads by himself. The user would then apply the 
splicing analysis of SANJUAN to her BAM files. The splicin 
analysis of SANJUAN relies amongst other things on the 
XS attribute field for junction reads in the BAM files. 
This attribute field is present by default when mapping 
with Tophat2. When using the STAR aligner the XS attribute 
for (canonical) junction reads can be generated for all 
types of libraries by setting the "--outSAMstrandField" 
switch to "intronMotif". The STAR aligner manual gives
more details.
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

2. The genome ids must have the prefix "chr".


7. OUTPUT OF SANJUAN
====================
SANJUAN applies a set of contraints on differential 
splicing junctions and reports only those which fulfill
the constraints.
The constraints are grouped according to different
levels of confidence on being differentially spliced.
Confidence levels:
1. VHC (very high confidence)
2.  HC (high confidence)
3.  MC (medium confidence)
4.  LC (low confidence)
5.  NC (no confidence)

When running SANJUAN, you may select on of these 
confidence levels. These confidence levels translate
into the following constraints on splicing junctions.

                                       VHC       HC       MC       LC
min. number of junction reads:           9        7        5        3
min. number of neighbor reads:         0.1     0.05     0.01    0.001
XXX                                  0.005    0.004    0.002    0.001
max. length of junction:            100000   100000   100000   200000
min length of junction:                 50       50       50       50
min. fold change of N_reads XXX       0.15      0.1      0.05   0.005
max p value of Hypergeometric test: 0.0001    0.001      0.01     0.3
XXX                                      1        1         1       1

The confidence level NC impies no constraints at all.

Depending on the user's choice of the confidence level,
differentially spliced junctions are reported in file
Annotated_Diff_Junctions.txt.
Differentially retained introns are reported in file
Annotated_Diff_Introns.txt.

Columns of Annotated_Diff_Junctions.txt are as follows.
Here, we assume that the label for the first group is
"cntr" and for the second group "ko".
INCL_COORDs XXX
Gene_Name(s)
High_Confidence_Junction
Competing_Junction
minJDist
COMPET_TYPE
HCJ_5'ss
HCJ_3'ss
HCJ_Junc
CompJ_5'ss
CompJ_3'ss
CompJ_Junc
HCJ_LR(cntr/ko)
CompJ_LogRatio(cntr/ko)
HCJ_Delta(cntr-ko)
CompJ_Delta(cntr-ko)
HCJ_Pval
CompJ_Pval
HCJ_Eff_ko
HCJ_Eff_cntr
HCJ_N_ko
HCJ_N_cntr
HCJ_PSI_ko
HCJ_PSI_cntr
CompJ_Eff_ko
CompJ_Eff_cntr
CompJ_N_ko
CompJ_N_cntr
CompJ_PSI_ko
CompJ_PSI_cntr

Columns in Annotated_Diff_Introns.txt are as follows.
INCL_COORDs
Gene_Name(s)
High_Confidence_Junction
COMPET_TYPE
HCJ_5'ss
HCJ_3'ss
HCJ_Junc
IRLR
PvalIR
HCJ_Delta
HCJ_Pval
HCJ_N_ko
HCJ_N_cntr
HCJ_PSI_ko
HCJ_PSI_cntr

IMPORTANT: XXX
SANJUAN reports fold changes for retained introns.
We plan to switch to PSI values for retained introns
in the future.

In addition, raw information are written into files
1. Diff_Junctions_<VHC,HC,MC,LC,NC>.txt
2. Diff_Junctions_LC.txt
3. Diff_Junctions_NC.txt
where the first file contains information for differential
splicing junctions according to the user's choice of
confidence level. The other two files always get created
for confidence levels LC and NC.

The pre-processing routine produces a TOPHAT
and a TRIM folder for each group / condition.


8. LICENSE
==========
SANJUAN is free software: you can redistribute it and 
modify under the terms of the 

GNU General Public License version 3 

or (at your option) any later version as published by the
Free Software Foundation.


9. CONTACT / BUG REPORTS
========================
Panagiotis Papasaikas started and developed SANJUAN: panagiotis.papasaikas@crg.eu
Andre Gohr: andre.gohr@crg.eu (Support)
