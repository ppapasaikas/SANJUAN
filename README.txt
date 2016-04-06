April 6th, 2016:
This file will updated very soon with
* explanations on where to get the genome gtfs from
* how to create the necessary index structures for tophat2 and bowtie
* how to use SANJUAN
* how to interprete results from SANJUAN


SANJUAN -- *S*plicing *AN*alysis & *JU*nction *AN*notation
version 1.0

SHORT INTRODUCTION
==================
SANJUAN determines differential inclusion levels of 
alternative splicing events between two conditions,
e.g., knock down versus wild type. 

It is not designed to determine inclusion levels
of alternative splicing events under only one 
condition.

SANJUAN is applicable to data from the 
following species.
* human (hg19)
* mouse (mm10)
* zebrafish (dr10)

As input SANJUAN takes either
* RNAseq data (only paired-end RNAseq) for each 
of the two conditions, or
* already mapped reads in form of a BAM file
for each of the two conditions  

The RNAseq data, after trimming and mapping,
is turned into two BAM files, one for each
condition. Both steps, trimming and mapping,
can be done by the user with other programs.
In these cases, the user might feed SANJUAN 
with already trimmed RNAseq data, or already 
mapped data in form of BAM files.

XXX name convetions and other constraints that
must be met when user does mapping by her-/himself.


Dependencies
============
SANJUAN is a Perl pipeline and was tested under
Perl v5.10.1 and v5.18.2.



Pre-processing of SANJUAN, i.e., trimming and mapping
of reads, depends on the following programs. 
* trim_galore
* cutadapt
* tophat2
These programs need to be installed and made executable.

Splicing analysis of SANJUAN does not depend on any
external programs. That means, if the user wants to
use the splicing analysis of SANJUAN only and 
omit pre-processing of RNAseq data, the user does not
have to install any further programs.

samtools, bedtools, overlapSelect


INSTALLATION
============
Checkout SANJUAN from its repository on GitHub.com

> git clone https://github.com/ppapasaikas/SANJUAN.git 

or download and un-compress the zip archive 
into any directory on your computer.

The SANJUAN (SANJUAN-master) directory should contain 
the following sub-directories:
* bin
* lib
* perllib
* db
* indexes

Directory bin contains the main script sanjuan.pl.
Directory lib contains helper scripts used by sanjuan.pl.
Directory perllib contains the two Perl modules Text and
Statistics which are used by SANJUAN.
Directory db contains AS annotations used 
by SANJUAN.

IMPORTANT: 
Pre-processing includes mapping using tophat2.
Tophat2, Bowtie2, and GTF files need to be made
available in directory indexes following the following
structure of subfolders.
/bowtie_dr10
/bowtie_hg19
/bowtie_mm10
/tophat_dr10
/tophat_hg19
/tophat_mm10

XXX description where to download the files and how to create them


If necessary make the file sanjuan.pl in the bin directory 
executable with

> chmod 700 sanjuan.pl

Now you can call sanjuan through its full
path

> /full/path/to/SANJUAN/bin/sanjuan.pl

You might define a link to sanjuan.pl in any of
your PATH directories to allow easier execution.

> ln -s /full/path/to/SANJUAN/bin/sanjuan.pl /one/of/my/directories/inPATH/sanjuan

Having set this link, you might call SANJUAN in any place by

> sanjuan


USE OUTSIDE OF CRG, BARCELONA
============
SANJUAN was developed at the CRG, Barcelona, and 
intially tailored for using local cluster resources 
at the CRG. For use without access to the cluster 
resources of the CRG, SANJUAN can be run with
the option -noqsub.


USAGE
=====
Once SANJUAN is installed, the user might call

sanjuan

to obtain a help message with further explanations.

SANJUAN allows the following two ways for definition
of parameters.
* by arguments of the command line call
* by a text files containing parameters

Both ways of defining parameters are almost identical
with one important exception regarding the way of
defining the input files.

If input files are given as arguments to the command 
line call, they have to be specified with full path 
and be given in a predefined order. The file names
do not have to follow any specific naming convention.

If input files are given through the SANJUAN parameters
file, they are not specified individually but given
by a link to a directory which contains all input files.
These files are then assigned to the two groups / conditions 
automatically by SANJUAN and, if they are RNAseq FASTQ files, 
they are ordered with respect to pairing of the 
paired-RNAseq data. For this, it is necessary that the files
follow a specific naming convention.

Details on all parameters, naming conventions, and the 
two different ways to bypass parameters to SANJUAN, 
are given within the help message of SANJUAN.


CONTACT
=======
Panagiotis Papasaikas: panagiotis.papasaikas@crg.eu
Andre Gohr:            andre.gohr@crg.eu
