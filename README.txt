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
* mouse (mm10) XXX we also have mm8 and mm9 ?
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

XXX if mapping is done by user, does the user 
has to guaranty that the mapping was done against
a certain version of a genome or do the chromosomes
have to fulfil any name conventions?


GENERAL NOTE
============
SANJUAN may be used only at the CRG as it is tightly
tailored to the cluster infra-structure of the CRG.
Access to the ant-cluster of the CRG is essential,
as SANJUAN will run all computations as cluster jobs. 


Dependencies
============
Preprocessing of SANJUAN, i.e., trimming and mapping
of reads, depends on the following programs. 
* trim_galore
* cutadapt
* tophat2
These programs should be executable under the user's
account on the CRG cluster.

Splicing analysis of SANJUAN does not depend on any
external programs. That means, if the user wants to
use the splicing analysis of SANJUAN only but not
the preprocessing of RNAseq data, the user does not
have to install any further programs.  


INSTALLATION
============
Checkout SANJUAN from its repository on GitHub.com
or download and un-compress the tar.gz-archive 
into any directory on your computer.

The directory should contain the following sub-directories:
* bin
* lib 
* perllib

Directory bin contains the main script sanjuan.pl.
Directory lib contains helper scripts used by sanjuan.pl.
Directory perllib contains the two Perl modules Text and
Statistics which are used by SANJUAN.

Make the main script sanjuan.pl in the bin directory 
executable with

> chmod 770 sanjuan.pl

Now you can call sanjuan through its full
path.

> /full/path/to/SANJUAN/bin/sanjuan.pl

You might define a link to sanjuan.pl in any of
your PATH directories to allow execution of 
sanjuan anywhere.

> ln -s /one/of/my/directories/inPATH/sanjuan /full/path/to/SANJUAN/bin/sanjuan.pl 
   
Having set this link, you might call SANJUAN in any place by

> sanjuan 


IMPORTANT: 
SANJUAN relies on predefined 
exon-exon-junctions for each species. These 
have to be downloaded from XXX.
It is recommended to download and un-compress 
the tar.gz-archive into a sub-directory called db 
of the SANJUAN installation.
The SANJUAN installation then finally contains
the following sub-directories.
* bin
* lib
* perllib
* db
and db should contain the two directories 
genomes and annotation_files.
SANJUAN will then use these sub-directories.

If you decide to put the sub-directory db 
to another place, you will have to tell SANJUAN
the place where you put it each time you use SANJUAN.  
See the help message of SANJUAN for more details
on how to specify the location of the db sub-directory.


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