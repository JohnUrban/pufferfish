pufferfish
==========

HMM-based approach(es) to analyzing genomics data.

Original application:

	- Copy Number segmentation 

	- e.g. to identify DNA *puff* sequences *for FISH* experiments in Sciara coprophila.

Other applications:

	- Chromatin modification ChIP-seq segmentation (e.g. H3K27me3)

	- Transcript factor binding ChIP-seq/Cut-and-Run (e.g. EcR)


Version Info:
============
- pufferfish - 0.1.20200925 (09/25/2020)
	- Many new features and updates including new normalization protocols, discrete HMMs, Viterbi Training, new utilities, segmented motif stuff, etc.

- pufferfish - 0.0.0b (06/01/2020)

	- Now its own Git Repo.

- pufferfish - 0.0.0 (02/15/2016)

	- early development was in subdir of sciaraTools.





Usage and Tools:
================


usage: pufferfish [-h] [-v]

                  {mapreads,getcov,findpuffs,dump,puffcn,summits,normalize,generate,filter,help}

                  ...

 PufferFish - HMM-based approach(es) to finding and analyzing developmentally regulated amplicons (genomic sites that are programmed to increase in copy number over time).

optional arguments:

  -h, --help            show this help message and exit

  -v, --version         Installed pufferfish version.

[sub-commands]:


  {mapreads,getcov,findpuffs,dump,puffcn,summits,normalize,generate,filter,help}

    mapreads            Depends on Bowtie2 and SAMtools.

                        Maps fasta/fastq files to genome (can provide bt2 index or fasta reference (which will first be converted to bt2 index).
                        Maps reads to genome and filters out unmapped reads before sorting and indexing.

    getcov               Depends on SAMtools and BEDtools.

    findpuffs           Take in getcov bedGraphs, do stuff.

    dump                 Take in pickled object containing MultiCovBed object where statepath has been found.

                        Output bedGraph of statepath and/or BED file containing coordinates of states.

    puffcn              Given a latest-stage sample (where all or most puffs have grown) and an optional earliest stage sample

                        (for additional Fold-enrichment normalization), define whether a region is best explained by cn=1,2,4,8,16,32,64.

                        Ideally, this can give an idea of where replication forks approximately reach from each firing.

    summits              Find summits...

    normalize           Given a latest-stage sample (where all or most puffs have grown) and an optional earliest stage sample

                        (for additional Fold-enrichment normalization), just return the late-stage sample with normalized values as specified by protocol options below.

    generate            Generate emitted_data and statepath bedGraphs.

    filter              Given a latest-stage sample (where all or most puffs have grown) and an optional earliest stage sample

                            (for additional Fold-enrichment normalization), just return the late-stage sample with normalized values as specified by protocol options below.

    help                 Gives more extensive guidance on using pufferfish.





Requirements
==========

R >= 3.0.0

Python >= 2.7

BEDTools
cython
pybedtools

##MAYBE
#rpy2 >= 2.4.2
#h5py >= 2.0
#pandas >= 0.14.1
#matplotlib >= 1.4.0
#edgeR (only for kmer 'differential expression' analysis)


INSTALL:
=======
download zip

cd pufferfish-master/

python setup.py install


Dependency installation (assuming Mac OS X):
-------------------------------------------
#pip install pandas
#pip install matplotlib



Expected updates:
================
- I wrote an HMM python library compatible with pufferfish as part of another suite of tools.
- - That will soon be incorporated here and R/rpy2 phased out.
- This will be updated to python3 at some point
- If time, I will clean up and refine code as it was written in earlier days.



Please cite pufferfish as:
-------------------------------------------
