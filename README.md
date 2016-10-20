##Wang, Stapleton, Klesmith, Hewlett, Whitehead, and Maynard
##Fine epitope mapping of two antibodies neutralizing the Bordetella adenylate cyclase toxin 

###Summary of the experiment:
Antigen genes were subjected to error-prone PCR.
The resulting libraries were displayed on the surface of yeast.
Yeast libraries were sorted by flow cytometry on the basis of
fluorescently conjugated antibody binding.
Plasmid DNA was isolated from sorted yeast.
Isolated antigen genes were fragmented and sequenced with
125-nt paired-end Illumina reads.

###Summary of the analysis goals:
We would like to count the abundance of each mutation in the
library before and after selection.
These counts will be used to calculate enrichment scores
and, thereby, the effect on antibody-antigen binding of
each mutation in the library.

Fowler et al. have published a software package called
Enrich (http://depts.washington.edu/sfields/software/enrich/)
that takes sequencing data from selected and unselected
populations and calculates enrichment statistics.
However, Enrich requires each read to begin and end at
the same location in the gene. This is not the case in
this randomly sheared sequencing library.

The code provided here takes sequencing data from the
randomly sheared library and creates a new .fastq
file that is compatible with Enrich.


###To run:
We have provided a makefile that automatically runs the
pipeline. The paths will need to be modified to make it
work on your computer.
```
make clean
make
```


###Algorithm:
1. Reads are trimmed with Trimmomatic to remove Illumina
adapter sequences and low-quality bases.
2. Read pairs that overlap are merged with FLASH.
3. Merged sequencing reads are aligned to a wild-type
reference with the EMBOSS water local alignment program.
4. The read could be sense or antisense to the wild-type.
Because water doesn't automatically align the reverse
complement, generate and align the reverse complement.
5. Parse the water output files to find the identity % of
each alignment and the wild-type positions where each
alignment begins and ends, and print these to new files.
6. Compare these two files (from the forward and rev-comp
orientations). Write the line with the higher identity %
to a new file.
7. Every line in the new file now corresponds to a
seqencing read and contains the start and stop indexes
of its best alignment to the wild-type sequence.
These indexes are used to build a full-length
sequence "read" consisting of wild-type sequence from
the start of the gene to the start of the read, then the
read, then wild-type sequence from the end of the read to
the end of the wild-type sequence. A corresponding
quality score (all "A") is also generated to create a
.fastq file that Enrich will accept.
8. Steps 3-7 above are repeated for reads that were
orphaned from their mates by Trimmomatic and read pairs
that were not merged by FLASH. In the latter case, the
read is wild-type from the start of the wild-type
to the start of the first read, then the first read,
then wild-type until the start of the second read, then
the second read, then wild-type until the end of the
wild-type sequence.
9. Run Enrich using the generated .fastq files.



###The script calls the following software:
- Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)
- FLASH (http://ccb.jhu.edu/software/FLASH/)
- EMBOSS water aligner (http://emboss.sourceforge.net/)
- Enrich (http://depts.washington.edu/sfields/software/enrich/)
