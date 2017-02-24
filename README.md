# aftermerge
aftermerge: A script for merging paired reads after mapping, based on the locations on the reference genome

aftermerge goes over the given sam file, checks if paired end reads are overlapping, and merges these overlapping reads to 1 fragment.
The output is a new sam file (standard format, with only the 11 mandatory fields, extra fields are removed). 

aftermerge gives basic statistics including the initial number of reads, the number of mapped reads, the number of merged reads and the number of overlapping reads which could not be merged.
An error file can be generated with the locations of reads that could not be merged (in bed format). These locations could be checked in the sam/bam file to see the cause of this. Most occuring are indels included in the mapping, which are introduced by a small variation on the reference (mis assembly, or forced mapping). This could be blacklisted by the given bed file.

##Licence

All parts of this tool is licenced under GPLv3.  
A copy of this licence is included under LICENSE.

##Help
aftermerge is really simple to use:
It has 2 mandatory arguments, the input file and output file. It is possible to use stdin and stdout. Progress and statistics are writen to the stderr.
There is 1 optional argument, this is the filename of the error file. If not supplied, this file is not generated.

Example of usage:
```bash

python aftermerge.py mymappedfile.sam outputfile.sam

```

Advised example in a pipeline (mapper and parameters are only an example):
```bash

SAMPLE="mysamplename"

bowtie2 -q -p 4 --seed 1 --end-to-end -x genome -1 $SAMPLE.R1.fastq -2 $SAMPLE.R2.fastq | python aftermerge.py -error $SAMPLE.aftermerge.error - - | elprep /dev/stdin $SAMPLE".bam" --replace-read-group "ID:$SAMPLE LB:$SAMPLE PL:illumina SM:$SAMPLE" --sorting-order coordinate

```
