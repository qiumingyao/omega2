# omega2

You've landed at the page for an overlap-layout-consensus (OLC) metagenome assembler - **omega2**. 

### Setup
1. Download the tarball including source code and executables compiled on Linux system at: [https://bitbucket.org/omicsbio/omega2/downloads](https://bitbucket.org/omicsbio/omega2/downloads). The code has been tested on both Linux and MacOS systems, but not under Windows.
2. To recompile, go to Debug folder, type `make`.
3. If compiled successfully, omega2 will be generated, use `./omega2 -h` for help information.


### Features

* Quick summary
**omega2** is a massively improved version of [omega](http://bioinformatics.oxfordjournals.org/content/early/2014/07/06/bioinformatics.btu395.short), its unique capabilities include:

    1. Modularization of contained and duplicate reads removal, and initial graph construction after transitive edge reduction step: Big data set can be processed in chunks so that memory limitation problem is solved.

    2. More effective and accurate graph reduction steps result in shorter running time and much higher assembly quality.

    3. We have taken out the module that uses paired end reads to build scaffolds since it seems to generate undesired misassemblies. We currently recommend using one of the specialized scaffolding tools. And we have plans to bring the module back with improvement. 

    4. New module where longer accurate reads threading (such as error corrected PacBio reads) can help resolve ambiguity in the overlap graph. Using this module requires a SAM file for input, where the alignment of non-redundant Illumina reads to the PacBio reads, **sorted** (which can be achieved using samtools). Currently the threading is done in an agressive version, i.e. as long as there is one PacBio read that supports connecting two edges in the overlap graph together, the action is taken.

### Current Version
* v1.4

### How to assemble metagenomic sequencing reads from raw data

#### Preprocessing of the Illumina data
##### Trimming, filtering, (merging), and eror correction.

Since omega2 works best with long reads without errors, preprocessing plays an important role in deciding the quality of the assembly results. We have tested Brian Bushnell's suite of tools [BBTools](http://sourceforge.net/projects/bbmap/files/) extensively on the Mock Community Illumina data and have obtained good results, better than other tools compared (FLASH, usearch). Suppose the Illumina reads data set is called `$reads`, the steps we recommend are following:

```
#!sh

# Use bbduk.sh to quality and length trim the Illumina reads and remove adapter sequences
# 1. ftm = 5, right trim read length to a multiple of 5
# 2. k = 23, Kmer length used for finding contaminants
# 3. ktrim=r, Trim reads to remove bases matching reference kmers to the right
# 4. mink=11, look for shorter kmers at read tips down to 11 bps
# 5. qhdist=1, hamming distance for query kmers
# 6. tbo, trim adapters based on where paired reads overlap
# 7. tpe, when kmer right-trimming, trim both reads to the minimum length of either
# 8. qtrim=r, trim read right ends to remove bases with low quality
# 9. trimq=10, regions with average quality below 10 will be trimmed.
# 10. minlength=70, reads shorter than 70bps after trimming will be discarded.
# 11. ref=$adapters, adapters shipped with bbnorm tools
# 12. –Xmx8g, use 8G memory
# 13. 1>trim.o 2>&1, redirect stderr to stdout, and save both to file *trim.o*
adapters= bbmap_dir/resources/adapters.fa
phiX_adapters= bbmap_dir/resources/phix174_ill.ref.fa.gz
bbduk.sh in=$reads out=trim.fq.gz ftm=5 k=23 ktrim=r mink=11 qhdist=1 tbo tpe qtrim=r trimq=10 minlength=70 ref=$adapters -Xmx8g 1>trim.o 2>&1
bbduk.sh in=trim.fq.gz out=filter.fq.gz ref=$phiX_adapters hdist=1 k=31 threads=48
```


##### Paired end reads merging for the tight insert libraries and for longer insert libraries (using gap extension)

```
#!sh
bbmerge.sh in=filter.fq.gz out=merged.fq.gz outu=unmerged.fq.gz &> bbmerge.out
bbmerge.sh in=filter.fq.gz out=merged.fq.gz outu=unmerged.fq.gz extend2=20 iterations=5 -Xmx100g &> bbmerge.out
```


##### Error correction with Tadpole.sh 

```
#!bash
# 1. mode=correct, use tadpole for correction
# 2. ecc=t, error correct via kmer counts
# 3. shave, remove dead ends in the kmer graph
# 4. rinse, remove bubbles
# 5. -Xmx120g, use 120G memory
# 6. markbadbases=2, change bases fully covered by less than 2 kmers to N
tadpole.sh in=merged.fq.gz out=merged_EC.fq.gz mode=correct markbadbases=2 ecc=t shave rinse -Xmx120g &> ecc.out
# 1. maxns=0, reads more than 0 Ns after trimming wil be discarded.
reformat.sh in=merged_EC.fq.gz out=merged_EC_reformat.fq.gz maxns=0 -Xmx120g 1>> ecc.out 2>&1
```


#### Contained reads removal

Now that the size of the metagenomic data is continuously increasing, Qiuming developed as part of [his tools](https://bitbucket.org/yaoornl/align_test/overview) a program to remove duplicated and contained reads. His program streams the dataset so we can split the big data set into smaller pieces that fit in the limited memory, and remove contained reads in parallel. His tool also supports multi-threading for each of the smaller jobs.


+ The following command split the big data file into pieces with about 10G bps each. For a piece of this size, the memory usage is about 50GB (5X of the Fasta file size). *splitReads.py* can be downloaded [here](https://bitbucket.org/jjchai/poopy/downloads). Note that *splitReads.py* needs python module *argparse* to run.

```
#!sh
python splitReads.py -i ${big_sequence_file} -o ${output_prefix} -count 10000000000 -bp
```

+ Remove duplicated and contained reads with Qiuming's tool. Please refer to the help message of `align_test` for the explanations of the parameters.

```
#!sh
align_test -i RemoveContainedReads -q ${reads} -s ${reads} -ht omega --ID 0 -t 8 -l ${ovl} -o ${unique_reads}
```

#### Construct overlap graph and remove transitive edges in parallel with Qiuming's tool. 

Again the new program has the ability to process pieces of the data in parallel in the distributed fashion. Sample commands for one graph construction on data set in one piece and in multiple pieces are following: (usually uses memory 7X of the Fasta file size).

```
#!bash
# -l minimum overlap length for the graph construction.
align_test -i ConstructOverlapGraph --TransitiveReduction -q ${unique_reads} -s ${unique_reads} -ht omega -t 8 -l ${ovl} -o ${align_out}

# commands to do one subset of reads against whole data set:
align_test -i ConstructOverlapGraph --TransitiveReduction -q ${unique_reads_piece_i} -s ${unique_reads} -ht omega -t 8 -l ${ovl} -o ${align_out}
```


#### Run omega2 for assembly of the Illumina reads

```
#!sh
omega2 -f ${unique_reads} -e ${align_out} -ovl ${min_overlap_len} -o ${omega_out_base}
```

For all the options of omega2, use `omega2 -h`

Higher coverage, longer reads can use longer overlap length. For our test, with 270bp reads 80bp minimum overlap length gave the best result. Since the graph construction is separated from assembly, the graph constructed with smaller minimum overlap can be reused to assemble with longer minimum overlap, we recommend trying different minimum overlap lengths to find the best assembly result.

### Dependencies

* C++11, gcc4.8+

### Questions?

* [JJ Crosskey](mailto:crosskey.jj@gmail.com)
* [Qiuming Yao](mailto:yao.ornl@gmail.com)
* [Ted Ahn](mailto:ahn.no1@gmail.com)
