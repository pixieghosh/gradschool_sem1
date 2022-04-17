## Homework 4 - CHIP-seq

#### We are going to reanalyze some of the data from [Cheng et al. 2018](https://academic.oup.com/plphys/article/178/2/824/6116647), "INDETERMINATE SPIKELET1 Recruits Histone Deacetylase and a Transcriptional Repression Complex to Regulate Rice Salt Tolerance"

#### Copy Homework4 folder into your scratch directory and change into the new directory


```bash
cp -r /scratch/work/courses/AppliedGenomicsSec3/students/homework_4 $SCRATCH/
```


```bash
cd /scratch/sg5096/homework_4
```

#### Explore folder contents:


```bash
ls -la
```

    total 39240964
    drwxrwsr-x.  4 sg5096 sg5096        4096 Apr  9 22:07  .
    drwx--S---. 10 sg5096 sg5096        4096 Apr  8 12:48  ..
    drwxr-sr-x.  2 sg5096 sg5096        4096 Apr  9 22:07  .ipynb_checkpoints
    -rw-r--r--.  1 sg5096 sg5096   378473961 Apr  9 22:06  GCA_001433935.1_IRGSP-1.0_genomic.fna
    -rw-r--r--.  1 sg5096 sg5096       38458 Apr  9 22:06  Homework_4.ipynb
    -rw-r--r--.  1 sg5096 sg5096  9015057924 Apr  9 22:06  Input_1.fastq
    -rw-r--r--.  1 sg5096 sg5096  8878485252 Apr  9 22:06  Rep_1_1.fastq
    -rw-r--r--.  1 sg5096 sg5096  9662905860 Apr  9 22:06  Rep_2_1.fastq
    -rw-r--r--.  1 sg5096 sg5096 12247293828 Apr  9 22:06  Rep_3_1.fastq
    drwxr-sr-x.  2 sg5096 sg5096        4096 Apr  9 21:14 'Untitled Folder'
    -rw-r--r--.  1 sg5096 sg5096         760 Apr  8 12:48  hw4.ipynb
    -rw-r--r--.  1 sg5096 sg5096         217 Apr  9 22:06  index_fasta.sh
    -rw-r--r--.  1 sg5096 sg5096         109 Apr  9 21:53  merge_bam.sh
    -rw-r--r--.  1 sg5096 sg5096         799 Apr  9 21:58  run_bowtie_alignment.sh


#### Input is the total input DNA
#### Rep_1, Rep_2, and Rep_3 are three replicates of CHIP-seq on the transcription factor IDS1

#### Write a script to index the genome fasta file using bowtie2, calling it index_fasta.sh; name the prefix rice_genome
#### Print out script and submit the job (10 pts)


```bash
module avail bowtie2
```

    
    --------------------------- /share/apps/modulefiles ----------------------------
       bowtie2/2.3.2    bowtie2/2.4.2    bowtie2/2.4.4 (L)
    
      Where:
       L:  Module is loaded



```bash
module load bowtie2/2.4.4
```


```bash
bowtie2
```

    No index, query, or output file specified!
    Bowtie 2 version 2.4.4 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)
    Usage: 
      bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>} [-S <sam>]
    
      <bt2-idx>  Index filename prefix (minus trailing .X.bt2).
                 NOTE: Bowtie 1 and Bowtie 2 indexes are not compatible.
      <m1>       Files with #1 mates, paired with files in <m2>.
                 Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
      <m2>       Files with #2 mates, paired with files in <m1>.
                 Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
      <r>        Files with unpaired reads.
                 Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
      <i>        Files with interleaved paired-end FASTQ/FASTA reads
                 Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
      <bam>      Files are unaligned BAM sorted by read name.
      <sam>      File for SAM output (default: stdout)
    
      <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be
      specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'.
    
    Options (defaults in parentheses):
    
     Input:
      -q                 query input files are FASTQ .fq/.fastq (default)
      --tab5             query input files are TAB5 .tab5
      --tab6             query input files are TAB6 .tab6
      --qseq             query input files are in Illumina's qseq format
      -f                 query input files are (multi-)FASTA .fa/.mfa
      -r                 query input files are raw one-sequence-per-line
      -F k:<int>,i:<int> query input files are continuous FASTA where reads
                         are substrings (k-mers) extracted from a FASTA file <s>
                         and aligned at offsets 1, 1+i, 1+2i ... end of reference
      -c                 <m1>, <m2>, <r> are sequences themselves, not files
      -s/--skip <int>    skip the first <int> reads/pairs in the input (none)
      -u/--upto <int>    stop after first <int> reads/pairs (no limit)
      -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)
      -3/--trim3 <int>   trim <int> bases from 3'/right end of reads (0)
      --trim-to [3:|5:]<int> trim reads exceeding <int> bases from either 3' or 5' end
                         If the read end is not specified then it defaults to 3 (0)
      --phred33          qualities are Phred+33 (default)
      --phred64          qualities are Phred+64
      --int-quals        qualities encoded as space-delimited integers
    
     Presets:                 Same as:
      For --end-to-end:
       --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50
       --fast                 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50
       --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)
       --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
    
      For --local:
       --very-fast-local      -D 5 -R 1 -N 0 -L 25 -i S,1,2.00
       --fast-local           -D 10 -R 2 -N 0 -L 22 -i S,1,1.75
       --sensitive-local      -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)
       --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
    
     Alignment:
      -N <int>           max # mismatches in seed alignment; can be 0 or 1 (0)
      -L <int>           length of seed substrings; must be >3, <32 (22)
      -i <func>          interval between seed substrings w/r/t read len (S,1,1.15)
      --n-ceil <func>    func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)
      --dpad <int>       include <int> extra ref chars on sides of DP table (15)
      --gbar <int>       disallow gaps within <int> nucs of read extremes (4)
      --ignore-quals     treat all quality values as 30 on Phred scale (off)
      --nofw             do not align forward (original) version of read (off)
      --norc             do not align reverse-complement version of read (off)
      --no-1mm-upfront   do not allow 1 mismatch alignments before attempting to
                         scan for the optimal seeded alignments
      --end-to-end       entire read must align; no clipping (on)
       OR
      --local            local alignment; ends might be soft clipped (off)
    
     Scoring:
      --ma <int>         match bonus (0 for --end-to-end, 2 for --local) 
      --mp <int>         max penalty for mismatch; lower qual = lower penalty (6)
      --np <int>         penalty for non-A/C/G/Ts in read/ref (1)
      --rdg <int>,<int>  read gap open, extend penalties (5,3)
      --rfg <int>,<int>  reference gap open, extend penalties (5,3)
      --score-min <func> min acceptable alignment score w/r/t read length
                         (G,20,8 for local, L,-0.6,-0.6 for end-to-end)
    
     Reporting:
      (default)          look for multiple alignments, report best, with MAPQ
       OR
      -k <int>           report up to <int> alns per read; MAPQ not meaningful
       OR
      -a/--all           report all alignments; very slow, MAPQ not meaningful
    
     Effort:
      -D <int>           give up extending after <int> failed extends in a row (15)
      -R <int>           for reads w/ repetitive seeds, try <int> sets of seeds (2)
    
     Paired-end:
      -I/--minins <int>  minimum fragment length (0)
      -X/--maxins <int>  maximum fragment length (500)
      --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)
      --no-mixed         suppress unpaired alignments for paired reads
      --no-discordant    suppress discordant alignments for paired reads
      --dovetail         concordant when mates extend past each other
      --no-contain       not concordant when one mate alignment contains other
      --no-overlap       not concordant when mates overlap at all
    
     BAM:
      --align-paired-reads
                         Bowtie2 will, by default, attempt to align unpaired BAM reads.
                         Use this option to align paired-end reads instead.
      --preserve-tags    Preserve tags from the original BAM record by
                         appending them to the end of the corresponding SAM output.
    
     Output:
      -t/--time          print wall-clock time taken by search phases
      --un <path>        write unpaired reads that didn't align to <path>
      --al <path>        write unpaired reads that aligned at least once to <path>
      --un-conc <path>   write pairs that didn't align concordantly to <path>
      --al-conc <path>   write pairs that aligned concordantly at least once to <path>
        (Note: for --un, --al, --un-conc, or --al-conc, add '-gz' to the option name, e.g.
        --un-gz <path>, to gzip compress output, or add '-bz2' to bzip2 compress output.)
      --quiet            print nothing to stderr except serious errors
      --met-file <path>  send metrics to file at <path> (off)
      --met-stderr       send metrics to stderr (off)
      --met <int>        report internal counters & metrics every <int> secs (1)
      --no-unal          suppress SAM records for unaligned reads
      --no-head          suppress header lines, i.e. lines starting with @
      --no-sq            suppress @SQ header lines
      --rg-id <text>     set read group id, reflected in @RG line and RG:Z: opt field
      --rg <text>        add <text> ("lab:value") to @RG line of SAM header.
                         Note: @RG line only printed when --rg-id is set.
      --omit-sec-seq     put '*' in SEQ and QUAL fields for secondary alignments.
      --sam-no-qname-trunc
                         Suppress standard behavior of truncating readname at first whitespace 
                         at the expense of generating non-standard SAM.
      --xeq              Use '='/'X', instead of 'M,' to specify matches/mismatches in SAM record.
      --soft-clipped-unmapped-tlen
                         Exclude soft-clipped bases when reporting TLEN
      --sam-append-comment
                         Append FASTA/FASTQ comment to SAM record
    
     Performance:
      -p/--threads <int> number of alignment threads to launch (1)
      --reorder          force SAM output order to match order of input reads
      --mm               use memory-mapped I/O for index; many 'bowtie's can share
    
     Other:
      --qc-filter        filter out reads that are bad according to QSEQ filter
      --seed <int>       seed for random number generator (0)
      --non-deterministic
                         seed rand. gen. arbitrarily instead of using read attributes
      --version          print version information and quit
      -h/--help          print this usage message
    (ERR): bowtie2-align exited with value 1





```bash
bowtie2 -- build
```

    No index, query, or output file specified!
    Bowtie 2 version 2.4.4 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)
    Usage: 
      bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>} [-S <sam>]
    
      <bt2-idx>  Index filename prefix (minus trailing .X.bt2).
                 NOTE: Bowtie 1 and Bowtie 2 indexes are not compatible.
      <m1>       Files with #1 mates, paired with files in <m2>.
                 Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
      <m2>       Files with #2 mates, paired with files in <m1>.
                 Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
      <r>        Files with unpaired reads.
                 Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
      <i>        Files with interleaved paired-end FASTQ/FASTA reads
                 Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
      <bam>      Files are unaligned BAM sorted by read name.
      <sam>      File for SAM output (default: stdout)
    
      <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be
      specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'.
    
    Options (defaults in parentheses):
    
     Input:
      -q                 query input files are FASTQ .fq/.fastq (default)
      --tab5             query input files are TAB5 .tab5
      --tab6             query input files are TAB6 .tab6
      --qseq             query input files are in Illumina's qseq format
      -f                 query input files are (multi-)FASTA .fa/.mfa
      -r                 query input files are raw one-sequence-per-line
      -F k:<int>,i:<int> query input files are continuous FASTA where reads
                         are substrings (k-mers) extracted from a FASTA file <s>
                         and aligned at offsets 1, 1+i, 1+2i ... end of reference
      -c                 <m1>, <m2>, <r> are sequences themselves, not files
      -s/--skip <int>    skip the first <int> reads/pairs in the input (none)
      -u/--upto <int>    stop after first <int> reads/pairs (no limit)
      -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)
      -3/--trim3 <int>   trim <int> bases from 3'/right end of reads (0)
      --trim-to [3:|5:]<int> trim reads exceeding <int> bases from either 3' or 5' end
                         If the read end is not specified then it defaults to 3 (0)
      --phred33          qualities are Phred+33 (default)
      --phred64          qualities are Phred+64
      --int-quals        qualities encoded as space-delimited integers
    
     Presets:                 Same as:
      For --end-to-end:
       --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50
       --fast                 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50
       --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)
       --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
    
      For --local:
       --very-fast-local      -D 5 -R 1 -N 0 -L 25 -i S,1,2.00
       --fast-local           -D 10 -R 2 -N 0 -L 22 -i S,1,1.75
       --sensitive-local      -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)
       --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
    
     Alignment:
      -N <int>           max # mismatches in seed alignment; can be 0 or 1 (0)
      -L <int>           length of seed substrings; must be >3, <32 (22)
      -i <func>          interval between seed substrings w/r/t read len (S,1,1.15)
      --n-ceil <func>    func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)
      --dpad <int>       include <int> extra ref chars on sides of DP table (15)
      --gbar <int>       disallow gaps within <int> nucs of read extremes (4)
      --ignore-quals     treat all quality values as 30 on Phred scale (off)
      --nofw             do not align forward (original) version of read (off)
      --norc             do not align reverse-complement version of read (off)
      --no-1mm-upfront   do not allow 1 mismatch alignments before attempting to
                         scan for the optimal seeded alignments
      --end-to-end       entire read must align; no clipping (on)
       OR
      --local            local alignment; ends might be soft clipped (off)
    
     Scoring:
      --ma <int>         match bonus (0 for --end-to-end, 2 for --local) 
      --mp <int>         max penalty for mismatch; lower qual = lower penalty (6)
      --np <int>         penalty for non-A/C/G/Ts in read/ref (1)
      --rdg <int>,<int>  read gap open, extend penalties (5,3)
      --rfg <int>,<int>  reference gap open, extend penalties (5,3)
      --score-min <func> min acceptable alignment score w/r/t read length
                         (G,20,8 for local, L,-0.6,-0.6 for end-to-end)
    
     Reporting:
      (default)          look for multiple alignments, report best, with MAPQ
       OR
      -k <int>           report up to <int> alns per read; MAPQ not meaningful
       OR
      -a/--all           report all alignments; very slow, MAPQ not meaningful
    
     Effort:
      -D <int>           give up extending after <int> failed extends in a row (15)
      -R <int>           for reads w/ repetitive seeds, try <int> sets of seeds (2)
    
     Paired-end:
      -I/--minins <int>  minimum fragment length (0)
      -X/--maxins <int>  maximum fragment length (500)
      --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)
      --no-mixed         suppress unpaired alignments for paired reads
      --no-discordant    suppress discordant alignments for paired reads
      --dovetail         concordant when mates extend past each other
      --no-contain       not concordant when one mate alignment contains other
      --no-overlap       not concordant when mates overlap at all
    
     BAM:
      --align-paired-reads
                         Bowtie2 will, by default, attempt to align unpaired BAM reads.
                         Use this option to align paired-end reads instead.
      --preserve-tags    Preserve tags from the original BAM record by
                         appending them to the end of the corresponding SAM output.
    
     Output:
      -t/--time          print wall-clock time taken by search phases
      --un <path>        write unpaired reads that didn't align to <path>
      --al <path>        write unpaired reads that aligned at least once to <path>
      --un-conc <path>   write pairs that didn't align concordantly to <path>
      --al-conc <path>   write pairs that aligned concordantly at least once to <path>
        (Note: for --un, --al, --un-conc, or --al-conc, add '-gz' to the option name, e.g.
        --un-gz <path>, to gzip compress output, or add '-bz2' to bzip2 compress output.)
      --quiet            print nothing to stderr except serious errors
      --met-file <path>  send metrics to file at <path> (off)
      --met-stderr       send metrics to stderr (off)
      --met <int>        report internal counters & metrics every <int> secs (1)
      --no-unal          suppress SAM records for unaligned reads
      --no-head          suppress header lines, i.e. lines starting with @
      --no-sq            suppress @SQ header lines
      --rg-id <text>     set read group id, reflected in @RG line and RG:Z: opt field
      --rg <text>        add <text> ("lab:value") to @RG line of SAM header.
                         Note: @RG line only printed when --rg-id is set.
      --omit-sec-seq     put '*' in SEQ and QUAL fields for secondary alignments.
      --sam-no-qname-trunc
                         Suppress standard behavior of truncating readname at first whitespace 
                         at the expense of generating non-standard SAM.
      --xeq              Use '='/'X', instead of 'M,' to specify matches/mismatches in SAM record.
      --soft-clipped-unmapped-tlen
                         Exclude soft-clipped bases when reporting TLEN
      --sam-append-comment
                         Append FASTA/FASTQ comment to SAM record
    
     Performance:
      -p/--threads <int> number of alignment threads to launch (1)
      --reorder          force SAM output order to match order of input reads
      --mm               use memory-mapped I/O for index; many 'bowtie's can share
    
     Other:
      --qc-filter        filter out reads that are bad according to QSEQ filter
      --seed <int>       seed for random number generator (0)
      --non-deterministic
                         seed rand. gen. arbitrarily instead of using read attributes
      --version          print version information and quit
      -h/--help          print this usage message
    (ERR): bowtie2-align exited with value 1





```bash
cat index_fasta.sh
```

    #!/bin/bash
    #SBATCH --cpus-per-task=6
    #SBATCH --time=5:00:00
    #SBATCH --mem=10GB
    #SBATCH --job-name=bowtie
    
    
    module purge
    module load bowtie2/2.4.4
    
    
    bowtie2-build -f GCA_001433935.1_IRGSP-1.0_genomic.fna rice_genome



```bash
sbatch index_fasta.sh 
```

    Submitted batch job 17587053


#### Write code to align the provided fastq files and produce sorted and indexed bam files within a new sbatch script called run_bowtie_alignment.sh

#### Retain the prefixes of the file names, so that you eventually end up with Input.sorted.bam, Rep_1.sorted.bam, Rep_2.sorted.bam, and Rep_3.sorted.bam , along with their indexes

#### Print out the contents of the sbatch script here (20 pts)


```bash
cat run_bowtie_alignment.sh
```

    #!/bin/bash
    #SBATCH --cpus-per-task=8
    #SBATCH --time=8:00:00
    #SBATCH --mem=64GB
    #SBATCH --job-name=bowtie
    
    
    module purge
    module load bowtie2/2.4.4
    
    bowtie2 --threads=8 -x rice_genome -U Input_1.fastq -S Input_1.sam
    bowtie2 --threads=8 -x rice_genome -U Rep_1_1.fastq -S Rep_1.sam
    bowtie2 --threads=8 -x rice_genome -U Rep_2_1.fastq -S Rep_2.sam
    bowtie2 --threads=8 -x rice_genome -U Rep_3_1.fastq -S Rep_3.sam
    
    module load samtools/intel/1.14
    
    samtools view -h -b Input_1.sam > Input_1.bam
    samtools view -h -b Rep_1.sam > Rep_1.bam
    samtools view -h -b Rep_2.sam > Rep_2.bam
    samtools view -h -b Rep_3.sam > Rep_3.bam
    
    samtools sort Input_1.bam -o Input_1.sorted.bam
    samtools sort Rep_1.bam -o Rep_1.sorted.bam
    samtools sort Rep_2.bam -o Rep_2.sorted.bam
    samtools sort Rep_3.bam -o Rep_3.sorted.bam
    
    samtools index Input_1.sorted.bam > Input_1.index
    samtools index Rep_1.sorted.bam > Rep_2.index
    samtools index Rep_2.sorted.bam > Rep_3.index
    samtools index Rep_3.sorted.bam > Rep_4.index


#### Run the job


```bash
sbatch run_bowtie_alignment.sh
```

    Submitted batch job 17587295



```bash

```

#### Observe the size of the files produced


```bash
ls -la
```

    total 102390668
    drwxrwsr-x.  4 sg5096 sg5096        8192 Apr 10 11:32  .
    drwx--S---. 10 sg5096 sg5096        4096 Apr  8 12:48  ..
    drwxr-sr-x.  2 sg5096 sg5096        4096 Apr 10 11:32  .ipynb_checkpoints
    -rw-r--r--.  1 sg5096 sg5096   378473961 Apr  9 22:06  GCA_001433935.1_IRGSP-1.0_genomic.fna
    -rw-r--r--.  1 sg5096 sg5096       38206 Apr 10 11:31  Homework_4.ipynb
    -rw-rw-r--.  1 sg5096 sg5096  2784303385 Apr  9 23:40  Input_1.bam
    -rw-r--r--.  1 sg5096 sg5096  9015057924 Apr  9 22:06  Input_1.fastq
    -rw-rw-r--.  1 sg5096 sg5096           0 Apr 10 11:29  Input_1.index
    -rw-rw-r--.  1 sg5096 sg5096  9972444752 Apr  9 22:41  Input_1.sam
    -rw-rw-r--.  1 sg5096 sg5096  1780836327 Apr 10 00:16  Input_1.sorted.bam
    -rw-rw-r--.  1 sg5096 sg5096      521368 Apr 10 11:30  Input_1.sorted.bam.bai
    -rw-rw-r--.  1 sg5096 sg5096  2789909783 Apr  9 23:48  Rep_1.bam
    -rw-rw-r--.  1 sg5096 sg5096  9760868085 Apr  9 22:56  Rep_1.sam
    -rw-rw-r--.  1 sg5096 sg5096  1817291217 Apr 10 00:25  Rep_1.sorted.bam
    -rw-rw-r--.  1 sg5096 sg5096      505776 Apr 10 11:30  Rep_1.sorted.bam.bai
    -rw-r--r--.  1 sg5096 sg5096  8878485252 Apr  9 22:06  Rep_1_1.fastq
    -rw-rw-r--.  1 sg5096 sg5096  2998328071 Apr  9 23:56  Rep_2.bam
    -rw-rw-r--.  1 sg5096 sg5096           0 Apr 10 11:30  Rep_2.index
    -rw-rw-r--.  1 sg5096 sg5096 10627262595 Apr  9 23:12  Rep_2.sam
    -rw-rw-r--.  1 sg5096 sg5096  1887355571 Apr 10 00:34  Rep_2.sorted.bam
    -rw-rw-r--.  1 sg5096 sg5096      572368 Apr 10 11:31  Rep_2.sorted.bam.bai
    -rw-r--r--.  1 sg5096 sg5096  9662905860 Apr  9 22:06  Rep_2_1.fastq
    -rw-rw-r--.  1 sg5096 sg5096  3833184539 Apr 10 00:07  Rep_3.bam
    -rw-rw-r--.  1 sg5096 sg5096           0 Apr 10 11:30  Rep_3.index
    -rw-rw-r--.  1 sg5096 sg5096 13468346163 Apr  9 23:32  Rep_3.sam
    -rw-rw-r--.  1 sg5096 sg5096  2403591294 Apr 10 00:47  Rep_3.sorted.bam
    -rw-rw-r--.  1 sg5096 sg5096      830632 Apr 10 11:31  Rep_3.sorted.bam.bai
    -rw-r--r--.  1 sg5096 sg5096 12247293828 Apr  9 22:06  Rep_3_1.fastq
    -rw-rw-r--.  1 sg5096 sg5096           0 Apr 10 11:31  Rep_4.index
    drwxr-sr-x.  2 sg5096 sg5096        4096 Apr  9 21:14 'Untitled Folder'
    -rw-r--r--.  1 sg5096 sg5096         760 Apr  8 12:48  hw4.ipynb
    -rw-r--r--.  1 sg5096 sg5096         215 Apr  9 22:37  index_fasta.sh
    -rw-r--r--.  1 sg5096 sg5096         109 Apr  9 21:53  merge_bam.sh
    -rw-rw-r--.  1 sg5096 sg5096   128771101 Apr  9 22:18  rice_genome.1.bt2
    -rw-rw-r--.  1 sg5096 sg5096    93419356 Apr  9 22:18  rice_genome.2.bt2
    -rw-rw-r--.  1 sg5096 sg5096        8702 Apr  9 22:13  rice_genome.3.bt2
    -rw-rw-r--.  1 sg5096 sg5096    93419351 Apr  9 22:13  rice_genome.4.bt2
    -rw-rw-r--.  1 sg5096 sg5096   128771101 Apr  9 22:23  rice_genome.rev.1.bt2
    -rw-rw-r--.  1 sg5096 sg5096    93419356 Apr  9 22:23  rice_genome.rev.2.bt2
    -rw-r--r--.  1 sg5096 sg5096         986 Apr  9 22:36  run_bowtie_alignment.sh
    -rw-rw-r--.  1 sg5096 sg5096       12968 Apr  9 22:23  slurm-17587053.out
    -rw-rw-r--.  1 sg5096 sg5096        1144 Apr 10 00:38  slurm-17587295.out
    -rw-rw-r--.  1 sg5096 sg5096           0 Apr 10 11:29  slurm-17622987.out


#### Merge the replicate BAM files using samtools merge and call the output IDS1.sorted.bam
#### Write an sbatch script called merge_bam.sh to do this; print script and submit job (this can take a while, I recommend asking for 8 threads using the --threads argument; make sure to also set cpus-per-task=8 if you do) (10 pts)


```bash
module avail SAMtools
```

    
    --------------------------- /share/apps/modulefiles ----------------------------
       samtools/intel/1.11    samtools/intel/1.12    samtools/intel/1.14 (L)
    
      Where:
       L:  Module is loaded



```bash
module load samtools/intel/1.14
samtools
```

    
    Program: samtools (Tools for alignments in the SAM format)
    Version: 1.14 (using htslib 1.14)
    
    Usage:   samtools <command> [options]
    
    Commands:
      -- Indexing
         dict           create a sequence dictionary file
         faidx          index/extract FASTA
         fqidx          index/extract FASTQ
         index          index alignment
    
      -- Editing
         calmd          recalculate MD/NM tags and '=' bases
         fixmate        fix mate information
         reheader       replace BAM header
         targetcut      cut fosmid regions (for fosmid pool only)
         addreplacerg   adds or replaces RG tags
         markdup        mark duplicates
         ampliconclip   clip oligos from the end of reads
    
      -- File operations
         collate        shuffle and group alignments by name
         cat            concatenate BAMs
         merge          merge sorted alignments
         mpileup        multi-way pileup
         sort           sort alignment file
         split          splits a file by read group
         quickcheck     quickly check if SAM/BAM/CRAM file appears intact
         fastq          converts a BAM to a FASTQ
         fasta          converts a BAM to a FASTA
         import         Converts FASTA or FASTQ files to SAM/BAM/CRAM
    
      -- Statistics
         bedcov         read depth per BED region
         coverage       alignment depth and percent coverage
         depth          compute the depth
         flagstat       simple stats
         idxstats       BAM index stats
         phase          phase heterozygotes
         stats          generate stats (former bamcheck)
         ampliconstats  generate amplicon specific stats
    
      -- Viewing
         flags          explain BAM flags
         tview          text alignment viewer
         view           SAM<->BAM<->CRAM conversion
         depad          convert padded BAM to unpadded BAM
         samples        list the samples in a set of SAM/BAM/CRAM files
    
      -- Misc
         help [cmd]     display this help message or help for [cmd]
         version        detailed version information
    





```bash
samtools merge
```

    Usage: samtools merge [options] -o <out.bam> [options] <in1.bam> ... <inN.bam>
       or: samtools merge [options] <out.bam> <in1.bam> ... <inN.bam>
    
    Options:
      -n         Input files are sorted by read name
      -t TAG     Input files are sorted by TAG value
      -r         Attach RG tag (inferred from file names)
      -u         Uncompressed BAM output
      -f         Overwrite the output BAM if exist
      -o FILE    Specify output file via option instead of <out.bam> argument
      -1         Compress level 1
      -l INT     Compression level, from 0 to 9 [-1]
      -R STR     Merge file in the specified region STR [all]
      -h FILE    Copy the header in FILE to <out.bam> [in1.bam]
      -c         Combine @RG headers with colliding IDs [alter IDs to be distinct]
      -p         Combine @PG headers with colliding IDs [alter IDs to be distinct]
      -s VALUE   Override random seed
      -b FILE    List of input BAM filenames, one per line [null]
      -X         Use customized index files
      -L FILE    Specify a BED file for multiple region filtering [null]
      --no-PG    do not add a PG line
          --input-fmt-option OPT[=VAL]
                   Specify a single input file format option in the form
                   of OPTION or OPTION=VALUE
      -O, --output-fmt FORMAT[,OPT[=VAL]]...
                   Specify output format (SAM, BAM, CRAM)
          --output-fmt-option OPT[=VAL]
                   Specify a single output file format option in the form
                   of OPTION or OPTION=VALUE
          --reference FILE
                   Reference sequence FASTA FILE [null]
      -@, --threads INT
                   Number of additional threads to use [0]
          --write-index
                   Automatically index the output files [off]
          --verbosity INT
                   Set level of verbosity



```bash
cat merge_bam.sh
```

    #!/bin/bash
    #SBATCH --cpus-per-task=8
    #SBATCH --time=5:00:00
    #SBATCH --mem=10GB
    #SBATCH --job-name=samtools
    
    module load samtools/intel/1.14
    
    samtools merge --threads=8 -o IDS1.sorted.bam Rep_1.sorted.bam Rep_2.sorted.bam Rep_3.sorted.bam 


```bash
sbatch merge_bam.sh
```

    Submitted batch job 17628299



```bash

```

#### Use MACS2 to call peaks using Input as control and the merged bam file as treatment

#### When would you want to look at broad peaks? (10 pts)


```bash
module avail macs2
```

    
    --------------------------- /share/apps/modulefiles ----------------------------
       macs2/2.1.1.20160309    macs2/intel/2.2.7.1 (L)
    
      Where:
       L:  Module is loaded



```bash
module load macs2/intel/2.2.7.1
```

I would want to look at broad peaks to find histone modifications that cover entire gene bodies. 

#### Print out the contents of the sbatch script and submit the job (10 pts)


```bash
cat macs.sh
```

    #!/bin/bash
    #SBATCH --cpus-per-task=8
    #SBATCH --time=5:00:00
    #SBATCH --mem=10GB
    #SBATCH --job-name=samtools
    
    module load macs2/intel/2.2.7.1
    
    macs2 callpeak \
    -t IDS1.sorted.bam \
    -c Input_1.bam \
    --format=BAM  \
    --gsize=380000000 \
    --cutoff-analysis  \
    --qvalue=0.05 \
    --outdir=macs2_out \
    --name IDS1



```bash
sbatch macs.sh
```

    Submitted batch job 17794925



```bash
seff 17794925
```

    Job ID: 17794925
    Cluster: greene
    User/Group: sg5096/sg5096
    State: COMPLETED (exit code 0)
    Nodes: 1
    Cores per node: 8
    CPU Utilized: 00:05:40
    CPU Efficiency: 12.35% of 00:45:52 core-walltime
    Job Wall-clock time: 00:05:44
    Memory Utilized: 707.25 MB
    Memory Efficiency: 6.91% of 10.00 GB


#### We will try running Homer's findMotifsGenome.pl command with no background sequences provided; as per the documentation:
[Documentation](http://homer.ucsd.edu/homer/ngs/peakMotifs.html)
#### "If the background sequences were not explicitly defined, HOMER will automatically select them for you.  If you are using genomic positions, sequences will be randomly selected from the genome, matched for GC% content."
#### Create a new folder called motif_analysis to dump the results from Homer
#### Write and print an sbatch script called run_homer.sh to do this and submit the job (10 pts)


```bash
module avail homer
```

    
    --------------------------- /share/apps/modulefiles ----------------------------
       homer/4.11



```bash
module load homer/4.11
```


```bash
cat run_homer.sh
```

    #!/bin/bash
    #SBATCH --cpus-per-task=8
    #SBATCH --time=8:00:00
    #SBATCH --mem=64GB
    #SBATCH --job-name=homer
    
    module purge
    module load homer/4.11
    
    findMotifsGenome.pl macs2_out/IDS1_peaks.narrowPeak GCA_001433935.1_IRGSP-1.0_genomic.fna motif_analysis/ -size 200 -len 8


```bash
sbatch run_homer.sh
```

    Submitted batch job 17796371



```bash
seff 17796371
```

    Job ID: 17796371
    Cluster: greene
    User/Group: sg5096/sg5096
    State: COMPLETED (exit code 0)
    Nodes: 1
    Cores per node: 8
    CPU Utilized: 00:26:38
    CPU Efficiency: 12.43% of 03:34:16 core-walltime
    Job Wall-clock time: 00:26:47
    Memory Utilized: 607.49 MB
    Memory Efficiency: 0.93% of 64.00 GB


#### Look at the known and de novo motifs detected by Homer 
#### Across the known and de novo motifs, what are the top three most enriched motifs (note down their names)? (15 pts)
For the known motifs, the most enriched motifs with the lowest p-values are: 
1. TF3A(C2H2)/col-TF3A-DAP-Seq(GSE60143)/Homer
2. WRKY55(WRKY)/col-WRKY55-DAP-Seq(GSE60143)/Homer
3. WRKY75(WRKY)/col-WRKY75-DAP-Seq(GSE60143)/Homer

For the de novo motifs, the most enriched motifs are:
1. PB0133.1_Hic1_2/Jaspar(0.829)
2. ZNF711(Zf)/SHSY5Y-ZNF711-ChIP-Seq(GSE20673)/Homer(0.794)
3. RAV1(1)(AP2/EREBP)/Arabidopsis thaliana/AthaMap(0.819)
#### Compare these motifs to the ones reported in Figure 3 in the original [paper](https://academic.oup.com/plphys/article/178/2/824/6116647). Comment on similarities/differences. (5 pts)
The motifs in figure 3 of the original paper are predominantly T,C and G, specifically TCCTCC and GCCGCC. The known motifs are primarily T, A, and G while the de novo motifs are primarily T, C and G. The known motifs have T and G in common while the de novo motifs have all three bases in common. 
#### Looking at the results from the paper (generated using MEME) and your results using Homer, how would you approach a similar problem in your own work in the future? This is not asking for any single correct answer; just think of ways you would trust your results more. (10 pts)
I think I would approach a similar problem in my own work in the future similarly to how I did in this assignment. I would also maybe do a functional analysis along with a sequence analysis to look for potential differences. 