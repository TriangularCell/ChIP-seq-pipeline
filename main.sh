# This pipeline is for ChIP-seq analysis, start from fastq raw sequencing files.
# 1. fastqc check data quality
# 1.2. Trimmomatic filter low quality data/cut adapters
# 2. bwa/bowtie map reads to reference genome
# 3. samtools convert sam file format to bam format, and check alignment quality
# 4. MACS2/HOMER call peaks
# 5. HOMER peak annotation, motif detection/enrichment
# 6. UCSC Visualization
# 7. Comparative analysis
# 8. Integrative analysis
# 9. Others

############################## fastqc ################################
for file in ./*;
do
        /software/sequencing/FastQC_0.10.1/fastqc $file -t 10 -o $outputFolder
done

# FastQC
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# SYNOPSIS
fastqc seqfile1 seqfile2 .. seqfileN
fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
           [-c contaminant file] seqfile1 .. seqfileN

# Summary statistics and total sequences can be extracted from output "fastqc_data.txt"

############################## seqtk ################################
# Convert fastq file to fasta file
seqtk seq -a $file > $fafile

# seqtk
# https://github.com/lh3/seqtk

samtools bam2fq $bamfile > $fqfile

# If the $bamfile is from the fastq file, then the converted back fastq file is the same with its original file except without index.

############################## blat ################################
# Align a set of sequences to a set of sequences
blat $database $query $output".psl"

# blat
# https://genome.ucsc.edu/FAQ/FAQblat.html
# type "blat" to see more options

############################## Trimmomatic ############################
java -jar /home/yiming/TBRS/Tools/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 12 -trimlog $input_1 $input_2 $paired_1 $unpaired_1 $paired_2 $unpaired_2 ILLUMINACLIP:/home/yiming/TBRS/Tools/Trimmomatic-0.32/adapters/adapter0817-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:20

# Trimmomatic
# http://www.usadellab.org/cms/?page=trimmomatic

# Paired End Mode:
java -jar <path to trimmomatic.jar> PE [-threads <threads] [-phred33 | -phred64] [-trimlog <logFile>] <input 1> <input 2> <paired output 1> <unpaired output 1> <paired output 2> <unpaired output 2> <step 1> ...

or
java -classpath <path to trimmomatic jar> org.usadellab.trimmomatic.TrimmomaticPE [-threads <threads>] [-phred33 | -phred64] [-trimlog <logFile>] <input 1> <input 2> <paired output 1> <unpaired output 1> <paired output 2> <unpaired output 2> <step 1> ...

# Single End Mode:
java -jar <path to trimmomatic jar> SE [-threads <threads>] [-phred33 | -phred64] [-trimlog <logFile>] <input> <output> <step 1> ...

or
java -classpath <path to trimmomatic jar> org.usadellab.trimmomatic.TrimmomaticSE [-threads <threads>] [-phred33 | -phred64] [-trimlog <logFile>] <input> <output> <step 1> ...

############################### Alignment_bwa ##########################
twoBitToFa mm10.2bit mm10.fa
bwa index mm10.fa

bwa mem -t 10 /home/groups/kathryncheah/yiming/Genomes/mm10_UCSC_BWA_index/mm10.fa $input_1 $input_2 > $output".sam"

# BWA
# http://bio-bwa.sourceforge.net/bwa.shtml#7

# It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp. 
# BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads.

-M      Mark shorter split hits as secondary (for Picard compatibility).
# short split hists are sets of short reads in "Chimeric alignment". 
# Reference for the following: 
# http://samtools.github.io/hts-specs/SAMv1.pdf

# A chimeric alignment is primarily caused by structural variations, gene fusions, misassemblies, RNA-seq or experimental protocols. It is more frequent given longer reads. For a chimeric alignment, the linear alignments consisting of the aligment are largely non-overlapping; each linear alignment may have high mapping quality and is informative in SNP/INDEL calling. 
# In contrast, multiple mappings are caused primarily by repeats. They are less frequent given longer reads. If a read has multiple mappings, all these mappings are almost entirely overlapping with each other; except the single-best optimal mapping, all the other mappings get mapping quality <Q3 and are ignored by most SNP/INDEL callers.

# Typically, one of the linear alignments in a chimeric alignment is considered the "representative" alignment, and the others are called "supplementary" and are distinguished by the supplementary alignment flag. 
# In multiple mapping, One of these alignments is considered "primary". All the other alignments have the "secondary" alignment flag set in the SAM records that represent them.

# In my data, I found a lot of "supplementary alignments" but no "secondary alignment". And answers from BioStars also said aligners generally don't report chimeric alignments.

# So when "-M" is set, "supplementary alignment" (flag 2048) will be flaged as "secondary alignment" ("not primary alignment", flag 256). As I tested, it doesn't affect the result of Picard MarkDuplicates. Seems MarkDuplicates skips all 2048 and 256 reads since the number of alignments with -f 2048 and -f 256 are exactly the same before and after remove duplicates.

# Explaination for SAM flags: 
# http://broadinstitute.github.io/picard/explain-flags.html

################################ samtools_view ##############################
/software/sequencing/samtools-0.1.18/samtools view -bhS $input_sam > $output_bam
/software/sequencing/samtools-0.1.18/samtools sort $output_bam $output_bam"_sort"
/software/sequencing/samtools-0.1.18/samtools index $output_bam"_sort.bam"
/software/sequencing/samtools-0.1.18/samtools idxstats $output_bam"_sort.bam"
# idxstats check the number of mapped reads and unmapped reads for each sequence name. Usually the sequence name is chromosome number. 

/software/sequencing/samtools-0.1.18/samtools flagstat $file
#1      15701907 + 0 in total (QC-passed reads + QC-failed reads)
#2      0 + 0 duplicates
#3      15288534 + 0 mapped (97.37%:-nan%)
#4      15701907 + 0 paired in sequencing
#5      7849559 + 0 read1
#6      7852348 + 0 read2
#7      14728567 + 0 properly paired (93.80%:-nan%)
#8      15199161 + 0 with itself and mate mapped
#9      89373 + 0 singletons (0.57%:-nan%)
#10     356537 + 0 with mate mapped to a different chr
#11     261805 + 0 with mate mapped to a different chr (mapQ>=5)

# Total number of reads (include chimeric reads and multiple mapped reads), so larger than total number of raw reads count by FastQC. Same as flagstat #1
/software/sequencing/samtools-0.1.18/samtools view -c $file
# Number of read1. Same as flagstat #5
/software/sequencing/samtools-0.1.18/samtools view -c -f 0x0040 $file
# Number of read2. Same as flagstat #6
/software/sequencing/samtools-0.1.18/samtools view -c -f 0x0080 $file
# Total number of mapped reads. Same as flagstat #3
/software/sequencing/samtools-0.1.18/samtools view -c -F 4 $file
# Total number of proper paired reads. Same as flagstat #7
/software/sequencing/samtools-0.1.18/samtools view -c -f 0x2 $file


# Samtools 
# http://www.htslib.org/doc/samtools.html

Usage:   samtools view [options] <in.bam>|<in.sam> [region1 [...]]

Options: -b       output BAM
         -h       print header for the SAM output
         -H       print header only (no alignments)
         -S       input is SAM
         -u       uncompressed BAM output (force -b)
         -1       fast compression (force -b)
         -x       output FLAG in HEX (samtools-C specific)
         -X       output FLAG in string (samtools-C specific)
         -c       print only the count of matching records
         -L FILE  output alignments overlapping the input BED FILE [null]
         -t FILE  list of reference names and lengths (force -S) [null]
         -T FILE  reference sequence file (force -S) [null]
         -o FILE  output file name [stdout]
         -R FILE  list of read groups to be outputted [null]
         -f INT   required flag, 0 for unset [0]
         -F INT   filtering flag, 0 for unset [0]
         -q INT   minimum mapping quality [0]
         -l STR   only output reads in library STR [null]
         -r STR   only output reads in read group STR [null]
         -s FLOAT fraction of templates to subsample; integer part as seed [-1]
         -?       longer help

################################## Picard_CollectAlignmentSummaryMetrics ###########
java -Xmx10g -jar /home/groups/kathryncheah/yiming/Tools/picard-tools-1.141/picard.jar CollectAlignmentSummaryMetrics R=$reference I=$input O=$output
# Input should be coordinates sorted

# Picard
# Summarize the alignment results from BAM/SAM files.

# Usage
# http://broadinstitute.github.io/picard/command-line-overview.html#CollectAlignmentSummaryMetrics

# Output explanation
# http://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics

# The results of items reported by Picard are all different from Samtools because they do not include chimeric mapping and multiple mappings.

################################## Picard_MarkDuplicates ###########################
java -Xmx10g -jar /home/groups/kathryncheah/yiming/Tools/picard-tools-1.141/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$input O=$removed_duplicates.bam M=marked_dup_metrics.txt ASSUME_SORTED=true

# Picard
# Mark/Remove duplicates in BAM/SAM files.

# Usage
# https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates
# ASSUME_SORTED (Boolean): If true, assume that the input file is coordinate sorted even if the header says otherwise. 

# ASSUME_SORT_ORDER (SortOrder): If not null, assume that the input file has this order even if the header says otherwise. Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate} Cannot be used in conjuction with option(s) ASSUME_SORTED (AS)
# Unfortunately, version 1.141 doesn't support this parameter...
# When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates. 
# When the input is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads.


# Output explanation
# http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics

################################## Picard_CollectInsertSizeMetrics ########################
java -Xmx10g -jar /home/groups/kathryncheah/yiming/Tools/picard-tools-1.141/picard.jar CollectInsertSizeMetrics I=$input O=$insert_size_metrics H=$insert_size_histogram

# Picard
# Show the distribution of insert sizes of paired end reads.

#Usage
# https://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics

# Output explanation
# https://broadinstitute.github.io/picard/picard-metric-definitions.html#InsertSizeMetrics

################################## Picard_QualityScoreDistribution ########################
java -Xmx10g -jar /home/groups/kathryncheah/yiming/Tools/picard-tools-1.141/picard.jar QualityScoreDistribution I=$input O=$qual_score_dist CHART=$qual_score_dist_chart

# Picard
# Show the distribution of per base mapping quality scores from BAM/SAM files.

# Usage
# https://broadinstitute.github.io/picard/command-line-overview.html#QualityScoreDistribution

################################## MACS2 peak calling ########################
macs2 callpeak -t $chipfile -c $controlfile -g mm -n $name -f BAMPE --broad -B --SPMR --trackline --outdir $outF$name"_broad"

# MACS2
# Usage
# https://github.com/taoliu/MACS

# macs2 [-h] [--version] {callpeak,filterdup,bdgpeakcall,bdgcmp,randsample,bdgdiff,bdgbroadcall}

# Output "NAME_peaks.xls"
1. chromosome name
2. start position of peak
3. end position of peak
4. length of peak region
5. absolute peak summit position
6. pileup height at peak summit, -log10(pvalue) for the peak summit (e.g. pvalue =1e-10, then this value should be 10)
7. fold enrichment for this peak summit against random Poisson distribution with local lambda, -log10(qvalue) at peak summit

# Output "NAME_peaks.narrowPeak". BED6+4 format
5. integer score for display
7. fold-change
8. -log10pvalue
9. -log10qvalue
10. relative summit position to peak start

################################## MACS2 bdgcmp ########################
macs2 bdgcmp -t $name"_treat_pileup.bdg" -c $name"_control_lambda.bdg" -o $name"_FE.bdg" -m FE

# MACS2
# Usage
# https://github.com/taoliu/MACS/wiki/Build-Signal-Track

Options:
-m FE means to calculate fold enrichment. Other options can be logLR for log likelihood, subtract for subtracting noise from treatment sample.
-p sets pseudocount. This number will be added to 'pileup per million reads' value. You don't need it while generating fold enrichment track because control lambda will always >0. But in order to avoid log(0) while calculating log likelihood, we'd add pseudocount. Because I set precision as 5 decimals, here I use 0.00001.


################################## HOMER makeTagDirectory ########################
makeTagDirectory $outF$name"_TagDir" $file

# HOMER
# Usage
# http://homer.ucsd.edu/homer/ngs/tagDir.html

# makeTagDirectory <Output Directory Name> [options] <alignment file1> [alignment file 2] ...

# Be attention when input file is BED format.

################################## HOMER findPeaks ########################
findPeaks $outF$name"_TagDir" -style factor -o $outF$name"/"$name"_peaks.txt"

# HOMER
# Usage
# http://homer.ucsd.edu/homer/ngs/peaks.html

# findPeaks <tag directory> -style <factor | histone> -o auto -i <control tag directory>

# Below the header information are the peaks, listed in each row.  Columns contain information about each peak:
Column 1: PeakID - a unique name for each peak (very important that peaks have unique names...)
Column 2: chr - chromosome where peak is located
Column 3: starting position of peak
Column 4: ending position of peak
Column 5: Strand (+/-)
Column 6: Normalized Tag Counts - number of tags found at the peak, normalized to 10 million total mapped tags (or defined by the user)
Column 7: (-style factor): Focus Ratio - fraction of tags found appropriately upstream and downstream of the peak center. (see below)
                 (-style histone/-style groseq): Region Size - length of enriched region
Column 8: Peak score (position adjusted reads from initial peak region - reads per position may be limited)
Columns 9+: Statistics and Data from filtering

################################## HOMER annotatePeaks ########################
annotatePeaks.pl $file mm10 > $outF$name.xls

# HOMER
# Usage
# http://homer.ucsd.edu/homer/ngs/annotation.html

# annotatePeaks.pl <peak/BED file> <genome> > <output file>

# Description of Columns:
1. Peak ID
2. Chromosome
3. Peak start position
4. Peak end position
5. Strand
6. Peak Score
7. FDR/Peak Focus Ratio/Region Size
8. Annotation (i.e. Exon, Intron, ...)
9. Detailed Annotation (Exon, Intron etc. + CpG Islands, repeats, etc.)
10. Distance to nearest RefSeq TSS
11. Nearest TSS: Native ID of annotation file
12. Nearest TSS: Entrez Gene ID
13. Nearest TSS: Unigene ID
14. Nearest TSS: RefSeq ID
15. Nearest TSS: Ensembl ID
16. Nearest TSS: Gene Symbol
17. Nearest TSS: Gene Aliases
18. Nearest TSS: Gene description
Additional columns depend on options selected when running the program.

################################## HOMER findMotifs ########################
findMotifsGenome.pl $file mm10 $outF$name"_motifs" -size given -p 2

# HOMER
# Usage
# http://homer.ucsd.edu/homer/ngs/peakMotifs.html

# findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]

################################## MAnorm differential ChIP-seq ################################
bash ./MAnorm.sh  $sample1_peaks $sample2_peaks $sample1_read $sample2_read $sample1_readshift_lentgh $sample2_readshift_length

# MAnorm
# http://bcb.dfci.harvard.edu/~gcyuan/MAnorm/R_tutorial.html
# https://github.com/ying-w/chipseq-compare/tree/master/MAnorm

# ./MAnorm.sh  sample1_peaks.bed  sample2_peaks.bed  sample1_read.bed  sample2_read.bed sample1_readshift_lentgh[INT] sample2_readshift_length[INT]
# The first 2 files have ONLY 3 columns: chromosome, start, end.
# The next 2 files should have 4 columns: chromosome, start, end, strand (+/-)


################################## Visualization ################################
1. IGV
2. UCSC genome Browser
3. Gviz

#2. Example custom tracks:
track type=bigWig visibility=full name="NPC2-K" description="NPC2-K_normed_coverage" bigDataUrl=ftp://ftpqj:hku_wanglabqj@147.8.193.64/NPC2/NPC2-K_trim_noDuplicate_sort_name_fixed_0x2_sameChr_sort_fs500_norm_CP100M_sort.bw itemRgb=On color=34,139,34

track type=bigWig visibility=full name="NPC2-I" description="NPC2-I_normed_coverage" bigDataUrl=ftp://ftpqj:hku_wanglabqj@147.8.193.64/NPC2/NPC2-I_noDuplicate_sort_name_fixed_0x2_sameChr_sort_fs500_norm_CP100M_sort.bw itemRgb=On color=34,139,34

track type=bigBed visibility=dense name="NPC2_peaks" description="NPC2_distal_peaks" bigDataUrl=ftp://ftpqj:hku_wanglabqj@147.8.193.64/NPC2/NPC2-NPC2_Homer_UniqTag_peaks_500bp_distal_CoSort.bb itemRgb=On color=34,139,34

#3. References for Gviz
# https://bioconductor.org/packages/release/bioc/html/Gviz.html
# https://link-springer-com.eproxy2.lib.hku.hk/content/pdf/10.1007%2F978-1-4939-3578-9_16.pdf
# https://www.bioconductor.org/help/course-materials/2015/CSAMA2015/lab/Epigenetics_and_Chip_seqLab.pdf

