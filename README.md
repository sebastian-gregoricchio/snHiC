[![Snakemake](https://img.shields.io/badge/snakemake-≥7.8.5-brightgreen.svg)](https://snakemake.github.io)
![release](https://img.shields.io/github/v/release/sebastian-gregoricchio/snHiC)
[![license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://sebastian-gregoricchio.github.io/snHiC/LICENSE.md/LICENSE)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/snHiC?style=social)](https://github.com/sebastian-gregoricchio/snHiC/fork)
<!-- ![update](https://badges.pufler.dev/updated/sebastian-gregoricchio/snHiC)
![visits](https://badges.pufler.dev/visits/sebastian-gregoricchio/snHiC)
![downloads](https://img.shields.io/github/downloads/sebastian-gregoricchio/snHiC/total.svg)--->

<div style="text-align: justify">

# snHiC [<img src="https://raw.githubusercontent.com/sebastian-gregoricchio/snHiC/main/resources/snHiC_logo.svg" align="right" height = 150/>](https://sebastian-gregoricchio.github.io/snHiC)
## Introduction
`snHiC` is a snakemake based end-to-end pipeline to analyze Hi-C data. The input files required to run the pipeline are Paired-End fastq files. The pipeline performs data quality control, normalization and correction. It also includes the possibility to perform grouped analyses (e.g, merging of replicates) besides TAD calling, loops detection and differential compartment analyses. Notability, the latter is performed using `dcHiC`, a recently published method ([A. Chakraborty, *et al.*, Nat. Comm. 2022](https://www.nature.com/articles/s41467-022-34626-6)) that enables more precise and high-resolution differential compartment analyses.

</div>

### Citation
If you use this package, please cite:



<div class="warning" style='padding:2.5%; background-color:#ffffee; color:#787878; margin-left:5%; margin-right:5%; border-radius:15px;'>
<span>
<font size="-0.5">

<div style="margin-left:2%; margin-right:2%; text-align: justify">
*--- No publication associated yet ---*
</div>
</font>

</span>
</div>

<br/><br/>

## Installation an dependencies
<div style="text-align: justify">
To install the pipeline it is required to download this repository and the installation of a conda environment is strongly recommended.
Follow the steps below for the installation:
* place yourself in the directory where the repository should be downloaded by typing `cd </target/folder>`
* download the GitHub repository with `git clone https://github.com/sebastian-gregoricchio/snakeATAC`, or click on *Code > Download ZIP* on the [GitHub page](https://github.com/sebastian-gregoricchio/snakeATAC)
* install the conda environment from the .yaml environment file contained in the repository:<br>
`conda env create -f </target/folder>/snHiC/workflow/envs/snHiC_conda_env_stable.yaml`
* activate the conda environment: `conda activate snHiC` (if the env is not activated the pipeline won't work properly)

If you want to run the differential compartment analyses by [`dcHiC`](https://www.nature.com/articles/s41467-022-34626-6), you need to install this tools first. Consider that all the pre-required packages are already included in the `snHiC` environment.
* place yourself in the directory where the repository should be downloaded by typing `cd </target/folder>`
* download the GitHub repository with `git clone https://github.com/ay-lab/dcHiC`, or click on *Code > Download ZIP* on the [GitHub page](https://github.com/ay-lab/dcHiC)
* activate, if not already done, the `snHiC` conda environment: `conda activate snHiC` (if the env is not activated the dcHiC functions will be installed in the wrong R environment)
* install the `dcHiC`'s functions available in the repository --> *dcHiC/packages/functionsdchic_1.0.tar.gz*:<br>
`${CONDA_PREFIX}/bin/R CMD INSTALL </target/folder>/dcHiC/packages/functionsdchic_1.0.tar.gz`
</div>

<br/><br/>

## How to run the pipeline
<div style="text-align: justify">
The basic `snHiC` pipeline requires only two files:
* `snHiC.snakefile`, containing all the rules that will be run;
* `configuration.yaml` file, in which the user can define and customize all the parameters for the different pipeline steps. <br>

Then, if grouped or differential compartment analyses are required by the user, a third file is required. This file will contain the sample names (without the read file suffix and extension, e.g *sample.A_R1.fastq.gz* --> *sample.A*) and the group to which they belong. Here follows an example (tab-delimited txt file):
</div>

| *sample* | *group* |
|:---------|:--------|
| Sample_A | Tumor   |
| Sample_B | Normal  |
| Sample_C | Normal  |
| Sample_D | Tumor   |
| Sample_E | Normal  |
| Sample_F | Tumor   |

<div style="text-align: justify">
To partially avoid unexpected errors during the execution of the pipeline, a so called 'dry-run' is strongly recommended. Indeed, adding a `-n` at the end of the snakemake running command will allow snakemake to check that all links and file/parameters dependencies are satisfied before to run the "real" processes. This command will therefore help the debugging process.
</div>

```shell
snakemake \
-s </target/folder>/snHiC/workflow/snHiC.snakefile \
--configfile </target/folder>/snHiC/config/snHiC_config.yaml \
--cores 12 \
-n
```
<div style="text-align: justify">
If no errors occur, the pipeline can be run with the same command line without the terminal `-n`:
</div>
```shell
snakemake \
-s </target/folder>/snHiC/workflow/snHiC.snakefile \
--configfile </target/folder>/snHiC/config/snHiC_config.yaml \
--cores 12
```
<div style="text-align: justify">
Notice that the absence of errors does not mean that the pipeline will run without any issues; the "dry-run" is only checking whether all the resources are available. <br>
Furthermore, there is the possibility to run the pipeline only partially. An example of usage could be if someone wants to have a look to the fast quality controls (fastQC and multiQC reports) before to perform the alignment. To do that, it is sufficient to run a dry-run (`-n` mode), then pick the name of the rule at which you want the pipeline to stop and lastly type the following command:
</div>

```shell
snakemake \
-s </target/folder>/snHiC/workflow/snHiC.snakefile \
--configfile </target/folder>/snHiC/config/snHiC_config.yaml \
--cores 12 \
--until <stop_rule_name>
```

<br/><br/>

### Snakefile
<div style="text-align: justify">
The snakefile are contained all the condition checks and rules (processing steps) that will be performed by the pipeline. In the following schematic mapping the main steps performed by the pipeline are depicted. <br>
Briefly, first of all a quality control (fastQC) of the raw fastq data is performed. Then, bwa-mem is used to align the paired-end sequences, but separately, onto the genome of reference. The aligned reads are filtered and used to generate the Hi-C matrices at the lowest resolution using [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/index.html) ([J. Wolff *et. al*, Nuc. Acids Res. 2020](https://doi.org/10.1093/nar/gkaa220)). <br> If necessary, the matrices' bins are merge to generate higher resolution matrices. Further, if required by the user, matrices belonging to the same *group* are summed and processed in parallel to the "single-sample" ones. <br>
These matrices are then normalized and corrected and used to generate quality control plots as well preform sample correlation analyses. <br>
Finally, TAD and Loops detection can be performed by [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/index.html), while differential compartment analyses can be run using [dcHiC](https://www.nature.com/articles/s41467-022-34626-6) on both individual and grouped samples.

More details on [parameters](#Configuration-file) and structure/meaning of the [results](#Results) can be found in the next paragraphs.

</div>

![snHiC workflow](https://sebastian-gregoricchio.github.io/snHiC/resources/snHiC_workflow.png)


### Configuration file
<div style="text-align: justify">
The configuration file is a yaml-formatted file containing all the parameters that are passed to different steps of the pipelines such as the directory with the input files, reference genome, threads of the tools, etc. <br>
The snHiC configuration file is divided in two sections:
* *experiment-specific*, with al the parameters that most likely are changing from experiment to experiment;
* *common*, with parameters that are quite stable independently of the experiments design. The latter should be changed only for very specific needs. <br>
Hereafter, the meaning of the different parameters is described.
</div>

<br/><br/>


#### Experiment-specific section
*Required parameters*

| **Parameter**   |   **Description**   |
|------------:|:----------------|
|*runs_directory*| The full path to the directory were the input fastq files are contained, e.g. `/home/user/HiC/00_runs/`. Importantly, the name of the files, deprived of the read suffix (e.g., _R1/_R2) and file extension (e.g., .fastq.gz) will be used as sample name.|
|*output_directory*| The full path to the folder in which the results should be stored, e.g. `"/home/user/HiC/"`. |
|*fastq_extension*| Defaul: `.fastq.gz`. The extension of the input fastq files. This string will be removed from the input file names to obtain the sample names. Examples: `".fastq.gz"`, `".fq.gz"`, `".fasta"`, `".fq"`. |
|*runs_suffix*| Default: `["_R1", "_R2"]`. A list (python format) with the two reads suffixes corresponding to the read1 and read2 for the paired-end sequencing. |
|*genome_fasta*| The full path to a fasta-formatted file containing the sequence of there reference genome into which the reads will be aligned. If the index (.fai) file is not present in the same folder, it will be generated by the pipeline during its execution (this process may take long time to be finalized). The reference genomes can be downloaded, for instance, from the [UCSC golden path](https://hgdownload.soe.ucsc.edu/goldenPath/). |
|*genome_assembly_name*| A string with the name of the genome (e.g., `"hg19"`). This will be used to download/generate some files.|
|*matrix_resolution*| A numeric list (e.g., [10,20,40,50,100]) indicated all the matrices resolutions that should be generated expressed in kilo-base-pairs (kb). Importantly, each resolution must be a multiple of the lowest one, which must be greater or equal to 1kb. |
|*generate_bam*| Default: `False`. A logical value (True/False) indicating whether during the generation of the lowest resolution matrices a bam file with the filtered reads should be written. Notably, this process may be significantly time consuming. |
|*restriction_enzyme*| A string indicating the enzyme used to generate the HiC libraries. This parameter is used to identify the restriction sites on the genome. Possible choices are: 'DpnII', 'MboI', 'NlaIII', 'Csp6I', 'CviQI', 'HindIII', 'EcoRI', 'BamHI', 'BglII'. Other enzymes can be added in the pipeline by modifying the `restriction_table` at the beginning of the *snHiC.Snakefile*.|



*Optional parameters*

| **Parameter**   |   **Description**   |
|------------:|:----------------|
|*perform_grouped_analyses*| A logical value (True/False) indicating whether grouped analyses should be performed. This function requires the parameter `sample_metadata`.|
|*sample_metadata*| The path to a tab-delimited txt file containing two columns: *sample* and the name of the *group* at which it belongs. This sample names must be indicated without the read file suffix and extension (e.g *sample.A_R1.fastq.gz* --> *sample.A*). |
|*detect_loops*| A logical value (True/False) indicating whether loops detection should be performed by [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicDetectLoops.html). |
|*detect_compartments*| A logical value (True/False) indicating whether differential compartment analyses should be performed by [dcHiC](https://www.nature.com/articles/s41467-022-34626-6). This function requires the parameter `sample_metadata`. |
|*dcHiC_repository_folder*| "<target/folder>/dcHiC/" |
|*minimal_resolution_compartments*| A number indicating the minimal resolution needed to run the compartment analyses, in kilo-base-pairs (bp). `dcHiC` works up to 5kb resolution. |
|*character_subsitution_dashes_and_points_sample_name*| Default `""` (nothing). `dcHiC` does not allow any `.` or `-` in the sample names. Thus, this parameter is used to define which character should be used to substitute dots and dashes in the sample name before preforming the compartment analyses. |
|*chr_filtering_string*| To optimize `dcHiC` compartment computation, low coverage chromosomes should be excluded. Provide a string with the name of the chromosomes to exclude separated by a `|`. By default: `"hap|gl|random|NC|hs|Y|M|mt"`. |
|*dcHiC_chr_threads*| Default `4`. Number of threads to be used for parallel chromosome processing per sample by `dcHiC`. |
|*dcHiC_PCA_threads*| Default `4`. Number of threads to be used for PCA calculation per chromosome per sample by `dcHiC`. |
|*dcHiC_analyses_type*| Default: `"cis"`. A string to indicate whether the compartment analyses will be performed on `"cis"` or `"trans"` interaction matrices. |

<br/><br/>

#### Common section

| **Parameter**   |  **Description**   |
|------------:|:----------------|
|*fastQC_threads*| Default: `2`. Number of CPUs to use for [fastQC](https://github.com/s-andrews/FastQC) (fastq quality control). |
|*bwa_threads*| Default: `12`. Number of CPUs to use for the mapping performed by [bwa-mem](http://bio-bwa.sourceforge.net/bwa.shtml). |
|*mapQ_cutoff*| Default: `15`. All reads with a mapping quality (MAPQ) score lower than this value will be filtered out from the bam files. |
|*SAMtools_threads*| Default: `8`. Number of CPUs used by [samtools](http://www.htslib.org/doc/samtools.html) for bam indexing and filtering. |
|*heatmap_color*| Default: `'RdBu'`. A string indicating the color gradient pattern to use for the correlation heatmaps. This value is passed to matplotlib. Therefore, available options (see [matplotlib page](https://matplotlib.org/stable/tutorials/colors/colormaps.html) for examples) are the following: 'Accent', 'Blues', 'BrBG', 'BuGn', 'BuPu', 'CMRmap', 'Dark2', 'GnBu', 'Greens', 'Greys', 'OrRd', 'Oranges', 'PRGn', 'Paired', 'Pastel1', 'Pastel2', 'PiYG', 'PuBu', 'PuBuGn', 'PuOr', 'PuRd', 'Purples', 'RdBu', 'RdGy', 'RdPu', 'RdYlBu', 'RdYlGn', 'Reds', 'Set1', 'Set2', 'Set3', 'Spectral', 'Wistia', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'afmhot', 'autumn', 'binary', 'bone', 'brg', 'bwr', 'cividis', 'cool', 'coolwarm', 'copper', 'cubehelix', 'flag', 'gist_earth', 'gist_gray', 'gist_heat', 'gist_ncar', 'gist_rainbow', 'gist_stern', 'gist_yarg', 'gnuplot', 'gnuplot2', 'gray', 'hot', 'hsv', 'icefire', 'inferno', 'jet', 'magma', 'mako', 'nipy_spectral', 'ocean', 'pink', 'plasma', 'prism', 'rainbow', 'rocket', 'seismic', 'spring', 'summer', 'tab10', 'tab20', 'tab20b', 'tab20c', 'terrain', 'twilight', 'twilight_shifted', 'viridis', 'vlag', 'winter'. |
|*correlation_method*| Default: `'pearson'`. Possible choices: `'pearson'`, `'spearman'`. Method to use for the sample correlation.
|*normalization_method*| Default: `'smallest'`. Possible choices: `'norm_range'`, `'smallest'`, `'multiplicative'`. Method to use for the matrices normalization. Values passed to [`hicNormalize`](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicNormalize.html). |
|*correction_method*| Default: `"ICE"`. Possible choices: `'ICE'`, `'KR'`. Method to use for the matrices correction. Values passed to [`hicCorrectMatrix`](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicCorrectMatrix.html). |
|*hicFindTADs_threads*| Default: `10`. Number of threads to be used to detect the TADs by [hicFindTADs](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicFindTADs.html). |
|*extra_findTAD_parameters*| Default: `'--thresholdComparisons 0.01'`. A string containing any additional parameter, separated by a space (e.g., '--paramA X --paramB Y'), to pass to [hicFindTADs](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicFindTADs.html). |  
|*hicDetectLoops_threads*| Default: `10`. Number of threads to be used to detect the TADs by [hicDetectLoops](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicDetectLoops.html). |
|*maxLoopDistance*| Default: `2000000`bp (2Mb). Maximum loop distance, in base-pairs (bp), to be used to detect loops by [hicDetectLoops](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicDetectLoops.html). |
|*loop_windowSize*| Default: `10`kb. The window size for the neighborhood region the peak is located in. All values from this region (exclude the values from the peak region) are tested against the peak region for significant difference. For more info see [hicDetectLoops](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicDetectLoops.html). |
|*loop_peakWidth*| Default: `6`kb. The width of the peak region in bins. The square around the peak will include (2 * peakWidth)^2 bins. For more info see [hicDetectLoops](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicDetectLoops.html). |
|*loop_pValuePreselection*| Default: `0.05`. Only candidates with p-values less the given threshold will be considered as candidates. For each genomic distance a negative binomial distribution is fitted and for each pixel a p-value given by the cumulative density function is given. This does NOT influence the p-value for the neighborhood testing. For more info see [hicDetectLoops](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicDetectLoops.html). |
|*loop_pValue*| Default: `0.05`. Rejection level for Anderson-Darling or Wilcoxon-rank sum test for H0. H0 is peak region and background have the same distribution. For more info see [hicDetectLoops](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicDetectLoops.html). |



<br/><br/>

## Results
The structure of the *output_folder* is the following:

<pre>
<b><em>output_folder</em></b>
├── <b>01_fastQC_raw</b>
│   ├── <em>sample</em>_fastqc.html
│   ├── <em>sample</em>_fastqc.zip
│   └── <b>multiQC_raw</b>
│       ├── <b>multiQC_report_fastqRaw_data</b>
│       │   ├── multiqc_citations.txt
│       │   ├── multiqc_data.json
│       │   ├── multiqc_fastqc.txt
│       │   ├── multiqc_general_stats.txt
│       │   ├── multiqc.log
│       │   └── multiqc_sources.txt
│       └── multiQC_report_fastqRaw.html
|
├── <b>02_BAM</b>
│   ├── <em>sample</em>_mapQ20_sorted_woMT.bam
│   ├── <em>sample</em>_mapQ20_sorted_woMT.bam.bai
│   └── <b>flagstat</b>
│       ├── <em>sample</em>_flagstat_filtered_bam_woMT.txt
│       └── <em>sample</em>_flagstat_UNfiltered_bam.txt
|
├── <b>03_BAM_dedup</b> (or 03_BAM_mdup, if duplicates are not removed)
│   ├── <em>sample</em>_mapQ20_woMT_dedup_shifted_sorted.bam
│   ├── <em>sample</em>_mapQ20_woMT_dedup_shifted_sorted.bam.bai
│   ├── <b>fastQC</b>
│   │   ├── <em>sample</em>_sorted_woMT_dedup_fastqc.html
│   │   ├── <em>sample</em>_sorted_woMT_dedup_fastqc.zip
│   │   └── <b>multiQC_dedup_bams</b>
│   │       ├── <b>multiQC_report_BAMs_dedup_data</b>
│   │       │   ├── multiqc_citations.txt
│   │       │   ├── multiqc_data.json
│   │       │   ├── multiqc_fastqc.txt
│   │       │   ├── multiqc_general_stats.txt
│   │       │   ├── multiqc.log
│   │       │   └── multiqc_sources.txt
│   │       └── multiQC_report_BAMs_dedup.html
│   ├── <b>flagstat</b>
│   │   ├── <em>sample</em>_flagstat_filtered_bam_woMT_dedup.txt
│   │   └── <em>sample</em>_flagstat_woMT_dedup_shifted_sorted.txt
│   ├── <b>fragmentSizeDistribution_plots</b>
|   |   ├── ALL.samples_fragment_size_distribution.pdf
│   │   └── <em>sample</em>_fragment_size_distribution.pdf
│   ├── <b>metrics</b>
│   │   └── <em>sample</em>_metrics_woMT_dedup_bam.txt
│   └── <b>unshifted_bams</b>
│       ├── <em>sample</em>_mapQ20_sorted_woMT_dedup.bam
│       └── <em>sample</em>_mapQ20_sorted_woMT_dedup.bai
|
├── <b>04_Normalization</b>
│   ├── <b>HMCan_output</b> ### --> only if HMCan correction is performed ###
│   └── <b>normalized_bigWigs</b>
│       └── <em>sample</em>_mapQ20_woMT_dedup_shifted_normalized_bs5.bw
│
├── <b>05_Peaks_MACS3</b> ### --> if HMCan correction is not performed ###
│   ├── <em>sample</em>_mapQ20_woMT_dedup_shifted_FDR0.01_peaks.narrowPeak
│   ├── <em>sample</em>_mapQ20_woMT_dedup_shifted_FDR0.01_peaks.xls
│   ├── <em>sample</em>_mapQ20_woMT_dedup_shifted_FDR0.01_summits.bed
│   └── <b>log</b>
│       └── <em>sample</em>_mapQ20_woMT_dedup_shifted_FDR0.01.log
|
├── <b>06_Overall_quality_and_info</b>
|   ├── Lorenz_curve_deeptools.plotFingreprint_allSamples.pdf
|   ├── <b>Counts</b>
|   │   ├── counts_summary.txt
|   │   └── <b>subread_featureCounts_output</b>
|   │       └── <b>sample</b>
|   │           ├── <em>sample</em>.readCountInPeaks
|   │           ├── <em>sample</em>.readCountInPeaks.log
|   │           └── <em>sample</em>.readCountInPeaks.summary
|   └── <b>Sample_comparisons</b>
|       ├── multiBigWigSummary_matrix_allSamples.npz
|       ├── PCA_on_BigWigs_wholeGenome.pdf
|       ├── <b>Peak_comparison</b>
|       │   ├── all_samples_peaks_concatenation_collapsed_sorted.bed
|       │   ├── peaks_score_matrix_all_samples_MACS3.npz
|       │   └── peaks_score_matrix_all_samples_table_MACS3.tsv
|       |   └── <b>Heatmaps</b>
|       |       ├── Heatmap_on_log1p.rawScores_for_MACS3.peaks_union_population.pdf
|       │       └── Heatmap_on_zScores_for_MACS3.peaks_union_population.pdf
|       └── <b>Sample_correlation</b>
|           ├── Correlation_heatmap_on_BigWigs_wholeGenome_pearsonMethod.pdf
|           ├── Correlation_heatmap_on_BigWigs_wholeGenome_spearmanMethod.pdf
|           ├── Correlation_scatterplot_on_BigWigs_wholeGenome_pearsonMethod.pdf
|           └── Correlation_scatterplot_on_BigWigs_wholeGenome_spearmanMethod.pdf
|
└── <b>07_Variant_calling</b>
    ├── all_samples_peaks_concatenation_collapsed_sorted.bed
    ├── all.samples_dedup_gatk-indel_filtered.DP20.QUAL20
    ├── all.samples_dedup_gatk-snp_filtered.DP20.QUAL20
    ├── all.samples_INDEL_counts_plot.pdf
    ├── all.samples_SNP_counts_plot.pdf
    └── <b><em>sample</em></b>
        ├── <em>sample</em>_dedup_gatk-indel_filtered.DP20.QUAL20.txt
        ├── <em>sample</em>_dedup_gatk-indel_filtered.DP20.QUAL20.vcf.gz
        ├── <em>sample</em>_dedup_gatk-indel_filtered.DP20.QUAL20.vcf.gz.tbi
        ├── <em>sample</em>_dedup_gatk-snp_filtered.DP20.QUAL20.txt
        ├── <em>sample</em>_dedup_gatk-snp_filtered.DP20.QUAL20.vcf.gz
        ├── <em>sample</em>_dedup_gatk-snp_filtered.DP20.QUAL20.vcf.gz.tbi
        ├── <em>sample</em>_dedup_gatk.vcf.gz
        ├── <em>sample</em>_dedup_gatk.vcf.gz.tbi
        ├── <em>sample</em>_mapQ20_sorted_woMT_dedup_bsqr.bai
        ├── <em>sample</em>_mapQ20_sorted_woMT_dedup_bsqr.bam
        ├── <em>sample</em>_mapQ20_sorted_woMT_dedup_bsqr.table
        └── <em>sample</em>_plotCoverage.pdf
</pre>

<br/><br/>

### 01_fastQC_raw
This folder contains a the fastq quality control (fastQC) reports for each fastq file and a summary report of multiQC.

<br/><br/>

### 02_BAM
When the reads are aligned onto the reference genome by bwa, the resulting SAM files are filtered for mapping quality (MAPQ) and the mithocondrial (suffix: woMT) reads are removed before sorting. Flagstat metrics is generated for each file and stored in the homonym folder.

<br/><br/>

### 03_BAM_dedup / 03_BAM_mdup
PICARD is used to remove (suffix: dedup) or mark (suffix: mdup) duplicates in the BAM files. The resulting BAMs are stored in the subfolder "unshifted_bams", while the PICARD metrics is stored in the "metrics" folder. A fastq quality control (fastQC) and relative multiQC report is performed on the unshifted bams. <br>
Then, the Tn5 nick reparation bias is corrected by shifting of the reads using [deeptools alignmentSieve](https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html) (suffix: shifted). Notice that, after shifting, the read sequence information is lost in the shifted BAMs. <br>
Flagstat metrics is generated for each unshifted and shifted bam file and stored in the "falgstat" folder.

Furthermore, in the "fragmentSizeDistribution_plots" folder the distribution of the fragment sizes for each sample (shifted BAMs) and a file collecting all the plots in a single file. Here after an example of a good (left) and a bad (right) fragment size distribution.

![fragment size distribution examples](https://sebastian-gregoricchio.github.io/snakeATAC/resources/fragmentSize_distribution_examples.svg)

An optimal fragment size distribution should be included within a range of 50-800bp, with a periodicity of ~150bp (corrsponding to mono-, di-, tri-, ... nucleosomes) with a lower intensity for larger fragments.

<br/><br/>


### 04_Normalization
Shifted signal is normalized on sequencing depth library upon copy number variation correction by [HMCan](https://academic.oup.com/bioinformatics/article/29/23/2979/246862?login=false) (if requested by the user). The bin size used is indicated in the resulting bigWig file name (suffix: bs#). <br>
However, these bigWig files can be normalized more precisely normalized in the case that you dispone of a corresponding RNA-seq data set using [CHIPIN](https://doi.org/10.1186/s12859-021-04320-3) (L. Polit *et.al*, BMC Bioinformatics, 2021). Examples of CHIPIN usage can be found at [S. Gregoricchio *et al.*, Nucleic Acids Research (2022)](https://doi.org/10.1093/nar/gkac613).

<br/><br/>

### 05_Peaks_MACS3 (when HMCan correction is not performed)
Peaks and summits (if required by the user) are called by MACS3 on shifted BAMs. The FDR (False Discovery Ratio) threshold used is indicated in the file name (suffix: FDR#). When HMCan correction is active, the peaks are called by HMCan itself.

<br/><br/>

### 06_Overall_quality_and_info
This folder contains multiple quality controls, feature counts and sample correlation plots:

*  `Lorenz_curve_deeptools.plotFingreprint_allSamples.pdf` is a plot showing the enrichment of the signal allover the genome. Indeed, if a sample does not show any enrichment the reads will equally distributed over the genome resulting in a diagonal line in the plot (left panel). When instead the signal is specific for the feature sought (e.g., open chromatin) it will be enriched only at specific location and the curve will be closer to the bottom-right corner of the plot (right panel).

![lorenz curve examples](https://sebastian-gregoricchio.github.io/snakeATAC/resources/lorenz_curve_examples.svg)

<br/><br/>

* `Counts`: contains the results of featureCounts (from subread) with the counts of reads and other statistics on called peaks for each sample. It is availble also tab-separated file containing a summary of the main features counts for each sample: <br><br>
**Summary counts table description**

| **Column**   |   **Description**   |
|------------:|:----------------|
| *Sample* | Sample name |
| *Reads_R1* | Number of reads in read.1 fastq file. |
| *Reads_R2* | Number of reads in read.2 fastq file. |
| *Reads_total* | Total number of reads (read.1 + read.2). |
| *unfiltered_BAM* | Total number of reads in the bam file after filtering by map quality (MAPQ). |
| *Percentage_MT* | Approximative percentage of reads represented by the mithocondrial DNA. Ideally lower than 10-20%. |
| *dedup_BAM* | Total number of reads left after BAM reads deduplication. |
| *duplicated_reads* | Number of duplicated reads. If the duplicates are not remove the value will be 0. |
| *shifted_BAM* | Number of reads in the shifted BAMs. |
| *loss_post_shifting* | Number of reads lost upon BAM shifting. Consider that reads falling in blacklisted regions are removed. |
| *n.peaks* | Total number of peaks called. |
| *FRiP.perc* | Frequency Reads in Peaks percentage, corresponds to the number of reads falling in peak regions divide by the total number of reads and multiplied by 100. |
| *FRiP.quality* | A label ("good" or "bad") to indicate whether the FRiP score is good or not for a given sample. The threshold can be changed in the config file by the user, by the default 20 (as suggested by the [ENCODE guidelines](https://www.encodeproject.org/atac-seq/)). |

<br/><br/>

* `Sample_comparisons`: the plots in this folder help the user to understand the variability of the samples.
  + `multiBigWigSummary_matrix_allSamples.npz`: result of deeptools multiBigWigSummary used to plot the PCA and correlation plots;
  + `PCA_on_BigWigs_wholeGenome.pdf`: Principal Component Analyses results of the signal allover the genome;
  + `Peak_comparison`:
    - `all_samples_peaks_concatenation_collapsed_sorted.bed`: the peaks called in all samples are merged and collapsed in this bed file;
    - `peaks_score_matrix_all_samples_MACS3.npz`: a matrix containing the average score at each peak (previous bed file) for each samples is generated;
    - `peaks_score_matrix_all_samples_table_MACS3.tsv`: same matrix as before, but in tab-separated format.
    - `Heatmaps`: the matrix generated on all peaks is used to cluster the samples and two heatmaps are plotted: one on the log1p of the raw scores, and one on the z-score (on rows)
  + `Sample_correlation`: scatter and heatmap correlation plots are generated based on the signal over the whole genome. Both Pearson and Spearman methods are used.

<br/><br/>




### 07_Variant_calling
If required by the user, the pipeline can call altered single-nucleotide polymorphism (SNP) and insertion/deletions (InDel). [PICARD CreateSequenceDictionary](https://gatk.broadinstitute.org/hc/en-us/articles/360037068312-CreateSequenceDictionary-Picard-) is used to create a genome dictionary in order to perform a Base Quality Score Recalibration ([BQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-)) of unshifted BAM files (after filling of the ReadGroups by [PICARD AddOrReplaceReadGroups](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-)) by [GATK BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator) and [GATK ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR). <br>
The recalibrated BAMs are used by [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) to generate a [GVCF](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format) (Genomic Variant Call Format) file containing the genomic variants individuated at the regions included in the file resulting by the merge of all the called peaks. This GVCF file is then recalibrated by [GATK GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037435831-GenotypeGVCFs) resulting in a VCF file. <br>
The VCF recalibrated file is selected by [GATK SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants) to obtain two separate VCF files corresponding to SNPs and InDels. <br>
Ultimately, these VCF are hard-filtered for sequencing depth (DP) and quality (QUAL) by [SnpSift Filter](https://pcingola.github.io/SnpEff/ss_filter/) and then exported in a txt table by [SnpSift Extract Fields](https://pcingola.github.io/SnpEff/ss_extractfields/). Notably, the genotype '0|0' (both alleles not mutated) is filtered out from the .txt table. <br>
The SNP, or InDel, .txt tables from all samples are merged in a unique one with the addition of a column corresponding to the sample name. These tables are used to generate a plot of the counts of variants found in each sample.
Further, a coverage plot at the merged peaks is generated by [deeptools plotCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/plotCoverage.html).


<br/><br/>

-----------------
## Package history and releases
A list of all releases and respective description of changes applied could be found [here](https://sebastian-gregoricchio.github.io/snakeATAC/NEWS).

## Contact
For any suggestion, bug fixing, commentary please report it in the [issues](https://github.com/sebastian-gregoricchio/snakeATAC/issues)/[request](https://github.com/sebastian-gregoricchio/snakeATAC/pulls) tab of this repository.

## License
This repository is under a [GNU General Public License (version 3)](https://sebastian-gregoricchio.github.io/Rseb/LICENSE.md/LICENSE).

<br/>

#### Contributors
[![contributors](https://badges.pufler.dev/contributors/sebastian-gregoricchio/Rseb?size=50&padding=5&bots=true)](https://sebastian-gregoricchio.github.io/)
