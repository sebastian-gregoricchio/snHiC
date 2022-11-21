[![Snakemake](https://img.shields.io/badge/snakemake-≥7.8.5-brightgreen.svg)](https://snakemake.github.io)
![release](https://img.shields.io/github/v/release/sebastian-gregoricchio/snHiC)
[![license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://sebastian-gregoricchio.github.io/snHiC/LICENSE.md/LICENSE)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/snHiC?style=social)](https://github.com/sebastian-gregoricchio/snHiC/fork)
<!-- ![update](https://badges.pufler.dev/updated/sebastian-gregoricchio/snHiC)
![visits](https://badges.pufler.dev/visits/sebastian-gregoricchio/snHiC)
![downloads](https://img.shields.io/github/downloads/sebastian-gregoricchio/snHiC/total.svg)--->



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

<br/><br/>

## How to run the pipeline
The basic `snHiC` pipeline requires only two files:
* `snHiC.snakefile`, containing all the rules that will be run;
* `configuration.yaml` file, in which the user can define and customize all the parameters for the different pipeline steps. <br>

Then, if grouped or differential compartment analyses are required by the user, a third file is required. This file will contain the sample names (without the read file suffix and extension, e.g *sample.A_R1.fastq.gz* --> *sample.A*) and the group to which they belong. Here follows an example (tab-delimited txt file):

| *sample* | *group* |
|:---------|:--------|
| Sample_A | Tumor   |
| Sample_B | Normal  |
| Sample_C | Normal  |
| Sample_D | Tumor   |
| Sample_E | Normal  |
| Sample_F | Tumor   |

To partially avoid unexpected errors during the execution of the pipeline, a so called 'dry-run' is strongly recommended. Indeed, adding a `-n` at the end of the snakemake running command will allow snakemake to check that all links and file/parameters dependencies are satisfied before to run the "real" processes. This command will therefore help the debugging process.
```shell
snakemake \
-s </target/folder>/snHiC/workflow/snHiC.snakefile \
--configfile </target/folder>/snHiC/config/snHiC_config.yaml \
--cores 12 \
-n
```

If no errors occur, the pipeline can be run with the same command line without the terminal `-n`:
```shell
snakemake \
-s </target/folder>/snHiC/workflow/snHiC.snakefile \
--configfile </target/folder>/snHiC/config/snHiC_config.yaml \
--cores 12
```

Notice that the absence of errors does not mean that the pipeline will run without any issues; the "dry-run" is only checking whether all the resources are available. <br>
Furthermore, there is the possibility to run the pipeline only partially. An example of usage could be if someone wants to have a look to the fast quality controls (fastQC and multiQC reports) before to perform the alignment. To do that, it is sufficient to run a dry-run (`-n` mode), then pick the name of the rule at which you want the pipeline to stop and lastly type the following command:
```shell
snakemake \
-s </target/folder>/snHiC/workflow/snHiC.snakefile \
--configfile </target/folder>/snHiC/config/snHiC_config.yaml \
--cores 12 \
--until <stop_rule_name>
```

<br/><br/>

### Snakefile
The snakefile are contained all the condition checks and rules (processing steps) that will be performed by the pipeline. In the following schematic mapping the main steps performed by the pipeline are depicted. <br>
Briefly, first of all a quality control (fastQC) of the raw fastq data is performed. Then, bwa-mem is used to align the paired-end sequences, but separately, onto the genome of reference. The aligned reads are filtered and used to generate the Hi-C matrices at the lowest resolution using [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/index.html) ([J. Wolff *et. al*, Nuc. Acids Res. 2020](https://doi.org/10.1093/nar/gkaa220)). <br> If necessary, the matrices' bins are merge to generate higher resolution matrices. Further, if required by the user, matrices belonging to the same *group* are summed and processed in parallel to the "single-sample" ones. <br>
These matrices are then normalized and corrected and used to generate quality control plots as well preform sample correlation analyses. <br>
Finally, TAD and Loops detection can be performed by [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/index.html), while differential compartment analyses can be run using [dcHiC](https://www.nature.com/articles/s41467-022-34626-6) on both individual and grouped samples.

More details on [parameters](#Configuration-file) and structure/meaning of the [results](#Results) can be found in the next paragraphs.

![snHiC workflow](https://sebastian-gregoricchio.github.io/snHiC/resources/snHiC_workflow.png)


All the possible rule names are listed below.
The example is based on 4 samples grouped in two groups (2 samples per group) for which two resolutions are generated:
```shell
---------------------------------------------------------------------  -------
Rule name                                                                Count
---------------------------------------------------------------------  -------
AAA_initialization                                                           1
A_fastQC_raw                                                                 8
B_multiQC_raw                                                                1
C_bwa_align                                                                  8
D_generate_restriction_file_and_get_chrSizes                                 1
E1_interaction_matrix_generation_at_smallest_resolution                      4
E2_multiQC_report_for_HiC_matrices                                           1
E3_merging_interaction_matrix_bins_for_all_resolutions                       4
F1_matrices_normalization                                                    2
F2_samples_correlation                                                       1
G1_matrices_correction__diagnosticPlot_and_MAD                               8
G2_matrices_correction__getting_threshold_values                             8
G3_matrices_correction__correction                                           8
H1_matrices_format_conversion__cool                                          1
H2_matrices_format_conversion__hicpro                                        1
I_call_TADs                                                                  8
J_plotting_intraChr_distances                                                2
L1_sum_matrices_by_group                                                     1
L2_merging_grouped_interaction_matrix_bins_for_all_resolutions               2
M_grouped_matrices_normalization                                             2
N1_summed_matrices_correction__diagnosticPlot_and_MAD                        4
N2_summed_matrices_correction__getting_threshold_values                      4
N3_summed_matrices_correction__correction                                    4
N4_summed_matrices_correction__cool_conversion                               1
N5_summed_matrices_correction__hicpro_conversion                             1
O_call_TADs_on_summed_matrices                                               4
P_detect_loops_singleSamples                                                 8
Q_detect_loops_groupedSamples                                                4
R1_detect_compartments_dcHiC_singleSamples__inputFile_all_vs_all             2
R2_detect_compartments_dcHiC_singleSamples__call_compartments                2
R3_detect_compartments_dcHiC_singleSamples__bedGraphToBigWig                 2
R4_detect_compartments_dcHiC_singleSamples__call_compartments_combos         2
R5_detect_compartments_dcHiC_singleSamples__bedGraphToBigWig_combos          2
S1_detect_compartments_dcHiC_groupedSamples___inputFile_all_vs_all           2
S2_detect_compartments_dcHiC_groupedSamples__call_compartments               2
S3_detect_compartments_dcHiC_groupedSamples__bedGraphToBigWig                2
S4_detect_compartments_dcHiC_groupedSamples__call_compartments_combos        2
S5_detect_compartments_dcHiC_groupedSamples__bedGraphToBigWig_combos         2
---------------------------------------------------------------------  -------
TOTAL                                                                      122
---------------------------------------------------------------------  -------
```


### Configuration file
The configuration file is a yaml-formatted file containing all the parameters that are passed to different steps of the pipelines such as the directory with the input files, reference genome, threads of the tools, etc. <br>
The snHiC configuration file is divided in two sections:
* *experiment-specific*, with al the parameters that most likely are changing from experiment to experiment;
* *common*, with parameters that are quite stable independently of the experiments design. The latter should be changed only for very specific needs. <br>
Hereafter, the meaning of the different parameters is described.

<br/><br/>


#### Experiment-specific section
*Required parameters*

| **Parameter**   |   **Description**   |
|------------:|:----------------|
|*runs_directory*| The full path to the directory were the input fastq files are contained, e.g. `/home/user/HiC/00_runs/`. Importantly, the name of the files, deprived of the read suffix (e.g., _R1/_R2) and file extension (e.g., .fastq.gz) will be used as sample name.|
|*output_directory*| The full path to the folder in which the results should be stored, e.g. `"/home/user/HiC_results/"`. If not already existing, it will be generated automatically. |
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
|*dcHiC_repository_folder*| Example `'<target/folder>/dcHiC/'`. A string with the full path to the directory of the `dcHiC` cloned repository. |
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
For an experiment including 4 samples (sampleA-B-C-D) assigned to two groups (Normal and Tumor) and for which 2 resolutions will be analyzed (500kb and 1000kb), the structure of the *output_folder* will be the following:

################## REMEMBER TO PUT THE OUTPUT BAM DIRECTORY ##########

<pre>
  <b><em>output_folder</em></b>
  <b>results_manual</b>
  ├── <b>01_fastQC_raw</b>
  │   ├── <b>multiQC_raw</b>
  │   │   ├── <b>multiQC_report_fastqRaw_data</b>
  │   │   │   ├── multiqc_citations.txt
  │   │   │   ├── multiqc_data.json
  │   │   │   ├── multiqc_fastqc.txt
  │   │   │   ├── multiqc_general_stats.txt
  │   │   │   ├── multiqc.log
  │   │   │   └── multiqc_sources.txt
  │   │   └── multiQC_report_fastqRaw.html
  │   ├── <em>sampleA</em>_R1_fastqc.html
  │   ├── <em>sampleA</em>_R1_fastqc.zip
  │   ├── <em>sampleA</em>_R2_fastqc.html
  │   ├── <em>sampleA</em>_R2_fastqc.zip
  │   └── ...
  |
  ├── <b>02_Alignements
  │   ├── <b>log</b>
  │   │   ├── <em>sampleA</em>_R1_bwa-mem.err
  │   │   ├── <em>sampleA</em>_R1_bwa-mem.out
  │   │   ├── <em>sampleA</em>_R2_bwa-mem.err
  │   │   ├── <em>sampleA</em>_R2_bwa-mem.out
  │   │   └── ...
  │   ├── <em>sampleA</em>_R1.bam
  │   ├── <em>sampleA</em>_R2.bam
  │   └── ...
  |
  ├── <b>03_BAM</b> / <b>03_BAM__not_generated</b>
  |   ├── <b>falgstat</b>
  │   │   ├── <em>sampleA</em>_flagstat_bam.txt
  │   │   ├── <em>sampleB</em>_flagstat_bam.txt
  │   |   └── ...
  |   ├── <em>sampleA</em>_mapQ15_sorted.bam
  |   ├── <em>sampleA</em>_mapQ15_sorted.bai
  │   └── ...
  |
  ├── <b>04_Interaction_matrices</b>
  │   ├── <b>log</b>
  │   │   ├── <em>sampleA</em>.<em>100kb</em>.hicBuildMatrix.err
  │   │   ├── <em>sampleA</em>.<em>100kb</em>.hicBuildMatrix.log
  │   │   ├── <em>sampleA</em>.<em>50kb</em>.hicBuildMatrix.err
  │   │   ├── <em>sampleA</em>.<em>50kb</em>.hicBuildMatrix.log
  │   │   └── ...
  │   ├── <b>QC_matrices</b>
  │   │   ├── <b>ALL_SAMPLES</b>
  │   │   │   ├── <b>HiQC_matrices</b>
  │   │   │   │   ├── discarded_table.txt
  │   │   │   │   ├── distance.png
  │   │   │   │   ├── distance_table.txt
  │   │   │   │   ├── hicQC.html
  │   │   │   │   ├── pairs_discarded.png
  │   │   │   │   ├── pairs_sequenced.png
  │   │   │   │   ├── QC_table.txt
  │   │   │   │   ├── read_orientation.png
  │   │   │   │   ├── read_orientation_table.txt
  │   │   │   │   ├── unmapable_table.txt
  │   │   │   │   └── unmappable_and_non_unique.png
  │   │   │   └── <b>multiQC_matrices</b>
  │   │   │       ├── <b>multiqc_data</b>
  │   │   │       │   ├── multiqc_citations.txt
  │   │   │       │   ├── multiqc_data.json
  │   │   │       │   ├── multiqc_general_stats.txt
  │   │   │       │   ├── multiqc_hicexplorer.txt
  │   │   │       │   ├── multiqc.log
  │   │   │       │   └── multiqc_sources.txt
  │   │   │       └── multiqc_report.html
  │   │   ├── <b><em>sampleA</em></b>
  │   │   │   ├── discarded_table.txt
  │   │   │   ├── distance.png
  │   │   │   ├── distance_table.txt
  │   │   │   ├── hicQC.html
  │   │   │   ├── pairs_discarded.png
  │   │   │   ├── pairs_sequenced.png
  │   │   │   ├── QC.log
  │   │   │   ├── QC_table.txt
  │   │   │   ├── read_orientation.png
  │   │   │   ├── read_orientation_table.txt
  │   │   │   ├── unmapable_table.txt
  │   │   │   └── unmappable_and_non_unique.png
  │   │   └── <b><em>sample...</em></b>
  │   │       └── ...
  │   ├── <em>sampleA</em>_mapQ15_<em>100kb</em>.h5
  │   ├── <em>sampleA</em>_mapQ15_<em>50kb</em>.cool
  │   ├── <em>sampleA</em>_mapQ15_<em>50kb</em>.h5
  │   └── ...
  !
  ├── <b>05_Interaction_matrices_<em>Normal</em>ized</b>
  │   ├── <em>sampleA</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized.h5
  │   ├── <em>sampleA</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized.h5
  │   ├── ...
  │   └── <b>sample_correlation</b>
  │       ├── heatmap_correlation.pdf
  │       └── scatter_correlation.pdf
  ├── <b>06_Interaction_matrices_<em>Normal</em>ized_and_corrected</b>
  │   ├── <b>corrected_matrices</b>
  │   │   ├── <b><em>sampleA</em></b>
  │   │   │   ├── <b>cool_format</b>
  │   │   │   │   ├── <em>sampleA</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized_corrected.cool
  │   │   │   │   └── <em>sampleA</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized_corrected.cool
  │   │   │   ├── <b>h5_format</b>
  │   │   │   │   ├── <em>sampleA</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized_corrected.h5
  │   │   │   │   └── <em>sampleA</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized_corrected.h5
  │   │   │   └── <b>hicpro_format</b>
  │   │   │       ├── <em>sampleA</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized_corrected.hicpro
  │   │   │       ├── <em>sampleA</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized_corrected_hicpro.bed
  │   │   │       ├── <em>sampleA</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized_corrected.hicpro
  │   │   │       └── <em>sampleA</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized_corrected_hicpro.bed
  │   │   └── <b><em>sample...</em></b>
  │   │       └── ...
  │   ├── <b>diagnostic_plots</b>
  │   │   ├── <em>sampleA</em>_<em>100kb</em>_<em>Normal</em>ized_diagnosticPlot.png
  │   │   ├── <em>sampleA</em>_<em>50kb</em>_<em>Normal</em>ized_diagnosticPlot.png
  │   │   └── ...
  │   └── <b>median_absolute_deviation</b>
  │       ├── <em>sampleA</em>_<em>100kb</em>_<em>Normal</em>ized_MedianAbsoluteDeviation.mad
  │       ├── <em>sampleA</em>_<em>50kb</em>_<em>Normal</em>ized_MedianAbsoluteDeviation.mad
  │       ├── ...
  │       └── <b>thresholds</b>
  │           ├── <em>sampleA</em>_<em>100kb</em>_<em>Normal</em>ized_thresholdValues.txt
  │           ├── <em>sampleA</em>_<em>50kb</em>_<em>Normal</em>ized_thresholdValues.txt
  │           └── ...
  |
  ├── <b>07_TADs_calling_HiCexplorer</b>
  │   ├── <b><em>sampleA</em></b>
  │   │   ├── <b><em>100kb</em>_resolution</b>
  │   │   │   ├── <b>log</b>
  │   │   │   │   ├── <em>sampleA</em>_mapQ15_<em>100kb</em>_TAD.calling.err
  │   │   │   │   └── <em>sampleA</em>_mapQ15_<em>100kb</em>_TAD.calling.out
  │   │   │   ├── <em>sampleA</em>_mapQ15_<em>100kb</em>_boundaries.bed
  │   │   │   ├── <em>sampleA</em>_mapQ15_<em>100kb</em>_boundaries.gff
  │   │   │   ├── <em>sampleA</em>_mapQ15_<em>100kb</em>_domains.bed
  │   │   │   ├── <em>sampleA</em>_mapQ15_<em>100kb</em>_score.bedgraph
  │   │   │   ├── <em>sampleA</em>_mapQ15_<em>100kb</em>_tad_score.bm
  │   │   │   └── <em>sampleA</em>_mapQ15_<em>100kb</em>_zscore_matrix.h5
  │   │   └── <b><em>50kb</em>_resolution</b>
  │   │       ├── <b>log</b>
  │   │       │   ├── <em>sampleA</em>_mapQ15_<em>50kb</em>_TAD.calling.err
  │   │       │   └── <em>sampleA</em>_mapQ15_<em>50kb</em>_TAD.calling.out
  │   │       ├── <em>sampleA</em>_mapQ15_<em>50kb</em>_boundaries.bed
  │   │       ├── <em>sampleA</em>_mapQ15_<em>50kb</em>_boundaries.gff
  │   │       ├── <em>sampleA</em>_mapQ15_<em>50kb</em>_domains.bed
  │   │       ├── <em>sampleA</em>_mapQ15_<em>50kb</em>_score.bedgraph
  │   │       ├── <em>sampleA</em>_mapQ15_<em>50kb</em>_tad_score.bm
  │   │       └── <em>sampleA</em>_mapQ15_<em>50kb</em>_zscore_matrix.h5
  │   ├── <b><em>sample...</em></b>
  │       └── ...
  |
  ├── <b>08_Interaction_distances</b>
  │   ├── intraChr_distances_all_samples_<em>100kb</em>_resolution_distancesValues.tsv
  │   ├── intraChr_distances_all_samples_<em>100kb</em>_resolution.pdf
  │   ├── intraChr_distances_all_samples_<em>50kb</em>_resolution_distancesValues.tsv
  │   ├── intraChr_distances_all_samples_<em>50kb</em>_resolution.pdf
  │   └── <b>log</b>
  │       ├── intraChr_distances_all_samples_<em>100kb</em>_resolution.err
  │       ├── intraChr_distances_all_samples_<em>100kb</em>_resolution.out
  │       ├── intraChr_distances_all_samples_<em>50kb</em>_resolution.err
  │       └── intraChr_distances_all_samples_<em>50kb</em>_resolution.out
  |
  ├── <b>09_Loop_detection_HiCexplorer</b>
  │   ├── <b><em>sampleA</em></b>
  │   │   ├── <b>log</b>
  │   │   │   ├── <em>sampleA</em>_mapQ15_<em>100kb</em>_loop.detectiong.err
  │   │   │   ├── <em>sampleA</em>_mapQ15_<em>100kb</em>_loop.detection.out
  │   │   │   ├── <em>sampleA</em>_mapQ15_<em>50kb</em>_loop.detectiong.err
  │   │   │   └── <em>sampleA</em>_mapQ15_<em>50kb</em>_loop.detection.out
  │   │   ├── <em>sampleA</em>_mapQ15_<em>100kb</em>_loops.bedpe
  │   │   └── <em>sampleA</em>_mapQ15_<em>50kb</em>_loops.bedpe
  │   └── <b><em>sample...</em></b>
  │       └── ...
  |
  ├── <b>10_Compartments_detection_dcHiC</b>
  │   ├── <b><em>100kb</em>_resolution</b>
  │   │   ├── dcHiC_input_file_individual_samples_<em>100kb</em>.txt
  │   │   ├── <b>DifferentialResult</b>
  │   │   │   └── <b>all_vs_all</b>
  │   │   │       ├── <b><em>Normal</em>_data</b>
  │   │   │       ├── <b>pcOri</b>
  │   │   │       ├── <b>pcQnm</b>
  │   │   │       ├── <b><em>Tumor</em>_data</b>
  │   │   │       └── <b>viz</b>
  │   │   │           └── <b>files</b>
  │   │   ├── <b>filtered_hicpro_beds</b>
  │   │   │   ├── <em>sampleA</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized_corrected_hicpro_FILTERED.bed
  │   │   │   └── ...
  │   │   ├── <b>hg19_1000000_goldenpathData</b>
  │   │   │   ├── cytoBand.txt.gz
  │   │   │   ├── hg19.binned.bed
  │   │   │   ├── hg19.chrom.sizes
  │   │   │   ├── hg19.fa
  │   │   │   ├── hg19.fa.fai
  │   │   │   ├── hg19.fa.gz
  │   │   │   ├── hg19.GCpt.bedGraph
  │   │   │   ├── hg19.GCpt.tss.bedGraph
  │   │   │   ├── hg19.refGene.gtf.gz
  │   │   │   └── hg19.tss.bed
  │   │   ├── <b><em>sampleA</em>_mapQ15_<em>100kb</em>_pca</b>
  │   │   │   └── <b>intra_pca</b>
  │   │   │       └── <b><em>sampleA</em>_mapQ15_<em>100kb</em>_mat</b>
  │   │   │           ├── chr1.bed
  │   │   │           ├── chr1.cmat.txt
  │   │   │           ├── chr1.distparam
  │   │   │           ├── chr1.pc.txt
  │   │   │           ├── chr1.precmat.txt
  │   │   │           ├── chr1.svd.rds
  │   │   │           ├── chr1.txt
  │   │   │           └── ...
  │   │   ├── <b><em>sample...</em>_mapQ15_<em>100kb</em>_pca</b>
  │   │   │   └── ...
  │   └── <b><em>50kb</em>_resolution</b>
  │       ├── dcHiC_input_file_individual_samples_<em>50kb</em>_<em>Normal</em>_vs_<em>Tumor</em>.txt
  │       ├── dcHiC_input_file_individual_samples_<em>50kb</em>.txt
  │       ├── <b>DifferentialResult</b>
  │       │   ├── <b>all_vs_all</b>
  │       │   │   ├── <b>fdr_result</b>
  │       │   │   │   ├── differential.intra_sample_chr1_combined.pcQnm.bedGraph
  │       │   │   │   ├── ...
  │       │   │   │   ├── differential.intra_sample_combined.Filtered.pcQnm.bedGraph
  │       │   │   │   ├── differential.intra_sample_combined.pcQnm.bedGraph
  │       │   │   │   ├── differential.intra_sample_group.Filtered.pcOri.bedGraph
  │       │   │   │   ├── differential.intra_sample_group.Filtered.pcQnm.bedGraph
  │       │   │   │   ├── differential.intra_sample_group.pcOri.bedGraph
  │       │   │   │   └── differential.intra_sample_group.pcQnm.bedGraph
  │       │   │   ├── <b><em>Normal</em>_data</b>
  │       │   │   │   ├── intra_chr1_combined.pcOri.bedGraph
  │       │   │   │   ├── intra_chr1_combined.pcQnm.bedGraph
  │       │   │   │   ├── ...
  │       │   │   │   ├── <em>sampleB</em>_mapQ15_<em>50kb</em>_intra_chr1.pc.bedGraph
  │       │   │   │   ├── ...
  │       │   │   │   ├── <em>sampleC</em>_mapQ15_<em>50kb</em>_intra_chr12.pc.bedGraph
  │       │   │   │   └── ...
  │       │   │   ├── <b>pcOri</b>
  │       │   │   │   ├── intra_sample_chr1_combined.pcOri.bedGraph
  │       │   │   │   └── ...
  │       │   │   ├── <b>pcQnm</b>
  │       │   │   │   ├── intra_sample_chr1_combined.pcQnm.bedGraph
  │       │   │   │   └── ...
  │       │   │   ├── <b><em>Tumor</em>_data</b>
  │       │   │   │   ├── intra_chr1_combined.pcOri.bedGraph
  │       │   │   │   └── ...
  │       │   │   │   ├── <em>sampleA</em>_mapQ15_<em>50kb</em>_intra_chr1.pc.bedGraph
  │       │   │   │   ├── <em>sample...</em>_mapQ15_<em>50kb</em>_intra_chr....pc.bedGraph
  │       │   │   │   └── ...
  │       │   │   └── <b>viz</b>
  │       │   │       ├── <b>files</b>
  │       │   │       │   ├── intra_compartment.bedGraph
  │       │   │       │   ├── intra_<em>sampleA</em>_mapQ15_<em>50kb</em>_PC_A.compartment.bedGraph
  │       │   │       │   ├── intra_<em>sampleA</em>_mapQ15_<em>50kb</em>_PC_B.compartment.bedGraph
  │       │   │       │   ├── intra_<em>sampleA</em>_mapQ15_<em>50kb</em>_PC.bedGraph
  │       │   │       │   └── ...
  │       │   │       ├── <b>files_bigWig</b>
  │       │   │       │   ├── intra_<em>sampleA</em>_mapQ15_<em>50kb</em>_PC.bw
  │       │   │       │   ├── intra_<em>sampleB</em>_mapQ15_<em>50kb</em>_PC.bw
  │       │   │       │   ├── intra_<em>sampleC</em>_mapQ15_<em>50kb</em>_PC.bw
  │       │   │       │   └── intra_<em>sampleD</em>_mapQ15_<em>50kb</em>_PC.bw
  │       │   │       ├── <b>files_compartment_beds</b>
  │       │   │       │   ├── intra_<em>sampleA</em>_mapQ15_<em>50kb</em>_PC_compartments_sorted.bed
  │       │   │       │   ├── intra_<em>sampleB</em>_mapQ15_<em>50kb</em>_PC_compartments_sorted.bed
  │       │   │       │   ├── intra_<em>sampleC</em>_mapQ15_<em>50kb</em>_PC_compartments_sorted.bed
  │       │   │       │   └── intra_<em>sampleD</em>_mapQ15_<em>50kb</em>_PC_compartments_sorted.bed
  │       │   │       └── <b>vizIGV_intra</b>
  │       │   │           ├── <b>data</b>
  │       │   │           │   ├── differential_compartment.log10Padj.bedGraph.gz
  │       │   │           │   ├── differential_compartment.Mahalanobis.bedGraph
  │       │   │           │   ├── differential_compartment.Mahalanobis.bedGraph.gz
  │       │   │           │   ├── <em>Normal</em>.PC.bedGraph.gz
  │       │   │           │   └── <em>Tumor</em>.PC.bedGraph.gz
  │       │   │           ├── <b>data_bigWig</b>
  │       │   │           │   ├── differential_compartment.log10Padj.bw
  │       │   │           │   ├── differential_compartment.Mahalanobis.bw
  │       │   │           │   ├── <em>Normal</em>.PC.bw
  │       │   │           │   └── <em>Tumor</em>.PC.bw
  │       │   │           ├── <b>data_compartment_beds</b>
  │       │   │           │   ├── <em>Normal</em>.PC_compartments_sorted.bed
  │       │   │           │   └── <em>Tumor</em>.PC_compartments_sorted.bed
  │       │   │           └── intra_igv_pcQnm.html
  │       │   └── <b><em>Normal</em>_vs_<em>Tumor</em></b>
  │       │       ├── <b>fdr_result</b>
  │       │       │   ├── differential.intra_sample_chr1_combined.pcQnm.bedGraph
  │       │       │   └── ...
  │       │       │   ├── differential.intra_sample_combined.Filtered.pcQnm.bedGraph
  │       │       │   ├── differential.intra_sample_combined.pcQnm.bedGraph
  │       │       │   ├── differential.intra_sample_group.Filtered.pcOri.bedGraph
  │       │       │   ├── differential.intra_sample_group.Filtered.pcQnm.bedGraph
  │       │       │   ├── differential.intra_sample_group.pcOri.bedGraph
  │       │       │   └── differential.intra_sample_group.pcQnm.bedGraph
  │       │       ├── <b><em>Normal</em>_data</b>
  │       │       │   ├── intra_chr1_combined.pcOri.bedGraph
  │       │       │   ├── intra_chr1_combined.pcQnm.bedGraph
  │       │       │   └── ...
  │       │       │   ├── <em>sample...</em>_mapQ15_<em>50kb</em>_intra_chr...pc.bedGraph
  │       │       │   └── ...
  │       │       ├── <b>pcOri</b>
  │       │       │   ├── intra_sample_chr1_combined.pcOri.bedGraph
  │       │       │   └── ...
  │       │       ├── <b>pcQnm</b>
  │       │       │   ├── intra_sample_chr1_combined.pcQnm.bedGraph
  │       │       │   └── ...
  │       │       ├── <b><em>Tumor</em>_data</b>
  │       │       │   ├── intra_chr1_combined.pcOri.bedGraph
  │       │       │   ├── intra_chr1_combined.pcQnm.bedGraph
  │       │       │   └── ...
  │       │       │   ├── <em>sampleA</em>_mapQ15_<em>50kb</em>_intra_chr1.pc.bedGraph
  │       │       │   ├── <em>sample...</em>_mapQ15_<em>50kb</em>_intra_chr....pc.bedGraph
  │       │       │   └── ...
  │       │       └── <b>viz</b>
  │       │           ├── <b>files</b>
  │       │           │   ├── intra_compartment.bedGraph
  │       │           │   ├── intra_<em>sampleA</em>_mapQ15_<em>50kb</em>_PC_A.compartment.bedGraph
  │       │           │   ├── intra_<em>sampleA</em>_mapQ15_<em>50kb</em>_PC_B.compartment.bedGraph
  │       │           │   ├── intra_<em>sampleA</em>_mapQ15_<em>50kb</em>_PC.bedGraph
  │       │           │   └── ...
  │       │           ├── <b>files_bigWig</b>
  │       │           │   ├── intra_<em>sampleA</em>_mapQ15_<em>50kb</em>_PC.bw
  │       │           │   └── ...
  │       │           ├── <b>files_compartment_beds</b>
  │       │           │   ├── intra_<em>sampleA</em>_mapQ15_<em>50kb</em>_PC_compartments_sorted.bed
  │       │           │   └── ...
  │       │           └── <b>vizIGV_intra</b>
  │       │               ├── <b>data</b>
  │       │               │   ├── differential_compartment.log10Padj.bedGraph.gz
  │       │               │   ├── differential_compartment.Mahalanobis.bedGraph
  │       │               │   ├── differential_compartment.Mahalanobis.bedGraph.gz
  │       │               │   ├── <em>Normal</em>.PC.bedGraph.gz
  │       │               │   └── <em>Tumor</em>.PC.bedGraph.gz
  │       │               ├── <b>data_bigWig</b>
  │       │               │   ├── differential_compartment.log10Padj.bw
  │       │               │   ├── differential_compartment.Mahalanobis.bw
  │       │               │   ├── <em>Normal</em>.PC.bw
  │       │               │   └── <em>Tumor</em>.PC.bw
  │       │               ├── <b>data_compartment_beds</b>
  │       │               │   ├── <em>Normal</em>.PC_compartments_sorted.bed
  │       │               │   └── <em>Tumor</em>.PC_compartments_sorted.bed
  │       │               └── intra_igv_pcQnm.html
  │       ├── <b>filtered_hicpro_beds</b>
  │       │   ├── <em>sampleA</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized_corrected_hicpro_FILTERED.bed
  │       │   └── ...
  │       ├── <b>hg19_500000_goldenpathData</b>
  │       │   ├── cytoBand.txt.gz
  │       │   ├── hg19.binned.bed
  │       │   ├── hg19.chrom.sizes
  │       │   ├── hg19.fa
  │       │   ├── hg19.fa.fai
  │       │   ├── hg19.fa.gz
  │       │   ├── hg19.GCpt.bedGraph
  │       │   ├── hg19.GCpt.tss.bedGraph
  │       │   ├── hg19.refGene.gtf.gz
  │       │   └── hg19.tss.bed
  │       ├── <em>Normal</em>_chr_pc_selected.txt
  │       ├── <em>Normal</em>_clus.txt
  │       ├── <em>Normal</em>_cor.txt
  │       ├── <em>Normal</em>_vals.txt
  │       ├── <b><em>sampleA</em>_mapQ15_<em>50kb</em>_pca</b>
  │       │   └── <b>intra_pca</b>
  │       │       └── <b><em>sampleA</em>_mapQ15_<em>50kb</em>_mat</b>
  │       │           ├── chr1.bed
  │       │           ├── chr1.cmat.txt
  │       │           ├── chr1.distparam
  │       │           ├── chr1.PC1.bedGraph
  │       │           ├── chr1.PC2.bedGraph
  │       │           ├── chr1.pc.bedGraph
  │       │           ├── chr1.pc.txt
  │       │           ├── chr1.precmat.txt
  │       │           ├── chr1.svd.rds
  │       │           ├── chr1.txt
  │       │           └── ...
  │       ├── <b><em>sample...</em>
  │       │   └── ...
  │       ├── <em>Tumor</em>_chr_pc_selected.txt
  │       ├── <em>Tumor</em>_clus.txt
  │       ├── <em>Tumor</em>_cor.txt
  │       └── <em>Tumor</em>_vals.txt
  |
  | ********************************************************************************************************
  |
  └── <b>11_Grouped_analyses</b>
      ├── <b>A_summed_matrices</b>
      │   ├── <b>log</b>
      │   │   ├── <em>Normal</em>.<em>100kb</em>.hicBuildMatrix.err
      │   │   ├── <em>Normal</em>.<em>100kb</em>.hicBuildMatrix.log
      │   │   ├── <em>Tumor</em>.<em>100kb</em>.hicBuildMatrix.err
      │   │   └── <em>Tumor</em>.<em>100kb</em>.hicBuildMatrix.log
      │   ├── <em>Normal</em>_mapQ15_<em>100kb</em>.h5
      │   ├── <em>Normal</em>_mapQ15_<em>50kb</em>.h5
      │   ├── <em>Tumor</em>_mapQ15_<em>100kb</em>.h5
      │   └── <em>Tumor</em>_mapQ15_<em>50kb</em>.h5
      |
      ├── <b>B_summed_matrices_<em>Normal</em>ized</b>
      │   ├── <em>Normal</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized.h5
      │   ├── <em>Normal</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized.h5
      │   ├── <em>Tumor</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized.h5
      │   └── <em>Tumor</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized.h5
      |
      ├── <b>C_summed_matrices_<em>Normal</em>ized_and_corrected</b>
      │   ├── <b>corrected_matrices</b>
      │   │   ├── <b><em>Normal</em></b>
      │   │   │   ├── <b>cool_format</b>
      │   │   │   │   ├── <em>Normal</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized_corrected.cool
      │   │   │   │   └── <em>Normal</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized_corrected.cool
      │   │   │   ├── <b>h5_format</b>
      │   │   │   │   ├── <em>Normal</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized_corrected.h5
      │   │   │   │   └── <em>Normal</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized_corrected.h5
      │   │   │   └── <b>hicpro_format</b>
      │   │   │       ├── <em>Normal</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized_corrected.hicpro
      │   │   │       ├── <em>Normal</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized_corrected_hicpro.bed
      │   │   │       ├── <em>Normal</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized_corrected.hicpro
      │   │   │       └── <em>Normal</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized_corrected_hicpro.bed
      │   │   └── <b><em>Tumor</em></b>
      │   │       ├── <b>cool_format</b>
      │   │       │   ├── <em>Tumor</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized_corrected.cool
      │   │       │   └── <em>Tumor</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized_corrected.cool
      │   │       ├── <b>h5_format</b>
      │   │       │   ├── <em>Tumor</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized_corrected.h5
      │   │       │   └── <em>Tumor</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized_corrected.h5
      │   │       └── <b>hicpro_format</b>
      │   │           ├── <em>Tumor</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized_corrected.hicpro
      │   │           ├── <em>Tumor</em>_mapQ15_<em>100kb</em>_<em>Normal</em>ized_corrected_hicpro.bed
      │   │           ├── <em>Tumor</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized_corrected.hicpro
      │   │           └── <em>Tumor</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized_corrected_hicpro.bed
      │   ├── <b>diagnostic_plots</b>
      │   │   ├── <em>Normal</em>_<em>100kb</em>_<em>Normal</em>ized_diagnosticPlot.png
      │   │   ├── <em>Normal</em>_<em>50kb</em>_<em>Normal</em>ized_diagnosticPlot.png
      │   │   ├── <em>Tumor</em>_<em>100kb</em>_<em>Normal</em>ized_diagnosticPlot.png
      │   │   └── <em>Tumor</em>_<em>50kb</em>_<em>Normal</em>ized_diagnosticPlot.png
      │   └── <b>median_absolute_deviation</b>
      │       ├── <em>Normal</em>_<em>100kb</em>_<em>Normal</em>ized_MedianAbsoluteDeviation.mad
      │       ├── <em>Normal</em>_<em>50kb</em>_<em>Normal</em>ized_MedianAbsoluteDeviation.mad
      │       ├── <b>thresholds</b>
      │       │   ├── <em>Normal</em>_<em>100kb</em>_<em>Normal</em>ized_thresholdValues.txt
      │       │   ├── <em>Normal</em>_<em>50kb</em>_<em>Normal</em>ized_thresholdValues.txt
      │       │   ├── <em>Tumor</em>_<em>100kb</em>_<em>Normal</em>ized_thresholdValues.txt
      │       │   └── <em>Tumor</em>_<em>50kb</em>_<em>Normal</em>ized_thresholdValues.txt
      │       ├── <em>Tumor</em>_<em>100kb</em>_<em>Normal</em>ized_MedianAbsoluteDeviation.mad
      │       └── <em>Tumor</em>_<em>50kb</em>_<em>Normal</em>ized_MedianAbsoluteDeviation.mad
      |
      ├── <b>D_TADs_calling_HiCexplorer</b>
      │   ├── <b><em>Normal</em></b>
      │   │   ├── <b><em>50kb</em>_resolution</b>
      │   │   |   ├── <b>log</b>
      │   │   |   │   ├── <em>Normal</em>_mapQ15_<em>50kb</em>_TAD.calling.err
      │   │   |   │   └── <em>Normal</em>_mapQ15_<em>50kb</em>_TAD.calling.out
      │   │   |   ├── <em>Normal</em>_mapQ15_<em>50kb</em>_boundaries.bed
      │   │   |   ├── <em>Normal</em>_mapQ15_<em>50kb</em>_boundaries.gff
      │   │   |   ├── <em>Normal</em>_mapQ15_<em>50kb</em>_domains.bed
      │   │   |   ├── <em>Normal</em>_mapQ15_<em>50kb</em>_score.bedgraph
      │   │   |   ├── <em>Normal</em>_mapQ15_<em>50kb</em>_tad_score.bm
      │   │   |   └── <em>Normal</em>_mapQ15_<em>50kb</em>_zscore_matrix.h5
      |   |   └── <b><em>100kb</em>_resolution</b>
      |   |       └── ...
      │   └── <b><em>Tumor</em></b>
      │       └── ...
      ├── <b>E_Loop_detection_HiCexplorer</b>
      │   ├── <b><em>Normal</em></b>
      │   │   ├── <b>log</b>
      │   │   │   ├── <em>Normal</em>_mapQ15_<em>100kb</em>_loop.detection.err
      │   │   │   ├── <em>Normal</em>_mapQ15_<em>100kb</em>_loop.detection.out
      │   │   │   ├── <em>Normal</em>_mapQ15_<em>50kb</em>_loop.detection.err
      │   │   │   └── <em>Normal</em>_mapQ15_<em>50kb</em>_loop.detection.out
      │   │   ├── <em>Normal</em>_mapQ15_<em>100kb</em>_loops.bedpe
      │   │   └── <em>Normal</em>_mapQ15_<em>50kb</em>_loops.bedpe
      │   └── <b><em>Tumor</em></b>
      │       └── ...
      |
      └── <b>F_Compartments_detection_dcHiC</b>
          ├── <b><em>50kb</em>_resolution</b>
          |   ├── dcHiC_input_file_grouped_samples_<em>50kb</em>.txt
          |   ├── dcHiC_input_file_individual_samples_<em>50kb</em>_<em>Normal</em>_vs_<em>Tumor</em>.txt
          |   ├── <b>DifferentialResult</b>
          |   │   ├── <b>all_vs_all</b>
          |   │   │   ├── <b>fdr_result</b>
          |   │   │   │   ├── differential.intra_sample_chr1_combined.pcQnm.bedGraph
          |   │   │   │   ├── ...
          |   │   │   │   ├── differential.intra_sample_combined.Filtered.pcQnm.bedGraph
          |   │   │   │   ├── differential.intra_sample_combined.pcQnm.bedGraph
          |   │   │   │   ├── differential.intra_sample_group.Filtered.pcOri.bedGraph
          |   │   │   │   ├── differential.intra_sample_group.Filtered.pcQnm.bedGraph
          |   │   │   │   ├── differential.intra_sample_group.pcOri.bedGraph
          |   │   │   │   └── differential.intra_sample_group.pcQnm.bedGraph
          |   │   │   ├── <b><em>Normal</em>_data</b>
          |   │   │   │   ├── intra_chr1_combined.pcOri.bedGraph
          |   │   │   │   ├── intra_chr1_combined.pcQnm.bedGraph
          |   │   │   │   └── ...
          |   │   │   │   ├── <em>Normal</em>_mapQ15_<em>50kb</em>_intra_chr1.pc.bedGraph
          |   │   │   │   └── ...
          |   │   │   ├── <b>pcOri</b>
          |   │   │   │   ├── intra_sample_chr1_combined.pcOri.bedGraph
          |   │   │   │   └── ...
          |   │   │   ├── <b>pcQnm</b>
          |   │   │   │   ├── intra_sample_chr1_combined.pcQnm.bedGraph
          |   │   │   │   └── ...
          |   │   │   ├── <b><em>Tumor</em>_data</b>
          |   │   │   │   ├── intra_chr1_combined.pcOri.bedGraph
          |   │   │   │   └── ...
          |   │   │   │   ├── <em>Tumor</em>_mapQ15_<em>50kb</em>_intra_chr11.pc.bedGraph
          |   │   │   │   └── ...
          |   │   │   └── <b>viz</b>
          |   │   │       ├── <b>files</b>
          |   │   │       │   ├── intra_compartment.bedGraph
          |   │   │       │   ├── intra_<em>Normal</em>_mapQ15_<em>50kb</em>_PC_A.compartment.bedGraph
          |   │   │       │   ├── intra_<em>Normal</em>_mapQ15_<em>50kb</em>_PC_B.compartment.bedGraph
          |   │   │       │   ├── intra_<em>Normal</em>_mapQ15_<em>50kb</em>_PC.bedGraph
          |   │   │       │   ├── intra_<em>Tumor</em>_mapQ15_<em>50kb</em>_PC_A.compartment.bedGraph
          |   │   │       │   ├── intra_<em>Tumor</em>_mapQ15_<em>50kb</em>_PC_B.compartment.bedGraph
          |   │   │       │   └── intra_<em>Tumor</em>_mapQ15_<em>50kb</em>_PC.bedGraph
          |   │   │       ├── <b>files_bigWig</b>
          |   │   │       │   ├── intra_<em>Normal</em>_mapQ15_<em>50kb</em>_PC.bw
          |   │   │       │   └── intra_<em>Tumor</em>_mapQ15_<em>50kb</em>_PC.bw
          |   │   │       ├── <b>files_compartment_beds</b>
          |   │   │       │   ├── intra_<em>Normal</em>_mapQ15_<em>50kb</em>_PC_compartments_sorted.bed
          |   │   │       │   └── intra_<em>Tumor</em>_mapQ15_<em>50kb</em>_PC_compartments_sorted.bed
          |   │   │       └── <b>vizIGV_intra</b>
          |   │   │           ├── <b>data</b>
          |   │   │           │   ├── differential_compartment.log10Padj.bedGraph.gz
          |   │   │           │   ├── differential_compartment.Mahalanobis.bedGraph
          |   │   │           │   ├── differential_compartment.Mahalanobis.bedGraph.gz
          |   │   │           │   ├── <em>Normal</em>.PC.bedGraph.gz
          |   │   │           │   └── <em>Tumor</em>.PC.bedGraph.gz
          |   │   │           ├── <b>data_bigWig</b>
          |   │   │           │   ├── differential_compartment.log10Padj.bw
          |   │   │           │   ├── differential_compartment.Mahalanobis.bw
          |   │   │           │   ├── <em>Normal</em>.PC.bw
          |   │   │           │   └── <em>Tumor</em>.PC.bw
          |   │   │           ├── <b>data_compartment_beds</b>
          |   │   │           │   ├── <em>Normal</em>.PC_compartments_sorted.bed
          |   │   │           │   └── <em>Tumor</em>.PC_compartments_sorted.bed
          |   │   │           └── intra_igv_pcQnm.html
          |   │   └── <b><em>Normal</em>_vs_<em>Tumor</em></b>
          |   │       ├── <b>fdr_result</b>
          |   │       │   ├── differential.intra_sample_chr1_combined.pcQnm.bedGraph
          |   │       │   ├── ...
          |   │       │   ├── differential.intra_sample_combined.Filtered.pcQnm.bedGraph
          |   │       │   ├── differential.intra_sample_combined.pcQnm.bedGraph
          |   │       │   ├── differential.intra_sample_group.Filtered.pcOri.bedGraph
          |   │       │   ├── differential.intra_sample_group.Filtered.pcQnm.bedGraph
          |   │       │   ├── differential.intra_sample_group.pcOri.bedGraph
          |   │       │   └── differential.intra_sample_group.pcQnm.bedGraph
          |   │       ├── <b><em>Normal</em>_data</b>
          |   │       │   ├── intra_chr1_combined.pcOri.bedGraph
          |   │       │   ├── intra_chr1_combined.pcQnm.bedGraph
          |   │       │   └── ...
          |   │       │   ├── <em>Normal</em>_mapQ15_<em>50kb</em>_intra_chr1.pc.bedGraph
          |   │       │   └── ...
          |   │       ├── <b>pcOri</b>
          |   │       │   ├── intra_sample_chr1_combined.pcOri.bedGraph
          |   │       │   └── ...
          |   │       ├── <b>pcQnm</b>
          |   │       │   ├── intra_sample_chr1_combined.pcQnm.bedGraph
          |   │       │   └── ...
          |   │       ├── <b><em>Tumor</em>_data</b>
          |   │       │   ├── intra_chr1_combined.pcOri.bedGraph
          |   │       │   └── ...
          |   │       └── <b>viz</b>
          |   │           ├── <b>files</b>
          |   │           │   ├── intra_compartment.bedGraph
          |   │           │   ├── intra_<em>Normal</em>_mapQ15_<em>50kb</em>_PC_A.compartment.bedGraph
          |   │           │   ├── intra_<em>Normal</em>_mapQ15_<em>50kb</em>_PC_B.compartment.bedGraph
          |   │           │   ├── intra_<em>Normal</em>_mapQ15_<em>50kb</em>_PC.bedGraph
          |   │           │   ├── intra_<em>Tumor</em>_mapQ15_<em>50kb</em>_PC_A.compartment.bedGraph
          |   │           │   ├── intra_<em>Tumor</em>_mapQ15_<em>50kb</em>_PC_B.compartment.bedGraph
          |   │           │   └── intra_<em>Tumor</em>_mapQ15_<em>50kb</em>_PC.bedGraph
          |   │           ├── <b>files_bigWig</b>
          |   │           │   ├── intra_<em>Normal</em>_mapQ15_<em>50kb</em>_PC.bw
          |   │           │   └── intra_<em>Tumor</em>_mapQ15_<em>50kb</em>_PC.bw
          |   │           ├── <b>files_compartment_beds</b>
          |   │           │   ├── intra_<em>Normal</em>_mapQ15_<em>50kb</em>_PC_compartments_sorted.bed
          |   │           │   └── intra_<em>Tumor</em>_mapQ15_<em>50kb</em>_PC_compartments_sorted.bed
          |   │           └── <b>vizIGV_intra</b>
          |   │               ├── <b>data</b>
          |   │               │   ├── differential_compartment.log10Padj.bedGraph.gz
          |   │               │   ├── differential_compartment.Mahalanobis.bedGraph
          |   │               │   ├── differential_compartment.Mahalanobis.bedGraph.gz
          |   │               │   ├── <em>Normal</em>.PC.bedGraph.gz
          |   │               │   └── <em>Tumor</em>.PC.bedGraph.gz
          |   │               ├── <b>data_bigWig</b>
          |   │               │   ├── differential_compartment.log10Padj.bw
          |   │               │   ├── differential_compartment.Mahalanobis.bw
          |   │               │   ├── <em>Normal</em>.PC.bw
          |   │               │   └── <em>Tumor</em>.PC.bw
          |   │               ├── <b>data_compartment_beds</b>
          |   │               │   ├── <em>Normal</em>.PC_compartments_sorted.bed
          |   │               │   └── <em>Tumor</em>.PC_compartments_sorted.bed
          |   │               └── intra_igv_pcQnm.html
          |   ├── <b>filtered_hicpro_beds</b>
          |   │   ├── <em>Normal</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized_corrected_hicpro_FILTERED.bed
          |   │   └── <em>Tumor</em>_mapQ15_<em>50kb</em>_<em>Normal</em>ized_corrected_hicpro_FILTERED.bed
          |   ├── <b>hg19_500000_goldenpathData</b>
          |   │   ├── cytoBand.txt.gz
          |   │   ├── hg19.binned.bed
          |   │   ├── hg19.chrom.sizes
          |   │   ├── hg19.fa
          |   │   ├── hg19.fa.fai
          |   │   ├── hg19.fa.gz
          |   │   ├── hg19.GCpt.bedGraph
          |   │   ├── hg19.GCpt.tss.bedGraph
          |   │   ├── hg19.refGene.gtf.gz
          |   │   └── hg19.tss.bed
          |   ├── <em>Normal</em>_chr_pc_selected.txt
          |   ├── <em>Normal</em>_clus.txt
          |   ├── <em>Normal</em>_cor.txt
          |   ├── <b><em>Normal</em>_mapQ15_<em>50kb</em>_pca</b>
          |   │   └── <b>intra_pca</b>
          |   │       └── <b><em>Normal</em>_mapQ15_<em>50kb</em>_mat</b>
          |   │           ├── chr1.bed
          |   │           ├── chr1.cmat.txt
          |   │           ├── chr1.distparam
          |   │           ├── chr1.PC1.bedGraph
          |   │           ├── chr1.PC2.bedGraph
          |   │           ├── chr1.pc.bedGraph
          |   │           ├── chr1.pc.txt
          |   │           ├── chr1.precmat.txt
          |   │           ├── chr1.svd.rds
          |   │           ├── chr1.txt
          |   │           └── ...
          |   ├── <em>Normal</em>_vals.txt
          |   ├── <em>Tumor</em>_chr_pc_selected.txt
          |   ├── <em>Tumor</em>_clus.txt
          |   ├── <em>Tumor</em>_cor.txt
          |   ├── <b><em>Tumor</em>_mapQ15_<em>50kb</em>_pca</b>
          |   │   └── <b>intra_pca</b>
          |   │       └── <b><em>Tumor</em>_mapQ15_<em>50kb</em>_mat</b>
          |   │           ├── chr1.bed
          |   │           ├── chr1.cmat.txt
          |   │           ├── chr1.distparam
          |   │           ├── chr1.PC1.bedGraph
          |   │           ├── chr1.PC2.bedGraph
          |   │           ├── chr1.pc.bedGraph
          |   │           ├── chr1.pc.txt
          |   │           ├── chr1.precmat.txt
          |   │           ├── chr1.svd.rds
          |   │           ├── chr1.txt
          |   │           └── ...
          |   └── <em>Tumor</em>_vals.txt
          └── <b><em>50kb</em>_resolution</b>
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
