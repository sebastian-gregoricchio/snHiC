[![Snakemake](https://img.shields.io/badge/snakemake-≥7.8.5-brightgreen.svg)](https://snakemake.github.io)
![release](https://img.shields.io/github/v/release/sebastian-gregoricchio/snHiC)
[![license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://sebastian-gregoricchio.github.io/snHiC/LICENSE.md/LICENSE)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/snHiC?style=social)](https://github.com/sebastian-gregoricchio/snHiC/fork)
<!-- ![update](https://badges.pufler.dev/updated/sebastian-gregoricchio/snHiC)
![visits](https://badges.pufler.dev/visits/sebastian-gregoricchio/snHiC)
![downloads](https://img.shields.io/github/downloads/sebastian-gregoricchio/snHiC/total.svg)--->

<div text-align: justify">

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
<div text-align: justify">
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
<div text-align: justify">
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

<div text-align: justify">
To partially avoid unexpected errors during the execution of the pipeline, a so called 'dry-run' is strongly recommended. Indeed, adding a `-n` at the end of the snakemake running command will allow snakemake to check that all links and file/parameters dependencies are satisfied before to run the "real" processes. This command will therefore help the debugging process.
</div>

```shell
snakemake \
-s </target/folder>/snHiC/workflow/snHiC.snakefile \
--configfile </target/folder>/snHiC/config/snHiC_config.yaml \
--cores 12 \
-n
```
<div text-align: justify">
If no errors occur, the pipeline can be run with the same command line without the terminal `-n`:
</div>
```shell
snakemake \
-s </target/folder>/snHiC/workflow/snHiC.snakefile \
--configfile </target/folder>/snHiC/config/snHiC_config.yaml \
--cores 12
```
<div text-align: justify">
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
<div text-align: justify">
The snakefile are contained all the condition checks and rules (processing steps) that will be performed by the pipeline. In the following schematic mapping the main steps performed by the pipeline are depicted. <br>
Briefly, first of all a quality control (fastQC) of the raw fastq data is performed. Then, bwa-mem is used to align the paired-end sequences, but separately, onto the genome of reference. The aligned reads are filtered and used to generate the Hi-C matrices at the lowest resolution using [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/index.html) ([J. Wolff *et. al*, Nuc. Acids Res. 2020](https://doi.org/10.1093/nar/gkaa220)). <br> If necessary, the matrices' bins are merge to generate higher resolution matrices. Further, if required by the user, matrices belonging to the same *group* are summed and processed in parallel to the "single-sample" ones. <br>
These matrices are then normalized and corrected and used to generate quality control plots as well preform sample correlation analyses. <br>
Finally, TAD and Loops detection can be performed by [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/index.html), while differential compartment analyses can be run using [dcHiC](https://www.nature.com/articles/s41467-022-34626-6) on both individual and grouped samples.

More details on [parameters](#Configuration-file) and structure/meaning of the [results](#Results) can be found in the next paragraphs.

</div>

![snHiC workflow](https://sebastian-gregoricchio.github.io/snHiC/resources/snHiC_workflow.png)


### Configuration file
<div text-align: justify">
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
