#########################
## Snakefile for snHiC ##
#########################

from typing import List
import pathlib
import re
import numpy
import os
import pandas as pd
import math
from itertools import combinations

### Function to handle the values for the wilcards
def constraint_to(values: List[str]) -> str:
    """
    From a list, return a regular expression allowing each
    value and not other.
    ex: ["a", "b", "v"] -> (a|b|v)
    """
    if isinstance(values, str):
            raise ValueError("constraint_to(): Expected a list, got str instead")
    return "({})".format("|".join(values))

# variable with any possibile message
messages = []

### working diirectory
home_dir = os.path.join(config["output_directory"],"")
shell('mkdir -p {home_dir}')
workdir: home_dir

### get the unique samples names and other variables
FILENAMES = next(os.walk(config["runs_directory"]))[2]
RUNNAMES = [re.sub(rf"{config['fastq_extension']}$", "", i) for i in FILENAMES]
SAMPLENAMES = numpy.unique([re.sub(rf"{config['runs_suffix'][0]}|{config['runs_suffix'][1]}.*$", "", i) for i in RUNNAMES])

# Matrices resolutions
# Detect whether the resolutions are multiple, and if not check that the variable is a list
if (type(config["matrix_resolution"]) == list):
    if (len(config["matrix_resolution"])>1):
        requested_resolutions = config["matrix_resolution"]
        requested_resolutions.sort()
    else:
        requested_resolutions = config["matrix_resolution"]
else:
    requested_resolutions = [config["matrix_resolution"]]

# Calculate the matrices to generate depending on the bin merging possibility (resolution depends on the lowest matrix resolution)
if len(requested_resolutions)>1:
    BINMERGE = ((numpy.array(requested_resolutions) / requested_resolutions[0]).astype(int))
    BINMERGE = BINMERGE[BINMERGE > 1]
    merging_resolutions = numpy.unique((numpy.array(BINMERGE) * requested_resolutions[0]).astype(int))
    NEW_RESOLUTIONS = numpy.concatenate(([requested_resolutions[0]], merging_resolutions))
    NEW_RESOLUTIONS = NEW_RESOLUTIONS.tolist()
    bins_table = pd.DataFrame({'binsToMerge': BINMERGE, 'newResolution': merging_resolutions})
else:
    BINMERGE = 1
    merging_resolutions = requested_resolutions
    NEW_RESOLUTIONS = requested_resolutions
    bins_table = pd.DataFrame({'binsToMerge': BINMERGE, 'newResolution': merging_resolutions})


# Define whether to performe compartmentalization analyses or not depending on the resolutions available
if (eval(str(config["compartments"]["detect_compartments"])) == True):
    if (NEW_RESOLUTIONS[0] >= config["compartments"]["minimal_resolution_compartments"]):
        perform_compartment_analyses = True
        MIN_RESOLUTIONS = NEW_RESOLUTIONS[0]
    else:
        perform_compartment_analyses = False
        comp_message = "printf '\033[0;95mThe minimum resolution required by the user is lower than the minimal resolution defined for the compartment analyses ("+str(config["compartments"]["minimal_resolution_compartments"])+"kb): compartment analyses will be skipped!\\n\033'"
        shell(comp_message)
        messages.append(comp_message)
        MIN_RESOLUTIONS = [0]
else:
    perform_compartment_analyses = False
    MIN_RESOLUTIONS = [0]


# Grouped analyses (if required): definition of the outputs to pass to iniztialization rule
if ((eval(str(config["groups"]["perform_grouped_analyses"])) == True) | (perform_compartment_analyses == True)):
    # read metadata files
    metadata_table = pd.read_csv(str(config["groups"]["sample_metadata"]),  sep='\t+', engine='python')
    groups = list(numpy.unique(list(metadata_table.iloc[:,1])))
    samples_by_group_list = []

    for i in groups:
        samples_by_group_list.append(list(numpy.unique(list(metadata_table[metadata_table.iloc[:,1]==i].iloc[:,0]))))
    if (config["TAD_caller"].lower() == "hicexplorer"):
        grouped_analyses_inputs = expand(os.path.join("12_Grouped_analyses/D_TADs_calling_HiCexplorer/{group}/{merged_res}kb_resolution/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_domains.bed"])), group = groups, merged_res = NEW_RESOLUTIONS)
    else:
        grouped_analyses_inputs = expand(os.path.join("12_Grouped_analyses/D_TADs_calling_GENOVA/{group}/{merged_res}kb_resolution/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_domains.bed"])), group = groups, merged_res = NEW_RESOLUTIONS)
else:
    # just use the normal matrices as inputs
    groups = "no_groups",
    grouped_analyses_inputs = expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/h5_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{resolution}kb_normalized_corrected.h5"])), sample = SAMPLENAMES, resolution=[str(x) for x in NEW_RESOLUTIONS])


# TAD calling
if (config["TAD_caller"].lower() == "hicexplorer"):
    domains_bed_single = expand(os.path.join("07_TADs_calling_HiCexplorer/{sample}/{resolution}kb_resolution/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{resolution}kb_domains.bed"])), sample = SAMPLENAMES, resolution=[str(x) for x in NEW_RESOLUTIONS])
else:
    domains_bed_single = expand(os.path.join("07_TADs_calling_GENOVA/{sample}/{resolution}kb_resolution/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{resolution}kb_domains.bed"])), sample = SAMPLENAMES, resolution=[str(x) for x in NEW_RESOLUTIONS])



# Loop detection (if required): definition of the outputs to pass to iniztialization rule
if (eval(str(config["loops"]["detect_loops"])) == True):
    # Define resolutions and understand wheter loop calling can be performed
    MAX_LOOPS_RESOLUTIONS = [x for x in NEW_RESOLUTIONS if x <= int(config["loops"]["max_resolution_loops"])]
    if (len(MAX_LOOPS_RESOLUTIONS) > 0):
        call_loops = True
    else:
        call_loops = False
else:
     call_loops = False

if (call_loops == True):
    if (config["loops"]["loop_caller"].lower() == "hicexplorer"):
        loops_single = expand(os.path.join("09_Loop_detection_HiCexplorer/{sample}/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{resolution}kb_loops.bedpe"])), sample = SAMPLENAMES, resolution=[str(x) for x in MAX_LOOPS_RESOLUTIONS])
        if (eval(str(config["groups"]["perform_grouped_analyses"])) == True):
            loops_group = expand(os.path.join("12_Grouped_analyses/E_Loop_detection_HiCexplorer/{group}/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_loops.bedpe"])), group = groups, merged_res = MAX_LOOPS_RESOLUTIONS)
        else:
            loops_group = []
    else:
        loops_single = expand(os.path.join("09_Loop_detection_mustache/{sample}/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{resolution}kb_loops.bedpe"])), sample = SAMPLENAMES, resolution=[str(x) for x in MAX_LOOPS_RESOLUTIONS])
        if (eval(str(config["groups"]["perform_grouped_analyses"])) == True):
            loops_group = expand(os.path.join("12_Grouped_analyses/E_Loop_detection_mustache/{group}/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_loops.bedpe"])), group = groups, merged_res = MAX_LOOPS_RESOLUTIONS)
        else:
            loops_group = []
else:
    loops_single = []
    loops_group = []



# Compartments detection (if required): definition of the outputs to pass to iniztialization rule
if ((perform_compartment_analyses == True) | (eval(str(config["perform_differential_contacts_analyses"])) == True)):
    # define compartmentalization analyses specific variables

    group_combinations = pd.DataFrame(combinations(groups, 2))
    comboA = [re.sub("\.|-", str(config["character_subsitution_dashes_and_points_sample_name"]), i) for i in list(group_combinations.iloc[:,0])]
    comboB = [re.sub("\.|-", str(config["character_subsitution_dashes_and_points_sample_name"]), i) for i in list(group_combinations.iloc[:,1])]

    combo_list = []
    for (A, B) in zip(comboA, comboB):
        combo_list.append('_vs_'.join([str(A), str(B)]))

    # Define differential compartments
    if (eval(str(config["perform_differential_contacts_analyses"])) == True & eval(str(config["groups"]["perform_grouped_analyses"])) == True):
        selfish_group = expand("12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH/{combo}/{combo}_{resolution}kb_SELFISH.txt", combo = combo_list, resolution=[str(x) for x in NEW_RESOLUTIONS])
    else:
        selfish_group = []

    # Outputs from compartimentalization analysis
    if (perform_compartment_analyses == True):
        compartments_single = expand(''.join(["10_Compartments_detection_dcHiC/{comp_res}kb_resolution/DifferentialResult/all_vs_all/viz/files_compartment_beds/intra_{sample}_mapQ", str(config["mapQ_cutoff"]), "_{comp_res}kb_PC_compartments_sorted.bed"]), sample = [re.sub("\.|-", str(config["character_subsitution_dashes_and_points_sample_name"]), i) for i in SAMPLENAMES], comp_res = MIN_RESOLUTIONS)
        compartments_single_combo = expand("10_Compartments_detection_dcHiC/{comp_res}kb_resolution/DifferentialResult/{combo}/viz/vizIGV_intra/data_bigWig/differential_compartment.Mahalanobis.bw", combo = combo_list, comp_res = MIN_RESOLUTIONS)

        if (eval(str(config["groups"]["perform_grouped_analyses"])) == True):
            compartments_group = expand(''.join(["12_Grouped_analyses/F_Compartments_detection_dcHiC/{comp_res}kb_resolution/DifferentialResult/all_vs_all/viz/files_compartment_beds/intra_{group}_mapQ", str(config["mapQ_cutoff"]), "_{comp_res}kb_PC_compartments_sorted.bed"]), group=[re.sub("\.|-", str(config["character_subsitution_dashes_and_points_sample_name"]), i) for i in groups], comp_res = MIN_RESOLUTIONS)
            compartments_group_combo = expand("12_Grouped_analyses/F_Compartments_detection_dcHiC/{comp_res}kb_resolution/DifferentialResult/{combo}/viz/vizIGV_intra/data_bigWig/differential_compartment.Mahalanobis.bw", combo = combo_list, comp_res = MIN_RESOLUTIONS)
        else:
            compartments_group = []
            compartments_group_combo = []
    else:
        compartments_single = []
        compartments_single_combo = []
        compartments_group = []
        compartments_group_combo = []
else:
    MIN_RESOLUTIONS = NEW_RESOLUTIONS
    compartments_single = []
    compartments_group = []
    combo_list = groups
    compartments_single_combo = []
    compartments_group_combo = []
    selfish_group = []


# Stripes analyses (STRIPENN)
if (eval(str(config["perform_stripes_analyses"])) == True):
    # Define resolutions and understand wheter stripes calling can be performed
    MAX_STRIPES_RESOLUTIONS = [x for x in NEW_RESOLUTIONS if x <= int(config["stripenn_params"]["max_stripes_resolution"])]

    if (len(MAX_STRIPES_RESOLUTIONS) > 0):
        stripes_single = expand("11_Stripes_analyses_STREPENN/{sample}/{sample}_{max_res}kb/result_filtered.tsv", sample = SAMPLENAMES, max_res = MAX_STRIPES_RESOLUTIONS[0])
        if (eval(str(config["groups"]["perform_grouped_analyses"])) == True):
            stripes_group = expand("12_Grouped_analyses/H_Stripes_analyses_STRIPPEN/{group}/{group}_{max_res}kb/result_filtered.tsv", group = groups, max_res = MAX_STRIPES_RESOLUTIONS[0])
        else:
            stripes_group = []
    else:
        stripes_single = []
        stripes_group = []
else:
    stripes_single = []
    stripes_group = []
    MAX_STRIPES_RESOLUTIONS = NEW_RESOLUTIONS



# Define the rescriction sites table
restriction_table = pd.DataFrame({'Enzyme': ['DpnII', 'MboI', 'NlaIII', 'Csp6I', 'CviQI', 'HindIII', 'EcoRI', 'BamHI', 'BglII', 'Sau3AI', 'Arima_4combo'],
                                  'RestrictionSeq': ['GATC', 'GATC', 'CATG', 'GTAC', 'GTAC', 'AAGCTT', 'GAATTC', 'GGATCC', 'AGATCT', 'GATC', 'GATC|GA.TC'],
                                  'DanglingSeq': ['GATC', 'GATC', 'CATG', 'TA', 'TA', 'AGCT', 'AATT', 'GATC', 'GATC', 'GATC', 'GATC|GA.TC']})

# Define other global variables
#script_dchicf_file = os.path.join(config["dcHiC_repository_folder"], "dchicf.r")
script_dchicf_file = "$CONDA_PREFIX/bin/dcHiC/dchicf.r"
stripeDiff_file = "$CONDA_PREFIX/bin/stripeDiff/src/stripeDiff.sh"

# generation of global wildcard_constraints
wildcard_constraints:
    RUNS = constraint_to(RUNNAMES),
    SAMPLES = constraint_to(SAMPLENAMES),
    READS = constraint_to(config["runs_suffix"]),
    MERGED_RESOLUTIONS = constraint_to([str(x) for x in merging_resolutions]), # for names of only merged matrices
    ALL_NEW_MATRIX_RESOLUTIONS = constraint_to([str(x) for x in NEW_RESOLUTIONS]), # to get all matrices resolutions for normalization/correction
    COMPARTMENT_RESOLUTIONS = MIN_RESOLUTIONS, # to filter all matrices resolutions for compartment calling
    STRIPES_RESOLUTIONS = MAX_STRIPES_RESOLUTIONS[0], # to filter all matrices resolutions for stripes calling
    LOOPS_RESOLUTIONS = constraint_to([str(x) for x in MAX_LOOPS_RESOLUTIONS]), # to filter all matrices resolutions for loops calling
    GROUPS = constraint_to(groups),
    COMBOS = constraint_to(combo_list)


# Defining inputs for ALL Rule
if len(requested_resolutions)>1:
    new_matrix = expand(os.path.join("04_Interaction_matrices/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb.h5"])), sample=SAMPLENAMES, merged_res = [str(x) for x in merging_resolutions])
else:
    new_matrix = expand(os.path.join("04_Interaction_matrices/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{resolution}kb.cool"])), sample=SAMPLENAMES, resolution=requested_resolutions[0])


ruleorder:
    A_fastQC_raw > B_multiQC_raw > C_bwa_align


# ========================================================================================
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================
# Function to run all funtions
rule AAA_initialization:
    input:
        fastQC_raw_zip = expand(os.path.join("01_fastQC_raw/", "{run}_fastqc.zip"), run=RUNNAMES),
        multiQC_raw_html = "01_fastQC_raw/multiQC_raw/multiQC_report_fastqRaw.html",
        #sam = expand(os.path.join("02_Alignements/{sample}{read}.bam"), sample=SAMPLENAMES, read=config["runs_suffix"]),
        interaction_matrices = expand(os.path.join("04_Interaction_matrices/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{resolution}kb.cool"])), sample=SAMPLENAMES, resolution=requested_resolutions[0]),
        new_matrices = new_matrix,
        multiQC_report_matrices = "04_Interaction_matrices/QC_matrices/ALL_SAMPLES/multiQC_matrices/multiqc_report.html",
        normalized_matrices = expand(os.path.join("05_Interaction_matrices_normalized/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{resolution}kb_normalized.h5"])), sample = SAMPLENAMES, resolution=[str(x) for x in NEW_RESOLUTIONS]),
        heatmap_plot = "05_Interaction_matrices_normalized/sample_correlation/heatmap_correlation.pdf",
        scatter_plot = "05_Interaction_matrices_normalized/sample_correlation/scatter_correlation.pdf",
        h5_matrix_corrected_single_sample = expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/h5_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{resolution}kb_normalized_corrected.h5"])), sample = SAMPLENAMES, resolution=[str(x) for x in NEW_RESOLUTIONS]),
        cool_matrix_corrected_single_sample = expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/cool_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{resolution}kb_normalized_corrected.cool"])), sample = SAMPLENAMES, resolution=[str(x) for x in NEW_RESOLUTIONS]),
        #hicpro_matrix_corrected_single_sample = expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/hicpro_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{resolution}kb_normalized_corrected.hicpro"])), sample = SAMPLENAMES, resolution=[str(x) for x in NEW_RESOLUTIONS]),
        domains_bed = domains_bed_single,
        plot_chr_dist = expand("08_Interaction_distances/intraChr_distances_all_samples_{resolution}kb_resolution.pdf", resolution=[str(x) for x in NEW_RESOLUTIONS]),
        grouped_analyses_inputs = grouped_analyses_inputs,
        loops_single = loops_single,
        loops_group = loops_group,
        compartments_single = compartments_single,
        compartments_group = compartments_group,
        compartments_single_combo = compartments_single_combo,
        compartments_group_combo = compartments_group_combo,
        selfish_group = selfish_group,
        stripes_single = stripes_single,
        stripes_group = stripes_group
    shell:
        """
        printf '\033[1;36mPipeline ended!\\n\033[0m'
        """
# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================




# ----------------------------------------------------------------------------------------
# Perform the FastQC on raw fastq.gz
rule A_fastQC_raw:
    input:
        fastq_gz = os.path.join(config["runs_directory"], ''.join(["{RUNS}", config['fastq_extension']]))
    output:
        html = os.path.join("01_fastQC_raw/","{RUNS}_fastqc.html"),
        zip =  os.path.join("01_fastQC_raw/","{RUNS}_fastqc.zip")
    params:
        build_fastqcDir = os.path.dirname("01_fastQC_raw/multiQC_raw/"),
        fastQC_raw_outdir = os.path.join(home_dir, "01_fastQC_raw/"),
        run = "{RUNS}"
    threads:
        math.floor(workflow.cores/len(RUNNAMES))
    benchmark:
        "benchmarks/A_fastQC_raw/A_fastQC_raw-{RUNS}.tsv"
    shell:
        """
        printf '\033[1;36m{params.run}: Performing fastQC on raw fastq...\\n\033[0m'

        mkdir -p {params.build_fastqcDir}
        fastqc -t {threads} --outdir {params.fastQC_raw_outdir} {input.fastq_gz}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform multiQC for raw fastq reports
rule B_multiQC_raw:
    input:
        fastqc_zip = expand(os.path.join("01_fastQC_raw/", "{run}_fastqc.zip"), run=RUNNAMES)
    output:
        multiqcReportRaw = "01_fastQC_raw/multiQC_raw/multiQC_report_fastqRaw.html"
    params:
        fastqc_Raw_reports = os.path.join("01_fastQC_raw/", "*.zip"),
        multiQC_raw_outdir = os.path.join(home_dir, "01_fastQC_raw/multiQC_raw/")
    benchmark:
        "benchmarks/B_multiQC_raw/B_multiQC_raw.tsv"
    shell:
        """
        printf '\033[1;36mGenerating multiQC report for fastq quality test...\\n\033[0m'
        multiqc -f --outdir {params.multiQC_raw_outdir} -n multiQC_report_fastqRaw.html {params.fastqc_Raw_reports}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
if not os.path.exists(''.join([re.sub(".gz", "", config["genome_fasta"], count=0, flags=0),".bwt.2bit.64"])):
    # ----------------------------------------------------------------------------------------
    # Reads alignement
    rule Cextra_generate_genome_index:
        input:
            genome = ancient(config["genome_fasta"])
        output:
            genome_fai = ''.join([re.sub(".gz", "", config["genome_fasta"], count=0, flags=0),".bwt.2bit.64"])
        threads: 1
        benchmark:
            "benchmarks/Cextra_generate_genome_index/Cextra_generate_genome_index.tsv"
        shell:
            """
            printf '\033[1;36mGenerating the genome index...\\n\033[0m'
            $CONDA_PREFIX/bin/bwa-mem2 index {input.genome}
            samtools faidx {input.genome}
            printf '\033[1;36mGenome index done.\\n\033[0m'
            """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Reads alignement
rule C_bwa_align:
    input:
        read = os.path.join(config["runs_directory"], ''.join(["{SAMPLES}", "{READS}", config['fastq_extension']])),
        genome_fai = ancient(''.join([re.sub(".gz", "", config["genome_fasta"], count=0, flags=0),".bwt.2bit.64"]))
    output:
        align = os.path.join("02_Alignements/{SAMPLES}{READS}.bam")
    params:
        build_align = "02_Alignements/log/",
        genome = re.sub(".gz", "", config["genome_fasta"], count=0, flags=0),
        sample = "{SAMPLES}",
        read = "{READS}"
    threads:
        #max([math.floor(workflow.cores-2), 1])
        math.floor(workflow.cores/2)
    log:
        out = os.path.join("02_Alignements/log/{SAMPLES}{READS}_bwa-mem2.out"),
        err = os.path.join("02_Alignements/log/{SAMPLES}{READS}_bwa-mem2.err")
    benchmark:
        "benchmarks/C_bwa_align/C_bwa_align-{SAMPLES}_{READS}.tsv"
    shell:
        """
        mkdir -p {params.build_align}

        printf '\033[1;36m{params.sample}{params.read}: alignment of the reads...\\n\033[0m'
        echo 'Mapping of {params.sample}{params.read}:' > {log.out}
        $CONDA_PREFIX/bin/bwa-mem2 mem -t {threads} -A1 -B4 -E50 -L0 {params.genome} {input.read} 2> {log.err} | samtools view -@ 1 -Shb - > {output.align} 2>> {log.err}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
rule D_generate_restriction_file_and_get_chrSizes:
    input:
        genome_fai = ''.join([re.sub(".gz", "", config["genome_fasta"], count=0, flags=0),".fai"]),
        genome_fasta = config["genome_fasta"]
    output:
        restriction_bed = temp("temp_resctriction_bed_file.bed"),
        chrSizes = temp("temp_chrSizes_file.txt")
    params:
        genome_assembly_name = config["genome_assembly_name"],
        restriction_sequence = (restriction_table[restriction_table['Enzyme']==config["restriction_enzyme"]].RestrictionSeq).iloc[0],
        enzyme = config["restriction_enzyme"]
    benchmark:
        "benchmarks/D_generate_restriction_file_and_get_chrSizes/D_generate_restriction_file_and_get_chrSizes.tsv"
    shell:
        """
        printf '\033[1;36mGeneration of the restriction sites locatation file for {params.enzyme}...\\n\033[0m'
        hicFindRestSite -f {input.genome_fasta} --searchPattern {params.restriction_sequence} -o {output.restriction_bed}

        printf '\033[1;36mGetting chromosome sizes for {params.genome_assembly_name} genome...\\n\033[0m'
        cut -f 1,2 {input.genome_fai} > {output.chrSizes}
        """
# ----------------------------------------------------------------------------------------


if (eval(str(config["generate_bam"])) == True):
    # ----------------------------------------------------------------------------------------
    # Matrix builing and bam generation
    rule E1_interaction_matrix_and_bam_generation_at_smallest_resolution:
        input:
            align_R1 = os.path.join("02_Alignements/", ''.join(["{SAMPLES}", config["runs_suffix"][0], ".bam"])),
            align_R2 = os.path.join("02_Alignements/", ''.join(["{SAMPLES}", config["runs_suffix"][1], ".bam"])),
            restriction_bed = "temp_resctriction_bed_file.bed",
            chrSizes = "temp_chrSizes_file.txt"
        output:
            interaction_matrix = os.path.join("04_Interaction_matrices/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_", str(requested_resolutions[0]), "kb.cool"])),
            interaction_matrix_h5 = os.path.join("04_Interaction_matrices/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_", str(requested_resolutions[0]), "kb.h5"])),
            QC_log_file = os.path.join("04_Interaction_matrices/QC_matrices/{SAMPLES}/QC.log")
        params:
            build_BAM_dir = os.path.dirname("03_BAM/flagstat/"),
            build_matrix_dir = os.path.dirname("04_Interaction_matrices/log/"),
            resolution = ''.join([str(requested_resolutions[0]), "000"]),
            resolution_kb = str(requested_resolutions[0]),
            restriction_sequence = (restriction_table[restriction_table['Enzyme']==config["restriction_enzyme"]].RestrictionSeq).iloc[0],
            dangling_sequence = (restriction_table[restriction_table['Enzyme']==config["restriction_enzyme"]].DanglingSeq).iloc[0],
            mapQ_cutoff = str(config["mapQ_cutoff"]),
            min_dist = config["min_distance"],
            max_dist = config["max_distance"],
            QCfolder = os.path.dirname("04_Interaction_matrices/QC_matrices/{SAMPLES}/"),
            genome_assembly_name = config["genome_assembly_name"],
            sample = "{SAMPLES}",
            filtBAM = os.path.join("03_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), ".bam"])),
            filtBAM_sorted = os.path.join("03_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted.bam"])),
            filtBAM_sorted_index = os.path.join("03_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted.bam.bai"])),
            flagstat_BAM = os.path.join("03_BAM/flagstat/", "{SAMPLES}_flagstat_bam.txt")
        log:
            out = os.path.join(''.join(["04_Interaction_matrices/log/{SAMPLES}.", str(requested_resolutions[0]), "kb.hicBuildMatrix.log"])),
            err = os.path.join(''.join(["04_Interaction_matrices/log/{SAMPLES}.", str(requested_resolutions[0]), "kb.hicBuildMatrix.err"]))
        threads:
            math.floor(workflow.cores/2)
        benchmark:
            "benchmarks/E1_interaction_matrix_and_bam_generation_at_smallest_resolution/E1_interaction_matrix_and_bam_generation_at_smallest_resolution-{SAMPLES}.tsv"
        shell:
            """
            mkdir -p {params.build_BAM_dir}
            mkdir -p {params.build_matrix_dir}
            mkdir -p {params.QCfolder}

            printf '\033[1;36m{params.sample}: building interaction matrix at {params.resolution_kb}kb resolution...\\n\033[0m'

            hicBuildMatrix \
            --samFiles {input.align_R1} {input.align_R2} \
            --restrictionCutFile {input.restriction_bed} \
            --binSize {params.resolution} \
            --restrictionSequence {params.restriction_sequence} \
            --danglingSequence {params.dangling_sequence} \
            --minDistance {params.min_dist} \
            --maxLibraryInsertSize {params.max_dist} \
            --QCfolder {params.QCfolder} \
            --minMappingQuality {params.mapQ_cutoff} \
            --outBam {params.filtBAM} \
            --chromosomeSizes {input.chrSizes} \
            --genomeAssembly {params.genome_assembly_name} \
            --threads {threads} \
            --outFileName {output.interaction_matrix} > {log.out} 2> {log.err}

            hicConvertFormat -m {output.interaction_matrix} -o {output.interaction_matrix_h5} --inputFormat cool --outputFormat h5

            printf '\033[1;36m{params.sample}: sorting BAM...\\n\033[0m'
            samtools sort -@ {threads} {params.filtBAM} -o {params.filtBAM_sorted}
            samtools index -@ {threads} -b {params.filtBAM_sorted} {params.filtBAM_sorted_index}

            printf '\033[1;36m{params.sample}: Getting flagstat from unfiltered BAM...\\n\033[0m'
            samtools flagstat {params.filtBAM_sorted} -@ {threads} > {params.flagstat_BAM}

            rm {params.filtBAM}
            """
    # ----------------------------------------------------------------------------------------
else:
    # ----------------------------------------------------------------------------------------
    # Matrix building (only)
    rule E1_interaction_matrix_generation_at_smallest_resolution:
        input:
            align_R1 = os.path.join("02_Alignements/", ''.join(["{SAMPLES}", config["runs_suffix"][0], ".bam"])),
            align_R2 = os.path.join("02_Alignements/", ''.join(["{SAMPLES}", config["runs_suffix"][1], ".bam"])),
            restriction_bed = "temp_resctriction_bed_file.bed",
            chrSizes = "temp_chrSizes_file.txt"
        output:
            interaction_matrix = os.path.join("04_Interaction_matrices/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_", str(requested_resolutions[0]), "kb.cool"])),
            interaction_matrix_h5 = os.path.join("04_Interaction_matrices/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_", str(requested_resolutions[0]), "kb.h5"])),
            QC_log_file = os.path.join("04_Interaction_matrices/QC_matrices/{SAMPLES}/QC.log")
        params:
            build_BAM_dir = os.path.dirname("03_BAM__not_generated/"),
            build_matrix_dir = os.path.dirname("04_Interaction_matrices/log/"),
            resolution = ''.join([str(requested_resolutions[0]), "000"]),
            resolution_kb = str(requested_resolutions[0]),
            restriction_sequence = (restriction_table[restriction_table['Enzyme']==config["restriction_enzyme"]].RestrictionSeq).iloc[0],
            dangling_sequence = (restriction_table[restriction_table['Enzyme']==config["restriction_enzyme"]].DanglingSeq).iloc[0],
            mapQ_cutoff = str(config["mapQ_cutoff"]),
            min_dist = config["min_distance"],
            max_dist = config["max_distance"],
            QCfolder = os.path.dirname("04_Interaction_matrices/QC_matrices/{SAMPLES}/"),
            genome_assembly_name = config["genome_assembly_name"],
            sample = "{SAMPLES}"
        log:
            out = os.path.join(''.join(["04_Interaction_matrices/log/{SAMPLES}.", str(requested_resolutions[0]), "kb.hicBuildMatrix.log"])),
            err = os.path.join(''.join(["04_Interaction_matrices/log/{SAMPLES}.", str(requested_resolutions[0]), "kb.hicBuildMatrix.err"]))
        benchmark:
            "benchmarks/E1_interaction_matrix_and_bam_generation_at_smallest_resolution/E1_interaction_matrix_and_bam_generation_at_smallest_resolution-{SAMPLES}.tsv"
        threads:
            math.floor(workflow.cores/2)
        shell:
            """
            mkdir -p {params.build_BAM_dir}
            mkdir -p {params.build_matrix_dir}
            mkdir -p {params.QCfolder}

            printf '\033[1;36m{params.sample}: building interaction matrix at {params.resolution_kb}kb resolution...\\n\033[0m'

            hicBuildMatrix \
            --samFiles {input.align_R1} {input.align_R2} \
            --restrictionCutFile {input.restriction_bed} \
            --binSize {params.resolution} \
            --restrictionSequence {params.restriction_sequence} \
            --danglingSequence {params.dangling_sequence} \
            --minDistance {params.min_dist} \
            --maxLibraryInsertSize {params.max_dist} \
            --QCfolder {params.QCfolder} \
            --minMappingQuality {params.mapQ_cutoff} \
            --chromosomeSizes {input.chrSizes} \
            --genomeAssembly {params.genome_assembly_name} \
            --threads {threads} \
            --outFileName {output.interaction_matrix} > {log.out} 2> {log.err}

            printf '\033[1;36m{params.sample}: conversion of {params.resolution_kb}kb resolution matrix from .h5 to .cool...\\n\033[0m'
            hicConvertFormat -m {output.interaction_matrix} -o {output.interaction_matrix_h5} --inputFormat cool --outputFormat h5
            """
        # ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Print the multiQC report
rule E2_multiQC_report_for_HiC_matrices:
    input:
        QC_log_files = expand(os.path.join("04_Interaction_matrices/QC_matrices/{sample}/QC.log"), sample = SAMPLENAMES)
    output:
        multiQC_report = os.path.join("04_Interaction_matrices/QC_matrices/ALL_SAMPLES/multiQC_matrices/multiqc_report.html")
    params:
        QC_logs_names = ' '.join(expand(os.path.join("04_Interaction_matrices/QC_matrices/{sample}/QC.log"), sample = SAMPLENAMES)),
        labels = ' '.join(SAMPLENAMES),
        QC_general_folder = os.path.dirname("04_Interaction_matrices/QC_matrices/"),
        multiQC_dir = os.path.dirname("04_Interaction_matrices/QC_matrices/ALL_SAMPLES/multiQC_matrices/"),
        hiQC_dir = os.path.dirname("04_Interaction_matrices/QC_matrices/ALL_SAMPLES/HiQC_matrices/")
    threads: 1
    benchmark:
        "benchmarks/E2_multiQC_report_for_HiC_matrices/E2_multiQC_report_for_HiC_matrices.tsv"
    shell:
        """
        mkdir -p {params.hiQC_dir}
        mkdir -p {params.multiQC_dir}

        printf '\033[1;36mGenerating multiQC report on matrices...\\n\033[0m'
        multiqc {params.QC_general_folder} -o {params.multiQC_dir}


        printf '\033[1;36mGenerating HiQC report on matrices...\\n\033[0m'
        hicQC \
        --logfiles {params.QC_logs_names} \
        --labels {params.labels} \
        -o {params.hiQC_dir}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Merging bins to make all the other resolution matrix (if multiple)
if len(requested_resolutions) > 1 :
    rule E3_merging_interaction_matrix_bins_for_all_resolutions:
        input:
            base_matrix = os.path.join("04_Interaction_matrices/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_", str(requested_resolutions[0]), "kb.cool"]))
        output:
            new_matrix = os.path.join("04_Interaction_matrices/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{MERGED_RESOLUTIONS}kb.h5"]))
        params:
            sample = "{SAMPLES}",
            base_matrix_h5 = os.path.join("04_Interaction_matrices/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_", str(requested_resolutions[0]), "kb.h5"])),
            base_resolution = requested_resolutions[0],
            merged_resolution = "{MERGED_RESOLUTIONS}"
        log:
            out = "04_Interaction_matrices/log/{SAMPLES}.{MERGED_RESOLUTIONS}kb.hicBuildMatrix.log",
            err = "04_Interaction_matrices/log/{SAMPLES}.{MERGED_RESOLUTIONS}kb.hicBuildMatrix.err"
        benchmark:
            "benchmarks/E3_merging_interaction_matrix_bins_for_all_resolutions/E3_merging_interaction_matrix_bins_for_all_resolutions-{SAMPLES}.{MERGED_RESOLUTIONS}kb.tsv"
        shell:
            """
            FINAL_RESOLUTION=$(({params.merged_resolution} / {params.base_resolution}))

            printf '\033[1;36m{params.sample}: merging %s adjacent bins in the original matrix to create a new one with a resolution of {params.merged_resolution}kb...\\n\033[0m' $FINAL_RESOLUTION
            hicMergeMatrixBins -m {params.base_matrix_h5} -nb $FINAL_RESOLUTION -o {output.new_matrix} >{log.out} 2>{log.err}
            """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Matrices normalization
rule F1_matrices_normalization:
    input:
        matrices_to_normalize = expand(os.path.join("04_Interaction_matrices/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb.h5"])), sample = SAMPLENAMES, allow_missing=True)
    output:
        normalized_matrices = expand(os.path.join("05_Interaction_matrices_normalized/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized.h5"])), sample = SAMPLENAMES, allow_missing=True)
    params:
        resolution = "{ALL_NEW_MATRIX_RESOLUTIONS}",
        normalization_method = config["normalization_method"]
    benchmark:
        "benchmarks/F1_matrices_normalization/F1_matrices_normalization-{ALL_NEW_MATRIX_RESOLUTIONS}kb.tsv"
    shell:
        """
        mkdir -p 05_Interaction_matrices_normalized/

        printf '\033[1;36mNormalization of {params.resolution}kb matrices \\n\033[0m'

        hicNormalize \
        --matrices {input.matrices_to_normalize} \
        --normalize {params.normalization_method} \
        --outFileName {output.normalized_matrices}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Correlation of the matrices
rule F2_samples_correlation:
    input:
        base_matrices = expand(os.path.join("05_Interaction_matrices_normalized/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_", str(requested_resolutions[0]), "kb_normalized.h5"])), sample=SAMPLENAMES)
    output:
        heatmap_plot = "05_Interaction_matrices_normalized/sample_correlation/heatmap_correlation.pdf",
        scatter_plot = "05_Interaction_matrices_normalized/sample_correlation/scatter_correlation.pdf"
    params:
        matrices_names = ' '.join(expand(os.path.join("05_Interaction_matrices_normalized/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_", str(requested_resolutions[0]), "kb_normalized.h5"])), sample=SAMPLENAMES)),
        labels = ' '.join(SAMPLENAMES),
        correlation_dir = os.path.dirname("05_Interaction_matrices_normalized/sample_correlation/"),
        correlation_method = config["correlation_method"],
        heatmap_color = config["heatmap_color"]
    benchmark:
        "benchmarks/F2_samples_correlation/F2_samples_correlation.tsv"
    shell:
        """
        mkdir -p {params.correlation_dir}

        printf '\033[1;36mGenerating correlation plots...\\n\033[0m'

        hicCorrelate \
        --matrices {params.matrices_names} \
        --method {params.correlation_method} \
        --colorMap {params.heatmap_color} \
        --labels {params.labels} \
        --plotNumbers \
        --outFileNameHeatmap {output.heatmap_plot} \
        --outFileNameScatter {output.scatter_plot}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Matrices correction: generate the thresholds
rule G1_matrices_correction__diagnosticPlot_and_MAD:
    input:
        matrix_to_correct = os.path.join("05_Interaction_matrices_normalized/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized.h5"]))
    output:
        diagnostic_plot = "06_Interaction_matrices_normalized_and_corrected/diagnostic_plots/{SAMPLES}_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_diagnosticPlot.png",
        mad = "06_Interaction_matrices_normalized_and_corrected/median_absolute_deviation/{SAMPLES}_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_MedianAbsoluteDeviation.mad"
    params:
        plots_dir = os.path.dirname("06_Interaction_matrices_normalized_and_corrected/diagnostic_plots/"),
        mad_dir = os.path.dirname("06_Interaction_matrices_normalized_and_corrected/median_absolute_deviation/"),
        sample = "{SAMPLES}",
        resolution = "{ALL_NEW_MATRIX_RESOLUTIONS}"
    benchmark:
        "benchmarks/G1_matrices_correction__diagnosticPlot_and_MAD/G1_matrices_correction__diagnosticPlot_and_MAD-{SAMPLES}_{ALL_NEW_MATRIX_RESOLUTIONS}kb.tsv"
    shell:
        """
        mkdir -p {params.plots_dir}
        mkdir -p {params.mad_dir}

        printf '\033[1;36m{params.sample}, matrix at {params.resolution}kb: diagnostic plot and MAD generation...\\n\033[0m'

        hicCorrectMatrix diagnostic_plot --matrix {input.matrix_to_correct} -o {output.diagnostic_plot} 2> {output.mad}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Matrices correction: get the thresholds
rule G2_matrices_correction__getting_threshold_values:
    input:
        mad = os.path.join("06_Interaction_matrices_normalized_and_corrected/median_absolute_deviation/{SAMPLES}_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_MedianAbsoluteDeviation.mad")
    output:
        threshold_file = os.path.join("06_Interaction_matrices_normalized_and_corrected/median_absolute_deviation/thresholds/{SAMPLES}_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_thresholdValues.txt")
    params:
        threshold_dir = os.path.dirname("06_Interaction_matrices_normalized_and_corrected/median_absolute_deviation/thresholds/"),
        sample = "{SAMPLES}",
        resolution = "{ALL_NEW_MATRIX_RESOLUTIONS}"
    benchmark:
        "benchmarks/G2_matrices_correction__getting_threshold_values/G2_matrices_correction__getting_threshold_values-{SAMPLES}_{ALL_NEW_MATRIX_RESOLUTIONS}kb.tsv"
    shell:
        """
        mkdir -p {params.threshold_dir}

        printf '\033[1;36m{params.sample}, matrix at {params.resolution}kb: getting threshold values for correction...\\n\033[0m'

        madscore=$(grep \"mad threshold \" {input.mad} | sed 's/INFO:hicexplorer.hicCorrectMatrix:mad threshold //g');
        upper=$(echo -3*$madscore | bc);
        echo $madscore \" \" $upper >> {output}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Matrices correction: get the thresholds
rule G3_matrices_correction__correction:
    input:
        matrix_to_correct = os.path.join("05_Interaction_matrices_normalized/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized.h5"])),
        threshold_file = os.path.join("06_Interaction_matrices_normalized_and_corrected/median_absolute_deviation/thresholds/{SAMPLES}_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_thresholdValues.txt")
    output:
        h5_matrix_corrected = os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{SAMPLES}/h5_format/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_corrected.h5"]))
    params:
        corrected_matrix_dir = os.path.dirname("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{SAMPLES}/h5_format/"),
        correction_method = config["correction_method"],
        sample = "{SAMPLES}",
        resolution = "{ALL_NEW_MATRIX_RESOLUTIONS}"
    benchmark:
        "benchmarks/G3_matrices_correction__correction/G3_matrices_correction__correction-{SAMPLES}_{ALL_NEW_MATRIX_RESOLUTIONS}kb.tsv"
    shell:
        """
        mkdir -p {params.corrected_matrix_dir}

        printf '\033[1;36m{params.sample}, matrix at {params.resolution}kb: matrix correction...\\n\033[0m'

        THRESHOLDS=$(head -n 1 {input.threshold_file} | cat);

        hicCorrectMatrix correct \
        --correctionMethod {params.correction_method} \
        --filterThreshold $THRESHOLDS \
        -m {input.matrix_to_correct} \
        -o {output.h5_matrix_corrected} >> {input.threshold_file}
        """
# ----------------------------------------------------------------------------------------


# # ----------------------------------------------------------------------------------------
# # Conversion of the corrected matrices format H5 -> cool
rule H1_matrices_format_conversion__cool:
    input:
        h5_matrix_corrected = expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/h5_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.h5"])), sample=SAMPLENAMES, merged_res = [str(x) for x in NEW_RESOLUTIONS])
    output:
        cool_matrix_corrected = expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/cool_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.cool"])), sample=SAMPLENAMES, merged_res = [str(x) for x in NEW_RESOLUTIONS])
    params:
        all_names_h5 = ' '.join(expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/h5_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.h5"])), sample=SAMPLENAMES, merged_res = [str(x) for x in NEW_RESOLUTIONS])),
        all_names_cool = ' '.join(expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/cool_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.cool"])), sample=SAMPLENAMES, merged_res = [str(x) for x in NEW_RESOLUTIONS])),
        cool_corrected_matrix_dir = ' '.join(expand(os.path.dirname("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/cool_format/"), sample=SAMPLENAMES)),
        resolutions = ' '.join([str(x*1000) for x in NEW_RESOLUTIONS]),
        correction_method = str(config["correction_method"])
    benchmark:
        "benchmarks/H1_matrices_format_conversion__cool/H1_matrices_format_conversion__cool.tsv"
    shell:
        """
        for DIR in {params.cool_corrected_matrix_dir}
        do
            mkdir -p $DIR
        done

        printf '\033[1;36mConversion of corrected normalized matrices from .h5 to .cool format ...\\n\033[0m'

        hicConvertFormat \
        -m {params.all_names_h5} \
        -o {params.all_names_cool} \
        --resolutions {params.resolutions} \
        --correction_name {params.correction_method} \
        --inputFormat h5 \
        --outputFormat cool
        """
# ----------------------------------------------------------------------------------------


# # ----------------------------------------------------------------------------------------
# # Conversion of the corrected matrices format H5 -> hicpro
rule H2_matrices_format_conversion__hicpro:
    input:
        h5_matrix_corrected = expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/h5_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.h5"])), sample=SAMPLENAMES, merged_res = [str(x) for x in NEW_RESOLUTIONS])
    output:
        hicpro_matrix_corrected = expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/hicpro_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.hicpro"])), sample=SAMPLENAMES, merged_res = [str(x) for x in NEW_RESOLUTIONS]),
        hicpro_matrix_corrected_bed = expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/hicpro_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected_hicpro.bed"])), sample=SAMPLENAMES, merged_res = [str(x) for x in NEW_RESOLUTIONS])
    params:
        all_names_h5 = ' '.join(expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/h5_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.h5"])), sample=SAMPLENAMES, merged_res = [str(x) for x in NEW_RESOLUTIONS])),
        all_names_hicpro = ' '.join(expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/hicpro_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.hicpro"])), sample=SAMPLENAMES, merged_res = [str(x) for x in NEW_RESOLUTIONS])),
        all_names_hicpro_bed = ' '.join(expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/hicpro_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected_hicpro.bed"])), sample=SAMPLENAMES, merged_res = [str(x) for x in NEW_RESOLUTIONS])),
        hicpro_corrected_matrix_dir = ' '.join(expand(os.path.dirname("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/hicpro_format/"), sample=SAMPLENAMES))
    benchmark:
        "benchmarks/H2_matrices_format_conversion__hicpro/H2_matrices_format_conversion__hicpro.tsv"
    shell:
        """
        for DIR in {params.hicpro_corrected_matrix_dir}
        do
            mkdir -p $DIR
        done

        printf '\033[1;36mConversion of corrected matrices from .h5 to .hicpro format ...\\n\033[0m'

        hicConvertFormat \
        -m {params.all_names_h5} \
        -o {params.all_names_hicpro} \
        --bedFileHicpro {params.all_names_hicpro_bed} \
        --inputFormat h5 \
        --outputFormat hicpro
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Calling TADs (HiCexplorer)
rule I_call_TADs_HiCexplorer:
    input:
        h5_matrix_corrected = os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{SAMPLES}/h5_format/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_corrected.h5"]))
    output:
        domains_bed = os.path.join("07_TADs_calling_HiCexplorer/{SAMPLES}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_domains.bed"]))
    params:
        TADs_dir = os.path.dirname("07_TADs_calling_HiCexplorer/{SAMPLES}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/log/"),
        prefix = str(os.path.join("07_TADs_calling_HiCexplorer/{SAMPLES}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb"]))),
        extra_parameters = config["extra_findTAD_parameters"],
        sample = "{SAMPLES}",
        resolution = "{ALL_NEW_MATRIX_RESOLUTIONS}"
    log:
        out = os.path.join("07_TADs_calling_HiCexplorer/{SAMPLES}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/log/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_TAD.calling.out"])),
        err = os.path.join("07_TADs_calling_HiCexplorer/{SAMPLES}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/log/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_TAD.calling.err"]))
    threads:
        math.floor(workflow.cores/4)
    benchmark:
        "benchmarks/I_call_TADs_HiCexplorer/I_call_TADs_HiCexplorer-{SAMPLES}_{ALL_NEW_MATRIX_RESOLUTIONS}.tsv"
    shell:
        """
        printf '\033[1;36m{params.sample}, matrix at {params.resolution}kb: calling TADs (HiCexplorer)...\\n\033[0m'

        mkdir -p {params.TADs_dir}

        hicFindTADs \
        -m {input.h5_matrix_corrected} \
        {params.extra_parameters} \
        --correctForMultipleTesting 'bonferroni' \
        -p {threads} \
        --outPrefix {params.prefix} > {log.out} 2> {log.err}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Calling TADs (GENOVA)
rule I_call_TADs_GENOVA:
    input:
        cool_matrix_corrected = os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{SAMPLES}/cool_format/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_corrected.cool"]))
    output:
        domains_bed = os.path.join("07_TADs_calling_GENOVA/{SAMPLES}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_domains.bed"]))
    params:
        TADs_dir = os.path.dirname("07_TADs_calling_GENOVA/{SAMPLES}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/log/"),
        prefix = str(os.path.join("07_TADs_calling_GENOVA/{SAMPLES}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb"]))),
        sample = "{SAMPLES}",
        resolution = "{ALL_NEW_MATRIX_RESOLUTIONS}",
        workdir = str(home_dir)
    log:
        out = os.path.join("07_TADs_calling_GENOVA/{SAMPLES}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/log/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_TAD.calling.log"]))
    threads:
        math.floor(workflow.cores/4)
    benchmark:
        "benchmarks/I_call_TADs_GENOVA/I_call_TADs_GENOVA-{SAMPLES}_{ALL_NEW_MATRIX_RESOLUTIONS}.tsv"
    shell:
        """
        printf '\033[1;36m{params.sample}, matrix at {params.resolution}kb: calling TADs (GENOVA)...\\n\033[0m'

        mkdir -p {params.TADs_dir}

        $CONDA_PREFIX/bin/Rscript \
        -e "matrix=GENOVA::load_contacts(signal_path='{params.workdir}{input.cool_matrix_corrected}',balancing=T)" \
        -e "insulation=GENOVA::insulation_score(explist=matrix,window=25)" \
        -e "saveRDS(insulation,file='{params.workdir}{params.prefix}_GENOVA.insulation.object.Rdata')" \
        -e "options(scipen=999)" \
        -e "write.table(insulation[[1]],file='{params.workdir}{params.prefix}_insulation.score.bedGraph',row.names=F,col.names=F,quote=F)" \
        -e "write.table(GENOVA::call_TAD_insulation(insulation),file='{params.workdir}{params.prefix}_domains.bed',row.names=F,col.names=F,quote=F)" &> {log.out}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Plotting the intra-chr distances among all samples
rule J_plotting_intraChr_distances:
    input:
        h5_matrix_corrected = expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/h5_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_corrected.h5"])), sample = SAMPLENAMES, allow_missing=True),
    output:
        plot_TAD_dist = "08_Interaction_distances/intraChr_distances_all_samples_{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution.pdf",
        plot_TAD_dist_tb = "08_Interaction_distances/intraChr_distances_all_samples_{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution_distancesValues.tsv"
    params:
        resolution = "{ALL_NEW_MATRIX_RESOLUTIONS}",
        matrices = ' '.join(expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/h5_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_corrected.h5"])), sample = SAMPLENAMES, allow_missing=True)),
        plot_dir = os.path.dirname("08_Interaction_distances/log/"),
        labels = ' '.join(SAMPLENAMES)
    log:
        out = "08_Interaction_distances/log/intraChr_distances_all_samples_{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution.out",
        err = "08_Interaction_distances/log/intraChr_distances_all_samples_{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution.err"
    benchmark:
        "benchmarks/J_plotting_intraChr_distances/J_plotting_intraChr_distances-{ALL_NEW_MATRIX_RESOLUTIONS}.tsv"
    shell:
        """
        printf '\033[1;36m{params.resolution}kb matrices: plotting intra-chromosome distances...\\n\033[0m'
        mkdir -p {params.plot_dir}

        hicPlotDistVsCounts \
        --matrices {params.matrices} \
        --labels {params.labels} \
        --outFileData {output.plot_TAD_dist_tb} \
        -o {output.plot_TAD_dist} > {log.out} 2> {log.err}
        """
# ----------------------------------------------------------------------------------------


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@  Grouped analyses (if required) @@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if (eval(str(config["groups"]["perform_grouped_analyses"])) == True):
    rule L1_sum_matrices_by_group:
        input:
            normalized_matrices = expand(os.path.join("04_Interaction_matrices/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{resolution}kb.h5"])), sample = SAMPLENAMES, resolution=str(requested_resolutions[0]))
        output:
            summed_matrices = expand(os.path.join("12_Grouped_analyses/A_summed_matrices/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{resolution}kb.h5"])), group = groups, resolution=str(requested_resolutions[0]))
        params:
            samples_group_list = samples_by_group_list,
            groups = groups,
            resolution = str(requested_resolutions[0]),
            grouped_analyses_dir = os.path.dirname("12_Grouped_analyses/A_summed_matrices/")
        benchmark:
            "benchmarks/L1_sum_matrices_by_group/L1_sum_matrices_by_group.tsv"
        run:
            shell("mkdir -p {params.grouped_analyses_dir}")

            for g in range(len(params.groups)):
                # Get a unique string with the names of the matrices to sum for the current group
                matrices_to_merge = ' '.join(expand(os.path.join("04_Interaction_matrices/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{resolution}kb.h5"])), sample = params.samples_group_list[g], resolution = params.resolution))
                samples_merged_names = ", ".join(params.samples_group_list[g])
                current_group = groups[g]

                # Define the output matrix name
                output_matrix = os.path.join("12_Grouped_analyses/A_summed_matrices/", ''.join([groups[g], "_mapQ", str(config["mapQ_cutoff"]), "_", params.resolution, "kb.h5"]))

                # Sum matrices
                shell("printf '\033[1;36mGroup {current_group}: summing matrices at {params.resolution}kb-resolution of samples {samples_merged_names}...\\n\033[0m'")
                shell("hicSumMatrices --matrices {matrices_to_merge} --outFileName {output_matrix}")
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Merging bins to make all the other resolution matrix 9if multiple)
    if len(requested_resolutions) > 1 :
        rule L2_merging_grouped_interaction_matrix_bins_for_all_resolutions:
            input:
                base_matrix_h5 = os.path.join("12_Grouped_analyses/A_summed_matrices/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_", str(requested_resolutions[0]), "kb.h5"]))
            output:
                new_matrix = os.path.join("12_Grouped_analyses/A_summed_matrices/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{MERGED_RESOLUTIONS}kb.h5"]))
            params:
                group = "{GROUPS}",
                base_resolution = requested_resolutions[0],
                merged_resolution = "{MERGED_RESOLUTIONS}",
                log_dir = os.path.dirname("12_Grouped_analyses/A_summed_matrices/log/")
            log:
                out = os.path.join(''.join(["12_Grouped_analyses/A_summed_matrices/log/{GROUPS}.{MERGED_RESOLUTIONS}kb.hicBuildMatrix.log"])),
                err = os.path.join(''.join(["12_Grouped_analyses/A_summed_matrices/log/{GROUPS}.{MERGED_RESOLUTIONS}kb.hicBuildMatrix.err"]))
            benchmark:
                "benchmarks/L2_merging_grouped_interaction_matrix_bins_for_all_resolutions/L2_merging_grouped_interaction_matrix_bins_for_all_resolutions-{GROUPS}_{MERGED_RESOLUTIONS}kb.tsv"
            shell:
                """
                mkdir -p {params.log_dir}

                FINAL_RESOLUTION=$(({params.merged_resolution} / {params.base_resolution}))

                printf '\033[1;36mGroup {params.group}: merging %s adjacent bins in the original matrix to create a new one with a resolution of {params.merged_resolution}kb...\\n\033[0m' $FINAL_RESOLUTION
                hicMergeMatrixBins -m {input.base_matrix_h5} -nb $FINAL_RESOLUTION -o {output.new_matrix} >{log.out} 2>{log.err}
                """
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Grouped Matrices normalization
    rule M_grouped_matrices_normalization:
        input:
            matrices_to_normalize = expand(os.path.join("12_Grouped_analyses/A_summed_matrices/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb.h5"])), group = groups, allow_missing = True)
        output:
            normalized_matrices = expand(os.path.join("12_Grouped_analyses/B_summed_matrices_normalized/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized.h5"])), group = groups, allow_missing = True)
        params:
            normalization_dir = os.path.dirname("12_Grouped_analyses/B_summed_matrices_normalized/"),
            normalization_method = config["normalization_method"]
        benchmark:
            "benchmarks/M_grouped_matrices_normalization/M_grouped_matrices_normalization-{ALL_NEW_MATRIX_RESOLUTIONS}kb.tsv"
        shell:
            """
            mkdir -p {params.normalization_dir}

            printf '\033[1;36mNormalization of the summed matrices...\\n\033[0m'

            hicNormalize \
            --matrices {input.matrices_to_normalize} \
            --normalize {params.normalization_method} \
            --outFileName {output.normalized_matrices}
            """
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Summed matrices correction: generate the thresholds
    rule N1_summed_matrices_correction__diagnosticPlot_and_MAD:
        input:
            matrix_to_correct = os.path.join("12_Grouped_analyses/B_summed_matrices_normalized/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized.h5"]))
        output:
            diagnostic_plot = os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/diagnostic_plots/{GROUPS}_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_diagnosticPlot.png"),
            mad = os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/median_absolute_deviation/{GROUPS}_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_MedianAbsoluteDeviation.mad")
        params:
            plots_dir = os.path.dirname("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/diagnostic_plots/"),
            mad_dir = os.path.dirname("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/median_absolute_deviation/"),
            group = "{GROUPS}",
            resolution = "{ALL_NEW_MATRIX_RESOLUTIONS}"
        benchmark:
            "benchmarks/N1_summed_matrices_correction__diagnosticPlot_and_MAD/N1_summed_matrices_correction__diagnosticPlot_and_MAD-{GROUPS}_{ALL_NEW_MATRIX_RESOLUTIONS}kb.tsv"
        shell:
            """
            mkdir -p {params.plots_dir}
            mkdir -p {params.mad_dir}

            printf '\033[1;36mGroup {params.group}, matrix at {params.resolution}kb: diagnostic plot and MAD generation...\\n\033[0m'

            hicCorrectMatrix diagnostic_plot --matrix {input.matrix_to_correct} -o {output.diagnostic_plot} 2> {output.mad}
            """
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Matrices correction: get the thresholds
    rule N2_summed_matrices_correction__getting_threshold_values:
        input:
            mad = os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/median_absolute_deviation/{GROUPS}_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_MedianAbsoluteDeviation.mad")
        output:
            threshold_file = os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/median_absolute_deviation/thresholds/{GROUPS}_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_thresholdValues.txt")
        params:
            threshold_dir = os.path.dirname("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/median_absolute_deviation/thresholds/"),
            group = "{GROUPS}",
            resolution = "{ALL_NEW_MATRIX_RESOLUTIONS}"
        benchmark:
            "benchmarks/N2_summed_matrices_correction__getting_threshold_values/N2_summed_matrices_correction__getting_threshold_values-{GROUPS}_{ALL_NEW_MATRIX_RESOLUTIONS}kb.tsv"
        shell:
            """
            mkdir -p {params.threshold_dir}

            printf '\033[1;36mGroup {params.group}, matrix at {params.resolution}kb: getting threshold values for correction...\\n\033[0m'

            madscore=$(grep \"mad threshold \" {input.mad} | sed 's/INFO:hicexplorer.hicCorrectMatrix:mad threshold //g');
            upper=$(echo -3*$madscore | bc);
            echo $madscore \" \" $upper >> {output}
            """
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Matrices correction: correction
    rule N3_summed_matrices_correction__correction:
        input:
            matrix_to_correct = os.path.join("12_Grouped_analyses/B_summed_matrices_normalized/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized.h5"])),
            threshold_file = os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/median_absolute_deviation/thresholds/{GROUPS}_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_thresholdValues.txt")
        output:
            h5_matrix_corrected = os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{GROUPS}/h5_format/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_corrected.h5"]))
        params:
            corrected_matrix_dir = os.path.dirname("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{GROUPS}/h5_format/"),
            correction_method = config["correction_method"],
            group = "{GROUPS}",
            resolution = "{ALL_NEW_MATRIX_RESOLUTIONS}"
        benchmark:
            "benchmarks/N3_summed_matrices_correction__correction/N3_summed_matrices_correction__correction-{GROUPS}_{ALL_NEW_MATRIX_RESOLUTIONS}kb.tsv"
        shell:
            """
            mkdir -p {params.corrected_matrix_dir}

            printf '\033[1;36mGroup {params.group}, matrix at {params.resolution}kb: matrix correction...\\n\033[0m'

            THRESHOLDS=$(head -n 1 {input.threshold_file} | cat);

            hicCorrectMatrix correct \
            --correctionMethod {params.correction_method} \
            --filterThreshold $THRESHOLDS \
            -m {input.matrix_to_correct} \
            -o {output.h5_matrix_corrected} >> {input.threshold_file}
            """
    # ----------------------------------------------------------------------------------------


    # # ----------------------------------------------------------------------------------------
    # # Matrices correction: conversion of the corrected matrices format H5 -> cool
    rule N4_summed_matrices_correction__cool_conversion:
        input:
            h5_matrix_corrected = expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/h5_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.h5"])), group=groups, merged_res = [str(x) for x in NEW_RESOLUTIONS])
        output:
            cool_matrix_corrected = expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/cool_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.cool"])), group=groups, merged_res = [str(x) for x in NEW_RESOLUTIONS])
        params:
            all_names_h5 = ' '.join(expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/h5_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.h5"])), group=groups, merged_res = [str(x) for x in NEW_RESOLUTIONS])),
            all_names_cool =  ' '.join(expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/cool_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.cool"])), group=groups, merged_res = [str(x) for x in NEW_RESOLUTIONS])),
            cool_corrected_matrix_dir = ' '.join(expand(os.path.dirname("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/cool_format/"), group=groups)),
            resolutions = ' '.join([str(x*1000) for x in NEW_RESOLUTIONS]),
            correction_method = str(config["correction_method"])
        benchmark:
            "benchmarks/N4_summed_matrices_correction__cool_conversion/N4_summed_matrices_correction__cool_conversion.tsv"
        shell:
            """
            for DIR in {params.cool_corrected_matrix_dir}
            do
                mkdir -p $DIR
            done

            printf '\033[1;36mConversion of corrected normalized matrices from .h5 to .cool format ...\\n\033[0m'

            hicConvertFormat \
            -m {params.all_names_h5} \
            -o {params.all_names_cool} \
            --correction_name {params.correction_method} \
            --resolutions {params.resolutions} \
            --inputFormat h5 \
            --outputFormat cool
            """


    # # ----------------------------------------------------------------------------------------
    # # Matrices correction: conversion of the corrected matrices format H5 -> HiC-pro
    rule N5_summed_matrices_correction__hicpro_conversion:
        input:
            h5_matrix_corrected = expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/h5_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.h5"])), group=groups, merged_res = [str(x) for x in NEW_RESOLUTIONS])
        output:
            hicpro_matrix_corrected = expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/hicpro_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.hicpro"])), group=groups, merged_res = [str(x) for x in NEW_RESOLUTIONS]),
            hicpro_matrix_corrected_bed = expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/hicpro_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected_hicpro.bed"])), group=groups, merged_res = [str(x) for x in NEW_RESOLUTIONS])
        params:
            all_names_h5 = ' '.join(expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/h5_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.h5"])), group=groups, merged_res = [str(x) for x in NEW_RESOLUTIONS])),
            all_names_hicpro = ' '.join(expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/hicpro_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.hicpro"])), group=groups, merged_res = [str(x) for x in NEW_RESOLUTIONS])),
            all_names_hicpro_bed = ' '.join(expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/hicpro_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected_hicpro.bed"])), group=groups, merged_res = [str(x) for x in NEW_RESOLUTIONS])),
            hicpro_corrected_matrix_dir = ' '.join(expand(os.path.dirname("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/hicpro_format/"), group=groups))
        benchmark:
            "benchmarks/N5_summed_matrices_correction__hicpro_conversion/N5_summed_matrices_correction__hicpro_conversion.tsv"
        shell:
            """
            for DIR in {params.hicpro_corrected_matrix_dir}
            do
                mkdir -p $DIR
            done

            printf '\033[1;36mConversion of corrected matrices from .h5 to .hicpro format ...\\n\033[0m'

            hicConvertFormat \
            -m {params.all_names_h5} \
            -o {params.all_names_hicpro} \
            --bedFileHicpro {params.all_names_hicpro_bed} \
            --inputFormat h5 \
            --outputFormat hicpro
            """
    # ----------------------------------------------------------------------------------------


#    # # ----------------------------------------------------------------------------------------
#    # # Matrices correction: conversion of the corrected matrices format H5 -> HiC-pro
#    rule N6_summed_matrices_correction__NxN_conversion:
#        input:
#            h5_matrix_corrected = expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/h5_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.h5"])), group=groups, merged_res = [str(x) for x in NEW_RESOLUTIONS])
#        output:
#            NxN_matrix_corrected = expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/NxN_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.txt.gz"])), group=groups, merged_res = [str(x) for x in NEW_RESOLUTIONS])
#        params:
#            all_names_h5 = ' '.join(expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/h5_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.h5"])), group=groups, merged_res = [str(x) for x in NEW_RESOLUTIONS])),
#            all_names_NxN = ' '.join(expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/NxN_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{merged_res}kb_normalized_corrected.txt.gz"])), group=groups, merged_res = [str(x) for x in NEW_RESOLUTIONS])),
#            NxN_corrected_matrix_dir = ' '.join(expand(os.path.dirname("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/NxN_format/"), group=groups))
#        benchmark:
#            "benchmarks/N6_summed_matrices_correction__NxN_conversion/N6_summed_matrices_correction__NxN_conversion.tsv"
#        shell:
#            """
#            for DIR in {params.NxN_corrected_matrix_dir}
#            do
#                mkdir -p $DIR
#            done
#
#            printf '\033[1;36mConversion of corrected matrices from .h5 to NxN.txt format ...\\n\033[0m'
#
#            hicConvertFormat \
#            -m {params.all_names_h5} \
#            -o {params.all_names_NxN} \
#            --inputFormat h5 \
#            --outputFormat homer
#            """
#    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Calling TADs summed matrices (HiCexplorer)
    rule O_call_TADs_on_summed_matrices_HiCexplorer:
        input:
            h5_matrix_corrected = os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{GROUPS}/h5_format/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_corrected.h5"]))
        output:
            domains_bed = os.path.join("12_Grouped_analyses/D_TADs_calling_HiCexplorer/{GROUPS}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_domains.bed"]))
        params:
            TADs_dir = os.path.dirname("12_Grouped_analyses/D_TADs_calling_HiCexplorer/{GROUPS}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/log/"),
            prefix = str(os.path.join("12_Grouped_analyses/D_TADs_calling_HiCexplorer/{GROUPS}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb"]))),
            extra_parameters = config["extra_findTAD_parameters"],
            group = "{GROUPS}",
            resolution = "{ALL_NEW_MATRIX_RESOLUTIONS}"
        log:
            out = os.path.join("12_Grouped_analyses/D_TADs_calling_HiCexplorer/{GROUPS}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/log/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_TAD.calling.out"])),
            err = os.path.join("12_Grouped_analyses/D_TADs_calling_HiCexplorer/{GROUPS}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/log/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_TAD.calling.err"]))
        threads:
            math.floor(workflow.cores/4)
        benchmark:
            "benchmarks/O_call_TADs_on_summed_matrices_HiCexplorer/O_call_TADs_on_summed_matrices_HiCexplorer-{GROUPS}_{ALL_NEW_MATRIX_RESOLUTIONS}kb.tsv"
        shell:
            """
            printf '\033[1;36mGroup {params.group}, matrix at {params.resolution}kb: calling TADs (HiCexplorer)...\\n\033[0m'

            mkdir -p {params.TADs_dir}

            hicFindTADs \
            -m {input.h5_matrix_corrected} \
            {params.extra_parameters} \
            --correctForMultipleTesting 'bonferroni' \
            -p {threads} \
            --outPrefix {params.prefix} > {log.out} 2> {log.err}
            """
    # ----------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------
    # Calling TADs (GENOVA)
    rule O_call_TADs_on_summed_matrices_GENOVA:
        input:
            cool_matrix_corrected = os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{GROUPS}/cool_format/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_corrected.cool"]))
        output:
            domains_bed = os.path.join("12_Grouped_analyses/D_TADs_calling_GENOVA/{GROUPS}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_domains.bed"]))
        params:
            TADs_dir = os.path.dirname("12_Grouped_analyses/D_TADs_calling_GENOVA/{GROUPS}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/log/"),
            prefix = str(os.path.join("12_Grouped_analyses/D_TADs_calling_GENOVA/{GROUPS}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb"]))),
            group = "{GROUPS}",
            resolution = "{ALL_NEW_MATRIX_RESOLUTIONS}",
            workdir = str(home_dir)
        log:
            out = os.path.join("12_Grouped_analyses/D_TADs_calling_GENOVA/{GROUPS}/{ALL_NEW_MATRIX_RESOLUTIONS}kb_resolution/log/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_TAD.calling.log"]))
        threads:
            math.floor(workflow.cores/4)
        benchmark:
            "benchmarks/O_call_TADs_on_summed_matrices_GENOVA/O_call_TADs_on_summed_matrices_GENOVA-{GROUPS}_{ALL_NEW_MATRIX_RESOLUTIONS}.tsv"
        shell:
            """
            printf '\033[1;36m{params.group}, matrix at {params.resolution}kb: calling TADs (GENOVA)...\\n\033[0m'

            mkdir -p {params.TADs_dir}

            $CONDA_PREFIX/bin/Rscript \
            -e "matrix=GENOVA::load_contacts(signal_path='{params.workdir}{input.cool_matrix_corrected}',balancing=T)" \
            -e "insulation=GENOVA::insulation_score(explist=matrix,window=25)" \
            -e "saveRDS(insulation,file='{params.workdir}{params.prefix}_GENOVA.insulation.object.Rdata')" \
            -e "options(scipen=999)" \
            -e "write.table(insulation[[1]],file='{params.workdir}{params.prefix}_insulation.score.bedGraph',row.names=F,col.names=F,quote=F)" \
            -e "write.table(GENOVA::call_TAD_insulation(insulation),file='{params.workdir}{params.prefix}_domains.bed',row.names=F,col.names=F,quote=F)" &> {log.out}
            """
    # ----------------------------------------------------------------------------------------


    #                                     **************************


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# **************************************************************************** OPTIONAL ************************************************************************************
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
if (config["loops"]["loop_caller"].lower() == "hicexplorer"):
    if (call_loops == True):
        shell("mkdir -p 09_Loop_detection_HiCexplorer_notPerformed")
        shell("rm -r 09_Loop_detection_HiCexplorer_notPerformed")
        if (eval(str(config["groups"]["perform_grouped_analyses"])) == True):
            shell("mkdir -p 12_Grouped_analyses/E_Loop_detection_HiCexplorer_notPerformed/")
            shell("rm -r 12_Grouped_analyses/E_Loop_detection_HiCexplorer_notPerformed/")

        # ----------------------------------------------------------------------------------------
        # Detect loops (single samples - HiCexplorer)
        rule P_detect_loops_singleSamples_HiCexplorer:
            input:
                h5_matrix_corrected = os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{SAMPLES}/h5_format/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_normalized_corrected.h5"]))
            output:
                loops_bedpe = os.path.join("09_Loop_detection_HiCexplorer/{SAMPLES}/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_loops.bedpe"]))
            params:
                loops_dir = os.path.dirname("09_Loop_detection_HiCexplorer/{SAMPLES}/log/"),
                sample = "{SAMPLES}",
                resolution = "{LOOPS_RESOLUTIONS}",
                maxLoopDistance = config["hicDetectLoops_params"]["maxLoopDistance"],
                loop_windowSize = config["hicDetectLoops_params"]["loop_windowSize"],
                loop_peakWidth = config["hicDetectLoops_params"]["loop_peakWidth"],
                loop_pValuePreselection = config["hicDetectLoops_params"]["loop_pValuePreselection"],
                loop_pValue = config["hicDetectLoops_params"]["loop_pValue"]
            log:
                out = os.path.join("09_Loop_detection_HiCexplorer/{SAMPLES}/log/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_loop.detection.out"])),
                err = os.path.join("09_Loop_detection_HiCexplorer/{SAMPLES}/log/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_loop.detectiong.err"]))
            threads:
                math.floor(workflow.cores/4)
            benchmark:
                "benchmarks/P_detect_loops_singleSamples_HiCexplorer/P_detect_loops_singleSamples_HiCexplorer-{SAMPLES}_{LOOPS_RESOLUTIONS}kb.tsv"
            shell:
                """
                printf '\033[1;36m{params.sample}, matrix at {params.resolution}kb: detecting loops (HiCexplorer)...\\n\033[0m'

                mkdir -p {params.loops_dir}

                hicDetectLoops \
                -m {input.h5_matrix_corrected} \
                -o {output.loops_bedpe} \
                --maxLoopDistance {params.maxLoopDistance} \
                --windowSize {params.loop_windowSize} \
                --peakWidth {params.loop_peakWidth} \
                --pValuePreselection {params.loop_pValuePreselection} \
                --pValue {params.loop_pValue} \
                -t {threads} > {log.out} 2> {log.err}

                if [ ! -f {output.loops_bedpe} ]
                then
                  > {output.loops_bedpe}
                fi
                """
                # ----------------------------------------------------------------------------------------



        # ----------------------------------------------------------------------------------------
        # Detect loops (grouped analyses - HiCexplorer)
        rule Q_detect_loops_groupedSamples_HiCexplorer:
            input:
                h5_matrix_corrected = os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{GROUPS}/h5_format/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_normalized_corrected.h5"]))
            output:
                loops_bedpe = os.path.join("12_Grouped_analyses/E_Loop_detection_HiCexplorer/{GROUPS}/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_loops.bedpe"]))
            params:
                loops_dir = os.path.dirname("12_Grouped_analyses/E_Loop_detection_HiCexplorer/{GROUPS}/log/"),
                sample = "{GROUPS}",
                resolution = "{LOOPS_RESOLUTIONS}",
                maxLoopDistance = config["hicDetectLoops_params"]["maxLoopDistance"],
                loop_windowSize = config["hicDetectLoops_params"]["loop_windowSize"],
                loop_peakWidth = config["hicDetectLoops_params"]["loop_peakWidth"],
                loop_pValuePreselection = config["hicDetectLoops_params"]["loop_pValuePreselection"],
                loop_pValue = config["hicDetectLoops_params"]["loop_pValue"]
            log:
                out = os.path.join("12_Grouped_analyses/E_Loop_detection_HiCexplorer/{GROUPS}/log/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_loop.detection.out"])),
                err = os.path.join("12_Grouped_analyses/E_Loop_detection_HiCexplorer/{GROUPS}/log/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_loop.detection.err"]))
            threads:
                math.floor(workflow.cores/4)
            benchmark:
                "benchmarks/Q_detect_loops_groupedSamples_HiCexplorer/Q_detect_loops_groupedSamples_HiCexplorer-{GROUPS}_{LOOPS_RESOLUTIONS}kb.tsv"
            shell:
                """
                printf '\033[1;36m{params.sample}, matrix at {params.resolution}kb: detecting loops (HiCexplorer)...\\n\033[0m'

                mkdir -p {params.loops_dir}

                hicDetectLoops \
                -m {input.h5_matrix_corrected} \
                -o {output.loops_bedpe} \
                --maxLoopDistance {params.maxLoopDistance} \
                --windowSize {params.loop_windowSize} \
                --peakWidth {params.loop_peakWidth} \
                --pValuePreselection {params.loop_pValuePreselection} \
                --pValue {params.loop_pValue} \
                -t {threads} > {log.out} 2> {log.err}

                if [ ! -f {output.loops_bedpe} ]
                then
                  > {output.loops_bedpe}
                fi
                """
        # ----------------------------------------------------------------------------------------

    else:
        shell("mkdir -p 09_Loop_detection_HiCexplorer_notPerformed/")
        if (eval(str(config["groups"]["perform_grouped_analyses"])) == True):
            shell("mkdir -p 12_Grouped_analyses/E_Loop_detection_HiCexplorer_notPerformed/")
else:
    if (call_loops == True):
        shell("mkdir -p 09_Loop_detection_Mustache_notPerformed")
        shell("rm -r 09_Loop_detection_Mustache_notPerformed")
        if (eval(str(config["groups"]["perform_grouped_analyses"])) == True):
            shell("mkdir -p 12_Grouped_analyses/E_Loop_detection_Mustache_notPerformed/")
            shell("rm -r 12_Grouped_analyses/E_Loop_detection_Mustache_notPerformed/")

        # ----------------------------------------------------------------------------------------
        # Detect loops (single samples - mustache)
        rule P_detect_loops_singleSamples_mustache:
            input:
                cool_matrix_corrected = os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{SAMPLES}/cool_format/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_normalized_corrected.cool"]))
            output:
                loops_bedpe = os.path.join("09_Loop_detection_mustache/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_loops.bedpe"]))
            params:
                loops_dir = os.path.dirname("09_Loop_detection_mustache/log/"),
                sample = "{SAMPLES}",
                resolution = "{LOOPS_RESOLUTIONS}",
                pThreshold = str(config["mustache_params"]["pThreshold"]),
                extra_params = str(config["mustache_params"]["extra_params"])
            log:
                out = os.path.join("09_Loop_detection_mustache/log/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_loop.detection.out"])),
                err = os.path.join("09_Loop_detection_mustache/log/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_loop.detectiong.err"]))
            threads:
                math.floor(workflow.cores/5)
            benchmark:
                "benchmarks/P_detect_loops_singleSamples_mustache/P_detect_loops_singleSamples_mustache-{SAMPLES}_{LOOPS_RESOLUTIONS}kb.tsv"
            shell:
                """
                printf '\033[1;36m{params.sample}, matrix at {params.resolution}kb: detecting loops (mustache)...\\n\033[0m'

                mkdir -p {params.loops_dir}

                mustache \
                -f {input.cool_matrix_corrected} \
                -r {params.resolution}kb \
                --pThreshold {params.pThreshold} \
                -o {output.loops_bedpe} \
                --processes {threads} {params.extra_params} > {log.out} 2> {log.err}

                if [ ! -f {output.loops_bedpe} ]
                then
                  > {output.loops_bedpe}
                fi
                """
                # ----------------------------------------------------------------------------------------


        # ----------------------------------------------------------------------------------------
        # Detect loops (grouped analyses - mustache)
        rule Q_detect_loops_groupedSamples_mustcahe:
            input:
                cool_matrix_corrected = os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{GROUPS}/cool_format/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_normalized_corrected.cool"]))
            output:
                loops_bedpe = os.path.join("12_Grouped_analyses/E_Loop_detection_mustache/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_loops.bedpe"]))
            params:
                loops_dir = os.path.dirname("12_Grouped_analyses/E_Loop_detection_mustache/log/"),
                sample = "{GROUPS}",
                resolution = "{LOOPS_RESOLUTIONS}",
                pThreshold = str(config["mustache_params"]["pThreshold"]),
                extra_params = str(config["mustache_params"]["extra_params"])
            log:
                out = os.path.join("12_Grouped_analyses/E_Loop_detection_mustache/log/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_loop.detection.out"])),
                err = os.path.join("12_Grouped_analyses/E_Loop_detection_mustache/log/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{LOOPS_RESOLUTIONS}kb_loop.detection.err"]))
            threads:
                math.floor(workflow.cores/5)
            benchmark:
                "benchmarks/Q_detect_loops_groupedSamples_mustcahe/Q_detect_loops_groupedSamples_mustcahe-{GROUPS}_{LOOPS_RESOLUTIONS}kb.tsv"
            shell:
                """
                printf '\033[1;36m{params.sample}, matrix at {params.resolution}kb: detecting loops (mustache)...\\n\033[0m'

                mkdir -p {params.loops_dir}

                mustache \
                -f {input.cool_matrix_corrected} \
                -r {params.resolution}kb \
                --pThreshold {params.pThreshold} \
                -o {output.loops_bedpe} \
                --processes {threads} {params.extra_params} > {log.out} 2> {log.err}

                if [ ! -f {output.loops_bedpe} ]
                then
                  > {output.loops_bedpe}
                fi
                """
        # ----------------------------------------------------------------------------------------
    else:
        shell("mkdir -p 09_Loop_detection_mustache_notPerformed/")
        if (eval(str(config["groups"]["perform_grouped_analyses"])) == True):
            shell("mkdir -p 12_Grouped_analyses/E_Loop_detection_mustache_notPerformed/")


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# COMPARTMENTS
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if (perform_compartment_analyses == True):
    shell("mkdir -p 10_Compartments_detection_dcHiC_notPerformed")
    shell("rm -r 10_Compartments_detection_dcHiC_notPerformed")
    if (eval(str(config["groups"]["perform_grouped_analyses"])) == True):
        shell("mkdir -p 12_Grouped_analyses/F_Compartments_detection_dcHiC_notPerformed/")
        shell("rm -r 12_Grouped_analyses/F_Compartments_detection_dcHiC_notPerformed/")

    # ----------------------------------------------------------------------------------------
    # Input file for compartments detections (single samples)
    rule R1_detect_compartments_dcHiC_singleSamples__inputFile_all_vs_all:
        input:
            hicpro_matrix_corrected = expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/hicpro_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{COMPARTMENT_RESOLUTIONS}kb_normalized_corrected.hicpro"])), sample=SAMPLENAMES, allow_missing=True),
            hicpro_matrix_corrected_bed = expand(os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{sample}/hicpro_format/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{COMPARTMENT_RESOLUTIONS}kb_normalized_corrected_hicpro.bed"])), sample=SAMPLENAMES, allow_missing=True)
        output:
            dcHiC_config_file = "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_input_file_individual_samples_{COMPARTMENT_RESOLUTIONS}kb.txt",
            filtered_beds = expand(os.path.join("10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/filtered_hicpro_beds/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_{COMPARTMENT_RESOLUTIONS}kb_normalized_corrected_hicpro_FILTERED.bed"])), sample=SAMPLENAMES, allow_missing=True)
        params:
            working_directory = os.path.dirname(home_dir),
            compartments_dir = os.path.dirname("10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/filtered_hicpro_beds/"),
            chr_filtering_string = config["chr_filtering_string"],
            sample = SAMPLENAMES,
            mapq = str(config["mapQ_cutoff"]),
            resolution = "{COMPARTMENT_RESOLUTIONS}",
            substitution_pattern = config["character_subsitution_dashes_and_points_sample_name"],
            metadata = metadata_table
        benchmark:
            "benchmarks/R1_detect_compartments_dcHiC_singleSamples__inputFile_all_vs_all/R1_detect_compartments_dcHiC_singleSamples__inputFile_all_vs_all-{COMPARTMENT_RESOLUTIONS}kb.tsv"
        run:
            shell("printf '\033[1;36mGeneration of the configuration input file for dcHiC (compartments calling)...\\n\033[0m'")
            shell("mkdir -p {params.compartments_dir}")

            import pandas as pd
            import re
            import os

            dir = os.path.join(params.working_directory, '')

            data_list = []
            for s in params.sample:
                shell("egrep -v -i '{params.chr_filtering_string}' 06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{s}/hicpro_format/{s}_mapQ{params.mapq}_{params.resolution}kb_normalized_corrected_hicpro.bed > 10_Compartments_detection_dcHiC/{params.resolution}kb_resolution/filtered_hicpro_beds/{s}_mapQ{params.mapq}_{params.resolution}kb_normalized_corrected_hicpro_FILTERED.bed")

                data_list.append([
                ''.join([dir, "06_Interaction_matrices_normalized_and_corrected/corrected_matrices/",s,"/hicpro_format/",s,"_mapQ",str(params.mapq),"_",str(params.resolution),"kb_normalized_corrected.hicpro"]),
                ''.join([dir, "10_Compartments_detection_dcHiC/", params.resolution, "kb_resolution/filtered_hicpro_beds/",s,"_mapQ",str(params.mapq),"_",str(params.resolution),"kb_normalized_corrected_hicpro_FILTERED.bed"]),
                ''.join([re.sub("\.|-", str(params.substitution_pattern), s, count=0, flags=0),"_mapQ",str(params.mapq),"_",str(params.resolution),"kb"]),
                re.sub("\.|-", str(params.substitution_pattern), params.metadata[params.metadata.iloc[:,0] == s].iloc[0,1], count=0, flags=0)])

            config_file = pd.DataFrame(data_list, columns = ["mat", "bed", "replicate_prefix", "experiment_prefix"])

            config_file.to_csv(output.dcHiC_config_file, index=False, header=False, sep="\t")
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Analyse compartments for single samples
    rule R2_detect_compartments_dcHiC_singleSamples__call_compartments:
        input:
            dcHiC_config_file = "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_input_file_individual_samples_{COMPARTMENT_RESOLUTIONS}kb.txt"
        output:
            IGV_report = "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/all_vs_all/viz/vizIGV_intra/intra_igv_pcQnm.html"
        params:
            dcHiC_config_file = os.path.join(home_dir, "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_input_file_individual_samples_{COMPARTMENT_RESOLUTIONS}kb.txt"),
            compartments_dir = os.path.dirname(os.path.join(home_dir, "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/")),
            resolution = "{COMPARTMENT_RESOLUTIONS}",
            script_dchicf_file = script_dchicf_file,
            dcHiC_analyses_type = config["compartments"]["dcHiC_analyses_type"],
            genome_name = config["genome_assembly_name"],
            differential_analyses_directory = "all_vs_all",
            home = os.path.join(home_dir,""),
            log_folder = os.path.join(home_dir, "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_logs/")
        threads:
            workflow.cores
        log:
            PCA_log = "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_logs/all_vs_all_single.sample_dcHiC.PCAs.log",
            select_log = "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_logs/all_vs_all_single.sample_dcHiC.select.log",
            analyze_log = "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_logs/all_vs_all_single.sample_dcHiC.analyze.log",
            viz_log = "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_logs/all_vs_all_single.sample_dcHiC.viz.log"
        benchmark:
            "benchmarks/R2_detect_compartments_dcHiC_singleSamples__call_compartments/R2_detect_compartments_dcHiC_singleSamples__call_compartments-{COMPARTMENT_RESOLUTIONS}kb.tsv"
        shell:
            """
            printf '\033[1;36mDetection compartments (single samples) for {params.resolution}kb resultion matrices (dcHiC)...\\n\033[0m'

            mkdir -p {params.log_folder}

            cd {params.compartments_dir}

            printf '\033[1;36m      - Generation of the {params.dcHiC_analyses_type} PCAs...\\n\033[0m'
            $CONDA_PREFIX/bin/Rscript {params.script_dchicf_file} \
            --file {params.dcHiC_config_file} \
            --pcatype {params.dcHiC_analyses_type} \
            --dirovwt T \
            --cthread {threads} \
            --pthread 1 >& {params.home}{log.PCA_log}


            printf '\033[1;36m      - Selection of the best PCA for each group...\\n\033[0m'
            $CONDA_PREFIX/bin/Rscript {params.script_dchicf_file} \
            --file {params.dcHiC_config_file} \
            --pcatype select \
            --dirovwt T \
            --genome {params.genome_name} >& {params.home}{log.select_log}


            printf '\033[1;36m      - Analyses of the compartments...\\n\033[0m'
            $CONDA_PREFIX/bin/Rscript {params.script_dchicf_file} \
            --file {params.dcHiC_config_file} \
            --pcatype analyze \
            --dirovwt T \
            --diffdir {params.differential_analyses_directory} >& {params.home}{log.analyze_log}


            printf '\033[1;36m      - Generate visualization tracks...\\n\033[0m'
            $CONDA_PREFIX/bin/Rscript {params.script_dchicf_file} \
            --file {params.dcHiC_config_file} \
            --pcatype viz \
            --genome {params.genome_name} \
            --diffdir {params.differential_analyses_directory} >& {params.home}{log.viz_log}
            """
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Convert bedGraphs to bigWig (single samples)
    rule R3_detect_compartments_dcHiC_singleSamples__bedGraphToBigWig:
        input:
            IGV_report = "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/all_vs_all/viz/vizIGV_intra/intra_igv_pcQnm.html",
            genome_fai = ''.join([re.sub(".gz", "", config["genome_fasta"], count=0, flags=0),".fai"])
        output:
            chrSizes = temp("10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/all_vs_all/viz/temp_chrSizes_file.txt"),
            bigWigs = expand(''.join(["10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/all_vs_all/viz/files_bigWig/intra_{sample}_mapQ", str(config["mapQ_cutoff"]), "_{COMPARTMENT_RESOLUTIONS}kb_PC.bw"]), sample = [re.sub("\.|-", str(config["character_subsitution_dashes_and_points_sample_name"]), i) for i in SAMPLENAMES], allow_missing=True),
            compartments_bed = expand(''.join(["10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/all_vs_all/viz/files_compartment_beds/intra_{sample}_mapQ", str(config["mapQ_cutoff"]), "_{COMPARTMENT_RESOLUTIONS}kb_PC_compartments_sorted.bed"]), sample = [re.sub("\.|-", str(config["character_subsitution_dashes_and_points_sample_name"]), i) for i in SAMPLENAMES], allow_missing=True)
        params:
            viz_dir = "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/all_vs_all/viz",
            basename_files = ' '.join(expand(''.join(["intra_{sample}_mapQ", str(config["mapQ_cutoff"]), "_{COMPARTMENT_RESOLUTIONS}kb_PC"]), sample = [re.sub("\.|-", str(config["character_subsitution_dashes_and_points_sample_name"]), i) for i in SAMPLENAMES], allow_missing=True)),
            genome_assembly_name = config["genome_assembly_name"]
        benchmark:
            "benchmarks/R3_detect_compartments_dcHiC_singleSamples__bedGraphToBigWig/R3_detect_compartments_dcHiC_singleSamples__bedGraphToBigWig-{COMPARTMENT_RESOLUTIONS}kb.tsv"
        shell:
            """
            printf '\033[1;36mConversion of compartment bedGraph in bigWigs...\\n\033[0m'

            cut -f 1,2 {input.genome_fai} > {output.chrSizes}
            mkdir -p {params.viz_dir}/files_bigWig/
            mkdir -p {params.viz_dir}/files_compartment_beds/

            for w in {params.basename_files}
            do
                $CONDA_PREFIX/bin/bedGraphToBigWig {params.viz_dir}/files/$w.bedGraph {output.chrSizes} {params.viz_dir}/files_bigWig/$w.bw

                VIZDIR={params.viz_dir}

                ### Selecting positive (A comp) and negative (B comp) sides of compartment scores
                grep -v '-' ${{VIZDIR}}/files/$w.bedGraph > ${{VIZDIR}}/files/${{w}}_A.compartment.bedGraph
                grep '-' ${{VIZDIR}}/files/$w.bedGraph > ${{VIZDIR}}/files/${{w}}_B.compartment.bedGraph

                ### Merging compartments bins by type
                # A compartments
                bedtools merge -c 4 -o mean -i ${{VIZDIR}}/files/${{w}}_A.compartment.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph
                cut -f 1-3 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed3
                cut -f 4 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.score
                cut -f 2,3 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.startEnd
                awk '{{print $0 "\\tA"}}' ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed3 > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed4
                paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed4 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.score | awk '{{print $0 "\\t+"}}' > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed6
                paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed6 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.startEnd | awk '{{print $0 "\\t234,100,0"}}' > ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed

                # B compartments
                bedtools merge -c 4 -o mean -i ${{VIZDIR}}/files/${{w}}_B.compartment.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph
                cut -f 1-3 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed3
                cut -f 4 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.score
                cut -f 2,3 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.startEnd
                awk '{{print $0 "\\tB"}}' ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed3 > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed4
                paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed4 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.score | awk '{{print $0 "\\t-"}}' > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed6
                paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed6 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.startEnd | awk '{{print $0 "\\t122,16,180"}}' >> ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed

                # Sorting compartments and clean files
                bedtools sort -i ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed | sed '1s/^/track type=bigBed itemRgb="On"\\n/' > ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments_sorted.bed
                rm ${{VIZDIR}}/files_compartment_beds/${{w}}_*.compartments.bed*
                rm ${{VIZDIR}}/files_compartment_beds/${{w}}_*.compartments.score
                rm ${{VIZDIR}}/files_compartment_beds/${{w}}_*.compartments.startEnd
                rm ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed
            done



            mkdir -p {params.viz_dir}/vizIGV_intra/data_bigWig/
            mkdir -p {params.viz_dir}/vizIGV_intra/data_compartment_beds/

            for w in $(cd {params.viz_dir}/vizIGV_intra/data/; ls *.bedGraph.gz | sed 's/.bedGraph.gz//')
            do
                ### unzipping bedGrpah and converting to bigWig
                VIZDIR={params.viz_dir}
                zcat ${{VIZDIR}}/vizIGV_intra/data/$w.bedGraph.gz > ${{VIZDIR}}/vizIGV_intra/data_bigWig/${{w}}.bedGraph
                $CONDA_PREFIX/bin/bedGraphToBigWig ${{VIZDIR}}/vizIGV_intra/data_bigWig/$w.bedGraph {output.chrSizes} ${{VIZDIR}}/vizIGV_intra/data_bigWig/${{w}}.bw

                if [ $w != "differential_compartment.log10Padj" ]; then
                    if [ $w != "differential_compartment.Mahalanobis" ]; then
                        ### Selecting positive (A comp) and negative (B comp) sides of compartment scores
                        grep -v '-' ${{VIZDIR}}/vizIGV_intra/data_bigWig/$w.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartment.bedGraph
                        grep '-' ${{VIZDIR}}/vizIGV_intra/data_bigWig/$w.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartment.bedGraph

                        ### Merging compartments bins by type
                        # A compartments
                        bedtools merge -c 4 -o mean -i ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartment.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph
                        cut -f 1-3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed3
                        cut -f 4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.score
                        cut -f 2,3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.startEnd
                        awk '{{print $0 "\\tA"}}' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed3 > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed4
                        paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.score | awk '{{print $0 "\\t+"}}' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed6
                        paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed6 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.startEnd | awk '{{print $0 "\\t234,100,0"}}' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed

                        # B compartments
                        bedtools merge -c 4 -o mean -i ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartment.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph
                        cut -f 1-3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed3
                        cut -f 4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.score
                        cut -f 2,3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.startEnd
                        awk '{{print $0 "\\tB"}}' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed3 > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed4
                        paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.score | awk '{{print $0 "\\t-"}}' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed6
                        paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed6 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.startEnd | awk '{{print $0 "\\t122,16,180"}}' >> ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed

                        # Sorting compartments and clean intermediary files
                        bedtools sort -i ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed | sed '1s/^/track type=bigBed itemRgb="On"\\n/' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments_sorted.bed
                        rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/*.bedGraph
                        rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_*.compartments.bed*
                        rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_*.compartments.score
                        rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_*.compartments.startEnd
                        rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed
                    fi
                fi
                rm ${{VIZDIR}}/vizIGV_intra/data_bigWig/${{w}}.bedGraph
            done
            """
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Analyse compartments for single samples -- COMBINATIONS
    rule R4_detect_compartments_dcHiC_singleSamples__call_compartments_combos:
        input:
            dcHiC_config_file = "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_input_file_individual_samples_{COMPARTMENT_RESOLUTIONS}kb.txt",
            IGV_report_all_vs_all = "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/all_vs_all/viz/vizIGV_intra/intra_igv_pcQnm.html"
        output:
            IGV_report_combo = expand("10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/{combo}/viz/vizIGV_intra/intra_igv_pcQnm.html", combo = combo_list, allow_missing=True)
        params:
            dcHiC_config_file = os.path.join(home_dir, "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_input_file_individual_samples_{COMPARTMENT_RESOLUTIONS}kb.txt"),
            comboA = ''.join(['"', '" "'.join(comboA), '"']),
            comboB = ''.join(['"', '" "'.join(comboB), '"']),
            all_combos = ''.join(['"', '" "'.join(combo_list), '"']),
            compartments_dir = os.path.dirname(os.path.join(home_dir, "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/")),
            genome_name = config["genome_assembly_name"],
            resolution = "{COMPARTMENT_RESOLUTIONS}",
            script_dchicf_file = script_dchicf_file,
            log_folder = os.path.join(home_dir, "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_logs/")
        threads:
            workflow.cores
        benchmark:
            "benchmarks/R4_detect_compartments_dcHiC_singleSamples__call_compartments_combos/R4_detect_compartments_dcHiC_singleSamples__call_compartments_combos-{COMPARTMENT_RESOLUTIONS}kb.tsv"
        shell:
            """
            printf '\033[1;36m{params.resolution}kb: Performig compartmentalization analyses for paired combinations (dcHiC)...\\n\033[0m'

            mkdir -p {params.log_folder}

            comboA=({params.comboA})
            comboB=({params.comboB})
            combolist=({params.all_combos})

            for i in $(seq 0 $[${{#combolist[@]}} - 1])
            do
                configFile={params.compartments_dir}/dcHiC_input_file_individual_samples_{params.resolution}kb_${{combolist[i]}}.txt
                grep -w ${{comboA[i]}} {params.dcHiC_config_file} > $configFile
                grep -w ${{comboB[i]}} {params.dcHiC_config_file} >> $configFile

                cd {params.compartments_dir}

                $CONDA_PREFIX/bin/Rscript {params.script_dchicf_file} \
                --file $configFile \
                --pcatype analyze \
                --dirovwt T \
                --diffdir ${{combolist[i]}} >& {params.log_folder}${{comboA[i]}}_${{comboB[i]}}_dcHiC.analyze.log

                $CONDA_PREFIX/bin/Rscript {params.script_dchicf_file} \
                --file $configFile \
                --pcatype viz \
                --genome {params.genome_name} \
                --diffdir ${{combolist[i]}} >& {params.log_folder}${{comboA[i]}}_${{comboB[i]}}_dcHiC.viz.log
            done
            """
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Analyse compartments for single samples -- COMBINATIONS
    rule R5_detect_compartments_dcHiC_singleSamples__bedGraphToBigWig_combos:
        input:
            IGV_report_combo = expand("10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/{combo}/viz/vizIGV_intra/intra_igv_pcQnm.html", combo = combo_list, allow_missing=True),
            genome_fai = ''.join([re.sub(".gz", "", config["genome_fasta"], count=0, flags=0),".fai"])
        output:
            bigWig_single_combos = expand("10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/{combo}/viz/vizIGV_intra/data_bigWig/differential_compartment.Mahalanobis.bw", combo = combo_list, allow_missing=True),
            chrSizes = temp("10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/temp_chrSizes_file.txt")
        params:
            comboA = ''.join(['"', '" "'.join(comboA), '"']),
            comboB = ''.join(['"', '" "'.join(comboB), '"']),
            all_combos = ''.join(['"', '" "'.join(combo_list), '"']),
            compartments_dir = os.path.dirname(os.path.join(home_dir, "10_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/")),
            genome_name = config["genome_assembly_name"],
            resolution = "{COMPARTMENT_RESOLUTIONS}"
        benchmark:
            "benchmarks/R5_detect_compartments_dcHiC_singleSamples__bedGraphToBigWig_combos/R5_detect_compartments_dcHiC_singleSamples__bedGraphToBigWig_combos-{COMPARTMENT_RESOLUTIONS}kb.tsv"
        shell:
            """
            printf '\033[1;36m{params.resolution}kb: Converting compartment scores and stats bedGraphs to BigWigs for paired combinations (dcHiC)...\\n\033[0m'

            comboA=({params.comboA})
            comboB=({params.comboB})
            combolist=({params.all_combos})

            cut -f 1,2 {input.genome_fai} > {output.chrSizes}

            for i in $(seq 0 $[${{#combolist[@]}} - 1])
            do
                VIZDIR={params.compartments_dir}/DifferentialResult/${{combolist[i]}}/viz

                mkdir -p ${{VIZDIR}}/files_bigWig/
                mkdir -p ${{VIZDIR}}/files_compartment_beds/
                for w in $(cd ${{VIZDIR}}/files/; ls *{params.resolution}kb_PC.bedGraph | sed 's/.bedGraph//')
                do
                    $CONDA_PREFIX/bin/bedGraphToBigWig ${{VIZDIR}}/files/$w.bedGraph {output.chrSizes} ${{VIZDIR}}/files_bigWig/$w.bw

                    ### Selecting positive (A comp) and negative (B comp) sides of compartment scores
                    grep -v '-' ${{VIZDIR}}/files/$w.bedGraph > ${{VIZDIR}}/files/${{w}}_A.compartment.bedGraph
                    grep '-' ${{VIZDIR}}/files/$w.bedGraph > ${{VIZDIR}}/files/${{w}}_B.compartment.bedGraph

                    ### Merging compartments bins by type
                    # A compartments
                    bedtools merge -c 4 -o mean -i ${{VIZDIR}}/files/${{w}}_A.compartment.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph
                    cut -f 1-3 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed3
                    cut -f 4 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.score
                    cut -f 2,3 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.startEnd
                    awk '{{print $0 "\\tA"}}' ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed3 > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed4
                    paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed4 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.score | awk '{{print $0 "\\t+"}}' > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed6
                    paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed6 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.startEnd | awk '{{print $0 "\\t234,100,0"}}' > ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed

                    # B compartments
                    bedtools merge -c 4 -o mean -i ${{VIZDIR}}/files/${{w}}_B.compartment.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph
                    cut -f 1-3 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed3
                    cut -f 4 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.score
                    cut -f 2,3 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.startEnd
                    awk '{{print $0 "\\tB"}}' ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed3 > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed4
                    paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed4 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.score | awk '{{print $0 "\\t-"}}' > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed6
                    paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed6 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.startEnd | awk '{{print $0 "\\t122,16,180"}}' >> ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed

                    # Sorting compartments and clean files
                    bedtools sort -i ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed | sed '1s/^/track type=bigBed itemRgb="On"\\n/' > ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments_sorted.bed
                    rm ${{VIZDIR}}/files_compartment_beds/${{w}}_*.compartments.bed*
                    rm ${{VIZDIR}}/files_compartment_beds/${{w}}_*.compartments.score
                    rm ${{VIZDIR}}/files_compartment_beds/${{w}}_*.compartments.startEnd
                    rm ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed
                done


                mkdir -p ${{VIZDIR}}/vizIGV_intra/data_bigWig/
                mkdir -p ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/
                for w in $(cd ${{VIZDIR}}/vizIGV_intra/data/; ls *.bedGraph.gz | sed 's/.bedGraph.gz//')
                do
                    ### unzipping bedGrpah and converting to bigWig
                    zcat ${{VIZDIR}}/vizIGV_intra/data/$w.bedGraph.gz > ${{VIZDIR}}/vizIGV_intra/data_bigWig/${{w}}.bedGraph
                    $CONDA_PREFIX/bin/bedGraphToBigWig ${{VIZDIR}}/vizIGV_intra/data_bigWig/$w.bedGraph {output.chrSizes} ${{VIZDIR}}/vizIGV_intra/data_bigWig/${{w}}.bw

                    if [ $w != "differential_compartment.log10Padj" ]; then
                        if [ $w != "differential_compartment.Mahalanobis" ]; then
                            ### Selecting positive (A comp) and negative (B comp) sides of compartment scores
                            grep -v '-' ${{VIZDIR}}/vizIGV_intra/data_bigWig/$w.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartment.bedGraph
                            grep '-' ${{VIZDIR}}/vizIGV_intra/data_bigWig/$w.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartment.bedGraph

                            ### Merging compartments bins by type
                            # A compartments
                            bedtools merge -c 4 -o mean -i ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartment.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph
                            cut -f 1-3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed3
                            cut -f 4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.score
                            cut -f 2,3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.startEnd
                            awk '{{print $0 "\\tA"}}' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed3 > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed4
                            paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.score | awk '{{print $0 "\\t+"}}' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed6
                            paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed6 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.startEnd | awk '{{print $0 "\\t234,100,0"}}' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed

                            # B compartments
                            bedtools merge -c 4 -o mean -i ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartment.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph
                            cut -f 1-3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed3
                            cut -f 4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.score
                            cut -f 2,3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.startEnd
                            awk '{{print $0 "\\tB"}}' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed3 > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed4
                            paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.score | awk '{{print $0 "\\t-"}}' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed6
                            paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed6 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.startEnd | awk '{{print $0 "\\t122,16,180"}}' >> ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed

                            # Sorting compartments and clean intermediary files
                            bedtools sort -i ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed | sed '1s/^/track type=bigBed itemRgb="On"\\n/' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments_sorted.bed
                            rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/*.bedGraph
                            rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_*.compartments.bed*
                            rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_*.compartments.score
                            rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_*.compartments.startEnd
                            rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed
                        fi
                    fi
                    rm ${{VIZDIR}}/vizIGV_intra/data_bigWig/${{w}}.bedGraph
                done

            done
            """
    # ----------------------------------------------------------------------------------------


    # @@@@@@@@@ GROUPED ANALYSES COMPARTMENTS @@@@@@@@@
    # ----------------------------------------------------------------------------------------
    # Input file for compartments detections (grouped samples)
    rule S1_detect_compartments_dcHiC_groupedSamples__inputFile_all_vs_all:
        input:
            hicpro_matrix_corrected = expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/hicpro_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{COMPARTMENT_RESOLUTIONS}kb_normalized_corrected.hicpro"])), group=groups, allow_missing=True),
            hicpro_matrix_corrected_bed = expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/hicpro_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{COMPARTMENT_RESOLUTIONS}kb_normalized_corrected_hicpro.bed"])), group=groups, allow_missing=True)
        output:
            dcHiC_config_file = "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_input_file_grouped_samples_{COMPARTMENT_RESOLUTIONS}kb.txt",
            filtered_beds = expand(os.path.join("12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/filtered_hicpro_beds/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{COMPARTMENT_RESOLUTIONS}kb_normalized_corrected_hicpro_FILTERED.bed"])), group=groups, allow_missing=True)
        params:
            working_directory = os.path.dirname(home_dir),
            compartments_dir = os.path.dirname("12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/filtered_hicpro_beds/"),
            chr_filtering_string = config["chr_filtering_string"],
            group = groups,
            mapq = str(config["mapQ_cutoff"]),
            resolution = "{COMPARTMENT_RESOLUTIONS}",
            substitution_pattern = config["character_subsitution_dashes_and_points_sample_name"],
            metadata = metadata_table
        benchmark:
            "benchmarks/S1_detect_compartments_dcHiC_groupedSamples__inputFile_all_vs_all/S1_detect_compartments_dcHiC_groupedSamples__inputFile_all_vs_all-{COMPARTMENT_RESOLUTIONS}kb.tsv"
        run:
            shell("printf '\033[1;36mGeneration of the configuration input file for dcHiC (compartments calling)...\\n\033[0m'")
            shell("mkdir -p {params.compartments_dir}")

            import pandas as pd
            import re
            import os

            dir = os.path.join(params.working_directory, '')

            data_list = []
            for s in params.group:
                shell("egrep -v -i '{params.chr_filtering_string}' 12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{s}/hicpro_format/{s}_mapQ{params.mapq}_{params.resolution}kb_normalized_corrected_hicpro.bed > 12_Grouped_analyses/F_Compartments_detection_dcHiC/{params.resolution}kb_resolution/filtered_hicpro_beds/{s}_mapQ{params.mapq}_{params.resolution}kb_normalized_corrected_hicpro_FILTERED.bed")

                data_list.append([
                ''.join([dir, "12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/",s,"/hicpro_format/",s,"_mapQ",str(params.mapq),"_",str(params.resolution),"kb_normalized_corrected.hicpro"]),
                ''.join([dir, "12_Grouped_analyses/F_Compartments_detection_dcHiC/", params.resolution, "kb_resolution/filtered_hicpro_beds/",s,"_mapQ",str(params.mapq),"_",str(params.resolution),"kb_normalized_corrected_hicpro_FILTERED.bed"]),
                ''.join([re.sub("\.|-", str(params.substitution_pattern), s, count=0, flags=0),"_mapQ",str(params.mapq),"_",str(params.resolution),"kb"]),
                re.sub("\.|-", str(params.substitution_pattern), s)])

            config_file = pd.DataFrame(data_list, columns = ["mat", "bed", "replicate_prefix", "experiment_prefix"])

            config_file.to_csv(output.dcHiC_config_file, index=False, header=False, sep="\t")
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Analyse compartments for grouped samples
    rule S2_detect_compartments_dcHiC_groupedSamples__call_compartments:
        input:
            dcHiC_config_file = "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_input_file_grouped_samples_{COMPARTMENT_RESOLUTIONS}kb.txt"
        output:
            IGV_report = "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/all_vs_all/viz/vizIGV_intra/intra_igv_pcQnm.html"
        params:
            dcHiC_config_file = os.path.join(home_dir, "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_input_file_grouped_samples_{COMPARTMENT_RESOLUTIONS}kb.txt"),
            compartments_dir = os.path.dirname(os.path.join(home_dir, "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/")),
            resolution = "{COMPARTMENT_RESOLUTIONS}",
            script_dchicf_file = script_dchicf_file,
            dcHiC_analyses_type = config["compartments"]["dcHiC_analyses_type"],
            genome_name = config["genome_assembly_name"],
            differential_analyses_directory = "all_vs_all",
            home = os.path.join(home_dir, ""),
            log_folder = os.path.join(home_dir, "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_logs/")
        log:
            PCA_log = "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_logs/all_vs_all_grouped__dcHiC.PCAs.log",
            select_log = "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_logs/all_vs_all_grouped_dcHiC.select.log",
            analyze_log = "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_logs/all_vs_all_grouped_dcHiC.analyze.log",
            viz_log = "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_logs/all_vs_all_grouped_dcHiC.viz.log"
        threads:
            workflow.cores
        benchmark:
            "benchmarks/S2_detect_compartments_dcHiC_groupedSamples__call_compartments/S2_detect_compartments_dcHiC_groupedSamples__call_compartments-{COMPARTMENT_RESOLUTIONS}kb.tsv"
        shell:
            """
            printf '\033[1;36mDetection compartments (grouped samples) for {params.resolution}kb resultion matrices (dcHiC)...\\n\033[0m'

            mkdir -p {params.log_folder}

            cd {params.compartments_dir}

            printf '\033[1;36m      - Generation of the {params.dcHiC_analyses_type} PCAs...\\n\033[0m'
            $CONDA_PREFIX/bin/Rscript {params.script_dchicf_file} \
            --file {params.dcHiC_config_file} \
            --pcatype {params.dcHiC_analyses_type} \
            --dirovwt T \
            --cthread {threads} \
            --pthread 1 >& {params.home}{log.PCA_log}


            printf '\033[1;36m      - Selection of the best PCA for each group...\\n\033[0m'
            $CONDA_PREFIX/bin/Rscript {params.script_dchicf_file} \
            --file {params.dcHiC_config_file} \
            --pcatype select \
            --dirovwt T \
            --genome {params.genome_name} >& {params.home}{log.select_log}


            printf '\033[1;36m      - Analyses of the compartments...\\n\033[0m'
            $CONDA_PREFIX/bin/Rscript {params.script_dchicf_file} \
            --file {params.dcHiC_config_file} \
            --pcatype analyze \
            --dirovwt T \
            --diffdir {params.differential_analyses_directory} >& {params.home}{log.analyze_log}


            printf '\033[1;36m      - Generate visualization tracks...\\n\033[0m'
            $CONDA_PREFIX/bin/Rscript {params.script_dchicf_file} \
            --file {params.dcHiC_config_file} \
            --pcatype viz \
            --genome {params.genome_name} \
            --diffdir {params.differential_analyses_directory} >& {params.home}{log.viz_log}
            """
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Convert bedGraphs to bigWig (grouped analyses)
    rule S3_detect_compartments_dcHiC_groupedSamples__bedGraphToBigWig:
        input:
            IGV_report = "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/all_vs_all/viz/vizIGV_intra/intra_igv_pcQnm.html",
            genome_fai = ''.join([re.sub(".gz", "", config["genome_fasta"], count=0, flags=0),".fai"])
        output:
            chrSizes = temp("12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/all_vs_all/viz/temp_chrSizes_file.txt"),
            bigWigs = expand("12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/all_vs_all/viz/vizIGV_intra/data_bigWig/{group}.PC.bw", group = [re.sub("\.|-", str(config["character_subsitution_dashes_and_points_sample_name"]), i) for i in groups], allow_missing=True),
            compartments_bed = expand(''.join(["12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/all_vs_all/viz/files_compartment_beds/intra_{group}_mapQ", str(config["mapQ_cutoff"]), "_{COMPARTMENT_RESOLUTIONS}kb_PC_compartments_sorted.bed"]), group = [re.sub("\.|-", str(config["character_subsitution_dashes_and_points_sample_name"]), i) for i in groups], allow_missing=True)
        params:
            viz_dir = "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/all_vs_all/viz",
            basename_files = ' '.join(expand(''.join(["intra_{group}_mapQ", str(config["mapQ_cutoff"]), "_{COMPARTMENT_RESOLUTIONS}kb_PC"]), group = [re.sub("\.|-", str(config["character_subsitution_dashes_and_points_sample_name"]), i) for i in groups], allow_missing=True)),
            genome_assembly_name = config["genome_assembly_name"]
        benchmark:
            "benchmarks/S3_detect_compartments_dcHiC_groupedSamples__bedGraphToBigWig/S3_detect_compartments_dcHiC_groupedSamples__bedGraphToBigWig-{COMPARTMENT_RESOLUTIONS}kb.tsv"
        shell:
            """
            printf '\033[1;36mConversion of compartment bedGraph in bigWigs...\\n\033[0m'

            cut -f 1,2 {input.genome_fai} > {output.chrSizes}
            mkdir -p {params.viz_dir}/files_bigWig/
            mkdir -p {params.viz_dir}/files_compartment_beds/

            for w in {params.basename_files}
            do
                $CONDA_PREFIX/bin/bedGraphToBigWig {params.viz_dir}/files/$w.bedGraph {output.chrSizes} {params.viz_dir}/files_bigWig/$w.bw

                VIZDIR={params.viz_dir}

                ### Selecting positive (A comp) and negative (B comp) sides of compartment scores
                grep -v '-' ${{VIZDIR}}/files/$w.bedGraph > ${{VIZDIR}}/files/${{w}}_A.compartment.bedGraph
                grep '-' ${{VIZDIR}}/files/$w.bedGraph > ${{VIZDIR}}/files/${{w}}_B.compartment.bedGraph

                ### Merging compartments bins by type
                # A compartments
                bedtools merge -c 4 -o mean -i ${{VIZDIR}}/files/${{w}}_A.compartment.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph
                cut -f 1-3 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed3
                cut -f 4 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.score
                cut -f 2,3 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.startEnd
                awk '{{print $0 "\\tA"}}' ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed3 > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed4
                paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed4 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.score | awk '{{print $0 "\\t+"}}' > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed6
                paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed6 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.startEnd | awk '{{print $0 "\\t234,100,0"}}' > ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed

                # B compartments
                bedtools merge -c 4 -o mean -i ${{VIZDIR}}/files/${{w}}_B.compartment.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph
                cut -f 1-3 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed3
                cut -f 4 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.score
                cut -f 2,3 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.startEnd
                awk '{{print $0 "\\tB"}}' ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed3 > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed4
                paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed4 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.score | awk '{{print $0 "\\t-"}}' > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed6
                paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed6 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.startEnd | awk '{{print $0 "\\t122,16,180"}}' >> ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed

                # Sorting compartments and clean files
                bedtools sort -i ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed | sed '1s/^/track type=bigBed itemRgb="On"\\n/' > ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments_sorted.bed
                rm ${{VIZDIR}}/files_compartment_beds/${{w}}_*.compartments.bed*
                rm ${{VIZDIR}}/files_compartment_beds/${{w}}_*.compartments.score
                rm ${{VIZDIR}}/files_compartment_beds/${{w}}_*.compartments.startEnd
                rm ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed
            done



            mkdir -p {params.viz_dir}/vizIGV_intra/data_bigWig/
            mkdir -p {params.viz_dir}/vizIGV_intra/data_compartment_beds/

            for w in $(cd {params.viz_dir}/vizIGV_intra/data/; ls *.bedGraph.gz | sed 's/.bedGraph.gz//')
            do
                ### unzipping bedGrpah and converting to bigWig
                VIZDIR={params.viz_dir}
                zcat ${{VIZDIR}}/vizIGV_intra/data/$w.bedGraph.gz > ${{VIZDIR}}/vizIGV_intra/data_bigWig/${{w}}.bedGraph
                $CONDA_PREFIX/bin/bedGraphToBigWig ${{VIZDIR}}/vizIGV_intra/data_bigWig/$w.bedGraph {output.chrSizes} ${{VIZDIR}}/vizIGV_intra/data_bigWig/${{w}}.bw

                if [ $w != "differential_compartment.log10Padj" ]; then
                    if [ $w != "differential_compartment.Mahalanobis" ]; then
                        ### Selecting positive (A comp) and negative (B comp) sides of compartment scores
                        grep -v '-' ${{VIZDIR}}/vizIGV_intra/data_bigWig/$w.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartment.bedGraph
                        grep '-' ${{VIZDIR}}/vizIGV_intra/data_bigWig/$w.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartment.bedGraph

                        ### Merging compartments bins by type
                        # A compartments
                        bedtools merge -c 4 -o mean -i ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartment.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph
                        cut -f 1-3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed3
                        cut -f 4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.score
                        cut -f 2,3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.startEnd
                        awk '{{print $0 "\\tA"}}' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed3 > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed4
                        paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.score | awk '{{print $0 "\\t+"}}' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed6
                        paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed6 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.startEnd | awk '{{print $0 "\\t234,100,0"}}' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed

                        # B compartments
                        bedtools merge -c 4 -o mean -i ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartment.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph
                        cut -f 1-3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed3
                        cut -f 4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.score
                        cut -f 2,3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.startEnd
                        awk '{{print $0 "\\tB"}}' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed3 > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed4
                        paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.score | awk '{{print $0 "\\t-"}}' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed6
                        paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed6 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.startEnd | awk '{{print $0 "\\t122,16,180"}}' >> ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed

                        # Sorting compartments and clean intermediary files
                        bedtools sort -i ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed | sed '1s/^/track type=bigBed itemRgb="On"\\n/' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments_sorted.bed
                        rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/*.bedGraph
                        rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_*.compartments.bed*
                        rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_*.compartments.score
                        rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_*.compartments.startEnd
                        rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed
                    fi
                fi
                rm ${{VIZDIR}}/vizIGV_intra/data_bigWig/${{w}}.bedGraph
            done
            """
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Analyse compartments for grouped samples -- COMBINATIONS
    rule S4_detect_compartments_dcHiC_groupedSamples__call_compartments_combos:
        input:
            dcHiC_config_file = "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_input_file_grouped_samples_{COMPARTMENT_RESOLUTIONS}kb.txt",
            IGV_report_all_vs_all = "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/all_vs_all/viz/vizIGV_intra/intra_igv_pcQnm.html"
        output:
            IGV_report_combo = expand("12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/{combo}/viz/vizIGV_intra/intra_igv_pcQnm.html", combo = combo_list, allow_missing=True)
        params:
            dcHiC_config_file = os.path.join(home_dir, "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_input_file_grouped_samples_{COMPARTMENT_RESOLUTIONS}kb.txt"),
            comboA = ''.join(['"', '" "'.join(comboA), '"']),
            comboB = ''.join(['"', '" "'.join(comboB), '"']),
            all_combos = ''.join(['"', '" "'.join(combo_list), '"']),
            compartments_dir = os.path.dirname(os.path.join(home_dir, "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/")),
            genome_name = config["genome_assembly_name"],
            resolution = "{COMPARTMENT_RESOLUTIONS}",
            script_dchicf_file = script_dchicf_file,
            log_folder = os.path.join(home_dir, "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/dcHiC_logs/")
        threads:
            workflow.cores
        benchmark:
            "benchmarks/S4_detect_compartments_dcHiC_groupedSamples__call_compartments_combos/S4_detect_compartments_dcHiC_groupedSamples__call_compartments_combos-{COMPARTMENT_RESOLUTIONS}kb.tsv"
        shell:
            """
            printf '\033[1;36m{params.resolution}kb: Performig compartmentalization analyses for paired combinations (dcHiC)...\\n\033[0m'

            comboA=({params.comboA})
            comboB=({params.comboB})
            combolist=({params.all_combos})

            mkdir -p {params.log_folder}

            for i in $(seq 0 $[${{#combolist[@]}} - 1])
            do
                configFile={params.compartments_dir}/dcHiC_input_file_individual_samples_{params.resolution}kb_${{combolist[i]}}.txt
                grep -w ${{comboA[i]}} {params.dcHiC_config_file} > $configFile
                grep -w ${{comboB[i]}} {params.dcHiC_config_file} >> $configFile

                cd {params.compartments_dir}

                $CONDA_PREFIX/bin/Rscript {params.script_dchicf_file} \
                --file $configFile \
                --pcatype analyze \
                --dirovwt T \
                --diffdir ${{combolist[i]}} >& {params.log_folder}${{comboA[i]}}_${{comboB[i]}}_grouped_dcHiC.analyze.log

                $CONDA_PREFIX/bin/Rscript {params.script_dchicf_file} \
                --file $configFile \
                --pcatype viz \
                --genome {params.genome_name} \
                --diffdir ${{combolist[i]}} >& {params.log_folder}${{comboA[i]}}_${{comboB[i]}}_grouped_dcHiC.viz.log
            done
            """
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Analyse compartments for single grouped -- COMBINATIONS
    rule S5_detect_compartments_dcHiC_groupedSamples__bedGraphToBigWig_combos:
        input:
            IGV_report_combo = expand("12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/{combo}/viz/vizIGV_intra/intra_igv_pcQnm.html", combo = combo_list, allow_missing=True),
            genome_fai = ''.join([re.sub(".gz", "", config["genome_fasta"], count=0, flags=0),".fai"])
        output:
            bigWig_single_combos = expand("12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/DifferentialResult/{combo}/viz/vizIGV_intra/data_bigWig/differential_compartment.Mahalanobis.bw", combo = combo_list, allow_missing=True),
            chrSizes = temp("12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/temp_chrSizes_file.txt")
        params:
            comboA = ''.join(['"', '" "'.join(comboA), '"']),
            comboB = ''.join(['"', '" "'.join(comboB), '"']),
            all_combos = ''.join(['"', '" "'.join(combo_list), '"']),
            compartments_dir = os.path.dirname(os.path.join(home_dir, "12_Grouped_analyses/F_Compartments_detection_dcHiC/{COMPARTMENT_RESOLUTIONS}kb_resolution/")),
            genome_name = config["genome_assembly_name"],
            resolution = "{COMPARTMENT_RESOLUTIONS}"
        benchmark:
            "benchmarks/S5_detect_compartments_dcHiC_groupedSamples__bedGraphToBigWig_combos/S5_detect_compartments_dcHiC_groupedSamples__bedGraphToBigWig_combos-{COMPARTMENT_RESOLUTIONS}kb.tsv"
        shell:
            """
            printf '\033[1;36m{params.resolution}kb: Converting compartment scores and stats bedGraphs to BigWigs for paired combinations (dcHiC) and\\n\033[0m'
            printf '\033[1;36mgeneration of compartment bed files...\\n\033[0m'

            comboA=({params.comboA})
            comboB=({params.comboB})
            combolist=({params.all_combos})

            cut -f 1,2 {input.genome_fai} > {output.chrSizes}

            for i in $(seq 0 $[${{#combolist[@]}} - 1])
            do
                VIZDIR={params.compartments_dir}/DifferentialResult/${{combolist[i]}}/viz

                mkdir -p ${{VIZDIR}}/files_bigWig/
                mkdir -p ${{VIZDIR}}/files_compartment_beds/
                for w in $(cd ${{VIZDIR}}/files/; ls *{params.resolution}kb_PC.bedGraph | sed 's/.bedGraph//')
                do
                    ### Coverting single sample compartment scores to bedGraph
                    $CONDA_PREFIX/bin/bedGraphToBigWig ${{VIZDIR}}/files/$w.bedGraph {output.chrSizes} ${{VIZDIR}}/files_bigWig/$w.bw

                    ### Selecting positive (A comp) and negative (B comp) sides of compartment scores
                    grep -v '-' ${{VIZDIR}}/files/$w.bedGraph > ${{VIZDIR}}/files/${{w}}_A.compartment.bedGraph
                    grep '-' ${{VIZDIR}}/files/$w.bedGraph > ${{VIZDIR}}/files/${{w}}_B.compartment.bedGraph

                    ### Merging compartments bins by type
                    # A compartments
                    bedtools merge -c 4 -o mean -i ${{VIZDIR}}/files/${{w}}_A.compartment.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph
                    cut -f 1-3 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed3
                    cut -f 4 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.score
                    cut -f 2,3 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.startEnd
                    awk '{{print $0 "\\tA"}}' ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed3 > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed4
                    paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed4 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.score | awk '{{print $0 "\\t+"}}' > ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed6
                    paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.bed6 ${{VIZDIR}}/files_compartment_beds/${{w}}_A.compartments.startEnd | awk '{{print $0 "\\t234,100,0"}}' > ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed

                    # B compartments
                    bedtools merge -c 4 -o mean -i ${{VIZDIR}}/files/${{w}}_B.compartment.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph
                    cut -f 1-3 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed3
                    cut -f 4 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.score
                    cut -f 2,3 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.startEnd
                    awk '{{print $0 "\\tB"}}' ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed3 > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed4
                    paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed4 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.score | awk '{{print $0 "\\t-"}}' > ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed6
                    paste -d '\\t' ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.bed6 ${{VIZDIR}}/files_compartment_beds/${{w}}_B.compartments.startEnd | awk '{{print $0 "\\t122,16,180"}}' >> ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed

                    # Sorting compartments and clean files
                    bedtools sort -i ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed | sed '1s/^/track type=bigBed itemRgb="On"\\n/' > ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments_sorted.bed
                    rm ${{VIZDIR}}/files_compartment_beds/${{w}}_*.compartments.bed*
                    rm ${{VIZDIR}}/files_compartment_beds/${{w}}_*.compartments.score
                    rm ${{VIZDIR}}/files_compartment_beds/${{w}}_*.compartments.startEnd
                    rm ${{VIZDIR}}/files_compartment_beds/${{w}}_compartments.bed
                done



                mkdir -p ${{VIZDIR}}/vizIGV_intra/data_bigWig/
                mkdir -p ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/

                for w in $(cd ${{VIZDIR}}/vizIGV_intra/data/; ls *.bedGraph.gz | sed 's/.bedGraph.gz//')
                do
                    ### unzipping bedGrpah and converting to bigWig
                    zcat ${{VIZDIR}}/vizIGV_intra/data/$w.bedGraph.gz > ${{VIZDIR}}/vizIGV_intra/data_bigWig/${{w}}.bedGraph
                    $CONDA_PREFIX/bin/bedGraphToBigWig ${{VIZDIR}}/vizIGV_intra/data_bigWig/$w.bedGraph {output.chrSizes} ${{VIZDIR}}/vizIGV_intra/data_bigWig/${{w}}.bw

                    if [ $w != "differential_compartment.log10Padj" ]; then
                        if [ $w != "differential_compartment.Mahalanobis" ]; then
                            ### Selecting positive (A comp) and negative (B comp) sides of compartment scores
                            grep -v '-' ${{VIZDIR}}/vizIGV_intra/data_bigWig/$w.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartment.bedGraph
                            grep '-' ${{VIZDIR}}/vizIGV_intra/data_bigWig/$w.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartment.bedGraph

                            ### Merging compartments bins by type
                            # A compartments
                            bedtools merge -c 4 -o mean -i ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartment.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph
                            cut -f 1-3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed3
                            cut -f 4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.score
                            cut -f 2,3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.startEnd
                            awk '{{print $0 "\\tA"}}' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed3 > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed4
                            paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.score | awk '{{print $0 "\\t+"}}' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed6
                            paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.bed6 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_A.compartments.startEnd | awk '{{print $0 "\\t234,100,0"}}' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed

                            # B compartments
                            bedtools merge -c 4 -o mean -i ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartment.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph
                            cut -f 1-3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed3
                            cut -f 4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.score
                            cut -f 2,3 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bedGraph > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.startEnd
                            awk '{{print $0 "\\tB"}}' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed3 > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed4
                            paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed4 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.score | awk '{{print $0 "\\t-"}}' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed6
                            paste -d '\\t' ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.bed6 ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_B.compartments.startEnd | awk '{{print $0 "\\t122,16,180"}}' >> ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed

                            # Sorting compartments and clean intermediary files
                            bedtools sort -i ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed | sed '1s/^/track type=bigBed itemRgb="On"\\n/' > ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments_sorted.bed
                            rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/*.bedGraph
                            rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_*.compartments.bed*
                            rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_*.compartments.score
                            rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_*.compartments.startEnd
                            rm ${{VIZDIR}}/vizIGV_intra/data_compartment_beds/${{w}}_compartments.bed
                        fi
                    fi
                    rm ${{VIZDIR}}/vizIGV_intra/data_bigWig/${{w}}.bedGraph
                done

            done
            """
    # ----------------------------------------------------------------------------------------

else:
    shell("mkdir -p 10_Compartments_detection_dcHiC_notPerformed")
    if (eval(str(config["groups"]["perform_grouped_analyses"])) == True):
        shell("mkdir -p 12_Grouped_analyses/F_Compartments_detection_dcHiC_notPerformed/")



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Differential Hi-C Contacts (SELFISH)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if (eval(str(config["perform_differential_contacts_analyses"])) == True & eval(str(config["groups"]["perform_grouped_analyses"])) == True):
    shell("mkdir -p 12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH_notPerformed/")
    shell("rm -r 12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH_notPerformed/")

    rule T_differential_contacts_SELFISH_groupedSamples:
        input:
            hicpro_matrix_corrected = expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/hicpro_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_corrected.hicpro"])), group = [re.sub("\.|-", str(config["character_subsitution_dashes_and_points_sample_name"]), i) for i in groups], allow_missing=True),
            hicpro_matrix_corrected_bed = expand(os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{group}/hicpro_format/", ''.join(["{group}_mapQ", str(config["mapQ_cutoff"]), "_{ALL_NEW_MATRIX_RESOLUTIONS}kb_normalized_corrected_hicpro.bed"])), group = [re.sub("\.|-", str(config["character_subsitution_dashes_and_points_sample_name"]), i) for i in groups], allow_missing=True)
        output:
            selfish_merge = expand("12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH/{combo}/{combo}_{ALL_NEW_MATRIX_RESOLUTIONS}kb_SELFISH.txt", combo = combo_list, resolution=[str(x) for x in NEW_RESOLUTIONS], allow_missing=True)
        params:
            comboA = ''.join(['"', '" "'.join(comboA), '"']),
            comboB = ''.join(['"', '" "'.join(comboB), '"']),
            mapQ = str(config["mapQ_cutoff"]),
            all_combos = ''.join(['"', '" "'.join(combo_list), '"']),
            selfish_extra_params = config["selfish_params"]["extra_params"],
            resolution = "{ALL_NEW_MATRIX_RESOLUTIONS}",
            chr_filtering_string = str(config["chr_filtering_string"]),
            qValue_threshold = str(config["selfish_params"]["qValue_threshold"])
        benchmark:
            "benchmarks/T_differential_contacts_SELFISH_groupedSamples/T_differential_contacts_SELFISH_groupedSamples-{ALL_NEW_MATRIX_RESOLUTIONS}kb.tsv"
        shell:
            """
            printf '\033[1;36m{params.resolution}kb: Performing grouped differential contacts analyses by SELFISH\\n\033[0m'

            comboA=({params.comboA})
            comboB=({params.comboB})
            combolist=({params.all_combos})

            for i in $(seq 0 $[${{#combolist[@]}} - 1])
            do
                mkdir -p 12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH/${{comboA[i]}}_vs_${{comboB[i]}}/log
                mkdir -p 12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH/${{comboA[i]}}_vs_${{comboB[i]}}/log

                # Get chromosome names
                BED="12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/${{comboA[i]}}/hicpro_format/${{comboA[i]}}_mapQ{params.mapQ}_{params.resolution}kb_normalized_corrected_hicpro.bed"
                CHROM=$(cut -f 1 $BED | uniq | sed 's/chr//' | grep -vE '{params.chr_filtering_string}')

                # Remove the merge file if exists
                file="{output.selfish_merge}"
                if [ -f "$file" ] ; then
                    rm "$file"
                fi

                # Run Selfish analyses
                for C in $CHROM
                do
                    selfish \
                    --resolution {params.resolution}kb \
                    --matrix1 12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/${{comboA[i]}}/hicpro_format/${{comboA[i]}}_mapQ{params.mapQ}_{params.resolution}kb_normalized_corrected.hicpro \
                    --bed1 12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/${{comboA[i]}}/hicpro_format/${{comboA[i]}}_mapQ{params.mapQ}_{params.resolution}kb_normalized_corrected_hicpro.bed \
                    --matrix2 12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/${{comboB[i]}}/hicpro_format/${{comboB[i]}}_mapQ{params.mapQ}_{params.resolution}kb_normalized_corrected.hicpro \
                    --bed2 12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/${{comboB[i]}}/hicpro_format/${{comboB[i]}}_mapQ{params.mapQ}_{params.resolution}kb_normalized_corrected_hicpro.bed \
                    --outfile 12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH/${{comboA[i]}}_vs_${{comboB[i]}}/${{comboA[i]}}_vs_${{comboB[i]}}_{params.resolution}kb_SELFISH.chr${{C}}.tsv \
                    --tsvout {params.qValue_threshold} \
                    --chromosome=${{C}} {params.selfish_extra_params} &> 12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH/${{comboA[i]}}_vs_${{comboB[i]}}/log/${{comboA[i]}}_vs_${{comboB[i]}}_{params.resolution}kb_SELFISH.chr${{C}}.log

                    # Combine all the outputs per resolution in one table
                    tail -n +2 12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH/${{comboA[i]}}_vs_${{comboB[i]}}/${{comboA[i]}}_vs_${{comboB[i]}}_{params.resolution}kb_SELFISH.chr${{C}}.tsv >> {output.selfish_merge}_temp
                done

                # Add the header to the merged table and remove per chr tables
                HEADER="CHR1\\tLOC1_start\\tLOC1_end\\tCHR2\\tLOC2_start\\tLOC2_end\\tQ_VAL\\tLOG_FOLD_CHANGE"
                echo -e $HEADER > {output.selfish_merge} && cat {output.selfish_merge}_temp >> {output.selfish_merge}
                rm {output.selfish_merge}_temp
                rm 12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH/${{comboA[i]}}_vs_${{comboB[i]}}/${{comboA[i]}}_vs_${{comboB[i]}}_{params.resolution}kb_SELFISH.chr*.tsv
            done
            """
else:
    if (eval(str(config["groups"]["perform_grouped_analyses"])) == True):
        shell("mkdir -p 12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH_notPerformed")



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Stripes calling (stripeDiff)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Stripes on single samples
if (eval(str(config["perform_stripes_analyses"])) == True):
    shell("mkdir -p 11_Stripes_analyses_STREPENN_notPerformed")
    shell("rm -r 11_Stripes_analyses_STREPENN_notPerformed")

    rule U1_stripe_detection_STRIPPEN_singleSamples:
        input:
            cool_matrix_corrected = os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{SAMPLES}/cool_format/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{STRIPES_RESOLUTIONS}kb_normalized_corrected.cool"])),
            hicpro_matrix_corrected_bed = os.path.join("06_Interaction_matrices_normalized_and_corrected/corrected_matrices/{SAMPLES}/hicpro_format/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_{STRIPES_RESOLUTIONS}kb_normalized_corrected_hicpro.bed"]))
        output:
            stripeDiff = "11_Stripes_analyses_STREPENN/{SAMPLES}/{SAMPLES}_{STRIPES_RESOLUTIONS}kb/result_filtered.tsv"
        params:
            sample = "{SAMPLES}",
            resolution = "{STRIPES_RESOLUTIONS}",
            out_dir = "11_Stripes_analyses_STREPENN/{SAMPLES}",
            extra_params = config["stripenn_params"]["extra_params"],
            chr_filtering_string = str(config["chr_filtering_string"]),
            pValue_threshold = str(config["stripenn_params"]["pValue_threshold"])
        benchmark:
            "benchmarks/U1_stripe_detection_STRIPPEN_singleSamples/U1_stripe_detection_STRIPPEN_singleSamples-{SAMPLES}_{STRIPES_RESOLUTIONS}kb.tsv"
        log:
            out = "11_Stripes_analyses_STREPENN/{SAMPLES}/logs/{SAMPLES}_{STRIPES_RESOLUTIONS}kb_stripenn_analyses.err"
        threads:
            workflow.cores
        shell:
            """
            printf '\033[1;36m{params.sample} ({params.resolution}kb): Performing stripes analyses by STRIPPEN\\n\033[0m'

            mkdir -p {params.out_dir}/{params.sample}_{params.resolution}kb
            rm -r {params.out_dir}/{params.sample}_{params.resolution}kb
            mkdir -p {params.out_dir}

            # Get chromosome list
            BED="{input.hicpro_matrix_corrected_bed}"
            CHROM=$(cut -f 1 $BED | uniq | grep -vE '{params.chr_filtering_string}')
            CHROMLIST=$(echo $CHROM | sed 's/ /,/g')

            # Run stripes detection
            stripenn compute \
            --cool {input.cool_matrix_corrected} \
            -o {params.out_dir}/{params.sample}_{params.resolution}kb \
            --numcores 1 \
            -k $CHROMLIST \
            --pvalue {params.pValue_threshold} {params.extra_params} &> {log.out}
            """
else:
    shell("mkdir -p 11_Stripes_analyses_STREPENN_notPerformed")


# Stripes on groups
if (eval(str(config["perform_stripes_analyses"])) == True & eval(str(config["groups"]["perform_grouped_analyses"])) == True):
    shell("mkdir -p 12_Grouped_analyses/H_Stripes_analyses_STRIPPEN_notPerformed/")
    shell("rm -r 12_Grouped_analyses/H_Stripes_analyses_STRIPPEN_notPerformed/")

    rule U2_stripe_detection_STRIPPEN_groupedSamples:
        input:
            cool_matrix_corrected = os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{GROUPS}/cool_format/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{STRIPES_RESOLUTIONS}kb_normalized_corrected.cool"])),
            hicpro_matrix_corrected_bed = os.path.join("12_Grouped_analyses/C_summed_matrices_normalized_and_corrected/corrected_matrices/{GROUPS}/hicpro_format/", ''.join(["{GROUPS}_mapQ", str(config["mapQ_cutoff"]), "_{STRIPES_RESOLUTIONS}kb_normalized_corrected_hicpro.bed"]))
        output:
            stripeDiff = "12_Grouped_analyses/H_Stripes_analyses_STRIPPEN/{GROUPS}/{GROUPS}_{STRIPES_RESOLUTIONS}kb/result_filtered.tsv"
        params:
            group = "{GROUPS}",
            resolution = "{STRIPES_RESOLUTIONS}",
            out_dir = "12_Grouped_analyses/H_Stripes_analyses_STRIPPEN/{GROUPS}",
            extra_params = config["stripenn_params"]["extra_params"],
            chr_filtering_string = str(config["chr_filtering_string"]),
            pValue_threshold = str(config["stripenn_params"]["pValue_threshold"])
        benchmark:
            "benchmarks/U2_stripe_detection_STRIPPEN_groupedSamples/U2_stripe_detection_STRIPPEN_groupedSamples-{GROUPS}_{STRIPES_RESOLUTIONS}kb.tsv"
        log:
            out = "12_Grouped_analyses/H_Stripes_analyses_STRIPPEN/{GROUPS}/logs/{GROUPS}_{STRIPES_RESOLUTIONS}kb_stripenn_analyses.err"
        threads:
            workflow.cores
        shell:
            """
            printf '\033[1;36m{params.group} ({params.resolution}kb): Performing grouped stripes analyses by STRIPPEN\\n\033[0m'

            mkdir -p {params.out_dir}/{params.group}_{params.resolution}kb
            rm -r {params.out_dir}/{params.group}_{params.resolution}kb
            mkdir -p {params.out_dir}

            # Get chromosome list
            BED="{input.hicpro_matrix_corrected_bed}"
            CHROM=$(cut -f 1 $BED | uniq | grep -vE '{params.chr_filtering_string}')
            CHROMLIST=$(echo $CHROM | sed 's/ /,/g')

            # Run stripes detection
            stripenn compute \
            --cool {input.cool_matrix_corrected} \
            -o {params.out_dir}/{params.group}_{params.resolution}kb \
            --numcores 1 \
            -k $CHROMLIST \
            --pvalue {params.pValue_threshold} {params.extra_params} &> {log.out}
            """
else:
    if (eval(str(config["groups"]["perform_grouped_analyses"])) == True):
        shell("mkdir -p 12_Grouped_analyses/H_Stripes_analyses_STRIPPEN_notPerformed")
