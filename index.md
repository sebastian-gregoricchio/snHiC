[![Snakemake](https://img.shields.io/badge/snakemake-≥7.8.5-brightgreen.svg)](https://snakemake.github.io)
![release](https://img.shields.io/github/v/release/sebastian-gregoricchio/snHiC)
[![license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://sebastian-gregoricchio.github.io/snHiC/LICENSE.md/LICENSE)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/snHiC?style=social)](https://github.com/sebastian-gregoricchio/snHiC/fork)
<!-- ![update](https://badges.pufler.dev/updated/sebastian-gregoricchio/snHiC)
![visits](https://badges.pufler.dev/visits/sebastian-gregoricchio/snHiC)
![downloads](https://img.shields.io/github/downloads/sebastian-gregoricchio/snHiC/total.svg)--->

<h1> snHiC </h1>

-------------------------

* TOC
{:toc}

--------------------------

<br/><br/>

# Introduction [<img src="https://raw.githubusercontent.com/sebastian-gregoricchio/snHiC/main/resources/snHiC_logo.svg" align="right" height = 150/>](https://sebastian-gregoricchio.github.io/snHiC)


`snHiC` is a snakemake based end-to-end pipeline to analyze Hi-C data. The input files required to run the pipeline are Paired-End fastq files. The pipeline performs data quality control, normalization and correction. It also includes the possibility to perform grouped analyses (e.g, merging of replicates) besides TAD, loops and stripes detection and differential contacts and compartment analyses. Notabily, the latter is performed using `dcHiC`, a recently published method ([A. Chakraborty, *et al.*, Nat. Comm. 2022](https://www.nature.com/articles/s41467-022-34626-6)) that enables more precise and high-resolution differential compartment analyses.


## Citation
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

# Documentation
Details on the [installation](https://github.com/sebastian-gregoricchio/snHiC/wiki/2.-Installation-and-dependencies) and [usage](https://github.com/sebastian-gregoricchio/snHiC/wiki/3.-Run-the-pipeline) of snHiC can be found at the dedicated [Wiki](https://github.com/sebastian-gregoricchio/snHiC/wiki/).

<br/><br/>

-----------------
# Package history and releases
A list of all releases and respective description of changes applied could be found [here](https://sebastian-gregoricchio.github.io/snHiC/NEWS).

# Contact
For any suggestion, bug fixing, commentary please report it in the [issues](https://github.com/sebastian-gregoricchio/snHiC/issues)/[request](https://github.com/sebastian-gregoricchio/snHiC/pulls) tab of this repository.

# License
This repository is under a [GNU General Public License (version 3)](https://sebastian-gregoricchio.github.io/snHiC/LICENSE.md/LICENSE).

<br/>

### Contributors
[![contributors](https://contrib.rocks/image?repo=sebastian-gregoricchio/snHiC)](https://sebastian-gregoricchio.github.io/)
