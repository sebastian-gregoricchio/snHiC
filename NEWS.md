---
title: "changeLog"
---

# snHiC [<img src="https://raw.githubusercontent.com/sebastian-gregoricchio/snHiC/main/resources/snHiC_logo.svg" align="right" height = 150/>](https://sebastian-gregoricchio.github.io/snHiC)
![update](https://badges.pufler.dev/updated/sebastian-gregoricchio/snHiC)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/snHiC?style=social)](https://github.com/sebastian-gregoricchio/snHiC/fork)


#### [v0.2.0](https://github.com/sebastian-gregoricchio/snHiC/releases/tag/0.2.0) - June 21<sup>st</sup> 2023
* New analyses implemented: Tad calling also by `GENOVA`, Loop detection also by `Mustache`, stripes detection by `STRIPENN`, differential contacts caling by `SELFISH`.
* Yaml environment updated
* Providing of a bash script for the installation, fixing and configuration of packages/softwares not included in the yaml file
* Reorganization of the manual in a GitHub-Wiki
* Implementation of the becnhmarking in the snakemake pipeline
* Re-formatting and update of the config file
* Providing of testa fastq files as well as resulting outputs
* Pipeline has been published on *Bioinformatics Advances* (DOI: <a href="https://doi.org/10.1093/bioadv/vbad080">10.1093/bioadv/vbad080</a>)


#### [v0.1.1](https://github.com/sebastian-gregoricchio/snHiC/releases/tag/0.1.1) - February 10<sup>th</sup> 2023
* Implemented the use of [bwa-mem2](https://ieeexplore.ieee.org/document/8820962) instead of bwa in order to reduce the time required for the mapping. It follows that both snakemake file and yaml-formatted environment files have been updated.
* Added a folder containing the `dcHiC` logs in the compartments folder.


#### [v0.1.0](https://github.com/sebastian-gregoricchio/snHiC/releases/tag/0.1.0) - December 29<sup>th</sup> 2022
* First release of the pipeline!
