#!/bin/bash

$1/bin/git clone https://github.com/ay-lab/dcHiC $1/bin/dcHiC
$1/bin/R CMD INSTALL $1/bin/dcHiC/packages/functionsdchic_1.0.tar.gz

$1/bin/Rscript -e 'if (!require("GENOVA")) devtools::install_github("robinweide/GENOVA@v1.0.0-TheCartographer")'

$1/bin/git clone https://github.com/akdemirlab/HiCPlotter $1/bin/HiCPlotter
cp $1/bin/HiCPlotter/HiCPlotter.py $1/bin/
chmod +x $1/bin/HiCPlotter.py

echo ''
echo '>>> dcHiC, GENOVA and HiCPlotter are installed <<<'
echo '>>> dcHiC folder is avilable at:' $1/bin/dcHiC/ '<<<'