#!/bin/bash

$1/bin/git clone https://github.com/ay-lab/dcHiC $1/bin/dcHiC
$1/bin/R CMD INSTALL $1/bin/dcHiC/packages/functionsdchic_1.0.tar.gz

echo ''
echo '>>> dcHiC installed in:' $1/bin/dcHiC/ '<<<'