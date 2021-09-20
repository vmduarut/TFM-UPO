#!/bin/bash

if [[ ! $# -eq 1 ]]
then
        echo "USAGE: $(basename $0) <inputFile>"
        exit 1
fi


set -x

rawInputFile=$1
outputDir="/mnt/d/Data/resultsTFM"

FREQ_CUTOFF="0.05"
POPULATION_PANEL="/mnt/d/Data/all_population.panel"
THREADS=$(nproc)

baseFileName=$(basename $rawInputFile)
outputSplitFile="${outputDir}/${baseFileName/.vcf.gz/.split.vcf.gz}"
outputAfFile=${outputSplitFile/.vcf.gz/.AF.vcf.gz}
outputFilteredFile="${outputAfFile/.vcf.gz/}_${FREQ_CUTOFF}_v2.vcf.gz"

# separamos multialélicas

bcftools norm -m - --threads ${THREADS}  ${rawInputFile} -O z -o ${outputSplitFile}

# Calculamos frecuencias de las 26 poblaciones

bcftools +fill-tags ${outputSplitFile} -Oz -o ${outputAfFile} -- -S ${POPULATION_PANEL} -t AF

# filtramos por frecuencia, según el filtro añadimos el resto de poblaciones o no

bcftools view --threads ${THREADS} -i 'AF_IBS>='${FREQ_CUTOFF}'' -O z -o ${outputFilteredFile} ${outputAfFile}

