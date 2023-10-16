#!/bin/bash

threads=120
WD="/fungen/funhome/Assembly/Krill/Metagenomic_assembly/final"
bins="${WD}/INITIAL_BINNING/final_bins"
NAME=KigumiPool
MinContigLength=5000

for (( i=1; i<=153; i++ )); do
 bin="bin.${i}"
 OUT="${WD}/Annotation/Prokka.v1.14.6/${bin}"
 mkdir -p "${OUT}"
 
 #################################################
 # de novo prediction : Prokka
 
 #conda activate prokka_env
 echo "$(date +"%d-%m-%Y %H:%M:%S"): Starting Prokka for bin.${i}" 
 prokka --outdir ${OUT} ${bins}/${bin}.fa --cpus=${threads} --metagenome --mincontiglen ${MinContigLength} --addgenes --prefix ${NAME}.${bin} --locustag ${bin} --centre IGB  --norrna --notrna --force 

done
