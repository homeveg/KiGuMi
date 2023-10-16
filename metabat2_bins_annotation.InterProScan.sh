#!/bin/bash

threads=60

WD="/fungen/funhome/Assembly/Krill/Metagenomic_assembly/final"
bins="${WD}/INITIAL_BINNING/final_bins"
NAME=KigumiPool
 
for (( i=1; i<=153; i++ )); do
 ###############################################
 # InterProScan annotation (under conda base)
 bin="bin.${i}"
 OUT="${WD}/Annotation/InterProScan/${bin}/"
 INPUT=${WD}/Annotation/Prokka.v1.14.6/${bin}/${NAME}.${bin}.faa
 
 # check if folder already exists
 if [ -d "$OUT" ]; then
  echo "folder $OUT exist"
 else
  mkdir -p $OUT
  echo "$(date +"%d-%m-%Y %H:%M:%S"): check if fasta file contains asterics in the amino-acid sequence bin.${i}"
  # Count the number of lines that contain asterisks
  NrHits=$(grep -c "\*" "$INPUT")
  if [[ $NrHits -gt 0 ]]; then
   echo "$(date +"%d-%m-%Y %H:%M:%S"): Removing $NrHits asterisks from $INPUT"
   perl -p -i -e 's/\*//g' "$INPUT"
  fi
	
  echo "$(date +"%d-%m-%Y %H:%M:%S"): Starting InterProScan for bin.${i}"
  echo "run InterPro on Prokka-annotated boins"
  interproscan.sh --output-file-base ${OUT} --cpu ${threads} --goterms --iprlookup --tempdir /mnt/ramdisk/tmp/ --input ${INPUT} -appl PANTHER,TIGRFAM
 fi
 
done
