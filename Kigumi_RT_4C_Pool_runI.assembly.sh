#!/bin/bash

#metagenomic hybrid assembly pipeline

Project="Kigumi_RT_4C_Pool_runI"
FCID="20230524_1449_MN15721_FAW86564_98fc40ab"
#SeqSummary="sequencing_summary_FAW86564_98fc40ab_9bbb0371.txt"
IlluminaSampleRT="kamo_699"
IlluminaSample4C="kamo_700"
BR=6

Basecaller="dorado"
ONTreads="Kigumi_RT_4C_Pool_runI.sup.fastq"

# set global varibles
threads=62 # cn5 
Qscore=15
MinLength=1000

# filter contigs after assembly
MinContigLength=10000
MinContigCoverage=10.0

# set PATH to working/temporary directories
TEMPDIR="/mnt/ramdisk/tmp"
WD="/fungen/funhome/Assembly/Krill/Metagenomic_assembly/final"
ONTProjectPath="/fungen/funhome/db05/OxNano_Sequencing_Projects/Project_Krill_UniOldenburg/Metagenomic_assembly"
IlluminaProject="/fungen/funhome/db02/Sequencing_Projects/Project_kigumi_dna"
ONTPath="$ONTProjectPath/$Project/$FCID/fastq_${Basecaller}/$ONTreads"
Spades_draft_genome=${WD}/Assembly/$Project/contigs.fasta
Filtered_draft_genome=${WD}/Assembly/$Project/$Project.contigs_L${MinLength}.fasta

# set alias to conda 
alias bbduk="/fungen/funhome/Software/miniconda3/envs/bbtools_env/bin/bbduk.sh"
alias spades="/fungen/funhome/Software/miniconda3/envs/spades_env/bin/spades.py"
alias merge_reads="python /fungen/funhome/Assembly/Krill/Metagenomic_assembly/final/NoSpace_fastqIDs.py"
alias quast="/fungen/funhome/Software/miniconda3/envs/quast_env/bin/quast"
alias samtools="/fungen/funhome/Software/miniconda3/envs/medaka_env/bin/samtools" 

# output path
mkdir -p "${WD}/QC/fastqc_out"
mkdir -p "${WD}/Assembly/$Project"
mkdir -p "${WD}/Illumina/${Project}"
mkdir -p "${WD}/ONTreads/${Project}"
mkdir -p "${WD}/Racon/${Project}"


# define functions
generate_yaml() {
  local WD=$1
  local Project=$2
  local TEMPDIR=$3
  local IlluminaSampleRT=$4
  local IlluminaSample4C=$5
  local BR=$6
  local Qscore=$7
  local MinLength=$8
  
  local outputFile="${WD}/Assembly/${Project}.input_data.yaml"
  local leftReads=()
  local rightReads=()
  
  # Generate left and right reads
  for ((b=1; b<=$BR; b++))
  do
    R1readsRT="${TEMPDIR}/reads/${IlluminaSampleRT}_br${b}_R1_AQ20.fastq"
    leftReads+=("\"${R1readsRT}\"")
  done
  
  for ((b=1; b<=$BR; b++))
  do
    R1reads4C="${TEMPDIR}/reads/${IlluminaSample4C}_br${b}_R1_AQ20.fastq"
    leftReads+=("\"${R1reads4C}\"")
  done
  
  for ((b=1; b<=$BR; b++))
  do
    R2readsRT="${TEMPDIR}/reads/${IlluminaSampleRT}_br${b}_R2_AQ20.fastq"
    rightReads+=("\"${R2readsRT}\"")
  done
  
  for ((b=1; b<=$BR; b++))
  do
    R2reads4C="${TEMPDIR}/reads/${IlluminaSample4C}_br${b}_R2_AQ20.fastq"
    rightReads+=("\"${R2reads4C}\"")
  done
  
  # Write to output file
  {
    echo "["
    echo "  {"
    echo '    type: "paired-end",'
    echo '    orientation: "fr",'
    echo '    left reads: ['
  
    echo "      ${leftReads[0]}"
    for read in "${leftReads[@]:1}"; do
      echo "      ${read},"
    done
  
    echo '    ],'
    echo '    right reads: ['
  
    echo "      ${rightReads[0]}"
    for read in "${rightReads[@]:1}"; do
      echo "      ${read},"
    done
  
    echo '    ]'
    echo "  },"
    echo "  {"
    echo '    type: "nanopore",'
    echo '    single reads: ['
    echo "      \"${TEMPDIR}/reads/${Project}.filtered.Q${Qscore}.L${MinLength}.fastq\""
    echo '    ]'
    echo "  }"
    echo "]"
  } > $outputFile
}

split_assembly() {
    # Define function parameters
    local parts=$1
    local Polisher=$2
    local WD=$3
    local Filtered_draft_genome=$4
    local Project=$5
    local MinLength=$6
	local OutDir=$7

    # split genome to 10 chunks and polish each separately
    echo "$(date +"%d-%m-%Y %H:%M:%S"):  split draft assembly into $parts chunks"
    local split_fasta=$WD/$Polisher/${Project}.contigs_L${MinLength}.fasta

    cp $Filtered_draft_genome $split_fasta

    # Calculate total sequence length
    local total_length=$(awk '
    BEGIN{RS=">"; ORS=""}
    $0 !~ /^$/{
        split($1, arr, "_");
        len = arr[4];
        total_length += len;
    }END{print total_length}' $split_fasta)
    echo "draft sequence length: $total_length bps"

    # Begin the awk command
    pv $Filtered_draft_genome | awk -v parts="$parts" -v total_length="$total_length" -v Project="$Project" -v OutDir="$OutDir" '
    BEGIN{RS=">"; ORS=""; file_num=1; current_length=0; max_length=total_length/parts}
    $0 !~ /^$/{
        split($1, arr, "_");
        len = arr[4];
        if (current_length + len > max_length && file_num < parts) {
            file_num++; current_length=0
        }
        print ">" $0 > sprintf("%s/%s.filtered_assembly.p%02d.fasta", OutDir, Project, file_num);
        current_length += len;
    }' 
}



# Get the current time
start_time=$(date +%s)
echo "$(date +"%d-%m-%Y %H:%M:%S"): Starting NanoPlot for raw $ONTSample"
NanoPlot --verbose -o ${WD}/QC/NanoPlot/${Project}.raw --fastq_rich ${ONTPath} --threads ${threads}

echo "$(date +"%d-%m-%Y %H:%M:%S"): filter low quality long reads"
pv ${ONTPath} | NanoFilt -q ${Qscore} -l ${MinLength} > ${WD}/ONTreads/${Project}.filtered.Q${Qscore}.L${MinLength}.fastq
mkdir -p ${TEMPDIR}/reads
pv ${WD}/ONTreads/${Project}.filtered.Q${Qscore}.L${MinLength}.fastq > ${TEMPDIR}/reads/${Project}.filtered.Q${Qscore}.L${MinLength}.fastq

echo "$(date +"%d-%m-%Y %H:%M:%S"): Starting NanoPlot for filtered $ONTSample"
NanoPlot --verbose -o ${WD}/QC/NanoPlot/${Project}.NanoFilt --fastq_rich ${WD}/ONTreads/${Project}.filtered.Q${Qscore}.L${MinLength}.fastq --threads ${threads}

for ((b=1; b<=$BR; b++))
do
  R1reads_RT=${IlluminaProject}/Sample_${IlluminaSampleRT}_br${b}/${IlluminaSampleRT}_br${b}_R1.fastq.gz
  R2reads_RT=${IlluminaProject}/Sample_${IlluminaSampleRT}_br${b}/${IlluminaSampleRT}_br${b}_R2.fastq.gz
  R1reads_4C=${IlluminaProject}/Sample_${IlluminaSample4C}_br${b}/${IlluminaSample4C}_br${b}_R1.fastq.gz
  R2reads_4C=${IlluminaProject}/Sample_${IlluminaSample4C}_br${b}/${IlluminaSample4C}_br${b}_R2.fastq.gz

  echo "$(date +"%d-%m-%Y %H:%M:%S"): QC raw Illumina reads"
  #fastqc -o ${WD}/QC/fastqc_out --noextract -t ${threads} ${R1reads_RT}
  #fastqc -o ${WD}/QC/fastqc_out --noextract -t ${threads} ${R2reads_RT}
  #fastqc -o ${WD}/QC/fastqc_out --noextract -t ${threads} ${R1reads_4C}
  #fastqc -o ${WD}/QC/fastqc_out --noextract -t ${threads} ${R2reads_4C}

  echo "$(date +"%d-%m-%Y %H:%M:%S"): filter low quality short reads"
  bbduk trimpolyg=10 hdist=1 mink=12 threads=${threads} maxns=10 minlen=50 trimq=20 qtrim=t ktrim=r k=28 ref=/fungen/funhome/db02/Reference_DB/adapterRef/adapter.fa in1=${R1reads_RT} in2=${R2reads_RT} out1=${WD}/Illumina/${Project}/${IlluminaSampleRT}_br${b}_R1_AQ20.fastq.gz out2=${WD}/Illumina/${Project}/${IlluminaSampleRT}_br${b}_R2_AQ20.fastq.gz
  
  bbduk trimpolyg=10 hdist=1 mink=12 threads=${threads} maxns=10 minlen=50 trimq=20 qtrim=t ktrim=r k=28 ref=/fungen/funhome/db02/Reference_DB/adapterRef/adapter.fa in1=${R1reads_4C} in2=${R2reads_4C} out1=${WD}/Illumina/${Project}/${IlluminaSample4C}_br${b}_R1_AQ20.fastq.gz out2=${WD}/Illumina/${Project}/${IlluminaSample4C}_br${b}_R2_AQ20.fastq.gz
  
  echo "$(date +"%d-%m-%Y %H:%M:%S"): decompress Illumina reads"
  pigz -d -p ${threads} -k ${WD}/Illumina/${Project}/${IlluminaSampleRT}_br${b}_R1_AQ20.fastq.gz
  pigz -d -p ${threads} -k ${WD}/Illumina/${Project}/${IlluminaSampleRT}_br${b}_R2_AQ20.fastq.gz
  pigz -d -p ${threads} -k ${WD}/Illumina/${Project}/${IlluminaSample4C}_br${b}_R1_AQ20.fastq.gz
  pigz -d -p ${threads} -k ${WD}/Illumina/${Project}/${IlluminaSample4C}_br${b}_R2_AQ20.fastq.gz

  echo "$(date +"%d-%m-%Y %H:%M:%S"): copy reads data to ramdisk"

  cp ${WD}/Illumina/${Project}/${IlluminaSampleRT}_br${b}_R1_AQ20.fastq ${TEMPDIR}/reads/
  cp ${WD}/Illumina/${Project}/${IlluminaSampleRT}_br${b}_R2_AQ20.fastq ${TEMPDIR}/reads/
  cp ${WD}/Illumina/${Project}/${IlluminaSample4C}_br${b}_R1_AQ20.fastq ${TEMPDIR}/reads/
  cp ${WD}/Illumina/${Project}/${IlluminaSample4C}_br${b}_R2_AQ20.fastq ${TEMPDIR}/reads/
done


echo "$(date +"%d-%m-%Y %H:%M:%S"): copy reads data to ramdisk"
for ((b=1; b<=$BR; b++))
do
  R1readsRT=${TEMPDIR}/reads/${IlluminaSampleRT}_br${b}_R1_AQ20.fastq
  R2readsRT=${TEMPDIR}/reads/${IlluminaSampleRT}_br${b}_R2_AQ20.fastq
  
  echo "$(date +"%d-%m-%Y %H:%M:%S"): copying ${IlluminaSampleRT}_br${b}_R1_AQ20.fastq"
  if [ -f "${R1readsRT}" ]; then
    echo "file exists"
  else
    pv ${WD}/Illumina/${Project}/${IlluminaSampleRT}_br${b}_R1_AQ20.fastq > ${R1readsRT}
  fi
  
  echo "$(date +"%d-%m-%Y %H:%M:%S"): copying ${IlluminaSampleRT}_br${b}_R2_AQ20.fastq"
  if [ -f "${R2readsRT}" ]; then
    echo "file exists"
  else
    pv ${WD}/Illumina/${Project}/${IlluminaSampleRT}_br${b}_R2_AQ20.fastq > ${R2readsRT}
  fi
done

for ((b=1; b<=$BR; b++))
do
  R1reads4C=${TEMPDIR}/reads/${IlluminaSample4C}_br${b}_R1_AQ20.fastq
  R2reads4C=${TEMPDIR}/reads/${IlluminaSample4C}_br${b}_R2_AQ20.fastq
  
  echo "$(date +"%d-%m-%Y %H:%M:%S"): copying ${IlluminaSampleRT}_br${b}_R1_AQ20.fastq"
  if [ -f "${R1reads4C}" ]; then
    echo "file exists"
  else
    pv ${WD}/Illumina/${Project}/${IlluminaSample4C}_br${b}_R1_AQ20.fastq > ${R1reads4C}
  fi
  
  echo "$(date +"%d-%m-%Y %H:%M:%S"): copying ${IlluminaSampleRT}_br${b}_R2_AQ20.fastq"
  if [ -f "${R2reads4C}" ]; then
    echo "file exists"
  else
    pv ${WD}/Illumina/${Project}/${IlluminaSample4C}_br${b}_R2_AQ20.fastq > ${R2reads4C}
  fi
done

# generate YAML
# please check for correctness of YAML manually!!
generate_yaml "${WD}" "${Project}" "${TEMPDIR}" "${IlluminaSampleRT}" "${IlluminaSample4C}" "$BR" "$Qscore" "$MinLength"
command="spades --meta --threads ${threads} --memory 1500 --tmp-dir ${TEMPDIR} --checkpoints all"
command+=" --dataset ${WD}/Assembly/${Project}.input_data.yaml -o ${WD}/Assembly/${Project}"

echo "$(date +"%d-%m-%Y %H:%M:%S"): Starting draft metagenome assembly"
echo $command
eval $command

echo "$(date +"%d-%m-%Y %H:%M:%S"): Quast QC draft assembly $Spades_draft_genome"
quast.py --output-dir ${WD}/QC/Quast --min-contig ${MinContigLength} --threads ${threads} -f -b --nanopore ${TEMPDIR}/reads/${Project}.filtered.Q${Qscore}.L${MinLength}.fastq -1 ${R1readsRT} -2 ${R2readsRT} ${Spades_draft_genome}


# Spades contig name: >NODE_17388_length_1000_cov_2.624339 
echo "$(date +"%d-%m-%Y %H:%M:%S"): remove contigs smaller than ${MinLength}bp with the coverage below $MinContigCoverage"
pv $Spades_draft_genome | awk -v min_length="$MinContigLength" -v min_cov="$MinContigCoverage" '
BEGIN{RS=">"; ORS=""}
$0 !~ /^$/{
    split($1, arr, "_");
    len = arr[4];
    cov = arr[6];
    if (len >= min_length && cov >= min_cov)
        print ">" $0
}' > $Filtered_draft_genome
cat $Filtered_draft_genome | grep ">" | wc -l

###############################
# merge paired-end reads
mkdir -p ${TEMPDIR}/reads

echo "$(date +"%d-%m-%Y %H:%M:%S"): merge paired-end reads"
for ((b=1; b<=$BR; b++))
do
 R1readsRT=${WD}/Illumina/${Project}/${IlluminaSampleRT}_br${b}_R1_AQ20.fastq
 R2readsRT=${WD}/Illumina/${Project}/${IlluminaSampleRT}_br${b}_R2_AQ20.fastq
 R1reads4C=${WD}/Illumina/${Project}/${IlluminaSample4C}_br${b}_R1_AQ20.fastq
 R2reads4C=${WD}/Illumina/${Project}/${IlluminaSample4C}_br${b}_R2_AQ20.fastq
 
 echo "$(date +"%d-%m-%Y %H:%M:%S"): extracting ${R1readsRT}"
 if [ -f "${R1readsRT}" ]; then
  echo "file exists"
 else
  pigz -d -p ${threads} -k ${R1readsRT}.gz -c | pv > ${R1readsRT}
 fi
 
 echo "$(date +"%d-%m-%Y %H:%M:%S"): extracting ${R2readsRT}"
 if [ -f "${R2readsRT}" ]; then
  echo "file exists"
 else
  pigz -d -p ${threads} -k ${R2readsRT}.gz -c | pv > ${R2readsRT}
 fi

 echo "$(date +"%d-%m-%Y %H:%M:%S"): extracting ${R1reads4C}"
 if [ -f "${R1reads4C}" ]; then
  echo "file exists"
 else
  pigz -d -p ${threads} -k ${R1reads4C}.gz -c | pv > ${R1reads4C}
 fi

 echo "$(date +"%d-%m-%Y %H:%M:%S"): extracting ${R2reads4C}"
 if [ -f "${R2reads4C}" ]; then
  echo "file exists"
 else
  pigz -d -p ${threads} -k ${R2reads4C}.gz -c | pv > ${R2reads4C}
 fi

 R1readsRT_size=$(stat -c%s "$R1readsRT")
 R2readsRT_size=$(stat -c%s "$R2readsRT")
 R1reads4C_size=$(stat -c%s "$R1reads4C")
 R2reads4C_size=$(stat -c%s "$R2reads4C")
 
 expected_size=$((R1readsRT_size+R2readsRT_size))
 echo "$(date +"%d-%m-%Y %H:%M:%S"): merge ${R1readsRT} ${R2readsRT} reads"
 if [ -f "${TEMPDIR}/reads/${IlluminaSampleRT}_br${b}.fastq.gz" ]; then
  echo "file exists"
 else
  merge_reads ${R1readsRT} ${R2readsRT} | pv -s $expected_size | pigz > ${TEMPDIR}/reads/${IlluminaSampleRT}_br${b}.fastq.gz
  rm ${R1readsRT} ${R2readsRT}
 fi

 expected_size=$((R1reads4C_size+R2reads4C_size))
 echo "$(date +"%d-%m-%Y %H:%M:%S"): merge ${R1reads4C} ${R2reads4C} reads"
  if [ -f "${TEMPDIR}/reads/${IlluminaSample4C}_br${b}.fastq.gz" ]; then
  echo "file exists"
 else
  merge_reads ${R1reads4C} ${R2reads4C} | pv -s $expected_size | pigz > ${TEMPDIR}/reads/${IlluminaSample4C}_br${b}.fastq.gz
  rm ${R1reads4C} ${R2reads4C}
 fi
done

#===
#pv ${TEMPDIR}/reads/${IlluminaSampleRT}_br* > ${TEMPDIR}/reads/${IlluminaSampleRT}_AQ20.ns.fastq.gz
#pv ${TEMPDIR}/reads/${IlluminaSample4C}_br* > ${TEMPDIR}/reads/${IlluminaSample4C}_AQ20.ns.fastq.gz
#rm ${TEMPDIR}/reads/kamo_70*br*

###############################
# polish with MinIon reads
# map extracted reads to the draft genome assembly

echo "$(date +"%d-%m-%Y %H:%M:%S"):  polishing draft assembly with long reads"
# split genome to 10 chunks and polish each separately
# Get the number of chunks to split the file into
parts=10
Polisher="Racon"

echo "$(date +"%d-%m-%Y %H:%M:%S"):  split draft assembly $Filtered_draft_genome into $parts chunks"
#split_assembly $parts $Polisher $WD $Filtered_draft_genome ${Project} ${MinLength} ${WD}/$Polisher/${Project}

echo "$(date +"%d-%m-%Y %H:%M:%S"): processing $i from $parts sequence chunks"
basecalls=${WD}/ONTreads/${Project}.filtered.Q${Qscore}.L${MinLength}.fastq
ONTalignments=${WD}/$Polisher/${Project}/${Project}.filtered_assembly
ONTpolished=${WD}/$Polisher/${Project}/${Project}.ONTpolished.fasta
 
#Align the ONT reads to the draft assembly
echo "$(date +"%d-%m-%Y %H:%M:%S"): map ONT reads to a polished draft genome"
minimap2 -ax map-ont -t ${threads} $Filtered_draft_genome ${basecalls} > ${ONTalignments}.sam

racon -t ${threads} -m 8 -x -6 -g -8 -w 500 \
${basecalls} \
${ONTalignments}.sam \
$Filtered_draft_genome > ${ONTpolished}

Racon_polished=${WD}/$Polisher/${Project}.ONT.polished.fasta
pv ${WD}/$Polisher/${Project}/${Project}.ONTpolished.fasta > ${Racon_polished}

for ((b=1; b<=$BR; b++))
 IlluminaRTalignments="${WD}/$Polisher/${Project}/${IlluminaSampleRT}_br${b}.accepted_hits.sam"
 Illumina_RT=${TEMPDIR}/reads/${IlluminaSampleRT}_br${b}.fastq.gz
 Illumina_RT_polished=$(printf "${WD}/$Polisher/${Project}/${Project}.Illumina_RT_polished.br_%02d.fasta" $b)

 echo "$(date +"%d-%m-%Y %H:%M:%S"): map merged reads to a ONT polished sequence chunk ${Racon_polished}"
 ngm -o ${IlluminaRTalignments} -r ${Racon_polished} -q ${Illumina_RT} -t ${threads} -g 0,1,2,3
 
 echo "$(date +"%d-%m-%Y %H:%M:%S"): Racon contigs poliching (consensus) with short reads"
 racon -t ${threads} -m 8 -x -6 -g -8 -w 500 \
 ${Illumina_RT} \
 ${IlluminaRTalignments} \
 ${Racon_polished} > ${Illumina_RT_polished}

 Illumina_4C=${TEMPDIR}/reads/${IlluminaSample4C}_br${b}.fastq.gz
 Illumina_4Calignments="${WD}/$Polisher/${Project}/${IlluminaSample4C}_br${b}.accepted_hits.sam"
 Illumina_RT_4C_polished=$(printf "${WD}/$Polisher/${Project}/${Project}.Illumina_RT_4C_polished.br_%02d.fasta" $b)

 echo "$(date +"%d-%m-%Y %H:%M:%S"): map merged reads to a Illumina_RT polished sequence chunk ${Illumina_RT_polished_chunk}"
 ngm -o ${Illumina_4Calignments} -r ${Illumina_RT_polished} -q ${Illumina_4C} -t ${threads} -g 0,1,2,3
 echo "$(date +"%d-%m-%Y %H:%M:%S"): Racon contigs poliching (consensus) with short reads"
 
 racon -t ${threads} -m 8 -x -6 -g -8 -w 500 \
 ${Illumina_4C} \
 ${Illumina_4Calignments} \
 ${Illumina_RT_polished} > ${Illumina_RT_4C_polished}

 Racon_polished=${Illumina_RT_4C_polished}
done

final_polished_assembly=${WD}/$Polisher/${Project}/${Project}.Illumina.polished.fasta
cp $Racon_polished $final_polished_assembly

# clean temporary
rm ${WD}/$Polisher/${Project}/${Project}.Illumina_RT_4C_polished.br_0*
rm ${WD}/$Polisher/${Project}/${Project}.Illumina_RT_polished.br_0*
rm ${WD}/$Polisher/${Project}/*.sam
rm ${TEMPDIR}/reads/*.gz

#clean temporary filies

R1readsRT=${WD}/Illumina/${Project}/${IlluminaSampleRT}_br1_R1_AQ20.fastq.gz
R2readsRT=${WD}/Illumina/${Project}/${IlluminaSampleRT}_br1_R2_AQ20.fastq.gz
ONTReads=${WD}/ONTreads/${Project}.filtered.Q${Qscore}.L${MinLength}.fastq

echo "$(date +"%d-%m-%Y %H:%M:%S"): Run Quast on the final polished assembly $final_polished_assembly"
quast.py --output-dir ${WD}/QC/Quast/${Project}/final --min-contig ${MinLength} --threads ${threads} -f -b --nanopore ${ONTReads} -1 ${R1readsRT} -2 ${R2readsRT} $final_polished_assembly

rm -rf ${TEMPDIR}/*


# binning
mkdir -p ${WD}/INITIAL_BINNING/${Project}
mkdir -p ${WD}/BIN_REFINEMENT/${Project}
mkdir -p ${WD}/BLOBOLOGY/${Project}

# rename file according metawrap format
dir="${WD}/Illumina/${Project}"

for file in "$dir"/*_AQ20.fastq; do
  # Remove the directory in file name
  base=$(basename "$file")

  # Extract required parts from filename
  part1=${base%%_br*}
  part2=${base#*_br}; part2=${part2%%_R*}
  r=${base##*_R}; r=${r%%_*}

  # Construct new filename
  # Construct new filename
  new_file="${dir}/${part1}_br${part2}_${r}.fastq"
  echo "renaming $base to ${part1}_br${part2}_${r}.fastq"
  
  # Rename the file
  mv -- "$file" "$new_file"
done

echo "$(date +"%d-%m-%Y %H:%M:%S"): MetaWRAP: INITIAL_BINNING"
metawrap binning -o ${WD}/INITIAL_BINNING/${Project} -t $threads -a $final_polished_assembly --metabat2 --maxbin2 --concoct ${WD}/Illumina/${Project}/*fastq

echo "$(date +"%d-%m-%Y %H:%M:%S"): MetaWRAP: BIN_REFINEMENT"
metawrap bin_refinement -o ${WD}/BIN_REFINEMENT/${Project} -t $threads -A ${WD}/INITIAL_BINNING/${Project}/metabat2_bins/ -B ${WD}/INITIAL_BINNING/${Project}/maxbin2_bins/ -C ${WD}/INITIAL_BINNING/${Project}/concoct_bins/ -c 50 -x 10

metawrap blobology -a $final_polished_assembly -t $threads -o ${WD}/BLOBOLOGY/${Project} --bins ${WD}/BIN_REFINEMENT/${Project}/metawrap_50_10_bins ${WD}/Illumina/${Project}/*fastq
