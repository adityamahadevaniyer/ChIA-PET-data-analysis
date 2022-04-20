#!/usr/bin/bash
#SBATCH -N 1 # number of nodes
#SBATCH -n 30  # number of threads
#SBATCH --mem 60Gb # memory pool for all cores
#SBATCH -t 0-72:00 # time (D-HH:MM)
#SBATCH --job-name LHG0147
#SBATCH --output="LHG0147.slurm.o"
#SBATCH --error="LHG0147.slurm.e"

NTHREAD="30"

## The help message:
function usage
{
    echo -e "usage: sbatch  chiapipe_sumner.sh -c config_file" 
}

## Parse the command-line argument (i.e., get the name of the config file)
while [ "$1" != "" ]; do
    case $1 in
        -c | --conf )           shift
                                conf=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

source ${conf}
hic_file="ChIA-PET_${genome}_${cell_type}_${ip_factor}_${run}_${run_type}_pairs.hic"
mkdir -p ${run}
cd ${run}

module load singularity

# input config
# The name of the sequencing run
LOGFILE="${run}.log"
addpath="/fastscratch/kimm/ChIA-PIPE_Sumner/src"
scom="singularity run $addpath/cpu0.0.1a-r2.sumner.sif"
pigz="$scom pigz"
cpuprog="$scom cpu"

#samt="singularity run $addpath/samtools_1.5.sif"
#samtools="$samt samtools"
scmd="singularity exec $addpath/samtools_1.5.sif samtools"
bcmd="singularity exec $addpath/bedtools_2.26.sif bedtools"
#bedt="singularity run $addpath/bedtools_2.26.sif"
#bedtools="$bedt bedtools"
kentUtil="singularity run $addpath/kentUtils_369.sif"
bedgraph2bigwig="$kentUtil bedGraphToBigWig"

bam2pair="/fastscratch/kimm/pairix_src/util/bam2pairs/bam2pairs"
#javatool="singularity run $addpath/java_1.7.0.sif"
#javarun="$javatool java"
juicer="$addpath/juicer_tools.1.7.5_linux_x64_jcuda.0.8.jar"

spp="$addpath/r-spp_1.13.sif"
macs2="$addpath/macs2_2.1.0.20151222.sif"

echo "--- ChIA-PIPE start ---" >> ${LOGFILE}
echo "`date`" >> ${LOGFILE}

$cpuprog stag -W -T 18 -t ${NTHREAD} -O ${run} ${r1_fastq} ${r2_fastq} 2>> ${LOGFILE}
echo "--- linker detection completed ---" >> ${LOGFILE}
echo "`date`" >> ${LOGFILE}

#Get the stat
$cpuprog stat -s -p -T 18 -t ${NTHREAD} ${run}.cpu 2>> ${LOGFILE} 1> ${run}.stat
echo "--- statistics done ---" >> ${LOGFILE}
echo "`date`" >> ${LOGFILE}


echo "--- pigziiping ---" >> ${LOGFILE}
$pigz -p ${NTHREAD} ${run}.singlelinker.paired.fastq 2>> ${LOGFILE}
$pigz -p ${NTHREAD} ${run}.singlelinker.single.fastq 2>> ${LOGFILE}
$pigz -p ${NTHREAD} ${run}.none.fastq 2>> ${LOGFILE}
$pigz -p ${NTHREAD} ${run}.conflict.fastq 2>> ${LOGFILE}
$pigz -p ${NTHREAD} ${run}.tied.fastq 2>> ${LOGFILE}

#------------------------------------------------------------------------------
#defined names
pairlabel="singlelinker.paired"
pair_map_qual="30"
pair_suffix="UU"
singlabel="singlelinker.single"
single_map_qual="10"
single_suffix="UxxU"
nonelabel="none"
none_map_qual="30"
none_suffix="UU"
selfbp="8000"
extbp="500"

#Mapping
echo "--- Mapping start ---" >> ${LOGFILE}

#Pair
#-- perform hybrid bwa men and bwa all mapping, de-duplication, span computation, and tag clustering --#
# mapping
echo  START  ${run} cpu memaln .. >> ${LOGFILE}
echo  Mapping paired tags .. >> ${LOGFILE}
$cpuprog memaln -T ${pair_map_qual} -t ${NTHREAD} $fasta ${run}.$pairlabel.fastq.gz 1> ${run}.$pairlabel.sam 2>> ${LOGFILE}

$pigz -p ${NTHREAD} ${run}.$pairlabel.sam >> ${LOGFILE}
echo  ENDED pair mapping >> ${LOGFILE}

#pairing
echo  STARTED ${run} cpu pair .. >> ${LOGFILE}
echo  Pairing paired tags .. >> ${LOGFILE}
$cpuprog pair -S -q ${pair_map_qual} -s $selfbp -t ${NTHREAD} ${run}.$pairlabel.sam.gz 1>${run}.$pairlabel.stat.xls 2>> ${LOGFILE}
echo  ENDED ${run} cpu pair .. >> ${LOGFILE}

# span
echo  STARTED ${run} cpu span .. >> ${LOGFILE}
echo  Computing span of paired tags .. >> ${LOGFILE}
$cpuprog span -g -t ${NTHREAD} -s $selfbp ${run}.$pairlabel.${pair_suffix}.bam 2>> ${LOGFILE} \
 1>${run}.$pairlabel.${pair_suffix}.span.xls
echo  ENDED ${run} span pair .. >> ${LOGFILE}

# deduplication
echo  STARTED ${run} cpu dedup .. >> ${LOGFILE}
echo  De-duplicating paired tags UU .. >> ${LOGFILE}
$cpuprog dedup -g -t ${NTHREAD} -s $selfbp ${run}.$pairlabel.${pair_suffix}.bam 1> ${run}.$pairlabel.${pair_suffix}.dedup.lc 2>> ${LOGFILE}
echo  ENDED ${run} cpu dedup .. >> ${LOGFILE}

# deduplicated span
echo  STARTED ${run} cpu dedup span.. >> ${LOGFILE}
echo  Computing span of paired tags UU nr .. >> ${LOGFILE}
$cpuprog span -t ${NTHREAD} -s ${selfbp} ${run}.$pairlabel.${pair_suffix}.nr.bam \
  2>> ${LOGFILE} 1>${run}.$pairlabel.${pair_suffix}.nr.span.xls
echo  ENDED ${run} cpu dedup span.. >>  ${LOGFILE}

# cluster tags
echo  STARTED ${run} cpu clustering.. >> ${LOGFILE}
$cpuprog cluster -m -s $selfbp -B 1000 -5 5,0 -3 3,$extbp -j -x -v 1 -g -t ${NTHREAD} -O ${run}.e$extbp \
   ${run}.$pairlabel.${pair_suffix}.nr.bam 1>> ${LOGFILE} 2>> ${LOGFILE}
echo  ENDED ${run} cpu clustering.. >> ${LOGFILE}


mv ${run}.e500.clusters.cis.chiasig.gz ${run}.e500.clusters.cis.gz
mv ${run}.e500.clusters.trans.chiasig.gz ${run}.e500.clusters.trans.gz
	
# Make subset file with intrachrom loops with PET_count >= 2
# for browsing in Excel
cis_file="${run}.e500.clusters.cis.gz"
be3_file="${run}.e500.clusters.cis.BE3"
be2_file="${run}.e500.clusters.cis.BE2"
zcat ${cis_file} | awk '{ if ( $7 >= 3 ) print }' > ${be3_file}
zcat ${cis_file} | awk '{ if ( $7 >= 2 ) print }' > ${be2_file}

#/projects/ruan-lab/USERS/kimm/chia_pet_tool_2/pairix_src/util/bam2pairs/bam2pairs -c ${chrom_sizes} ${run}.singlelinker.paired.UU.nr.bam ${run}
$bam2pair -c ${chrom_sizes} ${run}.singlelinker.paired.UU.nr.bam ${run}

java -Xmx16g -jar $juicer pre -r \
        2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 \
        ${run}.bsorted.pairs.gz ${hic_file} ${chrom_sizes} >> ${LOGFILE}

#None tag
# mapping
echo  STARTED ${run}.$nonelabel cpu memaln .. >> ${LOGFILE}
$cpuprog memaln -T ${none_map_qual} -t ${NTHREAD} $fasta ${run}.$nonelabel.fastq.gz 1> ${run}.$nonelabel.sam 2>> ${LOGFILE}
$pigz -p ${NTHREAD} ${run}.$nonelabel.sam >> ${LOGFILE}
echo  ENDED ${run} cpu memaln .. >> ${LOGFILE}

# pairing
echo Pairing $nonelabel tags .. >> ${LOGFILE}
$cpuprog pair -S -q ${none_map_qual} -t ${NTHREAD} -s ${selfbp} ${run}.$nonelabel.sam.gz 1>${run}.$nonelabel.stat.xls 2>> ${LOGFILE}
echo  ENDED ${run} cpu pair .. >> ${LOGFILE}

# span
echo  STARTED ${run}.$nonelabel cpu span .. >> ${LOGFILE}
$cpuprog span -g -t ${NTHREAD} -s ${selfbp} ${run}.$nonelabel.${none_suffix}.bam 2>> ${LOGFILE} 1>${run}.$nonelabel.${none_suffix}.span.xls
echo  ENDED ${run}.$nonelabel span pair .. >> ${LOGFILE}

# deduplication
echo  STARTED ${run}.$nonelabel cpu dedup .. >> ${LOGFILE}
$cpuprog dedup -g -t ${NTHREAD} -s ${selfbp} ${run}.$nonelabel.${none_suffix}.bam 1> ${run}.$nonelabel.${none_suffix}.dedup.lc 2>> ${LOGFILE}

# deduplicated span
$cpuprog span -t ${NTHREAD} -s ${selfbp} ${run}.$nonelabel.${none_suffix}.nr.bam \
  2>> ${LOGFILE} 1>${run}.$nonelabel.${none_suffix}.nr.span.xls
echo  ENDED ${run} cpu dedup span.. >>  ${LOGFILE}

#1tag
# mapping
echo STARTED ${run}.$singlabel cpu memaln .. >> ${LOGFILE}
$cpuprog memaln -T ${single_map_qual} -t ${NTHREAD} $fasta ${run}.$singlabel.fastq.gz 1> ${run}.$singlabel.sam 2>> ${LOGFILE}
$pigz -p ${NTHREAD} ${run}.$singlabel.sam >> ${LOGFILE}

# pairing
echo Pairing $singlabel tags .. >> ${LOGFILE}
$cpuprog pair -S -q ${single_map_qual} -t ${NTHREAD} -s ${selfbp} ${run}.$singlabel.sam.gz 1>${run}.$singlabel.stat.xls 2>> ${LOGFILE}

# span
$cpuprog span -g -t ${NTHREAD} -s ${selfbp} ${run}.$singlabel.${single_suffix}.bam \
  2>> ${LOGFILE} 1>${run}.$singlabel.${single_suffix}.span.xls

# deduplication
echo STARTED ${run}.$singlabel cpu dedup .. >> ${LOGFILE}
$cpuprog dedup -g -t ${NTHREAD} -s ${selfbp} ${run}.$singlabel.${single_suffix}.bam 1> ${run}.$singlabel.${single_suffix}.dedup.lc 2>> ${LOGFILE}
echo  ENDED ${run}.$singlabel cpu dedup .. >> ${LOGFILE}

# deduplicated span
$cpuprog span -t ${NTHREAD} -s ${selfbp} ${run}.$singlabel.${single_suffix}.nr.bam \
  2>> ${LOGFILE} 1>${run}.$singlabel.${single_suffix}.nr.span.xls
echo  ENDED ${run} cpu dedup span.. >>  ${LOGFILE}

# Convert to bedgraph
# Sort bam for samtools counting bases and for BASIC visualization
if [ ! -f ${run}.for.BROWSER.bam ]
then
    echo -e "`date` Converting file formats..\n" >> ${LOGFILE}
    $scmd view -F 2048 -h ${run}.singlelinker.paired.UU.nr.bam \
        | awk 'length($10) > 30 || $1 ~ /^@/' \
        | $scmd sort -@ 16 - -o ${run}.singlelinker.paired.UU.nr.sorted.bam
    $scmd view -F 2048 -h ${run}.singlelinker.single.UxxU.nr.bam \
        | awk 'length($10) > 30 || $1 ~ /^@/' \
        | $scmd sort -@ 16 - -o ${run}.singlelinker.single.UxxU.nr.sorted.bam
    $scmd view -F 2048 -h ${run}.none.UU.nr.bam \
        | awk 'length($10) > 30 || $1 ~ /^@/' \
        | $scmd sort -@ 16 - -o ${run}.none.UU.nr.sorted.bam

    $scmd merge ${run}.for.BROWSER.bam \
        ${run}.singlelinker.paired.UU.nr.sorted.bam \
        ${run}.singlelinker.single.UxxU.nr.sorted.bam \
        ${run}.none.UU.nr.sorted.bam
    
    # Create BAM index
    $scmd index ${run}.for.BROWSER.bam ${run}.for.BROWSER.bam.bai
    
    # Make bedgraph
    $bcmd genomecov -ibam ${run}.for.BROWSER.bam \
        -bg > ${run}.for.BROWSER.bedgraph

    if [ ${blacklist} != "none" ]
    then
        echo -e "`date` Removing blacklist..\n" >> ${LOGFILE}
        # Exclude blacklist regions
        # Re-name
        mv ${run}.for.BROWSER.bedgraph ${run}.for.BROWSER.orig.bedgraph
        # Subtract blacklist
        $bcmd subtract -a ${run}.for.BROWSER.orig.bedgraph \
            -b ${blacklist} > ${run}.for.BROWSER.bedgraph
        # Delete original file before subtraction
        #rm ${run}.for.BROWSER.orig.bedgraph
    fi

    # Sort bedgraph
    $scom bedSort ${run}.for.BROWSER.bedgraph ${run}.for.BROWSER.sorted.bedgraph
    
    # Make bigwig
    $bedgraph2bigwig ${run}.for.BROWSER.sorted.bedgraph ${chrom_sizes} ${run}.for.BROWSER.bigwig    
fi

if [ ${peak_caller} == "spp" ] || [ ${peak_caller} == "SPP" ]
then
    singularity exec $spp Rscript $addpath/spp_sumner.R ${run}.for.BROWSER.bam ${input_control} . 6 >> ${LOGFILE}
else
    # Call peaks using MACS2 without input control
    singularity exec $macs2 macs2 callpeak --keep-dup all --nomodel -t ${run}.for.BROWSER.bam \
            -f BAM -g hs -n ${run}.no_input_all 1>> ${LOGFILE} 2>> ${LOGFILE}
fi

## Set the output file
out_file=${run}.final_stats.tsv
rm -f ${out_file}

## Get library ID
echo -e "Library_ID\t"${run} >> ${out_file}

# Get library type
echo -e "Library_type\t"${run_type} >> ${out_file}

# Get reference genome
echo -e "Reference_genome\t"${genome} >> ${out_file}

# Get cell type
echo -e "Cell_type\t"${cell_type} >> ${out_file}

# Get IP factor
echo -e "Factor\t"${ip_factor} >> ${out_file}

# Create contact-map URL
hic_file="ChIA-PET_${genome}_${cell_type}_${ip_factor}_${run}_${run_type}_pairs.hic"
url="http://ctencode01.jax.org/chiapet/dev/${hic_file}"

echo -e "Contact-map_URL\t"${url} >> ${out_file}

## PET count
# Get PET count
n_read_pair=$( cat ${run}.stat | grep "Total pairs" | awk -F'[ \t]' '{print $3}' )

## Get linker statistics
read_pair_link=$( cat ${run}.stat | grep "Linker detected" | \
    awk -F '[ \t]' '{print $3}' )

frac_link=$( echo -e "${read_pair_link} / ${n_read_pair}" | bc -l | xargs printf "%.2f\n")

# Write PET count
n_read_pair=$( printf "%'.f\n" ${n_read_pair} )
echo -e "Total_read_pairs\t"${n_read_pair} >> ${out_file}

# Write linker statistics
read_pair_link=$( printf "%'.f\n" ${read_pair_link} )
echo -e "Read_pairs_with_linker\t"${read_pair_link} >> ${out_file}
echo -e "Fraction_read_pairs_with_linker\t"${frac_link} >> ${out_file}

# Write one tag vs two tag
one_tag=$( grep "Single Linker 1 tag (SL/ls)" ${run}.stat | cut -f2 )
two_tag=$( grep "Single Linker 2 tags (SL/ls)" ${run}.stat | cut -f2 )

one_tag=$( printf "%'.f\n" ${one_tag} )
two_tag=$( printf "%'.f\n" ${two_tag} )

echo -e "One_tag\t"${one_tag} >> ${out_file}
echo -e "PET\t"${two_tag} >> ${out_file}

## Mapping
# Get uniquely mapped PET count 
unique=$( cat ${run}.singlelinker.paired.UU.span.xls | grep "Total pairs" | \
    awk -F '[\t]' '{print $2}' )

# Get uniquely mapped and non-redundant PET count 
nr=$( cat ${run}.singlelinker.paired.UU.nr.span.xls | grep "Total pairs" | \
    awk -F '[\t]' '{print $2}' )

# Compute redundancy
redun=$( echo "(${unique} - ${nr}) / ${unique}" | bc -l )

# Write uniquely mapped PET count
unique=$( printf "%'.f" ${unique} )
echo -e "Uniquely_mapped_PET\t"${unique} >> ${out_file}

# Write unique mapped and non-redundant PET count
nr=$( printf "%'.f" ${nr} )
echo -e "Non-redundant_PET\t"${nr} >> ${out_file}

# Write redundancy
redun=$( printf %.2f ${redun} )
echo -e "Redundancy\t"${redun} >> ${out_file}

# Write non-redundant tags
nr_tag=$( samtools view -c ${run}.for.BROWSER.bam )
nr_tag=$( printf "%'.f" ${nr_tag} )
echo -e "Non-redundant_tag\t"${nr_tag} >> ${out_file}

## Get number of peaks
if [ ${peak_caller} == "spp" ] || [ ${peak_caller} == "SPP" ]
then
    n_peak=$( cat ${run}.for.BROWSER.spp.z6.broadPeak | wc -l )
else
    if [ ${input_control} == "none" ]
    then
        n_peak=$( cat ${run}.no_input_all_peaks.narrowPeak | wc -l )
    else
        n_peak=$( cat ${run}.all_peaks.narrowPeak | wc -l )
    fi
fi

n_peak=$( printf "%'.f" ${n_peak} )
echo -e "Peak\t"$n_peak >> ${out_file}


## Interaction types
# Get self-ligation PET count
self_lig=$( cat ${run}.singlelinker.paired.UU.nr.span.xls | \
    grep "second/best<0.95" -A5 | \
    awk -F '[\t]' '{if(NR==4)print $2}' )

self_lig=$( printf "%'.f" ${self_lig} )
echo -e "Self-ligation_PET\t"${self_lig} >> ${out_file}

# Get inter-ligation PET count (intra-chr)
intra_chr_pet=$( cat ${run}.singlelinker.paired.UU.nr.span.xls | \
    grep "second/best<0.95" -A5 | \
    awk -F '[\t]' '{if(NR==5)print $2}' )


# Get inter-ligation PET count (inter-chr)
inter_chr_pet=$( cat ${run}.singlelinker.paired.UU.nr.span.xls | \
    grep "second/best<0.95" -A5 | \
    awk -F '[\t]' '{if(NR==2)print $2}' )

# Compute ratio of intra-chr to inter-chr inter-ligation PETs
pet_ratio=$( echo "${intra_chr_pet} / ${inter_chr_pet}" | bc -l )

# Compute inter-ligation PET count (all)
inter_lig_all=$( echo "${intra_chr_pet} + ${inter_chr_pet}" | bc )

# Write inter-ligation PET count (all)
inter_lig_all=$( printf "%'.f" ${inter_lig_all} )
echo -e "Inter-ligation_PET\t"${inter_lig_all} >> ${out_file}

# Write inter-ligation PET count (intra-chr)
intra_chr_pet=$( printf "%'.f" ${intra_chr_pet} )
echo -e "Intra-chr_PET\t"${intra_chr_pet} >> ${out_file}

# Write inter-ligation PET count (inter-chr)
inter_chr_pet=$( printf "%'.f" ${inter_chr_pet} )
echo -e "Inter-chr_PET\t"${inter_chr_pet} >> ${out_file}

# Write ratio of intra-chr to inter-chr inter-ligation PETs
pet_ratio=$( printf %.2f ${pet_ratio} )
echo -e "ratio_of_intra/inter_PET\t"${pet_ratio} >> ${out_file}


## Singleton
# Get singleton PET count (all)
singleton=$(zcat *clusters*.gz | awk '$7==1{print}' | wc -l)
singleton=$( printf "%'.f" ${singleton} )
echo -e "Singleton\t"$singleton >> ${out_file}

# Get singleton PET count (intra-chr)
intra_singleton=$(zcat *cis.gz | awk '$7==1{print}' | wc -l)
intra_singleton=$( printf "%'.f" ${intra_singleton} )
echo -e "Intra-chr_singleton\t"$intra_singleton >> ${out_file}

# Get singleton PET count (inter-chr)
inter_singleton=$(zcat *trans.gz | awk '$7==1{print}' | wc -l)
inter_singleton=$( printf "%'.f" ${inter_singleton} )
echo -e "Inter-chr_singleton\t"$inter_singleton >> ${out_file}


## Clusters (overall)
# Get cluster count
total_cluster_number=$(zcat *clusters*.gz | awk '$7 != 1{print}' | wc -l)
total_cluster_number=$( printf "%'.f" ${total_cluster_number} )
echo -e "PET_cluster\t"${total_cluster_number} >> ${out_file}

# Get intra-chr cluster count
intra_cluster=$( zcat *cis.gz | awk '$7 >=2 {print}' | wc -l )

# Get inter-chr cluster count
inter_cluster=$( zcat *trans.gz | awk '$7 >=2 {print}' | wc -l)


# Compute ratio of intra-chr to inter-chr clusters
cluster_ratio=$( echo "${intra_cluster} / ${inter_cluster}" | bc -l )
cluster_ratio=$( printf %.2f ${cluster_ratio} )

# Write cluster ratio
echo -e "ratio_of_intra/inter_cluster\t"${cluster_ratio} >> ${out_file}


## Clusters (intra-chr)

# Write intra-chr cluster count
intra_cluster=$( printf "%'.f" ${intra_cluster} )
echo -e "Intra-chr_PET_cluster\t"${intra_cluster} >> ${out_file}

# Get intra-chr cluster count by number of PETs (1 - 10)
for i in $(seq 2 10)
do
	intra_pets_number=$(zcat *cis.gz | \
	    awk -v cutoff=${i} '$7 == cutoff {print}' | wc -l | \
	    xargs printf "%'.f")
	
	echo -e "pets_number_"${i}"\t"${intra_pets_number} >> ${out_file}
done

# Get intra-chr cluster count with > 10 PETs
echo -e "pets_number>10\t"$(zcat *cis.gz | awk '$7 >10 {print}' | \
    wc -l | xargs printf "%'.f") >> ${out_file}

# Write significant intra-chr cluster count (ChiaSig)
#chia_sig=$( cat *.BE2.sigf.interactions | grep -v '#' | wc -l )
#echo -e "ChiaSig\t"${chia_sig} >> ${out_file}

# Write loops with peak support
#if [ ${ip_factor} == 'CTCF' ]; then
#    loops_pk_supp=$( cat ${run}.e500.clusters.cis.BE3.peak_annot.E2 | wc -l )
#else
#    loops_pk_supp=$( cat ${run}.e500.clusters.cis.BE3.peak_annot.E2 | wc -l )
#fi

#echo -e "Loops_with_peak_support\t"${loops_pk_supp} >> ${out_file}

## Clusters (inter-chr)
# Write inter-chr cluster count
inter_cluster=$( printf "%'.f" ${inter_cluster} )
echo -e "Inter-chr_PET_cluster\t"${inter_cluster} >> ${out_file}


# Get inter-chr cluster count by number of PETs (1 - 10)
for i in $(seq 2 10)
do
	inter_pets_number=$(zcat *trans.gz | \
	    awk -v cutoff=${i} '$7 == cutoff {print}' | wc -l | \
	    xargs printf "%'.f")
	
	echo -e "pets_number_"${i}"\t"${inter_pets_number} >> ${out_file}
done

# Get inter-chr cluster count with > 10 PETs
echo -e "pets_number>10\t"$(zcat *trans.gz | \
    awk '$7 >10 {print}' | wc -l | xargs printf "%'.f") >> ${out_file}

# remove intermediated files
rm ${run}*for.BROWSER.bedgraph
rm ${run}*for.BROWSER.orig.bedgraph
rm ${run}*.sam.gz
rm ${run}*.*.fastq.gz
ls ${run}*.singlelinker.paired.*.bam | grep -v "sorted" | xargs rm -f
rm ${run}*cpu.dedup
rm ${run}*none.*.bam
rm ${run}*.singlelinker.single.*.bam 


