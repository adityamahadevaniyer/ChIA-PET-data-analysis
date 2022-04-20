### Config file for ChIA-PET Tool
### Sumner/Singularity version

# The name of the sequencing run
run="LHG0147" 

# The type of sequencing run:
#    "miseq" - around 30 million reads
#    "hiseq" - around 300 million reads
#    "pooled" - around 1 billion reads
run_type="qcseq"

# The factor for which the IP was performed
ip_factor="CTCF"

# Cell type
cell_type="GM12878"

# The directory containing the input FASTQ files
data_dir="/fastscratch/kimm/processing/fastq/${run}"
#data_dir="/projects/ruan-lab/ChIA-PIPE_Sumner/fastq/${run}"

if [ -f ${data_dir}/merged*_R1.fastq.gz ]
then
    r1_fastq="${data_dir}/merged*_R1.fastq.gz"
    r2_fastq="${data_dir}/merged*_R2.fastq.gz"
else
    r1_fastq="${data_dir}/*_R1*.fastq.gz"
    r2_fastq="${data_dir}/*_R2*.fastq.gz"
fi

# The name of the primary genome
# For example: "hg19", "hg38", "dm3", "mm9", "mm10"
genome="hg38"

# The reference genome FASTA file for aligning the reads
# (The same directory must also contain the BWA index files)
fasta="/fastscratch/kimm/processing/genomes/hg38/hg38.fa"

# The chrom.sizes file from UCSC Genome browser
# for the relevant genome build
chrom_sizes="/fastscratch/kimm/processing/genomes/hg38/hg38.chrom.sizes"

# The BAM file for the ChIP-seq input control
# (Required for spp; not required for macs2)
# If not available, set to "none"
#input_control="/projects/ruan-lab/processing/input_controls/hg38/HUVEC_input_hg38_CHH0003.bam"
input_control="/fastscratch/kimm/processing/input_controls/hg38/GM12878_input_hg38_CHG0002.bam"

# The peak-calling algorithm ("macs2" or "spp")
peak_caller="spp"

# The folder in BASIC browser to which to upload the tracks
basic_folder="Testing_kimm"

# ENCODE blacklist
# If not available, set to "none"
blacklist="/fastscratch/kimm/processing/genomes/hg38/hg38_blacklist.bed"
#blacklist="/projects/ch-lee-lab/processing/genomes/mm10/mm10_blacklist.bed"

if [ ${run_type} == "miseq" ] || [ ${run_type} == "novaseqsp" ] || [ ${run_type} == "qcseq" ]
then
    NTHREAD="20"
    hrs="12"
    GB="40"
elif [ ${run_type} == "hiseq" ] || [ ${run_type} == "nextseq" ] || [ ${run_type} == "novaseq" ]
then
    NTHREAD="30"
    hrs="72"
    GB="60"
else
    NTHREAD="40"
    hrs="72"
    GB="80"
fi
