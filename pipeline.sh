#########
# Make sure minimap2, samtools, and bedtools are installed
# Make sure DNAcopy and ggplot2 are installed in R
# See example code for installation
#########

export PATH="/PATH/TO/TOOLS:$PATH"

#########
# Reference genome name and its path
# NWR reference genome can be downloaded from:
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_021442165.1/
REF="CerSimCot1"
assemblyPath="/PATH/TO/REF_GENOME"
##################################
# Remember to index the reference genome, e.g.,
# samtools faidx $assemblyPath

# Path to your sequence file, usually in fastq.gz format
fastq_file="/PATH/TO/SEQUENCE_FILE"
###################################
# PATH to the alignment result
alignment_path="/PATH/TO/ALIGNMENT_RESULT"
##########################################
# PATH to the iPSC integrity result
result_path="/PATH/TO/GENOME_INTEGRITY_RESULT"
##########################################
# Set the number of threads used for alignment
thread=40
#########

# Set the bin size in Mbps
Bin_Size=3
##########

# Get the file name of sequence file
name="$(basename "${fastq_file}")"
name=${name%.fastq*}
###################

## running minimap2
mkdir -p $alignment_path
if [ ! -f $alignment_path/$name.bed ]
then
  minimap2 -ax map-ont -t $thread $assemblyPath $fastq_file | \
  samtools sort -o $alignment_path/$name.bam -@ $thread -T $alignment_path/$name.tmp

  ## running samtools and bedtools
  samtools view -F 2308 -bo $alignment_path/$name.sorted.bam $alignment_path/$name.bam -@ $thread
  samtools index -@ $thread $alignment_path/$name.sorted.bam
  bedtools bamtobed -i $alignment_path/$name.sorted.bam > $alignment_path/$name.bed
fi
## Initiate bin files based using reference genome
Coord_PATH="bin_data/${REF}"
mkdir -p $Coord_PATH
Bin_FILE="${Coord_PATH}/Bin_Coordi.${REF}.${Bin_Size}Mb.txt"
if [ ! -f $Bin_FILE ]
then
  Rscript Bin_Initial.R $assemblyPath,$Bin_Size,$Bin_FILE
fi

## Running CNV analysis pipeline to get genome integrity of the sequencing data
mkdir -p $result_path
Rscript genome_integrity.R $alignment_path/$name.bed,$result_path,$name,$Bin_FILE
Rscript plot_genome.R $result_path
