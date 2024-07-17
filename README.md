# Evaluating the genome integrity of iPSCs for NWR

- [Dependencies](#Dependencies)
- [Installation](#Installation)
- [Usage](#Usage)
- [Output](#Output)
- [Contact](#Contact)
- [Citation](#Citation)


## Dependencies

Make sure minimap2, samtools, and bedtools are installed

Make sure DNAcopy and ggplot2 are installed in R

See example code for dependencies installation in install directory

Download NWR reference genome from: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_021442165.1/


## Installation

No need to install, just download the repository and run pipeline.sh

## Usage
Add or modify variables: assemblyPath, fastq_file, alignment_path, result_path, thread and Bin_Size in pipeline.sh

assemblyPath: the path to the NWR assembly file in fasta format. The fasta file needs to be indexed. 

fastq_file: the path to the sequencing file. Usually in fastq.gz format. 

alignment_path: the path to store alignment results during processing. 

result_path: the path to store all the results.

thread: number of threads used for alignment.

Bin_Size: the size of the bin for genome integrity analysis (in Mbps).


## Output

You can find three files in result_path:

Read_Count.Rdata: The read counts of each bin and the bin file is stored in Rdata format. 

bin_result.txt: The copy number of each bin. 

segment.txt: The segmented copy number results.

genome_fig.pdf: The visualization of genome integrity. 


## Contact

Questions about the paper can be sent to correspondence authors: 

Franz-Josef Müller: fjmuellr@molgen.mpg.de

Marisa L. Korody: MKorody@sdzwa.org

Jeanne F. Loring: jloring@scripps.edu


## Citation

Chromosome-level genome assembly of the functionally extinct northern white rhinoceros (Ceratotherium simum cottoni). Gaojianyong Wang, Björn Brändl, Christian Rohrandt, Karl Hong, Andy Pang, Joyce Lee, Harris A. Lewin, Giovanna Migliorelli, Mario Stanke, Remy Schwab, Sarah Ford, Iris Pollmann, Bernhard M. Schuldt, Marlys Houck, Oliver A. Ryder, Alexander Meissner, Jeanne F. Loring, Franz-Josef Müller, Marisa L. Korody
bioRxiv 2021.12.11.472206; doi: https://doi.org/10.1101/2021.12.11.472206

