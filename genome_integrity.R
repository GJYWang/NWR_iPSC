vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))

library(DNAcopy)
library(stringr)

bed_file <- split.vars[1]
output_CNV_folder <- split.vars[2]
sample_name <- split.vars[3]
bin_file <- split.vars[4]


options(stringsAsFactors = F)

Bin_Coordinate = read.table(bin_file,sep = '\t',header = T)
Bin_Size = Bin_Coordinate[,3]-Bin_Coordinate[,2]
chrList = unique(Bin_Coordinate[,1])
Bin_Coordinate_List = list()
for (onrChr in chrList) {
  Bin_Coordinate_List[[onrChr]] = Bin_Coordinate[Bin_Coordinate[,1]==onrChr,]
}

Read_Bed = tryCatch(
  read.table(bed_file,
             sep = '\t',
             header = F, 
             quote = ""), 
  error=function(e) NULL)

if (is.null(Read_Bed)) {
  stop("Stopping script execution because no line in bed file.")
}

cat("Counting read count in each bin ... \n")
Read_Bed = Read_Bed[Read_Bed[,5]>=10,]
Read_Bed_List = list()
for (onrChr in chrList) {
  temp = Read_Bed[Read_Bed[,1]==onrChr,]
  temp = temp[order(temp[,2]),]
  Read_Bed_List[[onrChr]] = temp
}

Read_In_Bin = NULL
for (onrChr in chrList) {
  Bin_Coordinate_oneChr = Bin_Coordinate_List[[onrChr]]
  Read_Bed_oneChr = Read_Bed_List[[onrChr]]
  
  Read_In_Bin_oneChr = rep(NA,nrow(Bin_Coordinate_oneChr))
  
  for (i in 1:length(Read_In_Bin_oneChr)) {
    binStart = Bin_Coordinate_oneChr[i,2];binEnd = Bin_Coordinate_oneChr[i,3]
    
    noOverlap = (Read_Bed_oneChr[,3] < binStart) | (binEnd < Read_Bed_oneChr[,2])
    
    Read_In_Bin_oneChr[i] = sum(!noOverlap)
  }
  
  Read_In_Bin = c(Read_In_Bin,Read_In_Bin_oneChr)
}  
Read_In_Bin[Read_In_Bin==0]=1

Count_File = paste(output_CNV_folder,'Read_Count.Rdata',sep = '/')
save(Read_In_Bin,Bin_Coordinate,Bin_Size,file = Count_File)

cat("Segmentation ... \n")
expected_Read_In_Bin = sum(Read_In_Bin) * (Bin_Coordinate$end-Bin_Coordinate$start) / sum(Bin_Coordinate$end-Bin_Coordinate$start)
log_read_Norm = log2(Read_In_Bin/expected_Read_In_Bin)
chr_Num = match(Bin_Coordinate$chr,unique(Bin_Coordinate$chr))
chr_Names = str_pad(chr_Num,ceiling(log10(41)), pad = "0")

CNA.object <- CNA(cbind(log_read_Norm),
                  chr_Names,(Bin_Coordinate$start+Bin_Coordinate$end)/2,
                  data.type="logratio",sampleid=sample_name)
segment.smoothed.CNA.object <- segment(CNA.object, verbose=1)

CNA_result = data.frame(chrom=CNA.object[,1],
                    start=Bin_Coordinate$start,
                    end=Bin_Coordinate$end,
                    value=CNA.object[,3])
CNA_file = paste(output_CNV_folder,'bin_result.txt',sep = '/')
write.table(CNA_result,
            file = CNA_file,
            quote = F,sep = '\t',
            col.names = T, row.names = F)

CNV_segment_tmp = segment.smoothed.CNA.object$output
for (j in 1:nrow(CNV_segment_tmp)) {
  segStart = CNV_segment_tmp[j,3];segEnd = CNV_segment_tmp[j,4]
  which_CNV_loc = (! ( (CNA_result[,3] < segStart) | (segEnd < CNA_result[,2]) )) & (CNA_result[,1]==CNV_segment_tmp[j,2])
  CNV_segment_tmp[j,3] = min(CNA_result$start[which_CNV_loc==T])
  CNV_segment_tmp[j,4] = max(CNA_result$end[which_CNV_loc==T])
}

Segment_File = paste(output_CNV_folder,'segment.txt',sep = '/')
write.table(CNV_segment_tmp,
            file = Segment_File,
            quote = F,sep = '\t',
            col.names = T, row.names = F)