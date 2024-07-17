vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))


##  Setting input paths for normalized read count and experimental design ###
assembly_file <- split.vars[1]
Bin_size <- as.numeric(split.vars[2])
Bin_file = split.vars[3]

Bin_size = Bin_size*10^6
options(stringsAsFactors = F)

library(Biostrings)
assembly = readDNAStringSet(assembly_file)

Genome_Coordinate = data.frame(chr = c(paste('scaffold_',c(1:40),sep = ''),'scaffold_X'),
                               start = 1, end = width(assembly)[1:41])

chrList = c(paste('scaffold_',c(1:40),sep = ''),'scaffold_X')

Bin_Coordinate = NULL
for (oneChr in chrList) {
  Segment_Coordinate = data.frame(chr=rep(oneChr,1),
                                  start = 1,
                                  end = max(Genome_Coordinate[which(Genome_Coordinate[,1]==oneChr),3]))
  
  Bin_Coordinate_Start_End=NULL
  
  for (i in 1:nrow(Segment_Coordinate)) {
    numberOFBins = round((Segment_Coordinate[i,3]-Segment_Coordinate[i,2]+1)/Bin_size,digits = 0)
    tempSize = floor((Segment_Coordinate[i,3]-Segment_Coordinate[i,2]+1)/numberOFBins)
    if ((numberOFBins==0)||(numberOFBins==1)) {
      Bin_Coordinate_Start_End = rbind(Bin_Coordinate_Start_End,Segment_Coordinate[i,2:3])
    } else {
      Bin_Coordinate_Start_End_temp = cbind(start = Segment_Coordinate[i,2]+c(0:(numberOFBins-1))*tempSize,
                                            end = Segment_Coordinate[i,2]+c(1:numberOFBins)*tempSize-1)
      Bin_Coordinate_Start_End_temp[numberOFBins,2] = Segment_Coordinate[i,3]
      Bin_Coordinate_Start_End = rbind(Bin_Coordinate_Start_End,Bin_Coordinate_Start_End_temp)
    }
  }
  Bin_Coordinate_Start_End = Bin_Coordinate_Start_End[order(Bin_Coordinate_Start_End[,1]),]
  Bin_Coordinate = rbind(Bin_Coordinate,data.frame(chr = rep(oneChr),
                                                   start = Bin_Coordinate_Start_End[,1],
                                                   end = Bin_Coordinate_Start_End[,2]
                                                   )
                         )
}

Bin_Coordinate = Bin_Coordinate[Bin_Coordinate[,3]-Bin_Coordinate[,2]>1,]
write.table(Bin_Coordinate,file = Bin_file,sep = '\t',quote = F,row.names = F)


