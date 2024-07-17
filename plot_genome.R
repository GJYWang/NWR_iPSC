vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars,","))

##  Setting input paths for normalized read count and experimental design ###
output_CNV_folder <- split.vars[1]

Bin_File = paste(output_CNV_folder,'Read_Count.Rdata',sep = '/')
Segment_File = paste(output_CNV_folder,'segment.txt',sep = '/')
CNA_file = paste(output_CNV_folder,'bin_result.txt',sep = '/')
figure_file = paste(output_CNV_folder,'genome_fig.pdf',sep = '/')

options(stringsAsFactors = F)

library(ggplot2)

chr_to_plot = 41

#Bin_Coordinate = read.table(Bin_File,sep = '\t',header = T)
load(Bin_File)
CNV_segment_all = read.table(Segment_File,sep = '\t',header = T)
CNV_data = read.table(CNA_file,sep = '\t',header = T)
accumulate_start = rep(0,chr_to_plot)
for (i in 2:chr_to_plot) {
  previous_chr = i-1
  accumulate_start[i] = accumulate_start[i-1]+max(Bin_Coordinate[Bin_Coordinate[,1]==paste0('scaffold_',previous_chr),3])+10
}
max_coor = accumulate_start[chr_to_plot] + max(Bin_Coordinate[Bin_Coordinate[,1]==paste0('scaffold_','X'),3])

tick_pos = NULL
tick_label = paste('',c(1:chr_to_plot),sep = '')
tick_label[which(tick_label=='41')]='X'
CNV_segment_tmp = CNV_segment_all
CNV_segment_chr_tmp = CNV_segment_tmp[,2]
CNV_segment_tmp = data.frame(CNV_segment_tmp,loc.start.accumu = CNV_segment_tmp$loc.start+accumulate_start[CNV_segment_chr_tmp])
CNV_segment_tmp = data.frame(CNV_segment_tmp,loc.end.accumu = CNV_segment_tmp$loc.end+accumulate_start[CNV_segment_chr_tmp])
for (i in 1:chr_to_plot) {
  tmp = CNV_segment_tmp[CNV_segment_chr_tmp==i,]
  tick_pos = c(tick_pos,(min(tmp$loc.start.accumu)+max(tmp$loc.end.accumu))/2)
}
tick_pos = tick_pos[ (c(1:length(tick_pos)) %%2) ==1 ]
tick_label = tick_label[ (c(1:length(tick_label)) %%2) ==1 ]
##


##
chr_tmp_CNA = CNV_data[,1]
CNV_data = data.frame(CNV_data,loc.start.accumu = CNV_data$start+accumulate_start[chr_tmp_CNA])
CNV_data = data.frame(CNV_data,loc.end.accumu = CNV_data$end+accumulate_start[chr_tmp_CNA])


fontSize = 3
lineSize = 0.1
lineColor = 'blue'


titleStr=paste("Genome Integrity of",CNV_segment_all[1,1])

bin_str = round((Bin_Coordinate[1,3]-Bin_Coordinate[1,2])/1000000,digits = 0)
titleStr = paste(titleStr,paste0('Bin Size: ',bin_str,"Mb"), sep = '  ')

plotList = ggplot() +
  geom_segment(data = CNV_segment_tmp, aes(x = loc.start.accumu,
                                           y = seg.mean, 
                                           xend = loc.end.accumu,
                                           yend = seg.mean),linewidth = lineSize*3, color='red') +
  geom_segment(data = CNV_data, aes(x = loc.start.accumu, 
                                    y = value, 
                                    xend = loc.end.accumu, 
                                    yend = value),linewidth = lineSize*1, color='darkgrey') +
  ylab('Copy number') + labs(title = titleStr) + 
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=fontSize*2),
        axis.title.y = element_text(size=fontSize*3),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(size=fontSize*2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) + 
  scale_x_continuous(breaks = tick_pos, 
                     labels = tick_label
  ) +
  coord_cartesian(xlim = c((max_coor-1)*0.04+1, (max_coor-1)*0.96+1),ylim=c(-1.5,1.5))

plotList = plotList + geom_hline(yintercept = 0, linetype="dotdash", 
                                 color = 'black', linewidth=lineSize/3)

plotList = plotList + geom_vline(xintercept = 1, linetype="longdash", 
                                 color = lineColor, linewidth=lineSize/4)

for (i in 1:chr_to_plot) {
  oneChr = paste0('scaffold_',i)
  if (i==41) {oneChr='scaffold_X'}
  tmp = Bin_Coordinate[Bin_Coordinate[,1]==oneChr,]
  chrEnd = max(tmp[,3])
  plotList = plotList + geom_vline(xintercept = chrEnd+accumulate_start[i], linetype="longdash", 
                                   color = lineColor, linewidth=lineSize/4)
  
}

ggsave(figure_file, plot = plotList,
       height = 1.2, width = 6.5, units = 'in')
