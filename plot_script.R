
library(pafr)
library(tidyverse)
library(gggenes)


args = commandArgs(TRUE)
#               1       2       3      4       5        6      7       8       9       10     11     12    13      14     15       16      17      18     19      20 
colnames = c("qname","qlen","qstart","qend","strand","tname","tlen","tstart","tend","nmatch","aln","mapq","tag1","tag2", "tag3", "tag4", "tag5", "tag6", "tag7", "tag8")

fn = args[1]
outdir= args[2]
hap = args[3]


t_paf = read.table(fn, header=FALSE, col.names=colnames, sep="\t", quote="",comment="") %>% 
  separate(tag1,sep=":",into=c("gi1","gi2","id")) %>%
  mutate(id=as.numeric(id))

plot_contig = function(alns, contig, min_aln_per_contig=0){

  data_subset = alns %>% filter(tname==contig) %>%
              group_by(qname) %>%
              mutate(total_alned = sum(tend-tstart)) %>%
              arrange(total_alned) %>%
              filter(total_alned>min_aln_per_contig)

  title=paste("aligned to ",contig,"\nfiltered for query contig >",min_aln_per_contig,"aligned",sep="")
  
  g=ggplot(data_subset)
  g=g+geom_segment(aes(y=qname,yend=qname,x=tstart,xend=tend,color=id),size=2.5)+
    geom_gene_arrow(aes(y=qname,
                        xmin=tstart,
                          xmax=tend,
                          forward=strand=="+"),
                      arrowhead_height = unit(1, "mm"), 
                      arrowhead_width = unit(1, "mm"),
                      arrow_body_height = unit(.5, "mm"),
                      position=position_nudge(x=0,y=-.2),
                      size=.1)+
    ggtitle(title)+
    scale_color_viridis_c()+
    scale_size_continuous(range = c(0.01,2))+
    theme_bw(base_size=6)
    
    return(list(g=g, data_subset=data_subset))
}


filtered_paf_list = list()

i = 1
for (contig in unique(t_paf$tname)){
#for (contig in unique(t_paf$tname)[1:3]){  
  print(contig)
  ret = plot_contig(t_paf,contig,0)
  g = ret$g

  fn_out = paste(outdir,"/pdfs/",contig,".filt_0_aln",".pdf",sep="")
  pdf(fn_out,width=4,height=12)
  print(g)
  dev.off()
  ret = plot_contig(t_paf,contig,10e6)
  fn_out = paste(outdir,"/pdfs/",contig,".filt_10e6_aln",".pdf",sep="")
    
  g = ret$g
  pdf(fn_out,width=4,height=2)
  print(g)
  dev.off()
    
  filtered_paf_list[[i]] = ret$data_subset
  i=i+1
    
}

df_out = do.call("rbind", filtered_paf_list) %>%
            dplyr::select(-c("gi1","gi2","id"))

fn_df_out = paste(outdir,"/",hap,"_filt_10e6_aln",".paf",sep="")
write.table(df_out, fn_df_out, sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)

#fn=paste(paste(outdir,"/dummy.txt",sep=""))
#write.table(data.frame(),fn)



 
