#seqdata <- file.path(snakemake@input[["fastp_json"]])
seqdata <- c("/home/mostafa.abdel/aProjects/gitProjects/ClosTyper/123/2_pipeline_results/ST8-RT002-2-IP-e-D19MD028/ST8-RT002-2-IP-e-D19MD028.fastp.json", "/home/mostafa.abdel/aProjects/gitProjects/ClosTyper/123/2_pipeline_results/ST8-RT002-2-IP-e-D19MD036/ST8-RT002-2-IP-e-D19MD036.fastp.json", "/home/mostafa.abdel/aProjects/gitProjects/ClosTyper/123/2_pipeline_results/ST8-RT002-2-IP-e-D19MD370/ST8-RT002-2-IP-e-D19MD370.fastp.json")

#ref_bases=read from the reference genome size TODO
ref_bases=as.numeric(2500000)
percent <- function(x, digits = 2, format = "f") {
  paste0(formatC(x * 100, format = format, digits = digits), "%")
}

fastp_results=function(file,max=9){
  Sample = basename(dirname(file))
  total_bases=as.numeric(fromJSON(file=file)$summary$before_filtering$total_bases)
  total_reads=as.numeric(fromJSON(file=file)$summary$before_filtering$total_reads)
  gc_content=percent(fromJSON(file=file)$summary$before_filtering$gc_content)
  q30_rate=percent(fromJSON(file=file)$summary$before_filtering$q30_rate)
  q20_rate=percent(fromJSON(file=file)$summary$before_filtering$q20_rate)
  read1_mean_length=as.numeric(fromJSON(file=file)$summary$before_filtering$read1_mean_length)
  read2_mean_length=as.numeric(fromJSON(file=file)$summary$before_filtering$read2_mean_length)
  duplication_rate=percent(fromJSON(file=file)$duplication$rate)
  depth=paste0(formatC(total_bases/ref_bases), "X")
  return(c(Sample,total_bases,total_reads,gc_content,q30_rate,q20_rate,read1_mean_length,read2_mean_length,duplication_rate, depth))
}

dataList<- lapply(seqdata,fastp_results)
res<- do.call("rbind",dataList)
colnames(res)=c("Sample","Total.bases", "Total.reads", "GC.content", "%Q30", "%Q20", "R1.mean_length", "R2_mean_length", "%Duplication","Depth")
