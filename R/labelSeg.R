#function to give label as to what side of centromere

# ##use hg19 build for coordinates of centromeres:
# hg19chromosomes<-read.table("~/Documents/Sequencing/CancerTiming/chromosomeRegions.txt",header=TRUE)
# # hg19chromosomes<-read.table("~/chromosomeRegions.txt",header=TRUE)
# names(hg19chromosomes)[c(1:3)]<-c("chr","start","end")
# hg19chromosomes$start[hg19chromosomes$label=="q"]<-hg19chromosomes$start[hg19chromosomes$label=="q"]+1
# hg19chromosomes$end[hg19chromosomes$label=="p"]<-hg19chromosomes$end[hg19chromosomes$label=="p"]-1
# 
# ##My file, the p is *first* not the *smallest*, so need to fix that
# ##But in fact, this is equivalent for how the chromosomes are ordered
# hg19chromosomes$width<-hg19chromosomes$end-hg19chromosomes$start+1
# hg19chromosomes<-do.call("rbind",by(hg19chromosomes,list(hg19chromosomes$chr),function(x){
# 	whp<-which(x$label=="p")
# 	whq<-which(x$label=="q")
# 	if(x$width[whp]>x$width[whq]){ #switch the labels
# 		x$label[whp]<-"q"
# 		x$label[whq]<-"p"
# 	}
# 	return(x)
# 	
# }))
# row.names(hg19chromosomes)<-NULL
# save(hg19chromosomes,file="~/Documents/RfunctionsGenerally/InternalRPackages/cancerTiming/data/hg19chromosomes.rdata")

labelSeg<-function(chr,start,end,pctOv=0.10){
	require(GenomicRanges)
	data(hg19chromosomes,package="cancerTiming")
	centGr<-GRanges(hg19chromosomes$chr,IRanges(hg19chromosomes$start,hg19chromosomes$end),label=hg19chromosomes$label)
	segGr<-GRanges(paste("chr",chr,sep=""),IRanges(start,end))
	ov<-findOverlaps(subject=segGr,query=centGr)
	subsetByOverlaps(query=segGr,subject=centGr)
	shits <- segGr[subjectHits(ov)]
	chits <- centGr[queryHits(ov)]
	mint <- pintersect(shits, chits)
	spercent <- width(mint) / width(shits)	
	df<-data.frame(index=queryHits(ov),pct=spercent,val=values(centGr)$label[queryHits(ov)])
	lab<-by(df,list(subjectHits(ov)),
		function(x){
			wh<-which(x$pct>=pctOv)
			paste(sort(values(centGr)$label[x$index[wh]]),collapse="",sep="")
		})
	labInd<-as.numeric(names(lab)) #just to be safe.
	lab<-as.vector(lab)
	crossCent<-intersect(intersect(grep("centromere",lab),grep("q",lab)),grep("p",lab))
	lab[crossCent]<-"pq"
	whNotCross<-setdiff(grep("centromere",lab),which(lab=="centromere")) #take out those entirely in the centromere
	lab[whNotCross]<-gsub("centromere","",lab[whNotCross])
	lab[lab=="centromere"]<-"c"
	outdf<-data.frame(chr=chr,start=start,end=end)
	outdf$label<-NA
	outdf$label[labInd]<-lab
	return(outdf$label)
}
numChromosome<-function(chr){
	chr[chr=="X"]<-"23"
	chr[chr=="Y"]<-"24"
	chrN<-as.numeric(chr)
	return(chrN)
}