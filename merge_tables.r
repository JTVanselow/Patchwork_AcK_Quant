rm(list=ls(all=TRUE))

###simple pasting path vectors
pastePath	<- function(inVec){
	return(paste(inVec,collapse="\\"))
}

##arguments
superMainFolder	<- "D:\\R-data-analysis\\Rasha"

timestamp		<- format(Sys.time(), format="%Y%m%d.%H%M")
mainresultsFolder	<- pastePath(c(superMainFolder,"_results"))

tablesFolder	<- pastePath(c(superMainFolder,"_tables"))

tablesFolder	<- sort(list.files(tablesFolder,full.names = TRUE),decreasing = TRUE)[1]
mergedFolder	<- pastePath(c(superMainFolder,"_merged_tables"))
unlink(mergedFolder,recursive = TRUE)
if (!file.exists(mergedFolder)){
	    dir.create(mergedFolder)
}

tobe.merged.folders	<- list.files(tablesFolder,full.names = TRUE)
tobe.merged.names		<- list.files(tablesFolder,full.names = FALSE)

for(i in 1:length(tobe.merged.folders)){

	merge.me<-tobe.merged.folders[i]
	all.files			<- list.files(merge.me,full.names = TRUE)
	Combined_Acetylation	<- all.files[grepl("Combined_Acetylation",all.files)]
	pepMat_noFilt		<- all.files[grepl("pepMat_noFilt",all.files)]

	###merged site quant table
	outFilename		<- paste(c(mergedFolder,"\\",tobe.merged.names[i],"_merged_Combined_Acetylation.txt"),collapse="")
	unlink(outFilename)
	for(ii in 1:length(Combined_Acetylation)){
		tempMat	<- read.delim(Combined_Acetylation[ii], sep="\t",header=TRUE,stringsAsFactors=FALSE,na.strings = c("NA","n. def."))
		tempMat	<- tempMat[grepl("^agms|histone|Histone",tempMat$accession)|grepl("histone|Histone",tempMat$prot_desc),]
		tempMat	<- tempMat[!is.na(tempMat$filt.ambig2),]
		prefix	<- gsub(".*[/]|_Combined_Acetylation.*|_pepMat_noFilt.*","",Combined_Acetylation[ii])
		tempMat$experiment	<- prefix
		exclude.colnames	<- which(colnames(tempMat) %in% unlist(strsplit("selcol,filt.uni1,pep_local_mod_pos,pep_var_mod.2",",")))
		exclude.colnames	<- c(exclude.colnames,which(colnames(tempMat)=="selcol"):which(colnames(tempMat)=="ms2mz"),
				which(colnames(tempMat)=="z"):which(colnames(tempMat)=="Ac.mods.pos"),
				which(colnames(tempMat)=="C"):which(colnames(tempMat)=="Cx"),
				which(colnames(tempMat)=="minus1"):which(colnames(tempMat)=="intens"))
		tempMat	<- tempMat[,-exclude.colnames]
		if(ii==1){
			write.table(tempMat, file = outFilename,quote = FALSE, sep = "\t",na = "NA", dec = ".", col.names=TRUE,row.names = FALSE)
		}else{
			write.table(tempMat, file = outFilename,append=TRUE,quote = FALSE, sep = "\t",na = "NA", dec = ".", col.names=FALSE,row.names = FALSE)
		}
	}
	
	mergedMat2	<- read.delim(outFilename, sep="\t",header=TRUE,stringsAsFactors=FALSE,na.strings = c("NA","n. def."))
	unlink(outFilename)
	selcolnames					<- c("accession","AcK","pep_var_mod.4")
	mergedMat2$combselcol			<- do.call(paste, c(mergedMat2[,selcolnames], sep=";;;"))
	temptable					<- data.frame(table(mergedMat2$combselcol[mergedMat2$filt2==3]))
	colnames(temptable)[2]			<- "n.filt2"
	mergedMat2	<- merge(mergedMat2,temptable,by.x="combselcol",by.y="Var1")
	mergedMat2	<- mergedMat2[,!colnames(mergedMat2)%in%"combselcol"]
	mergedMat2	<- mergedMat2[order(-mergedMat2$Freq,mergedMat2$accession,as.numeric(gsub("N","-",mergedMat2$AcK)),-mergedMat2$n.filt2,-mergedMat2$filt2,mergedMat2$pep_var_mod.4,mergedMat2$pep_start,mergedMat2$pep_end),]
	write.table(mergedMat2, file = outFilename, quote = FALSE, sep = "\t",na = "", dec = ".", row.names = FALSE)


	###merged peptides table
	outFilename		<- paste(c(mergedFolder,"\\",tobe.merged.names[i],"_pepMat_noFilt.txt"),collapse="")
	unlink(outFilename)
	for(ii in 1:length(pepMat_noFilt)){
		tempMat	<- read.delim(pepMat_noFilt[ii], sep="\t",header=TRUE,stringsAsFactors=FALSE,na.strings = c("NA","n. def."))
		tempMat	<- tempMat[grepl("^agms|histone",tempMat$accession)|grepl("histone|Histone",tempMat$prot_desc),]

		prefix	<- gsub(".*[/]|_Combined_Acetylation.*|_pepMat_noFilt.*","",pepMat_noFilt[ii])
		tempMat$experiment	<- prefix
		exclude.colnames	<- which(colnames(tempMat) %in% unlist(strsplit("selcol",",")))
		
		tempMat	<- tempMat[,-exclude.colnames]
		if(ii==1){
			write.table(tempMat, file = outFilename,quote = FALSE, sep = "\t",na = "NA", dec = ".", col.names=TRUE,row.names = FALSE)
		}else{
			write.table(tempMat, file = outFilename,append=TRUE,quote = FALSE, sep = "\t",na = "NA", dec = ".", col.names=FALSE,row.names = FALSE)
		}
	}
	
	mergedMat2	<- read.delim(outFilename, sep="\t",header=TRUE,stringsAsFactors=FALSE,na.strings = c("NA","n. def."))
	unlink(outFilename)
	selcolnames					<- c("accession","pep_var_mod.3")
	mergedMat2$combselcol			<- do.call(paste, c(mergedMat2[,selcolnames], sep=";;;"))
	temptable					<- data.frame(table(mergedMat2$combselcol))
	colnames(temptable)[2]			<- "n.pep_vm3"
	mergedMat2	<- merge(mergedMat2,temptable,by.x="combselcol",by.y="Var1")
	mergedMat2	<- mergedMat2[,!colnames(mergedMat2)%in%"combselcol"]
	mergedMat2	<- mergedMat2[order(-mergedMat2$n.pep_vm3,mergedMat2$accession,mergedMat2$pep_var_mod.4,mergedMat2$pep_start,mergedMat2$pep_end),]
	write.table(mergedMat2, file = outFilename, quote = FALSE, sep = "\t",na = "", dec = ".", row.names = FALSE)
}

