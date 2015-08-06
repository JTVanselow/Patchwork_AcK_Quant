#args<-passargs
library("mzR")
library("plyr")
library("XML")

###parse mzData file from Distiller (optimized for quantitation) to obtain peak list
parse.mzData	<- function(dataLocation){
	mzDataFilename 			<- list.files(dataLocation)[grepl(".mzData$",list.files(dataLocation))]
	mzDataFile 				<- paste(c(dataLocation,mzDataFilename),collapse="\\")
	con  					<- file(mzDataFile , open = "r")
		x					<- readLines(mzDataFile)
		x					<- x[!grepl("<int>1</int>",x)]
		x					<- x[!grepl("supDesc|supDataDesc|Charge State|supDataArray",x)]
	close(con)
	mzDataFilename 			<- gsub(".mzData$",".mzData.slim",mzDataFilename)
	mzDataFile 				<- paste(c(dataLocation,mzDataFilename),collapse="\\")
	con  					<- file(mzDataFile , open = "w")
		writeLines(x,con=con)
	close(con)
	mzxml					<- mzDataFile
	aa 					<- openMSfile(mzxml)
	pl 					<- peaks(aa)
	head					<- header(aa)
	doc 					<- xmlTreeParse(mzxml, useInternalNodes=TRUE,ignoreBlanks=FALSE,trim=TRUE)
	r 					<- xmlRoot(doc)
	rootnames				<- as.vector(names(r))
	description				<- xmlToList(r[[which(rootnames %in% c("description"))]])
	sourceInfo				<- as.matrix(unlist(description$admin))
	instrument				<- as.matrix(unlist(description$instrument))
	software				<- data.frame(t(unlist(description$dataProcessing$software)))
	processingMethod			<- unlist(description$dataProcessing$processingMethod)
	names(processingMethod)		<- c()
	processingMethod2			<- data.frame(t(processingMethod[seq(from=2,to=length(processingMethod),by=2)]))
	names(processingMethod2)	<- processingMethod[seq(from=1,to=length(processingMethod)-1,by=2)]
	processingMethod			<- matrix(processingMethod,ncol=2,nrow=length(processingMethod)/2)
	selnodes 				<- getNodeSet(doc, "//spectrumDesc")
	nscans				<- length(selnodes)
	reslist				<- replicate(nscans, list()) 
	for(i in which(head$msLevel>1)){
		reslist[[i]]	<- parseXML(selnodes,i)
	}
	reslist				<- reslist[unlist(lapply(reslist,length))>0]
	dfout					<- ldply(reslist,stringsAsFactors =FALSE)
	spectrumID				<- as.numeric(ldply(xpathApply(doc, "//spectrum", xmlAttrs))[head$msLevel>1,])
	dfout					<- cbind(spectrumID,dfout)
	dfout					<- data.frame(apply(dfout[,-c(3:4)],2,as.numeric))
	close(aa)
	unlink(mzDataFile)
	gc()
	return(list(cbind(dfout,head),pl,mzDataFilename))
}

###manual parse XML
parseXML	<- function(selnodes,i){
	templist				<- xmlToList(selnodes[[i]])
	spectrummslevel			<- templist$spectrumSettings$spectrumInstrument$.attrs[1]
	spectrumInfoMat		<- matrix(unlist(templist$spectrumSettings$spectrumInstrument[-4]),4,3)[3:4,]
	spectrumInfo		<- spectrumInfoMat[2,]
	names(spectrumInfo)	<- spectrumInfoMat[1,]
	spectrumInfo		<- as.data.frame(t(c(spectrumInfo,spectrummslevel)),stringsAsFactors =FALSE)
	acqNumber			<- as.numeric(templist$spectrumSettings$acqSpecification$acquisition)
	precursormslevel		<- templist$precursorList$precursor$.attrs[1]
	precursorInfoMat		<- matrix(unlist(templist$precursorList$precursor$ionSelection),4,3)[3:4,]
	precursorInfo		<- as.numeric(precursorInfoMat[2,])
	names(precursorInfo)	<- precursorInfoMat[1,]
	precursorInfo		<- as.data.frame(t(c(precursorInfo,precursor=precursormslevel)),stringsAsFactors =FALSE)
	scaninfo			<- cbind(acqNumber,spectrumInfo,precursorInfo)
	rownames(scaninfo)	<-c()
	return(scaninfo)
}
###
args 			<- commandArgs(trailingOnly = TRUE)
dataLocation	<- args[1]
outDir		<- args[2]

tempList		<- parse.mzData(dataLocation)
save(tempList,file=paste(c(outDir,"mzData.RData"),collapse="\\"))
