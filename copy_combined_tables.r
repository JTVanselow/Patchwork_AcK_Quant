rm(list=ls(all=TRUE))

###simple pasting path vectors
pastePath	<- function(inVec){
	return(paste(inVec,collapse="\\"))
}

args 			<- commandArgs(trailingOnly = TRUE)
#args			<- passargs

superMainFolder		<- args[1]

timestamp		<- format(Sys.time(), format="%Y%m%d.%H%M")
mainresultsFolder	<- pastePath(c(superMainFolder,"_results"))

tablesFolder	<- pastePath(c(superMainFolder,"_tables"))
if (!file.exists(tablesFolder)){
	    dir.create(tablesFolder)
}
tablesFolder	<- pastePath(c(tablesFolder,timestamp))
if (!file.exists(tablesFolder)){
	    dir.create(tablesFolder)
}

###rename and copy
mainFolder.collection	<- list.files(mainresultsFolder,full.names=TRUE)
mainFolder.collection	<- mainFolder.collection[file_test("-d", mainFolder.collection)]
mainFolder.collection	<- gsub("[/]","\\\\",mainFolder.collection)
mainFolder.collection	<- mainFolder.collection[!grepl("_progress$|_results$",mainFolder.collection)]

for(mainFolder in mainFolder.collection){
	mainName	<- unlist(strsplit(mainFolder,"\\\\"))
	mainName	<- mainName[length(mainName)]
	tempFolders		<- list.files(mainFolder)
	tempFolders		<- tempFolders[!file_test("-d", tempFolders)&!tempFolders %in% "_results"]
	for(projectName in tempFolders){
		projectFolder	<- pastePath(c(mainFolder ,projectName ,"combined"))
		projectFolder	<- pastePath(c(projectFolder ,sort(list.files(projectFolder),decreasing = TRUE)[1]))
		txt.files		<- list.files(projectFolder)
		txt.files		<- txt.files[grepl(".txt$",txt.files)&!grepl(".slim_",txt.files)&!grepl("condensMat_Acetylation_quant",txt.files)]
		oldnames	 	 <- txt.files
		prefix		<- paste(c(mainName,projectName,""),collapse="_")		
		newnames		<- gsub("^",prefix,oldnames)
		for(i in 1:length(oldnames)){
			oldnames[i]	<- pastePath(c(projectFolder,oldnames[i]))
			newnames[i]	<- pastePath(c(tablesFolder,newnames[i]))
		}
		file.copy(oldnames,newnames)
	}
}

###end rename and copy


