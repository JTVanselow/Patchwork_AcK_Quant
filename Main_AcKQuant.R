
rm(list=ls(all=TRUE))

server	<- TRUE

library("doParallel")

###simple pasting path vectors
pastePath	<- function(inVec){
	return(paste(inVec,collapse="\\"))
}
	
############
###MAIN####
###########

###arguments###
###you need to change following lines according to your needs
reagent.purity	<- 0.99	## 0.99; enter value from LOT certificate
intens.cutoff	<- 0		## 0; intensity cutoff for selecetion of fragment ions (sum of isotope peaks)
corr.cutoff		<- 0.99	## 0.99; min correlation between theoretical and calculated isotope patterns
diff.cutoff		<- 20		## 20; max. tol. deviation from theoretical rel. isotope peak intensity
errortolPPM		<- 7E-6 	## 7E-6; error tolerance (in ppm*E6) for fragment ions
cores			<- 4		## number of CPU cores for parallel processing (requires 1-2 GB of RAM per core, dependent on data complexity)
mainloopcores	<- 2		## 



if(server){
	ScriptFolder	<- "C:\\apps\\R scripts\\_Acetyl.K.Quant\\AcK_Quant"
	superMainFolder	<- "D:\\R-data-analysis\\Rasha"
	cores			<- 6
	mainloopcores	<- 6

}else{
	ScriptFolder	<- "C:\\Users\\jev98xq\\CloudStation\\_Projects\\Postdoc_RVZ\\R scripts\\_Acetyl.K.Quant\\AcK_Quant"
	superMainFolder	<- "C:\\R-data-analysis\\RE_Acetylation\\fusion"

}


timestamp		<- format(Sys.time(), format="%Y%m%d.%H%M")
logFolder		<- pastePath(c(superMainFolder,"_progress"))
if (!file.exists(logFolder)){
		    dir.create(logFolder)
	}
logFolder		<- pastePath(c(logFolder,timestamp))
if (!file.exists(logFolder)){
		    dir.create(logFolder)
	}

mainresultsFolder	<- pastePath(c(superMainFolder,"_results"))
if (!file.exists(mainresultsFolder)){
	    dir.create(mainresultsFolder)
}

##fixed arguments
atomcols		<- c("C","H","N","O","S","Cx")
mainFolder.collection	<- list.files(superMainFolder,full.names=TRUE)
mainFolder.collection	<- mainFolder.collection[file_test("-d", mainFolder.collection)]
mainFolder.collection	<- gsub("[/]","\\\\",mainFolder.collection)
mainFolder.collection	<- mainFolder.collection[!grepl("_progress$|_results$|_tables$",mainFolder.collection)]
todo.df	<- c()

for(mainFolder in mainFolder.collection){
	mainName	<- unlist(strsplit(mainFolder,"\\\\"))
	mainName	<- mainName[length(mainName)]
	tempFolders		<- list.files(mainFolder)
	tempFolders		<- tempFolders[!file_test("-d", tempFolders)&!tempFolders %in% "_results"]
	todo.df	<- rbind(todo.df,data.frame(cbind(mainFolder,tempFolders),stringsAsFactors=FALSE))
}

cl1 <- makeCluster(mainloopcores)
registerDoParallel(cl1)
foreach(i = 1:nrow(todo.df)) %dopar% {
	mainFolder	<- todo.df$mainFolder[i]
	mainName	<- unlist(strsplit(mainFolder,"\\\\"))
	mainName	<- mainName[length(mainName)]
	resultsFolder	<- pastePath(c(mainresultsFolder,mainName))
	if (!file.exists(resultsFolder)){
		    dir.create(resultsFolder)
	}
	projectfolder	<- pastePath(c(mainFolder,todo.df$tempFolders[i]))
	passargs	<- c(paste(atomcols,collapse=";"),reagent.purity,projectfolder,resultsFolder,todo.df$tempFolders[i],errortolPPM,ScriptFolder,intens.cutoff,corr.cutoff,diff.cutoff,cores,logFolder,mainName)
	system(paste(c("RScript.exe --vanilla",shQuote(paste(c(ScriptFolder,"\\execute_acetyl.quant.R"),collapse="")),shQuote(passargs)),collapse=" "))
	sink()
	sink(type = "message")
}
stopCluster(cl1)
print(paste(" done",mainFolder))
closeAllConnections()

print(paste(" merging tables"))
passargs	<- c(superMainFolder)
system(paste(c("RScript.exe --vanilla",shQuote(paste(c(ScriptFolder,"\\copy_combined_tables.R"),collapse="")),shQuote(passargs)),collapse=" "))
	

print(paste(" all done!"))