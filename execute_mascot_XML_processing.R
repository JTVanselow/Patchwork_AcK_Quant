library("XML")
library("plyr")
###Read Mascot results (save as: XML, include peptide start and end values and masses), processed with optimized method for Mascot ID
parse.Mascot.XML	<- function(dataLocation,atomcols){
	XMLFilename		<- list.files(dataLocation)[grepl(".xml$",list.files(dataLocation))]	
	XMLFile 		<- paste(c(dataLocation,XMLFilename),collapse="\\")
	mascotXML		<- XMLFile
	doc 			<- xmlTreeParse(mascotXML, useInternalNodes=TRUE,ignoreBlanks=FALSE,trim=TRUE)
	r 			<- xmlRoot(doc)
	rootnames		<- as.vector(names(r))
	hits			<- getNodeSet(doc,"//r:hit",namespaces = c(r="http://www.matrixscience.com/xmlns/schema/mascot_search_results_2"))
	hits.number		<- as.numeric(sapply(hits, xmlGetAttr, "number"))
	pepMat		<- c()
	for(hit.i in hits.number){
		test			<- xmlToList(hits[[hit.i]])
		protein.number	<- which(names(test) %in% "protein")
		for(prot.i in protein.number){
				test2			<- test[[prot.i]]
				peptide.number	<- which(names(test2) %in% "peptide")
				protdesc.number	<- 1:min(peptide.number-1)
				pepinf		<- do.call(rbind, test2[peptide.number])
				rownames(pepinf)	<- c()
				pepattr		<- unlist(do.call(rbind,pepinf[,".attrs"]))
				pepinf		<- pepinf[,!(colnames(pepinf)==".attrs")]
				pepinf		<- data.frame(rbind(pepinf),stringsAsFactors=FALSE)
				pepinf		<- cbind(pepinf,pepattr)
				empty.prot		<- sapply(test2[protdesc.number],function(x) length(x)==0)
				test2[protdesc.number][empty.prot]	<- NA
				protinf		<- unlist(test2[protdesc.number])
				protattr		<- test2[[".attrs"]]
				if(!"member" %in% names(protattr) ){
					protattr	<- c(protattr,member=NA)
				}
				protinf		<- c(protattr,protinf)
				pepMat		<- rbind(pepMat,merge(rbind(protinf),pepinf))
		}
	}
	pepMat				<- data.frame(apply(pepMat,2,function(x) gsub("NULL",NA,x)),stringsAsFactors=FALSE)
	pepMat				<- cbind(pepMat[,-which(colnames(pepMat)=="pep_scan_title")],ldply(apply(rbind(pepMat$pep_scan_title),1,function(x) strsplit(x," "))[[1]]))
	colnames(pepMat)[(length(colnames(pepMat))-4):length(colnames(pepMat))]	<- c("spectrumID","rem","acqNumber","rt","raw")
	pepMat				<- pepMat[,-which(colnames(pepMat)=="rem")]
	pepMat$spectrumID	<- as.numeric(gsub(":","",pepMat$spectrumID))
	pepMat$rt			<- as.numeric(gsub("\\(rt=|\\)","",pepMat$rt))
	pepMat$acqNumber	<- as.numeric(pepMat$acqNumber)
	pepMat$query		<- as.numeric(pepMat$query)
	pepMat$rank			<- as.numeric(pepMat$rank)
	namespaces			<- c(r="http://www.matrixscience.com/xmlns/schema/mascot_search_results_2")
	var.mod.id				<- unlist(xpathApply(doc, "//r:variable_mods/r:modification", xmlAttrs,namespaces = namespaces))
	var.mod.names			<- unlist(xpathApply(doc, "//r:variable_mods/r:modification/r:name", xmlValue,namespaces = namespaces))
	var.mod.delta			<- unlist(xpathApply(doc, "//r:variable_mods/r:modification/r:delta", xmlValue,namespaces = namespaces))
	var.mod					<- data.frame(var.mod.id,var.mod.names,var.mod.delta,stringsAsFactors = FALSE)

	########!!!!!!!!!!!!!!!#############
	##Enter here atomic composition of variable modifications!!!
	var.mod.comp			<- matrix(NA,nrow=nrow(var.mod),ncol=length(atomcols)) ##change to use atomcols
	colnames(var.mod.comp)	<- atomcols
	for(i in grep("^Acetyl ",var.mod[,2])){
		var.mod.comp[i,]	<- c(2,2,0,1,0,0)
	}
	var.mod.comp[var.mod[,2] ==  "Gln->pyro-Glu (N-term Q)",]			<- c(0,-3,-1,0,0,0)
	var.mod.comp[var.mod[,2] ==  "Oxidation (M)",]						<- c(0,0,0,1,0,0)
	var.mod.comp[var.mod[,2] ==  "Methyl (KR)",]						<- c(1,2,0,0,0,0)
	var.mod.comp[var.mod[,2] ==  "Dimethyl (KR)",]						<- c(2,4,0,0,0,0)
	var.mod.comp[var.mod[,2] ==  "Methyl (K)",]							<- c(1,2,0,0,0,0)
	var.mod.comp[var.mod[,2] ==  "Dimethyl (K)",]						<- c(2,4,0,0,0,0)
	var.mod.comp[var.mod[,2] ==  "Trimethyl (K)",]						<- c(3,6,0,0,0,0)
	var.mod.comp[var.mod[,2] ==  "Acetyl 13C1 me (K)",]					<- c(2,4,0,1,0,1) ##special modification, artefact
	var.mod.comp[var.mod[,2] ==  "Acetyl 13C1 me (Protein N-term)",]	<- c(2,4,0,1,0,1) ##special modification, artefact
	var.mod.comp[var.mod[,2] ==  "AcQnt Acetyl 12C1 (O)",]				<- c(1,0,0,0,0,-1) ##special modification, artefact
	var.mod.comp[var.mod[,2] ==  "AcQnt Acetyl 12C1 (Protein N-term)",]	<- c(1,0,0,0,0,-1) ##special modification, artefact
	var.mod.comp[var.mod[,2] ==  "AcQnt di-methyl (O)",]				<- c(1,2,0,-1,0,-1) ##special modification, artefact
	var.mod.comp[var.mod[,2] ==  "AcQnt di-methyl (Protein N-term)",]	<- c(1,2,0,-1,0,-1) ##special modification, artefact
	var.mod.comp[var.mod[,2] ==  "AcQnt methyl (O)",]					<- c(1,2,0,0,0,0) ##special modification, artefact
	var.mod.comp[var.mod[,2] ==  "AcQnt methyl (Protein N-term)",]		<- c(1,2,0,0,0,0) ##special modification, artefact
	var.mod.comp[var.mod[,2] ==  "AcQnt tri-methyl (O)",]				<- c(2,4,0,-1,0,-1) ##special modification, artefact
	var.mod.comp[var.mod[,2] ==  "AcQnt tri-methyl (Protein N-term)",]	<- c(2,4,0,-1,0,-1) ##special modification, artefact
 
	var.mod					<- cbind(var.mod,var.mod.comp)
	var.mod$var.mod.names	<- gsub("Acetyl 13C1 me","Monomethyl",var.mod$var.mod.names)
	pepMat$pep_var_mod		<- gsub("Acetyl 13C1 me","Monomethyl",pepMat$pep_var_mod)


	fix.mod.id				<- unlist(xpathApply(doc, "//r:fixed_mods/r:modification", xmlAttrs,namespaces = namespaces))
	fix.mod.names			<- unlist(xpathApply(doc, "//r:fixed_mods/r:modification/r:name", xmlValue,namespaces = namespaces))
	fix.mod.delta			<- unlist(xpathApply(doc, "//r:fixed_mods/r:modification/r:delta", xmlValue,namespaces = namespaces))
	fix.mod					<- data.frame(fix.mod.id,fix.mod.names,fix.mod.delta)
	masses.names			<- unlist(xpathApply(doc, "//r:masses/r:mass", xmlAttrs,namespaces = namespaces))
	masses					<- as.numeric(unlist(xpathApply(doc, "//r:masses/r:mass", xmlValue,namespaces = namespaces)))
	names(masses)			<- masses.names
	pepMat$pep_var_mod.2	<- apply(cbind(gsub("\\.","",pepMat$pep_var_mod_pos)),1,function(x)	{
			y	<- var.mod[as.numeric(unlist(strsplit(x,""))),"var.mod.names"]
			y	<- sort(gsub("^Acetyl.*[(]","Acetyl (",y,perl=TRUE))
			return(paste(y,collapse=";"))})
	pepMat$sumIons			<- NA
	pepMat					<- pepMat[,!grepl("^.attrs|sumIons|pep_summed_mod_pos",colnames(pepMat))]
	pepMat					<- pepMat[pepMat$rank==1,]
	
	query			<- getNodeSet(doc,"//r:query",namespaces = c(r="http://www.matrixscience.com/xmlns/schema/mascot_search_results_2"))
	query.number	<- ldply(xpathApply(doc, "//r:query/r:q_peptide", xmlAttrs,namespaces = namespaces))
	query.number	<- table(as.numeric(query.number$query))
	query.number	<- as.numeric(names(query.number[query.number>1]))
	query.number	<- query.number[query.number%in% pepMat$query]
	queryMat		<- c()
	for(query.i in query.number){
		test			<- xmlToList(query[[query.i]])
		q_peptide.number	<- which(names(test) %in% "q_peptide")
		if(length(q_peptide.number)>0){
			tempMat	<- c()
			for(q_pep.i in q_peptide.number){
				test2	<- test[[q_pep.i]]
				test2	<- lapply(test2,function(x) {
					x[length(x)==0]	<- NA
					return(x)})
				tempMat	<- rbind(tempMat,unlist(test2)[!names(unlist(test2)) %in% c("pep_summed_mod_pos","pep_local_mod_pos","pep_var_mod_conf")])
			}
			tempMat		<- data.frame(tempMat,stringsAsFactors=FALSE)
			tempMat$pep_score	<- as.numeric(tempMat$pep_score)
			tempMat$pep_score_delta	<- tempMat$pep_score[1]-tempMat$pep_score
			tempMat		<- tempMat[tempMat$pep_score_delta < 10 & tempMat$pep_score>15 & tempMat$pep_seq == tempMat$pep_seq[1],]
			if(nrow(tempMat)>1){
				queryMat		<- rbind(queryMat,tempMat)
			}
		}
	}
	pepMat$pep_score_delta	<- NA
	if(length(queryMat)>0){
		queryMat					<- data.frame(apply(queryMat,2,function(x) gsub("NULL",NA,x)),stringsAsFactors=FALSE)
		queryMat					<- cbind(queryMat[,-which(colnames(queryMat)=="pep_scan_title")],ldply(apply(rbind(queryMat$pep_scan_title),1,function(x) strsplit(x," "))[[1]]))
		colnames(queryMat)[(length(colnames(queryMat))-4):length(colnames(queryMat))]	<- c("spectrumID","rem","acqNumber","rt","raw")
		queryMat					<- queryMat[-which(colnames(queryMat)=="rem")]
		queryMat$spectrumID			<- as.numeric(gsub(":","",queryMat$spectrumID))
		queryMat$rt					<- as.numeric(gsub("\\(rt=|\\)","",queryMat$rt))
		queryMat$acqNumber			<- as.numeric(queryMat$acqNumber)
		queryMat$query				<- as.numeric(queryMat$.attrs.query)
		queryMat$rank				<- as.numeric(queryMat$.attrs.rank)
		queryMat					<- queryMat[,!grepl("^.attrs|pep_summed_mod_pos",colnames(queryMat))]
		queryMat$ambiguous			<- 1
		chrToNum.cols				<- unlist(strsplit("pep_exp_mz,pep_exp_mr,pep_exp_z,pep_calc_mr,pep_delta,pep_start,pep_end,pep_miss,pep_score,pep_expect",","))
		pepMat[,chrToNum.cols]		<- apply(pepMat[,chrToNum.cols],2,as.numeric)
		chrToNum.cols				<- unlist(strsplit("pep_exp_mz,pep_exp_mr,pep_exp_z,pep_calc_mr,pep_delta,pep_miss,pep_score,pep_expect",","))
		queryMat[,chrToNum.cols]	<- apply(queryMat[,chrToNum.cols],2,as.numeric)
		mergeCols	<- c("query","spectrumID","acqNumber","rt","raw","pep_exp_mz","pep_exp_mr","pep_exp_z","pep_seq")
		mergedMat	<- merge(pepMat,queryMat,by=mergeCols,all.x=T)
		notmerged	<- colnames(queryMat) [!colnames(queryMat) %in% mergeCols]
		notmerged	<- notmerged[!notmerged %in% "ambiguous"]
		for(selcol in notmerged){
			tempcol.y	<- paste(c(selcol,".y"),collapse="")
			tempcol.x	<- paste(c(selcol,".x"),collapse="")
			selrow	<- which(is.na(mergedMat[,tempcol.y]))
			mergedMat[selrow,tempcol.y]	<- mergedMat[selrow,tempcol.x]	
		}
		mergedMat			<- mergedMat[,!grepl("[.]x$",colnames(mergedMat))]
		colnames(mergedMat)	<- gsub("[.]y$","",colnames(mergedMat))
		mergedMat$pep_var_mod.2	<- apply(cbind(gsub("\\.","",mergedMat$pep_var_mod_pos)),1,function(x)	{
				y	<- var.mod[as.numeric(unlist(strsplit(x,""))),"var.mod.names"]
				y	<- sort(gsub("^Acetyl.*[(]","Acetyl (",y,perl=TRUE))
				return(paste(y,collapse=";"))})
	}else{
		pepMat$ambiguous		<- NA
		mergedMat			<- pepMat
	}
	rm(doc)
	return(list(mergedMat,masses,var.mod,XMLFilename))
}
###
args 			<- commandArgs(trailingOnly = TRUE)
dataLocation	<- args[1]
outDir		<- args[2]
atomcols		<- unlist(strsplit(args[3],";"))

tempList		<- parse.Mascot.XML(dataLocation,atomcols)
save(tempList,file=paste(c(outDir,"mascotXML.RData"),collapse="\\"))
