###simple pasting path vectors
pastePath	<- function(inVec){
	return(paste(inVec,collapse="\\"))
}

args 			<- commandArgs(trailingOnly = TRUE)

#args			<- passargs
atomcols		<- unlist(strsplit(args[1],";"))
reagent.purity	<- as.numeric(args[2])
IDfolder		<- args[3]
quanfolder		<- args[3]
resultsFolder	<- args[4]
project			<- args[5]
errortolPPM		<- as.numeric(args[6])
ScriptFolder	<- args[7]
cutoff.intens	<- as.numeric(args[8]) #0
cutoff.corr		<- as.numeric(args[9]) #0.99
cutoff.maxdiff	<- as.numeric(args[10])  #20
cores			<- as.numeric(args[11])  #4
logFolder		<- args[12]
mainName		<- args[13]

logFolder2		<- pastePath(c(logFolder,paste(c(mainName,project),collapse="_")))
logfilename		<- paste(c(logFolder2,".log.txt"),collapse="")
error.logfilename	<- paste(c(logFolder2,".errorlog.txt"),collapse="")

logfilename	<- paste(c(logFolder2,".log.txt"),collapse="")
writeLines(c(""), logfilename)
error.logfilename	<- paste(c(logFolder2,".errorlog.txt"),collapse="")
error.logfile.con	<- file(error.logfilename,open="wt")
sink(logfilename,append=TRUE,type=c("output"))
sink(error.logfile.con,append=TRUE,type=c("message"))
x	<- file(paste(c(logFolder,"\\starting_",mainName,"_",project,".txt"),collapse=""),"w")
close(x)

cat(paste0(Sys.time(),paste(c("   starting project: ",project),collapse=""),"\n"))
flush.console()


###
library("mzR")
library("plyr")
library("XML")
library("OrgMassSpecR")
library("Rdisop")
library("gplots")
library("doParallel")
library("RColorBrewer")


###make output directory
make.outDir	<- function(dataLocation,subDir){
	subDir		<- "analysis"
	timestamp		<- format(Sys.time(), format="%Y%m%d.%H%M")
	outDir		<- paste(c(dataLocation, paste(c(timestamp,subDir),collapse="_")),collapse="\\")
	if (!file.exists(outDir)){
	    dir.create(outDir)
	}
	return(outDir)
}
##calulate theoretical fragment ion spectrum and correct for N- or C-terminal modifications

calc.fragmentions	<- function(pepMat,dfout2,reagent.purity,atomcols,masses,var.mod,pl,outDir,mzDataFilename,XMLFilename,errortolPPM,elements,cores){
	library("OrgMassSpecR")
	library("plyr")
	library("Rdisop")
	AcModKeys	<- var.mod$var.mod.id[grepl("^Acetyl \\(K)",var.mod$var.mod.names)]
	NtermAcModKeys	<- var.mod$var.mod.id[grepl(c("Acetyl \\(Protein N-term)"),var.mod$var.mod.names)]
	outList	<- vector("list",nrow(pepMat))
	cl <- makeCluster(cores)
	registerDoParallel(cl)
	outList	<- c()
	outList	<- foreach(line.i = 1:nrow(pepMat)) %dopar% {
	#for(line.i in 1:nrow(pepMat)){
		library("mzR")
		library("OrgMassSpecR")
		library("Rdisop")
		library("plyr")
		pepLine		<- pepMat[line.i,]
		pep.len		<- nchar(pepLine$pep_seq)
		pep.mods		<- pepLine$pep_var_mod_pos
			pep.mods			<- unlist(strsplit(pep.mods,"\\."))
			pep.mods.N			<- pep.mods[1]
			pep.mods.C			<- pep.mods[3]
			pep.mods			<- unlist(strsplit(pep.mods[2],""))
			pep.mods.mat	<- data.frame(c(pep.mods.N,pep.mods,pep.mods.C),mod.name=NA,mod.pep.pos=c(0:(length(pep.mods)+1)),stringsAsFactors=FALSE)
			for(i in 1:length(pep.mods)){
				if(pep.mods[i]!="0"){
					pep.mods.mat[pep.mods.mat$mod.pep.pos==i,"mod.name"]	<- var.mod[var.mod$var.mod.id == pep.mods[i],"var.mod.names"]
				}
			}
			if(pep.mods.N!="0"){
				pep.mods.mat[pep.mods.mat$mod.pep.pos==0,"mod.name"]	<- var.mod[var.mod$var.mod.id == pep.mods.N,"var.mod.names"]
			}
			if(pep.mods.C!="0"){
				pep.mods.mat[pep.mods.mat$mod.pep.pos==(length(pep.mods)+1),"mod.name"]	<- var.mod[var.mod$var.mod.id == pep.mods.C,"var.mod.names"]
			}
		calc.fragments		<- FragmentPeptide(pepLine$pep_seq)
		tempi				<- sapply(calc.fragments, is.factor)
		calc.fragments[,tempi] 	<- lapply(calc.fragments[,tempi], as.character)
		ion.type	<- ldply(apply(rbind(calc.fragments$ms2type),1,function(x) (strsplit(x,"\\[|\\]") ))[[1]])[,2:3]
		colnames(ion.type)	<- c("ion","z")
		ion.type$z	<- as.numeric(gsub("\\+","",ion.type$z))
		ion.type$by	<- gsub("[0-9]","",ion.type$ion)
		ion.type$n	<- as.numeric(gsub("[a-z]","",ion.type$ion))
		calc.fragments	<- cbind(calc.fragments,ion.type)
		calc.fragments	<- calc.fragments[calc.fragments$z==1,]
		calc.fragments	<- cbind(calc.fragments,mods=NA,mods.pos=NA,mods.name=NA,Ac.mods.pos=NA, corr.mass=NA,C=NA,H=NA,N=NA,O=NA,S=NA,Cx=NA,AcK=NA,AcK2a=NA,AcK2b=NA,AcK3a=NA,AcK3b=NA,AcK3c=NA)
		calc.fragments[,atomcols]	<- cbind(t(apply(calc.fragments,1,function(x) unlist(ConvertPeptide(as.character(x["ms2seq"]),IAA = TRUE)))),Cx=0)
		for(frag.i in 1:nrow(calc.fragments)){
			if(calc.fragments[frag.i,"by"]=="b"){
				calc.fragments[frag.i,"H"]	<- calc.fragments[frag.i,"H"]-1
				calc.fragments[frag.i,"O"]	<- calc.fragments[frag.i,"O"]-1
				frag.pos	<- 1:as.numeric(calc.fragments[frag.i,"n"])
				frag.mod	<- c(pep.mods.N,pep.mods[frag.pos])
				frag.mod	<- frag.mod[!frag.mod %in% "0"]
				if(sum(c(frag.mod %in% AcModKeys,frag.mod %in% NtermAcModKeys))>0){
					temp.atomcols	<- c()
					temp.mods		<- c()
					for(x in frag.mod){
						temp.atomcols	<- rbind(temp.atomcols,var.mod[var.mod$var.mod.id == x,atomcols])
						temp.mods		<- c(temp.mods,var.mod[var.mod$var.mod.id == x,"var.mod.names"])
					}
					calc.fragments[frag.i,atomcols]	<- calc.fragments[frag.i,atomcols] + apply(temp.atomcols,2,sum)
					calc.fragments[frag.i,"mods"]		<- paste(temp.mods,collapse=";")
					calc.fragments[frag.i,"Ac.mods.pos"]<- paste(c(if(pep.mods.N %in% NtermAcModKeys){0},which(pep.mods[frag.pos] %in% AcModKeys)),collapse=";")
					temp.mods	<- pep.mods.mat[!is.na(pep.mods.mat$mod.name)&pep.mods.mat$mod.pep.pos %in% c(0,1:pep.len),]

					calc.fragments[frag.i,"mods.pos"]	<- paste(temp.mods$mod.pep.pos,collapse=";")
					calc.fragments[frag.i,"mods.name"]	<- paste(temp.mods$mod.name,collapse=";")
				}
			}
			if(calc.fragments[frag.i,"by"]=="y"){
				calc.fragments[frag.i,"H"]	<- calc.fragments[frag.i,"H"]+1
				calc.fragments[frag.i,"O"]	<- calc.fragments[frag.i,"O"]
				frag.pos	<- (pep.len+1-as.numeric(calc.fragments[frag.i,"n"])):pep.len
				frag.mod	<- c(pep.mods[frag.pos],pep.mods.C)
				frag.mod	<- frag.mod[!frag.mod %in% "0"]
				if(sum(frag.mod %in% AcModKeys)>0){
					temp.atomcols	<- c()
					temp.mods		<- c()
					for(x in frag.mod){
						temp.atomcols	<- rbind(temp.atomcols,var.mod[var.mod$var.mod.id == x,atomcols])
						temp.mods		<- c(temp.mods,var.mod[var.mod$var.mod.id == x,"var.mod.names"])
					}
					calc.fragments[frag.i,atomcols]	<- calc.fragments[frag.i,atomcols] + apply(temp.atomcols,2,sum)
					calc.fragments[frag.i,"mods"]		<- paste(temp.mods,collapse=";")
					calc.fragments[frag.i,"Ac.mods.pos"]<- paste(frag.pos[which(pep.mods[frag.pos] %in% AcModKeys)],collapse=";")
					temp.mods	<- pep.mods.mat[!is.na(pep.mods.mat$mod.name)&pep.mods.mat$mod.pep.pos %in% c(1:pep.len,pep.len+1),]

					calc.fragments[frag.i,"mods.pos"]	<- paste(temp.mods$mod.pep.pos,collapse=";")
					calc.fragments[frag.i,"mods.name"]	<- paste(temp.mods$mod.name,collapse=";")
				}
			}
			calc.fragments[frag.i,"corr.mass"]	<- MonoisotopicMass(formula = calc.fragments[frag.i,atomcols[!atomcols %in% "Cx"]])-masses["Electron"]+calc.fragments[frag.i,"Cx"]*12
		}
		if(sum(!is.na(calc.fragments$Ac.mods.pos))>0){
			calc.fragments$AcetylMods			<- sapply(strsplit(calc.fragments$mods," "),function(x) sum(grepl("Acetyl",x)))
			calc.fragments					<- calc.fragments[calc.fragments$AcetylMods>0,]
		}else{
			calc.fragments	<- c()
		}
		if(length(calc.fragments)>0){
		 	###match to experimental peak list
			peaklist.index		<- which(dfout2$acqNumber==pepLine$acqNumber)
			if(length(peaklist.index)>1){
				peaklist			<- unique(ldply(pl[peaklist.index]))
			}else{
				peaklist			<- pl[[peaklist.index]]
			}
			pepLine$sumIons		<- sum(peaklist[,2])

			diffmat			<- abs(outer(peaklist[,1],calc.fragments$corr.mass,"-"))
			errorvec			<- calc.fragments$corr.mass*errortolPPM
			fragment.match		<- sapply(1:ncol(diffmat),function(i) list(which(diffmat[,i]<errorvec[i])))
			if(length(fragment.match)==0){
				fragment.match	<- rep(list(0),ncol(diffmat))
			}
			diffmat.minus1		<- abs(outer(peaklist[,1],calc.fragments$corr.mass-1.00335483,"-"))
			errorvec			<- (calc.fragments$corr.mass-1.00335483)*errortolPPM
			fragment.match.minus1	<- sapply(1:ncol(diffmat.minus1),function(i) list(which(diffmat.minus1[,i]<errorvec[i])))
			if(length(fragment.match.minus1)==0){
				fragment.match.minus1	<- rep(list(0),ncol(diffmat.minus1))
			}
			diffmat.plus1		<- abs(outer(peaklist[,1],calc.fragments$corr.mass+1.00335483,"-"))
			errorvec			<- (calc.fragments$corr.mass+1.00335483)*errortolPPM
			fragment.match.plus1	<- sapply(1:ncol(diffmat.plus1),function(i) list(which(diffmat.plus1[,i]<errorvec[i])))
			if(length(fragment.match.plus1)==0){
				fragment.match.plus1	<- rep(list(0),ncol(diffmat.plus1))
			}
			diffmat.plus2		<- abs(outer(peaklist[,1],calc.fragments$corr.mass+2*1.00335483,"-"))
			errorvec			<- (calc.fragments$corr.mass+2*1.00335483)*errortolPPM
			fragment.match.plus2	<- sapply(1:ncol(diffmat.plus2),function(i) list(which(diffmat.plus2[,i]<errorvec[i])))
			if(length(fragment.match.plus2)==0){
				fragment.match.plus2	<- rep(list(0),ncol(diffmat.plus2))
			}
			diffmat.plus3		<- abs(outer(peaklist[,1],calc.fragments$corr.mass+3*1.00335483,"-"))
			errorvec			<- (calc.fragments$corr.mass+3*1.00335483)*errortolPPM
			fragment.match.plus3	<- sapply(1:ncol(diffmat.plus3),function(i) list(which(diffmat.plus3[,i]<errorvec[i])))
			if(length(fragment.match.plus3)==0){
				fragment.match.plus3	<- rep(list(0),ncol(diffmat.plus3))
			}
	
			diffmat.plus4		<- abs(outer(peaklist[,1],calc.fragments$corr.mass+4*1.00335483,"-"))
			errorvec			<- (calc.fragments$corr.mass+4*1.00335483)*errortolPPM
			fragment.match.plus4	<- sapply(1:ncol(diffmat.plus4),function(i) list(which(diffmat.plus4[,i]<errorvec[i])))
			if(length(fragment.match.plus4)==0){
				fragment.match.plus4	<- rep(list(0),ncol(diffmat.plus4))
			}
			diffmat.plus5		<- abs(outer(peaklist[,1],calc.fragments$corr.mass+5*1.00335483,"-"))
			errorvec			<- (calc.fragments$corr.mass+5*1.00335483)*errortolPPM
			fragment.match.plus5	<- sapply(1:ncol(diffmat.plus5),function(i) list(which(diffmat.plus5[,i]<errorvec[i])))
			if(length(fragment.match.plus5)==0){
				fragment.match.plus5	<- rep(list(0),ncol(diffmat.plus5))
			}
			diffmat.plus6		<- abs(outer(peaklist[,1],calc.fragments$corr.mass+6*1.00335483,"-"))
			errorvec			<- (calc.fragments$corr.mass+6*1.00335483)*errortolPPM
			fragment.match.plus6	<- sapply(1:ncol(diffmat.plus6),function(i) list(which(diffmat.plus6[,i]<errorvec[i])))
			if(length(fragment.match.plus6)==0){
				fragment.match.plus6	<- rep(list(0),ncol(diffmat.plus6))
			}
			intensities4		<- data.frame(	minus1=sapply(fragment.match.minus1,function(x) sum(peaklist[x,2])),
										mono=sapply(fragment.match,function(x) sum(peaklist[x,2])),
										plus1=sapply(fragment.match.plus1,function(x) sum(peaklist[x,2])),
										plus2=sapply(fragment.match.plus2,function(x) sum(peaklist[x,2])),
										plus3=sapply(fragment.match.plus3,function(x) sum(peaklist[x,2])),
										plus4=sapply(fragment.match.plus4,function(x) sum(peaklist[x,2])),
										plus5=sapply(fragment.match.plus5,function(x) sum(peaklist[x,2])),
										plus6=sapply(fragment.match.plus6,function(x) sum(peaklist[x,2])))
	
			errorMat			<- data.frame(	delta.mono.Da=sapply(1:length(fragment.match),function(x) median(abs(diffmat[fragment.match[[x]],x]))),
										delta.plus1.Da=sapply(1:length(fragment.match.plus1),function(x) median(abs(diffmat.plus1[fragment.match.plus1[[x]],x]))),
										delta.plus2.Da=sapply(1:length(fragment.match.plus2),function(x) median(abs(diffmat.plus2[fragment.match.plus2[[x]],x]))),
										delta.plus3.Da=sapply(1:length(fragment.match.plus3),function(x) median(abs(diffmat.plus3[fragment.match.plus3[[x]],x]))),
										delta.plus4.Da=sapply(1:length(fragment.match.plus4),function(x) median(abs(diffmat.plus4[fragment.match.plus4[[x]],x]))),
										delta.plus5.Da=sapply(1:length(fragment.match.plus5),function(x) median(abs(diffmat.plus5[fragment.match.plus5[[x]],x]))),
										delta.plus6.Da=sapply(1:length(fragment.match.plus6),function(x) median(abs(diffmat.plus6[fragment.match.plus6[[x]],x]))))
			errorPPM			<- errorMat/calc.fragments$corr.mass*1E6
			colnames(errorPPM)	<- gsub("Da","ppm",colnames(errorPPM))
			calc.fragments		<- cbind(calc.fragments,intensities4,errorMat,errorPPM)
			calc.fragments$intens	<- log10(calc.fragments$mono+calc.fragments$plus1+calc.fragments$plus2)
			calc.fragments$intens2	<- log10(calc.fragments$mono+calc.fragments$plus1+calc.fragments$plus2+calc.fragments$plus3+calc.fragments$plus4)
			calc.fragments$intens3	<- log10(calc.fragments$mono+calc.fragments$plus1+calc.fragments$plus2+calc.fragments$plus3+calc.fragments$plus4+calc.fragments$plus5+calc.fragments$plus6)
			calc.fragments		<- calc.fragments[calc.fragments$AcetylMods %in% c(1,2,3) & calc.fragments$intens>0,]
			if(nrow(calc.fragments)>0){
				selectFrag	<- calc.fragments$AcetylMods==1 & calc.fragments$by=="b" & grepl(c("Acetyl \\(K)"),calc.fragments$mods)
				if(sum(selectFrag)>0){
					calc.fragments[selectFrag,"AcK"] <- as.numeric(pepLine$pep_start) + as.numeric(calc.fragments$Ac.mods.pos[selectFrag]) -1
				}
	
				selectFrag	<- calc.fragments$AcetylMods==1 & calc.fragments$by=="b" & grepl(c("Acetyl \\(Protein N-term)"),calc.fragments$mods)
				if(sum(selectFrag)>0){
					calc.fragments[selectFrag,"AcK"] <- paste(c("N",as.numeric(pepLine$pep_start)),collapse="")
				}
	
				selectFrag	<- calc.fragments$AcetylMods==1 & calc.fragments$by=="y"
				if(sum(selectFrag)>0){
					calc.fragments[selectFrag,"AcK"] <- as.numeric(pepLine$pep_start) + as.numeric(calc.fragments$Ac.mods.pos[selectFrag]) -1
				}
	
				selectFrag	<- calc.fragments$AcetylMods==2 & calc.fragments$by=="b" & grepl(c("Acetyl \\(Protein N-term)"),calc.fragments$mods)
				if(sum(selectFrag)>0){
					temp <-  as.numeric(pepLine$pep_start) + ldply(lapply(strsplit(calc.fragments$Ac.mods.pos[selectFrag],";"), as.numeric))  -1
					temp[,1]	<- gsub("^","N",temp[,1]+1)
					calc.fragments[selectFrag,c("AcK2a","AcK2b")]	<- temp 
				}
	
				selectFrag	<- calc.fragments$AcetylMods==2 & calc.fragments$by=="b" & !grepl(c("Acetyl \\(Protein N-term)"),calc.fragments$mods)
				if(sum(selectFrag)>0){
					calc.fragments[selectFrag,c("AcK2a","AcK2b")] <- as.numeric(pepLine$pep_start) + ldply(lapply(strsplit(calc.fragments$Ac.mods.pos[selectFrag],";"), as.numeric))-1
				}
	
				selectFrag	<- calc.fragments$AcetylMods==2 & calc.fragments$by=="y"
				if(sum(selectFrag)>0){
					calc.fragments[selectFrag,c("AcK2a","AcK2b")] <- as.numeric(pepLine$pep_start) + ldply(lapply(strsplit(calc.fragments$Ac.mods.pos[selectFrag],";"), as.numeric))-1
				}

				selectFrag	<- calc.fragments$AcetylMods==3 & calc.fragments$by=="b" & grepl(c("Acetyl \\(Protein N-term)"),calc.fragments$mods)
				if(sum(selectFrag)>0){
					temp <-  as.numeric(pepLine$pep_start) + ldply(lapply(strsplit(calc.fragments$Ac.mods.pos[selectFrag],";"), as.numeric))  -1
					temp[,1]	<- gsub("^","N",temp[,1]+1)
					calc.fragments[selectFrag,c("AcK3a","AcK3b","AcK3c")]	<- temp 
				}
	
				selectFrag	<- calc.fragments$AcetylMods==3 & calc.fragments$by=="b" & !grepl(c("Acetyl \\(Protein N-term)"),calc.fragments$mods)
				if(sum(selectFrag)>0){
					calc.fragments[selectFrag,c("AcK3a","AcK3b","AcK3c")] <- as.numeric(pepLine$pep_start) + ldply(lapply(strsplit(calc.fragments$Ac.mods.pos[selectFrag],";"), as.numeric))-1
				}
	
				selectFrag	<- calc.fragments$AcetylMods==3 & calc.fragments$by=="y"
				if(sum(selectFrag)>0){
					calc.fragments[selectFrag,c("AcK3a","AcK3b","AcK3c")] <- as.numeric(pepLine$pep_start) + ldply(lapply(strsplit(calc.fragments$Ac.mods.pos[selectFrag],";"), as.numeric))-1
				}

				occup3.results	<- vector("list",nrow(calc.fragments))
				for(i in 1:length(occup3.results)){
					x	<- calc.fragments[i,]
					if(!is.na(x["AcK"])){
						atoms			<- as.numeric(x[atomcols])
						names(atoms)	<- names(x[atomcols])
						atoms.light		<- paste0(names(atoms)[atoms>0],atoms[atoms>0],collapse="")
						atoms["C"]		<- atoms["C"]-1
						atoms["Cx"]		<- atoms["Cx"]+1
						atoms.heavy		<- paste0(names(atoms)[atoms>0],atoms[atoms>0],collapse="")	
	
						isodist.light	<- getMolecule(atoms.light,elements=elements,maxisotopes=10)$isotopes[[1]]
						isodist.light[2,]	<- isodist.light[2,]/sum(isodist.light[2,])*100
						isodist.heavy	<- getMolecule(atoms.heavy,elements=elements,maxisotopes=10)$isotopes[[1]]
						isodist.heavy[2,]	<- isodist.heavy[2,]/sum(isodist.heavy[2,])*100
	
						peaks.cols		<- c("minus1","mono","plus1","plus2","plus3","plus4","plus5","plus6")
						peaks.sel		<- which(peaks.cols %in% "mono"):(which(peaks.cols %in% "plus2")+atoms["Cx"]-1)		
						peaks.sel		<- peaks.sel[peaks.sel<=length(peaks.cols)]
						peaks3		<- as.numeric(c(x[peaks.cols[peaks.sel]]))
						peaks3		<- peaks3/max(peaks3)*100
	
						peaks.sel		<- which(peaks.cols %in% "minus1"):(which(peaks.cols %in% "plus2")+atoms["Cx"]-1)		
						peaks.sel		<- peaks.sel[peaks.sel<=length(peaks.cols)]
						peaks4		<- as.numeric(c(x[peaks.cols[peaks.sel]]))
						peaks4		<- peaks4/max(peaks4)*100

						theoratios			<- seq(0,100,by=.5)
						isodist.sum			<- t(outer(isodist.light[2,1:length(peaks3)],theoratios,"*") + outer(isodist.heavy[2,1:length(peaks3)],(100-theoratios),"*"))
						theopeaks			<- isodist.sum/apply(isodist.sum,1,max)*100
						diffpattern			<- apply(theopeaks,1,function(x) sum((peaks3-x)^2))


						which.min			<- which(diffpattern==min(diffpattern))
						occupancy			<- median(theoratios[which.min])
						theoratios			<- seq(occupancy-2.5,occupancy+2.5,by=0.01)
						theoratios			<- theoratios[theoratios>=0 & theoratios<=100]

						isodist.sum			<- t(outer(isodist.light[2,1:length(peaks3)],theoratios,"*") + outer(isodist.heavy[2,1:length(peaks3)],(100-theoratios),"*"))
						theopeaks			<- isodist.sum/apply(isodist.sum,1,max)*100
						diffpattern			<- apply(theopeaks,1,function(x) sum((peaks3-x)^2))
			
						which.min			<- which(diffpattern==min(diffpattern))
						occupancy			<- median(theoratios[which.min])
	
						theopeaks4			<- c(0,isodist.sum[which.min,]/max(isodist.sum[which.min,])*100)
						peaks4.ori			<- peaks4
						theopeaks4.ori		<- theopeaks4
						corr.coff			<- cor(peaks4,theopeaks4)
						sumdiff			<- diffpattern[which.min]
						maxdiff			<- max(abs(theopeaks4-peaks4))
						occup3.results[[i]]	<- (data.frame(occupancy,corr.coff,maxdiff,sumdiff,paste(peaks4.ori,collapse=";"),paste(theopeaks4.ori,collapse=";"),stringsAsFactors=FALSE))
					}else{
						occup3.results[[i]]	<- (data.frame(NA,NA,NA,NA,NA,NA,stringsAsFactors=FALSE))
					}}
				occup3.results				<- ldply(occup3.results)
				calc.fragments$occup3			<- occup3.results[,1]
				calc.fragments$occup3.corr.coff	<- occup3.results[,2]
				calc.fragments$occup3.maxdiff		<- occup3.results[,3]
				calc.fragments$occup3.sumdiff		<- occup3.results[,4]
				calc.fragments$peaks4			<- occup3.results[,5]
				calc.fragments$theopeaks4		<- occup3.results[,6]
				outMatparallel	<- merge(pepLine,calc.fragments)
				outMatparallel	<- outMatparallel[!outMatparallel$ion %in% "b1" | (outMatparallel$ion %in% "b1" & outMatparallel$pep_start<2),]
				outMatparallel
			}
			
		}
	}
	stopCluster(cl)
	outMat		<- ldply(outList)
	outMat		<- outMat[!is.na(outMat$occup3)|outMat$AcetylMods==2|outMat$AcetylMods==3,]
	mod.temp	<- apply(outMat,1,function(x){
	temp1	<- cbind(unlist(strsplit(x["mods.pos"],";")),unlist(strsplit(x["mods.name"],";")))
		temp2	<- apply(temp1,1,function(y){
			if(y[1]=="0"){
				pos.temp	<- gsub("^","N",as.numeric(x["pep_start"]))
			}else{
				pos.temp	<- sum(as.numeric(c(y[1],x["pep_start"])))-1
			}
			return(paste(c(gsub(" \\(.*)| .*","",y[2]),"(",pos.temp,")"),collapse=""))
		})
		temp2	<- temp2[!grepl("pyro-Glu|Oxidat",temp2)]
		return(paste(temp2,collapse=";"))
	})
	outMat$other.mods	<- mod.temp
	outMat		<- outMat[order(-as.numeric(outMat$prot_score),outMat$accession,outMat$AcK,outMat$other.mods,outMat$intens),]
	outFilename		<- paste(c(outDir,"\\",mzDataFilename,"_",XMLFilename,"_Acetylation_quant.txt"),collapse="")
	write.table(outMat, file = outFilename, quote = FALSE, sep = "\t",na = "NA", dec = ".", row.names = FALSE)
	return(outFilename)
}
#########combine result tables
combineResults	<- function(ProcData,outDir){
	outMat	<- c()
	for(inMatname in ProcData){
		inMat		<- read.delim(inMatname, sep="\t",header=TRUE,stringsAsFactors=FALSE)
		inMat$other.mods[is.na(inMat$other.mods)]	<- ""
		outMat	<- rbind(outMat,inMat)
	}
	outMat	<- merge(outMat,data.frame(table(outMat$accession)),by.x="accession",by.y="Var1")
	outMat	<- outMat[order(-outMat$Freq,outMat$accession,outMat$AcK,outMat$other.mods,outMat$intens),]
	outMat$id	<- 1:nrow(outMat)
	return(outMat)
}
###############
###condense results and create graphics
create.condens.1	<- function(inMat,outDir){
	condensMat	<- c()
	selcols	<- c("accession","AcK","other.mods","pep_var_mod.4")
	inMat$selcol2	<- do.call(paste,c(inMat[,selcols],"NA",sep=";;;"))
	tempsel	<- unique(inMat$selcol2)
	for(i in 1:length(tempsel)){
		flush.console()
		selected	<- tempsel[i]
		temp		<- inMat[which(inMat$selcol2 %in% selected),]
		if(nrow(temp)>0){
			filt			<- temp$occup3 >= boxplot(temp$occup3,plot=F)$stats[1] & temp$occup3 <= boxplot(temp$occup3,plot=F)$stats[5]
			temp.filt		<- temp[filt,]
			selected		<- ldply(strsplit(selected,";;;"),stringsAsFactors=FALSE)[,1:length(selcols)]
			colnames(selected)	<- selcols
			filt.result	<- data.frame(accession=selected$accession,pos=selected$AcK,pep_var_mod.4=selected$pep_var_mod.4,n.filt=nrow(temp.filt),
				med.filt=median(temp.filt$occup3,na.rm=T),mad.filt=mad(temp.filt$occup3,na.rm=T),
				mean.filt=mean(temp.filt$occup3,na.rm=T),sd.filt=sd(temp.filt$occup3,na.rm=T),
				wmean.filt=weighted.mean(temp.filt$occup3,temp.filt$intens,na.rm=T),stringsAsFactors=FALSE)
			condensMat	<- rbind(condensMat,cbind(filt.result))
		}
	}
	condensMat	<- condensMat[!is.na(condensMat$med.filt)&condensMat$n.filt>1,]
	return(condensMat)
}
###################
create.condens.2	<- function(outMat,outDir,cutoff.intens,timestamp,project,rawNames){
	outMat$filt	<- 1
	outMat$id	<- 1:nrow(outMat)
	condensMat	<- c()
	for(accession in unique(outMat$accession)){
		filtMat	<- c()
		AcKVec	<- unique(outMat[outMat$accession %in% accession,"AcK"])
		AcKorder	<- order(as.numeric(gsub("N","0.",AcKVec)))
		for(AcK in AcKVec[AcKorder]){
			for(other.mods in unique(outMat[outMat$accession %in% accession,"pep_var_mod.4"])){
				temp	<- rbind(outMat[outMat$accession %in% accession & 
						outMat$AcK == AcK & 
						outMat$pep_var_mod.4 %in% other.mods,])
				temp	<- temp[!is.na(temp$occup3),]
				temp0	<- temp
				temp.filt	<- rbind()
				if(nrow(temp)>0){
					temp	<- temp0[temp0$filt0>1,] ##only correlation survivors
					if(nrow(temp)>0){
						outMat$filt[outMat$id %in% temp$id]	<- 2 #survived intensity and correlation filter
						bxpstat	<- boxplot(temp$occup3,plot=FALSE,outline=FALSE,at=ceiling(max(temp$intens,na.rm=T)+3)-.1)$stats
						filt		<- temp$occup3 >= bxpstat[1] & temp$occup3 <= bxpstat[5]
						temp.filt	<- temp[filt,]
						if(nrow(temp.filt)>0){
							outMat$filt[outMat$id %in% temp.filt$id]	<- 3 #survived outlier test, red
							filt.result	<- data.frame(accession=accession,pos=AcK,pep_var_mod.4=other.mods,n.filt=nrow(temp.filt),
								med.filt=median(temp.filt$occup3,na.rm=T),mad.filt=mad(temp.filt$occup3,na.rm=T),
								mean.filt=mean(temp.filt$occup3,na.rm=T),sd.filt=sd(temp.filt$occup3,na.rm=T),
								wmean.filt=weighted.mean(temp.filt$occup3,temp.filt$intens,na.rm=T),stringsAsFactors=F)
							condensMat	<- rbind(condensMat,filt.result)
							filtMat	<- rbind(filtMat,temp.filt)
						}
						temp<-temp.filt
					}
				}
			}#other.mods
		}#AcK
	}
	return(list(condensMat,outMat))
}

create.condens.3	<- function(outMat,outDir,cutoff.intens,timestamp,project,rawNames){
	outMat$filt2	<- 1
	outMat$id	<- 1:nrow(outMat)
	condensMat	<- c()
	outDirSite	<- paste(c(outDir,"sites_details"),collapse="\\")
	if (!file.exists(outDirSite)){
	    dir.create(outDirSite)
	}	
	outDirCombi	<- paste(c(outDir,"sites_combined"),collapse="\\")
	if (!file.exists(outDirCombi)){
	    dir.create(outDirCombi)
	}

	for(accession in unique(outMat$accession)){
		GraphOutName	<- paste(c(outDirSite,"\\details_",gsub("^.*[|:]","",accession),".pdf"), collapse="")
		unlist(GraphOutName)
		pdf(GraphOutName)
		filtMat	<- c()
		AcKVec	<- unique(outMat[outMat$accession %in% accession,"AcK"])
		AcKorder	<- order(as.numeric(gsub("N","0.",AcKVec)))
		for(AcK in AcKVec[AcKorder]){
			for(other.mods in unique(outMat[outMat$accession %in% accession,"pep_var_mod.4"])){
				temp	<- rbind(outMat[outMat$accession %in% accession & 
						outMat$AcK == AcK & 
						outMat$pep_var_mod.4 %in% other.mods,])
				temp	<- temp[!is.na(temp$occup3),]
				temp0	<- temp
				temp.filt	<- rbind()
				if(nrow(temp)>0){
					par(mar=c(4,4,4,4))
					plot(temp0$intens,temp0$occup3,xlim=c(cutoff.intens,ceiling(max(temp$intens2,na.rm=T)+3)),ylim=c(0,100),xlab="",,ylab="",type="n")
					abline(h=50)
					rawfamily	<- as.factor(temp0$raw)
					col.rainbow <- brewer.pal(9,"Set1")
					palette(col.rainbow)
					#plot only intensity survivors
					points(temp0$intens[temp0$AcetylMods==1&temp0$filt<3],temp0$occup3[temp0$AcetylMods==1&temp0$filt<3],pch=1,cex=.8,col=rawfamily[temp0$AcetylMods==1&temp0$filt<3])
					points(temp0$intens2[temp0$AcetylMods==2&temp0$filt<3],temp0$occup3[temp0$AcetylMods==2&temp0$filt<3],pch=2,cex=.8,col=rawfamily[temp0$AcetylMods==2&temp0$filt<3])
					points(temp0$intens2[temp0$AcetylMods==3&temp0$filt0==1],temp0$occup3[temp0$AcetylMods==3&temp0$filt0==1],pch=5,cex=.8,col=rawfamily[temp0$AcetylMods==3&temp0$filt0==1])

					title(xlab="log10 intensity",line=2.5)
					title(ylab="occupancy (%)",line=2.5)
					temp	<- temp0[(temp0$AcetylMods==3&temp0$filt0>1)|(temp0$AcetylMods<3&temp0$filt==3),] ##only correlation survivors and boxplot survivors from 2K
					if(nrow(temp)>0){
						outMat$filt2[outMat$id %in% temp$id]	<- 2 #survived intensity and correlation filter
						bxpstat	<- boxplot(temp$occup3,add=T,outline=FALSE,at=ceiling(max(temp$intens,na.rm=T)+3)-.1)$stats
						filt		<- temp$occup3 >= bxpstat[1] & temp$occup3 <= bxpstat[5]
						temp.filt	<- temp[filt,]
						points(temp$intens[temp$AcetylMods==1&!filt],temp$occup3[temp$AcetylMods==1&!filt],pch=1,cex=.8,col=rawfamily[temp0$AcetylMods==1&temp0$filt0==1])
						points(temp$intens[temp$AcetylMods==1&!filt],temp$occup3[temp$AcetylMods==1&!filt],pch=3,cex=.8,col=rawfamily[temp0$AcetylMods==1&temp0$filt0==1])

						points(temp$intens2[temp$AcetylMods==2&!filt],temp$occup3[temp$AcetylMods==2&!filt],pch=2,cex=.8,col=rawfamily[temp0$AcetylMods==2&temp0$filt0==1])
						points(temp$intens2[temp$AcetylMods==2&!filt],temp$occup3[temp$AcetylMods==2&!filt],pch=3,cex=.8,col=rawfamily[temp0$AcetylMods==2&temp0$filt0==1])

						points(temp$intens2[temp$AcetylMods==3&!filt],temp$occup3[temp$AcetylMods==3&!filt],pch=5,cex=.8,col=rawfamily[temp0$AcetylMods==3&temp0$filt0==1])
						points(temp$intens2[temp$AcetylMods==3&!filt],temp$occup3[temp$AcetylMods==3&!filt],pch=3,cex=.8,col=rawfamily[temp0$AcetylMods==3&temp0$filt0==1])


						plotmed	<- bxpstat[3]
						abline(h=plotmed,col="blue",lwd=2)
						axis(side=4,at=plotmed,sprintf("%.2f", plotmed),las=2)
						if(nrow(temp.filt)>0){
							outMat$filt2[outMat$id %in% temp.filt$id]	<- 3 #survived outlier test, red
							points(temp.filt$intens[temp.filt$AcetylMods==1],temp.filt$occup3[temp.filt$AcetylMods==1],cex=.8,pch=21,bg=rawfamily,col=rawfamily)
							points(temp.filt$intens2[temp.filt$AcetylMods==2],temp.filt$occup3[temp.filt$AcetylMods==2],cex=.8,pch=24,bg=rawfamily,col=rawfamily)
							points(temp.filt$intens2[temp.filt$AcetylMods==3],temp.filt$occup3[temp.filt$AcetylMods==3],cex=.8,pch=23,bg=rawfamily,col=rawfamily)

							filt.result	<- data.frame(accession=accession,pos=AcK,mod=other.mods,n.filt=nrow(temp.filt),
								med.filt=median(temp.filt$occup3,na.rm=T),mad.filt=mad(temp.filt$occup3,na.rm=T),
								mean.filt=mean(temp.filt$occup3,na.rm=T),sd.filt=sd(temp.filt$occup3,na.rm=T),
								wmean.filt=weighted.mean(temp.filt$occup3,temp.filt$intens,na.rm=T))
							condensMat	<- rbind(condensMat,filt.result)
							filtMat	<- rbind(filtMat,temp.filt)
						}
						legend(ceiling(max(temp$intens,na.rm=T))+2.6, -3, gsub("^.*[\\]|.raw","",levels(rawfamily)), pch = 15, col = col.rainbow[1:length(levels(rawfamily))], cex = .8,bg=NULL,xjust=1,yjust=0)
						caption	<- paste(c("Ac position: ",AcK,"; other mod.: ",other.mods,", n=",length(temp.filt$occup3)," (",length(temp0$occup3),")"),collapse="")
						title(main=paste(c(project,caption),collapse="; "),cex.main=0.8,line=3)
						title(main=rawNames,cex.main=0.6,line=2,font.main = 1)
					
						temp<-temp.filt
						if(nrow(temp)>0){
							ionsorted.list	<- by(cbind(temp$intens,temp$occup3),temp$ms2seq,function(x) x)
							par(mar=c(10,4,4,4))
							plot(0:(length(ionsorted.list)+1),rep(0,length(ionsorted.list)+2),ylim=c(0,max(100,ceiling(max(temp$occup3,na.rm=T)))),ylab="occupancy (%)",type="n",xlab="",xaxt="n")
							abline(v=1:length(ionsorted.list),lty=2)
							qqq	<- sapply(1:length(ionsorted.list),function(i) points(rep(i,nrow(ionsorted.list[[i]])),unlist(ionsorted.list[[i]][2]),pch=16,col="cornflowerblue"))
							abline(h=seq(0,100,10),col="grey")
							abline(h=plotmed,col="red",lwd=2)
							title(main=AcK,line=3)
							axis(3,at=1:length(ionsorted.list),sapply(ionsorted.list,function(x) nrow(x)),line=0,las=2,cex.axis 	=0.8)
							axis(1,at=1:length(ionsorted.list),  names(ionsorted.list),line=0,las=2,cex.axis 	=0.8)
							boxplot(lapply(ionsorted.list,function(x) unlist(x[2])),las=2,ylim=c(0,max(100,ceiling(max(temp$occup3,na.rm=T)))),add=T,xaxt="n",yaxt="n")
						}
					}else{
						caption	<- paste(c("Ac position: ",AcK,"; other mod.: ",other.mods,", n=",length(temp.filt$occup3)," (",length(temp0$occup3),")"),collapse="")
						title(main=caption)
					}
				}
			}#other.mods
		}#AcK
		dev.off()

		if(length(filtMat)>0){
			filtMat$moduniAcK	<- apply(filtMat,1,function(x) paste(c(x["AcK"],x["pep_var_mod.4"]),collapse="_"))
			plotList	<- split(filtMat$occup3,filtMat$moduniAcK)
			plotList.length	<- unlist(lapply(plotList,length))
			namevecPlotlist	<- gsub("_.*","",names(plotList))
			namevecPlotlist	<- as.numeric(gsub("N","0.",namevecPlotlist))
			sortVec.Plotlist.names	<- order(namevecPlotlist,-plotList.length)
			plotList	<- plotList[names(plotList)[sortVec.Plotlist.names]]
			plotList	<- plotList[unlist(lapply(plotList,length))>2]
			if(length(plotList)>0){
				GraphOutName	<- paste(c(outDirCombi,"\\combi_",gsub("^.*[|:]","",accession),".pdf"), collapse="")
				unlist(GraphOutName)
				pdf(GraphOutName)
				xaxstxt	<-  as.data.frame(strsplit(names(plotList),"_"),stringsAsFactors=F)[1,]
				xaxstxt2	<- as.character(as.data.frame(strsplit(gsub("_","_ ",names(plotList)),"_"),stringsAsFactors=F)[2,])
				xaxstxt2	<- gsub("^ ","",xaxstxt2)
				xaxstxt2	<- gsub("NA","",xaxstxt2)
				xaxstxt2	<- gsub("Protein N-term","Pr.N",xaxstxt2)
				xaxstxt2	<- gsub("Acetyl","AcK",xaxstxt2)
				xaxstxt2	<- gsub("Monomethyl","Monomet.",xaxstxt2)
				xaxstxt2	<- gsub("Dimethyl","Dimet.",xaxstxt2)
				xaxstxt2	<- gsub("Trimethyl","Trimet.",xaxstxt2)
				xaxstxt2	<- gsub("\\(","_",xaxstxt2)
				xaxstxt2	<- gsub(")","",xaxstxt2)
				par(mar=c(10,5,4,2),cex=0.8)

				plot(c(0.5,length(plotList)+0.5),c(0,100),type="n",axes=F,xlab="",ylab="")
				abline(h=50,col="blue")
				abline(h=seq(5,95,5),col="lightgrey")
				abline(h=seq(0,100,10),col="grey")
				stripchart(plotList,vertical=T,cex=1,pch="-",add=T,col="darkgrey")
				median	<- unlist(lapply(plotList,function(x) median(x)))
				mad		<- unlist(lapply(plotList,function(x) mad(x)))
				sd		<- unlist(lapply(plotList,function(x) sd(x)))
				points(1:length(plotList),median,pch=16,col="black")
				arrows(1:length(plotList), median-mad-0.025,1:length(plotList), median+mad+0.025, length=0.05, angle=90, code=3,col="black")
				box("plot")
				axis(side=2,pretty(c(0,100)),las=1)
				axis(side=1,mgp = c(3, 3, 0),at=c(1:length(plotList)),labels=NA,las=2,line=0,cex.axis 	=0.8)
				axis(side=1,mgp = c(3, 3, 0),at=c(1:length(plotList)),labels=xaxstxt,las=2,line=-2,lty=0,cex.axis 	=0.8)
				axis(side=1,mgp = c(3, 3, 0),at=c(1:length(plotList)),labels=xaxstxt2,las=2,line=0,lty=0,cex.axis 	=0.6)
				axis(side=3,at=0:(length(plotList)),c("n=",unlist(lapply(plotList,function(x) length(x)))),las=2,line=-.5,tck=0,lty=0,cex.axis 	=0.8)
				title(ylab="Acetyl-Lysine occupancy (%)",line=2.5)
				title(main=paste(c(project,accession),collapse="; "),cex.main=1,line=3)
				title(main=rawNames,cex.main=0.6,line=2,font.main = 1)
				mtext(timestamp,side=1,adj=1,line=7,cex=0.8)
				dev.off()
			}
		}
	}
	write.table(condensMat, file = paste(c(outDir,"\\condensMat_Acetylation_quant.txt"),collapse=""), quote = FALSE, sep = "\t",na = "NA", dec = ".", row.names = FALSE)
	outFilename	<- paste(c(outDir,"\\Combined_Acetylation_quant.txt"),collapse="")
	write.table(outMat, file = outFilename, quote = FALSE, sep = "\t",na = "NA", dec = ".", row.names = FALSE)
}
##############################

write.elements.list	<- function(atomcols,reagent.purity){
	fracval	<- (reagent.purity-0.011)+0.011
	elements	<- initializeElements(atomcols[!atomcols %in% "Cx"])
	elements[[length(elements)+1]]	<- elements[[which(sapply(elements,function(x) x["name"] %in% "C"))]]
	elements[[length(elements)]]$name	<- "Cx"
	elements[[length(elements)]]$isotope$mass		<- elements[[which(sapply(elements,function(x) x["name"]%in%"C"))]]$isotope$mass[1:2]
	elements[[length(elements)]]$isotope$abundance	<- c(1-fracval,fracval)
	return(elements)
}
write.elements.list2	<- function(atomcols,reagent.purity,occup1){
	fracval	<- (reagent.purity-0.011)+0.011
	elements	<- initializeElements(atomcols)
	elements[[length(elements)+1]]	<- elements[[which(sapply(elements,function(x) x["name"] %in% "C"))]]
	elements[[length(elements)]]$name	<- "Cx"
	elements[[length(elements)]]$isotope$mass		<- elements[[which(sapply(elements,function(x) x["name"]%in%"C"))]]$isotope$mass[1:2]
	elements[[length(elements)]]$isotope$abundance	<- c(1-fracval,fracval)

	fracval	<- (reagent.purity-0.011)*occup1/100+0.011
	elements[[length(elements)+1]]	<- elements[[which(sapply(elements,function(x) x["name"] %in% "C"))]]
	elements[[length(elements)]]$name	<- "Cz"
	elements[[length(elements)]]$isotope$mass		<- elements[[which(sapply(elements,function(x) x["name"]%in%"C"))]]$isotope$mass[1:2]
	elements[[length(elements)]]$isotope$abundance	<- c(fracval,1-fracval)

	return(elements)
}

write.elements.list3	<- function(atomcols,reagent.purity,occup1){
	fracval	<- (reagent.purity-0.011)+0.011
	elements	<- initializeElements(atomcols)
	elements[[length(elements)+1]]	<- elements[[which(sapply(elements,function(x) x["name"] %in% "C"))]]
	elements[[length(elements)]]$name	<- "Cx"
	elements[[length(elements)]]$isotope$mass		<- elements[[which(sapply(elements,function(x) x["name"]%in%"C"))]]$isotope$mass[1:2]
	elements[[length(elements)]]$isotope$abundance	<- c(1-fracval,fracval)

	fracval	<- (reagent.purity-0.011)*occup1[1]/100+0.011
	elements[[length(elements)+1]]	<- elements[[which(sapply(elements,function(x) x["name"] %in% "C"))]]
	elements[[length(elements)]]$name	<- "Cz"
	elements[[length(elements)]]$isotope$mass		<- elements[[which(sapply(elements,function(x) x["name"]%in%"C"))]]$isotope$mass[1:2]
	elements[[length(elements)]]$isotope$abundance	<- c(fracval,1-fracval)

	fracval	<- (reagent.purity-0.011)*occup1[2]/100+0.011
	elements[[length(elements)+1]]	<- elements[[which(sapply(elements,function(x) x["name"] %in% "C"))]]
	elements[[length(elements)]]$name	<- "Cy"
	elements[[length(elements)]]$isotope$mass		<- elements[[which(sapply(elements,function(x) x["name"]%in%"C"))]]$isotope$mass[1:2]
	elements[[length(elements)]]$isotope$abundance	<- c(fracval,1-fracval)

	return(elements)
}

############################################################################
###MAIN#####################################################################
############################################################################
elements		<- write.elements.list(atomcols,reagent.purity)
subFolders		<- list.files(IDfolder)
subFolders		<- subFolders[!grepl("combined",subFolders)]
if (!file.exists(pastePath(c(resultsFolder,project)))){
    dir.create(pastePath(c(resultsFolder,project)))
}
rawNames	<- paste(subFolders,collapse=", ")
IDdataLocations	<- c()
quandataLocations	<- c()
resultsLocations	<- c()
for(temp in subFolders){
	IDdataLocations	<- c(IDdataLocations,pastePath(c(IDfolder,temp)))
	quandataLocations	<- c(quandataLocations,pastePath(c(quanfolder,temp)))
	resultsLocations	<- c(resultsLocations,pastePath(c(resultsFolder,project,temp)))
	if (!file.exists(pastePath(c(resultsFolder,project,temp)))){
	    dir.create(pastePath(c(resultsFolder,project,temp)))
	}
}
tempsplit	<- unlist(strsplit(resultsLocations[1],"\\\\"))
timestamp	<- format(Sys.time(), format="%Y%m%d.%H%M")
combineDir	<- paste(c(tempsplit[1:(length(tempsplit)-1)],"combined"),collapse="\\")
if (!file.exists(combineDir)){
	    dir.create(combineDir)
}
combineDir	<- paste(c(combineDir,timestamp),collapse="\\")
if (!file.exists(combineDir)){
	    dir.create(combineDir)
}
for(subi in  1:length(subFolders)){
	
	cat(paste0(Sys.time(),paste(c("   start: ",resultsLocations[subi]),collapse=""),"\n"))
	flush.console()
	outDir		<- make.outDir(dataLocation=resultsLocations[subi],subDir="analysis") 	###make output directory
	cat(paste0(Sys.time()," starting mzdata"),"\n")
	flush.console()
	passargs1	<- c(quandataLocations[subi],outDir)
	system(paste(c("RScript.exe --vanilla",shQuote(paste(c(ScriptFolder,"\\execute_mzData_processing.R"),collapse="")),shQuote(passargs1)),collapse=" "))
	load(pastePath(c(outDir,"mzData.RData")))
	dfout2		<- tempList[[1]]
	pl			<- tempList[[2]]
	mzDataFilename	<- tempList[[3]]
	unlink(pastePath(c(outDir,"mzData.RData")))

	cat(paste0(Sys.time(),"  mzdata - done"),"\n")
	flush.console()

	passargs2	<- c(IDdataLocations[subi],outDir,paste(atomcols,collapse=";"))
	system(paste(c("RScript.exe --vanilla",shQuote(paste(c(ScriptFolder,"\\execute_mascot_XML_processing.R"),collapse="")),shQuote(passargs2)),collapse=" "))
	load(pastePath(c(outDir,"mascotXML.RData")))
	pepMat		<- tempList[[1]]
	masses		<- tempList[[2]]
	var.mod		<- tempList[[3]]
	XMLFilename		<- tempList[[4]]
	unlink(pastePath(c(outDir,"mascotXML.RData")))

	cat(paste0(Sys.time(),"  xml - done"),"\n")
	flush.console()
	pepMat$pep_start	<- as.numeric(pepMat$pep_start)
	pepMat$pep_end	<- as.numeric(pepMat$pep_end)

	pepMat	<- pepMat[!(grepl("Protein N-term",pepMat$pep_var_mod)&!grepl("^J",pepMat$pep_seq)),]
	pepMat	<- pepMat[!pepMat$pep_start==2,]

###replace O and J in peptide sequences and modifications
	pepMat.temp	<- pepMat
	pepMat.temp$pep_var_mod_pos_add	<- NA

	for(i in 1:nrow(pepMat.temp)){
		x	<- pepMat.temp[i,]
		if(is.na(x$pep_var_mod_pos)){
			x$pep_var_mod_pos	<- paste(c("0.",rep(0,nchar(x$pep_seq)),".0"),collapse="")
		}
		temp_addvarmod	<- rep(0,nchar(x$pep_seq))
		if(TRUE%in%grepl("O",x$pep_seq)){
			temp_addvarmod[as.vector(gregexpr("O",x$pep_seq)[[1]])]	<- "A"
		}
		x$pep_var_mod_pos_add	<- paste(c("0.",temp_addvarmod,".0"),collapse="")
		if(!grepl("noMet$|_x$|x_",x$accession)) {
			x$pep_start 	<- x$pep_start -1
			x$pep_end	<- x$pep_end-1
		}
		if(grepl("^J",x$pep_seq)){
			x$pep_start 			<- x$pep_start
			x$pep_end			<- x$pep_end-1
			x$pep_seq			<- gsub("^J","",x$pep_seq)
			temp_addvarmod		<- temp_addvarmod[2:length(temp_addvarmod)]
			x$pep_var_mod_pos_add	<- paste(c("B.",temp_addvarmod,".0"),collapse="")
			temp.varmod			<- unlist(strsplit(x$pep_var_mod_pos,"\\."))
			temp.varmod2			<- paste(unlist(strsplit(temp.varmod[2],""))[2:nchar(temp.varmod[2])],collapse="")
			x$pep_var_mod_pos		<- paste(c(temp.varmod[1],temp.varmod2,temp.varmod[3]),collapse=".")
		}else{
			if(grepl("^agms:|histone-H|Histone-H|histone H|Histone H",x$accession)|grepl("histone-H|Histone-H|histone H|Histone H",x$prot_desc)){
				x$pep_start 	<- x$pep_start-1
				x$pep_end	<- x$pep_end-1
			}
		}
		x$pep_seq	<- gsub("O","K",x$pep_seq)
		temp.pep_var_mod_pos	<- unlist(strsplit(x$pep_var_mod_pos,""))
		temp.pep_var_mod_pos[temp.pep_var_mod_pos %in% 0 & unlist(strsplit(x$pep_var_mod_pos_add,"")) %in% "A"]	<- "A"
		temp.pep_var_mod_pos[temp.pep_var_mod_pos %in% 0 & unlist(strsplit(x$pep_var_mod_pos_add,"")) %in% "B"]	<- "B"

		x$pep_var_mod_pos	<- paste(temp.pep_var_mod_pos,collapse="")	

		pepMat.temp[i,]	<- x
	}
	pepMat	<- pepMat.temp

	pepMat$pep_calc_mr	<- as.numeric(pepMat$pep_calc_mr)
	pepMat$pep_delta		<- as.numeric(pepMat$pep_delta)
	pepMat$pep_delta_ppm	<- 1E6*pepMat$pep_delta/pepMat$pep_calc_mr
	precursor.cutoffs		<- boxplot(pepMat$pep_delta_ppm[is.na(pepMat$ambiguous)],range=3,plot=F)$stats[c(1,5)]
	pepMat			<- pepMat[!pepMat$pep_delta_ppm<precursor.cutoffs[1]&!pepMat$pep_delta_ppm>precursor.cutoffs[2],]#robust filter ppm outliers
	pepMat$pep_var_mod.3	<- ""

	var.mod	<- rbind(var.mod,data.frame(var.mod.id="A",var.mod.names="Acetyl 13C1 (K)",var.mod.delta=43.013920,C=1,H=2,N=0,O=1,S=0,Cx=1))
	var.mod	<- rbind(var.mod,data.frame(var.mod.id="B",var.mod.names="Acetyl 13C1 (Protein N-term)",var.mod.delta=43.013920,C=1,H=2,N=0,O=1,S=0,Cx=1))

	var.mod.temp	<- var.mod
	var.mod.temp[grepl("AcQnt",var.mod$var.mod.names),atomcols]	<-  ldply(apply(var.mod.temp[grepl("AcQnt",var.mod$var.mod.names),atomcols],1,function(x) x+var.mod.temp[var.mod.temp$var.mod.id =="A",atomcols]))[,atomcols]
	var.mod.temp$var.mod.delta	<- as.numeric(var.mod$var.mod.delta)
	var.mod.temp$var.mod.delta[grepl("AcQnt",var.mod.temp$var.mod.names)]	<-  apply(cbind(var.mod.temp$var.mod.delta[grepl("AcQnt",var.mod.temp$var.mod.names)]),1,function(x) x+var.mod.temp[var.mod.temp$var.mod.id =="A","var.mod.delta"])
	var.mod.temp$var.mod.names	<- gsub("AcQnt | 12C1| 13C1","",var.mod.temp$var.mod.names)
	var.mod.temp$var.mod.names	<- gsub("\\(O)","\\(K)",var.mod.temp$var.mod.names)
	var.mod.temp$var.mod.names	<- gsub("^Acetyl.*[(]","Acetyl (",var.mod.temp$var.mod.names,perl=TRUE)

	Ato	<- var.mod.temp$var.mod.id[var.mod.temp$var.mod.names %in% "Acetyl (K)" & !var.mod.temp$var.mod.id %in% c("A","B")]
	Bto	<- var.mod.temp$var.mod.id[var.mod.temp$var.mod.names %in% "Acetyl (Protein N-term)" & !var.mod.temp$var.mod.id %in% c("A","B")]
	var.mod.temp	<- var.mod.temp[!var.mod.temp$var.mod.id %in% c("A","B"),]
	var.mod	<- var.mod.temp
	pepMat$pep_var_mod_pos	<- gsub("A",Ato,pepMat$pep_var_mod_pos)
	pepMat$pep_var_mod_pos	<- gsub("B",Bto,pepMat$pep_var_mod_pos)

	for(line.i in 1:nrow(pepMat)){
		temp.line		<- pepMat[line.i,]
		test_pep_pos	<- (as.numeric(temp.line$pep_start)):(as.numeric(temp.line$pep_end))
		test_pep_mod	<- unlist(strsplit(unlist(strsplit(temp.line$pep_var_mod_pos,"[.]"))[2],""))
		tempMat		<- data.frame(pos=c("N",test_pep_pos,"C"),mod=NA,pos2=NA,stringsAsFactors=FALSE)

		test_N_mod		<- unlist(strsplit(temp.line$pep_var_mod_pos,"[.]"))[1]
		test_N_pos		<- "N"
		if(test_N_mod!="0"){
			if(grepl("Protein N-term",var.mod$var.mod.names[var.mod$var.mod.id==test_N_mod])){
				test_N_pos	<- gsub("^","N",test_pep_pos[1])
			}
		}
		test_N_mod		<- var.mod[var.mod$var.mod.id == test_N_mod,"var.mod.names"]
		if(length(test_N_mod)==0){test_N_mod<-""}
		tempMat[tempMat[,1]=="N",2:3]		<- c(test_N_mod,test_N_pos)

		test_C_mod		<- unlist(strsplit(temp.line$pep_var_mod_pos,"[.]"))[3]
		test_C_pos		<- "C"
		test_C_mod		<- paste(c(var.mod[var.mod$var.mod.id == test_C_mod,"var.mod.names"]))
		if(length(test_C_mod)==0){test_C_mod<-""}
		tempMat[tempMat[,1]=="C",2:3]		<- c(test_C_mod,test_C_pos)

		for(pos in 1:length(test_pep_pos)){
			test_pos		<- test_pep_pos[pos]
			test_mod		<- paste(c(var.mod[var.mod$var.mod.id == test_pep_mod[pos],"var.mod.names"]))
			if(length(test_mod)==0){test_mod<-""}
			tempMat[tempMat[,1]==test_pos,2:3]		<- c(test_mod,test_pos)
		}
		tempMat	<- tempMat[tempMat$mod!=""&!grepl("pyro-Glu|Oxidat",tempMat$mod),2:3]
		if(nrow(tempMat)>0){
			pepMat[line.i,]$pep_var_mod.3	<- paste((apply(tempMat,1,function(x) paste(c(gsub(" \\(.*)| .*","",x["mod"]),"(",x["pos2"],")"),collapse=""))),collapse=";")
		}
	}
	###duplicate peptides with different combinations of 12C and 13C acetyl are eliminated
	pepMat			<- pepMat[order(pepMat$rank),]
	selcolnames			<- c("query","accession","pep_seq","raw","pep_var_mod_pos","pep_var_mod_pos_add")
	pepMat$selcol		<- apply(pepMat,1,function(x) paste(c(x[selcolnames]),collapse=";;;"))
	selected			<- table(pepMat$selcol)
	pepMat$filt.uni1		<- NA
	duplicate.queries		<- unique(names(selected)[selected>1])
	unique.queries		<- names(selected)[selected==1]

	if(length(unique.queries)>0){
		pepMat$filt.uni1[pepMat$selcol%in%unique.queries]	<- 0
	}
	if(length(duplicate.queries)>0){
		for(selquery  in duplicate.queries){
			pepMat$filt.uni1[pepMat$selcol%in%selquery][1]	<- 1
		}
	}

	pepMat		<- pepMat[!is.na(pepMat$filt.uni1),]
	pepMat.nofilt	<- pepMat
	pepMat 		<- pepMat[grepl("Acetyl",pepMat$pep_var_mod.3),]
	ProcData	<- calc.fragmentions(pepMat,dfout2,reagent.purity,atomcols,masses,var.mod,pl,outDir,mzDataFilename,XMLFilename,errortolPPM,elements,cores)
	save(ProcData,file=pastePath(c(combineDir,paste0(subFolders[subi],".procData"))))

	outFilename		<- paste(c(combineDir,"\\",mzDataFilename,"_",XMLFilename,"_pepMat_noFilt.txt"),collapse="")
	write.table(pepMat.nofilt, file = outFilename, quote = FALSE, sep = "\t",na = "NA", dec = ".", row.names = FALSE)
	cat(paste0(Sys.time(),paste(c("   done: ",resultsLocations[subi]),collapse=""),"\n"))
	cat(paste0(Sys.time(),"  process - done"),"\n")
	flush.console()
}


save.image(pastePath(c(combineDir,"tempimage.RData")))

###open combineDir, load all ProcData Objects
ProcData2		<- c()
FileList		<- list.files(combineDir)

PepMat.files	<- FileList[grepl("._pepMat_noFilt.txt$",FileList)]
PepMat.files	<- apply(cbind(PepMat.files),1,function(x) pastePath(c(combineDir,x)))
PepMat2		<- c()
for(inMatname in PepMat.files){
	inMat<- read.delim(inMatname, sep="\t",header=TRUE,stringsAsFactors=FALSE)
	PepMat2<- rbind(PepMat2,inMat)
}

PepMat2$pep_var_mod.4 <- apply(PepMat2,1,function(x) {
	temp	<- unlist(strsplit(x["pep_var_mod.3"],";"))
	temp	<- temp[!grepl("^Acetyl",temp)]
	return(paste(temp,collapse=";"))})
PepMat2	<- PepMat2[order(-PepMat2$prot_matches,PepMat2$pep_start,PepMat2$pep_end,-PepMat2$pep_score,PepMat2$pep_var_mod.4),]
outFilename		<- paste(c(combineDir,"\\","pepMat_noFilt.txt"),collapse="")
write.table(PepMat2, file = outFilename, quote = FALSE, sep = "\t",na = "NA", dec = ".", row.names = FALSE)

ProcfileList	<- FileList[grepl(".procData$",FileList)]
for(ProcFile in ProcfileList){
	load(pastePath(c(combineDir,ProcFile)))
	ProcData2	<- c(ProcData2,ProcData)
}

combinedResults	<- combineResults(ProcData2,combineDir)
cat(paste0(Sys.time(),"  done: combine results"),"\n")
flush.console()
combinedResults$other.mods	<- unlist(lapply(strsplit(combinedResults$other.mods,";"),function(x) paste(x[!grepl("^Acetyl",x)],collapse=";")))
combinedResults$pep_var_mod.4	<- unlist(lapply(strsplit(combinedResults$pep_var_mod.3,";"),function(x) paste(x[!grepl("^Acetyl",x)],collapse=";")))

cat(paste0(Sys.time(),"  done: writing mod columns in combined results"),"\n")

#####new for 2 K residues in fragment

############################################
#+++++++++++++++++++++++++++++++++++++++++++
############################################
combinedResults$filt0			<- 0
combinedResults$filt0[combinedResults$AcetylMods==1 & combinedResults$intens>cutoff.intens]	<- 1
combinedResults$filt0[combinedResults$AcetylMods==2 & combinedResults$intens2>cutoff.intens]	<- 1
combinedResults$filt0[combinedResults$AcetylMods==3 & combinedResults$intens3>cutoff.intens]	<- 1

#!!!# here I filter/assign not resolved PTM isoforms
selcolnames			<- c("query","accession","pep_seq","raw","pep_var_mod.4")
combinedResults$selcol	<- do.call(paste, c(combinedResults[,selcolnames], sep=";;;"))
selected			<- unique(combinedResults$selcol)
selected			<- table(unlist(lapply(strsplit(selected,";;;"),function(x) paste(c(x[c(1,2,4)]),collapse=";;;"))))
combinedResults$selcol	<- do.call(paste, c(combinedResults[,selcolnames[c(1,2,4)]], sep=";;;"))
combinedResults$filt.ambig2	<- NA
ambigous.queries		<- unique(names(selected)[selected>1])
unique.queries		<- names(selected)[selected==1]

if(length(unique.queries)>0){
	combinedResults$filt.ambig2[combinedResults$selcol%in%unique.queries]	<- 0
}
if(length(ambigous.queries)>0){
	for(acc in unique(combinedResults$accession[combinedResults$selcol%in%ambigous.queries])){
		sure.mods		<- unique(unlist(strsplit(combinedResults[combinedResults$accession==acc & combinedResults$filt.ambig2==0,]$pep_var_mod.4,";")))
		sure.mods.no.filt	<- unique(unlist(strsplit(PepMat2[PepMat2$accession == acc&is.na(PepMat2$ambiguous),]$pep_var_mod.4,";")))
		sure.mods		<- unique(c(sure.mods,sure.mods.no.filt))
		ambigous.queries.temp	<- ambigous.queries[unlist(lapply(strsplit(ambigous.queries,";;;"),function(x) x[2] %in% acc))]
		if(length(sure.mods)>0){
			for(i in 1:length(ambigous.queries.temp)){
				temp			<- combinedResults[combinedResults$accession==acc & combinedResults$selcol%in%ambigous.queries.temp[i],c("id","pep_var_mod.4")]	
				isoform.mods		<- unique(temp$pep_var_mod.4)
				isoform.mods		<- isoform.mods[unlist(lapply(strsplit(isoform.mods,";"),function(x) sum(!x %in% sure.mods)==0))]
				if(length(isoform.mods)==0){
					combinedResults$filt.ambig2[combinedResults$id %in% temp$id]			<- 1
				}
				if(length(isoform.mods)>1){
					combinedResults$filt.ambig2[combinedResults$id %in% temp$id & combinedResults$pep_var_mod.4 %in% isoform.mods]	<- 1
				}
				if(length(isoform.mods)==1){
					combinedResults$filt.ambig2[combinedResults$id %in% temp$id & combinedResults$pep_var_mod.4 %in% isoform.mods]	<- 0
				}
			}
		}
	}
}
i<-2
combinedResults	<- combinedResults[!is.na(combinedResults$filt.ambig2),]

###get 2K fragment occupancy
combinedResults1	<- combinedResults[combinedResults$AcetylMods==1&combinedResults$filt0==1,]
combinedResults1$filt0[combinedResults1$occup3.maxdiff<cutoff.maxdiff&combinedResults1$occup3.corr.coff>cutoff.corr]	<- 2

condensMat1		<- create.condens.1(combinedResults1[combinedResults1$filt0==2,],combineDir)
twoKmat			<- combinedResults[combinedResults$AcetylMods==2&combinedResults$filt0==1,]

merge1			<- merge(twoKmat,condensMat1,by.x=c("accession","AcK2a","pep_var_mod.4"),by.y=c("accession","pos","pep_var_mod.4"))
merge1$AcK		<- merge1$AcK2b
merge2			<- merge(twoKmat,condensMat1,by.x=c("accession","AcK2b","pep_var_mod.4"),by.y=c("accession","pos","pep_var_mod.4"))
merge2$AcK		<- merge2$AcK2a

colnames.merge	<- unique(c(colnames(twoKmat),colnames(condensMat1)))
colnames.merge	<- colnames.merge[colnames.merge %in% colnames(merge1)]

merge1			<- merge1[,colnames.merge]
merge2			<- merge2[,colnames.merge]
merge12			<- rbind(merge1,merge2)
merge12			<- merge12[!is.na(merge12$med.filt),]

secOccup.calc	<- function(x){
		median.occ			<- as.numeric(x["med.filt"])
		elements			<- write.elements.list2(atomcols[!atomcols %in% "Cx"],reagent.purity,median.occ)
		atoms				<- as.numeric(x[atomcols])
		names(atoms)			<- names(x[atomcols])
		atoms["C"]			<- atoms["C"]-1
		atoms				<- c(atoms,Cz=1)
		atoms.light			<- paste0(names(atoms)[atoms>0],atoms[atoms>0],collapse="")
		atoms["C"]			<- atoms["C"]-1
		atoms["Cx"]			<- atoms["Cx"]+1
		atoms.heavy			<- paste0(names(atoms)[atoms>0],atoms[atoms>0],collapse="")	

		isodist.light		<- getMolecule(atoms.light,elements=elements,maxisotopes=10)$isotopes[[1]]
		isodist.light[2,]		<- isodist.light[2,]/sum(isodist.light[2,])*100
		isodist.heavy		<- getMolecule(atoms.heavy,elements=elements,maxisotopes=10)$isotopes[[1]]
		isodist.heavy[2,]		<- isodist.heavy[2,]/sum(isodist.heavy[2,])*100

		peaks.cols			<- c("minus1","mono","plus1","plus2","plus3","plus4","plus5","plus6")
		peaks.sel			<- which(peaks.cols %in% "mono"):(which(peaks.cols %in% "plus4")+atoms["Cx"]-1)	
		peaks.sel			<- peaks.sel[peaks.sel<=length(peaks.cols)]

		peaks3				<- as.numeric(c(x[peaks.cols[peaks.sel]]))
		peaks3				<- peaks3/max(peaks3)*100

		peaks.sel			<- which(peaks.cols %in% "minus1"):(which(peaks.cols %in% "plus4")+atoms["Cx"]-1)	
		peaks.sel			<- peaks.sel[peaks.sel<=length(peaks.cols)]

	
		peaks4				<- as.numeric(c(x[peaks.cols[peaks.sel]]))
		peaks4				<- peaks4/max(peaks4)*100

		theoratios			<- seq(0,100,by=.5)
		isodist.sum			<- t(outer(isodist.light[2,1:length(peaks3)],theoratios,"*") + outer(isodist.heavy[2,1:length(peaks3)],(100-theoratios),"*"))
		theopeaks			<- isodist.sum/apply(isodist.sum,1,max)*100
		diffpattern			<- apply(theopeaks,1,function(x) sum((peaks3-x)^2))


		which.min			<- which(diffpattern==min(diffpattern))
		occupancy			<- median(theoratios[which.min])
		theoratios			<- seq(occupancy-2.5,occupancy+2.5,by=0.01)
		theoratios			<- theoratios[theoratios>=0 & theoratios<=100]

		isodist.sum			<- t(outer(isodist.light[2,1:length(peaks3)],theoratios,"*") + outer(isodist.heavy[2,1:length(peaks3)],(100-theoratios),"*"))
		theopeaks			<- isodist.sum/apply(isodist.sum,1,max)*100
		diffpattern			<- apply(theopeaks,1,function(x) sum((peaks3-x)^2))

		which.min			<- which(diffpattern==min(diffpattern))
		occupancy			<- median(theoratios[which.min])

		theopeaks4			<- c(0,isodist.sum[which.min,]/max(isodist.sum[which.min,])*100)
		peaks4.ori			<- peaks4
		theopeaks4.ori		<- theopeaks4
		corr.coff			<- cor(peaks4,theopeaks4)
		sumdiff				<- diffpattern[which.min]
		maxdiff				<- max(abs(peaks4-theopeaks4),na.rm=T)
		return(data.frame(occupancy,corr.coff,maxdiff,sumdiff,peaks4=paste(peaks4.ori,collapse=";"),theopeaks4=paste(theopeaks4.ori,collapse=";"),stringsAsFactors=FALSE))
}

occup2.results	<- vector("list", nrow(merge12))
for(i in 1:nrow(merge12)){
		occup2.results[[i]]	<- secOccup.calc(merge12[i,])
}

cat(paste0(Sys.time(),"  second AcK calculation done","\n"))
flush.console()

occup2.results		<- ldply(occup2.results)
merge12[,c("occup3","occup3.corr.coff","occup3.maxdiff","occup3.sumdiff","peaks4","theopeaks4")]	<- occup2.results
merge12$filt0[merge12$occup3.maxdiff<cutoff.maxdiff&merge12$occup3.corr.coff>cutoff.corr]	<- 2

merge12		<- merge12[,colnames(combinedResults1)]
combinedResults2	<- rbind(combinedResults1,merge12)


cat(paste0(Sys.time(),"  merging - done"),"\n")
flush.console()

tempList			<- create.condens.2(combinedResults2,combineDir,cutoff.intens,timestamp,project,rawNames)
condensMat2	<- tempList[[1]]
combined12	<- tempList[[2]]


#############################################
###Calculation for 3 lysine fragments (O)~(O)
#############################################
threeKmat		<- combinedResults[combinedResults$AcetylMods==3&combinedResults$filt0==1,]

inMat	<- threeKmat
outMat	<- c()
for(i in 1:nrow(threeKmat)){
	tempRow	<- inMat[i,]
	outRow	<- rbind(tempRow,tempRow,tempRow)
	all.AcK3	<-  as.character(c(tempRow$AcK3a,tempRow$AcK3b,tempRow$AcK3c))
	outRow$AcK	<- all.AcK3
	selcondens	<- condensMat2[condensMat2$accession %in% tempRow$accession,]
	selcondens	<- selcondens[selcondens$pep_var_mod.4 %in% tempRow$pep_var_mod.4,]
	selcondens$pos	<- as.character(selcondens$pos)

	if(sum(outRow$AcK %in% selcondens$pos)>=2){
		for(combi in 1:3){
			if(sum(selcondens$pos %in% all.AcK3[!all.AcK3 %in% outRow$AcK[combi]])==2){
				median.occ	<- selcondens$med.filt[selcondens$pos %in% all.AcK3[!all.AcK3 %in% outRow$AcK[combi]]]
				elements			<- write.elements.list3(atomcols[!atomcols %in% "Cx"],reagent.purity,median.occ)
				atoms				<- as.numeric(outRow[combi,atomcols])
				names(atoms)			<- names(outRow[combi,atomcols])
				atoms["C"]			<- atoms["C"]-2
				atoms				<- c(atoms,Cz=1)
				atoms				<- c(atoms,Cy=1)
				atoms.light			<- paste0(names(atoms)[atoms>0],atoms[atoms>0],collapse="")

				atoms["C"]			<- atoms["C"]-1
				atoms["Cx"]			<- atoms["Cx"]+1
				atoms.heavy			<- paste0(names(atoms)[atoms>0],atoms[atoms>0],collapse="")	

				isodist.light		<- getMolecule(atoms.light,elements=elements,maxisotopes=10)$isotopes[[1]]
				isodist.light[2,]		<- isodist.light[2,]/sum(isodist.light[2,])*100
				isodist.heavy		<- getMolecule(atoms.heavy,elements=elements,maxisotopes=10)$isotopes[[1]]
				isodist.heavy[2,]		<- isodist.heavy[2,]/sum(isodist.heavy[2,])*100

				peaks.cols			<- c("minus1","mono","plus1","plus2","plus3","plus4","plus5","plus6")
				peaks.sel			<- which(peaks.cols %in% "mono"):(which(peaks.cols %in% "plus4")+atoms["Cx"])	
				peaks.sel			<- peaks.sel[peaks.sel<=length(peaks.cols)]

				peaks3				<- as.numeric(c(outRow[combi,peaks.cols[peaks.sel]]))
				peaks3				<- peaks3/max(peaks3)*100	

				peaks.sel			<- which(peaks.cols %in% "minus1"):(which(peaks.cols %in% "plus4")+atoms["Cx"])	

				peaks.sel			<- peaks.sel[peaks.sel<=length(peaks.cols)]

	
				peaks4				<- as.numeric(c(outRow[combi,peaks.cols[peaks.sel]]))
				peaks4				<- peaks4/max(peaks4)*100

				theoratios			<- seq(0,100,by=.5)
				isodist.sum			<- t(outer(isodist.light[2,1:length(peaks3)],theoratios,"*") + outer(isodist.heavy[2,1:length(peaks3)],(100-theoratios),"*"))
				theopeaks			<- isodist.sum/apply(isodist.sum,1,max)*100
				diffpattern			<- apply(theopeaks,1,function(x) sum((peaks3-x)^2))


				which.min			<- which(diffpattern==min(diffpattern))
				occupancy			<- median(theoratios[which.min])
				theoratios			<- seq(occupancy-2.5,occupancy+2.5,by=0.01)
				theoratios			<- theoratios[theoratios>=0 & theoratios<=100]
	
				isodist.sum			<- t(outer(isodist.light[2,1:length(peaks3)],theoratios,"*") + outer(isodist.heavy[2,1:length(peaks3)],(100-theoratios),"*"))
				theopeaks			<- isodist.sum/apply(isodist.sum,1,max)*100
				diffpattern			<- apply(theopeaks,1,function(x) sum((peaks3-x)^2))
	
				which.min			<- which(diffpattern==min(diffpattern))
				occupancy			<- median(theoratios[which.min])
	
				theopeaks4			<- c(0,isodist.sum[which.min,]/max(isodist.sum[which.min,])*100)
				peaks4.ori			<- peaks4
				theopeaks4.ori		<- theopeaks4
				corr.coff			<- cor(peaks4,theopeaks4)
				sumdiff			<- diffpattern[which.min]
				maxdiff			<- max(abs(peaks4-theopeaks4),na.rm=T)
				tempdf	<- data.frame(occup3=occupancy,occup3.corr.coff=corr.coff,occup3.maxdiff=maxdiff,occup3.sumdiff=sumdiff,peaks4=paste(peaks4.ori,collapse=";"),theopeaks4=paste(theopeaks4.ori,collapse=";"),stringsAsFactors=FALSE)
				outRow[combi,c("occup3","occup3.corr.coff","occup3.maxdiff","occup3.sumdiff","peaks4","theopeaks4")]	<- tempdf
				outMat	<- rbind(outMat,outRow[combi,])
			}
		}

	}
}

occup3.results<-outMat

cat(paste0(Sys.time(),"  third AcK calculation done","\n"))
flush.console()

occup3.results$filt0[occup3.results$occup3.maxdiff<cutoff.maxdiff&occup3.results$occup3.corr.coff>cutoff.corr]	<- 2
occup3.results$filt	<- NA
combinedResults3	<- rbind(combined12,occup3.results)

plot_corr	<- FALSE
if(plot_corr == TRUE){
	cl <- makeCluster(cores)
	registerDoParallel(cl)
	GraphOutName	<- paste(c(combineDir,"\\correlation_",".pdf"), collapse="")
	unlist(GraphOutName)
	pdf(GraphOutName)

	for(row in 1:nrow(combinedResults2)){
		x		<- combinedResults2[row,]
		peaks		<- as.numeric(unlist(strsplit(x$peaks4,";")))
		theopeaks	<- as.numeric(unlist(strsplit(x$theopeaks4,";")))
		r		<- round((x$occup3.corr.coff),4)
		max.diff	<- round(max(abs(peaks-theopeaks)),2)
		sum.diff	<- round(x$occup3.sumdiff,1)
		plot(-1:(length(theopeaks)-2),theopeaks,xaxp=c(-1,length(theopeaks)-2,length(theopeaks)-1),type="h",xlab="",xlim=c(-1,(length(theopeaks)-2))+0.1,lwd=2)
		lines(-1:(length(theopeaks)-2)+0.05,peaks,xaxp=c(-1,length(theopeaks)-2,length(theopeaks)-1),type="h",col="red",lwd=2)
		title(main=paste(c("corr (r): ",r,", max.diff: ",max.diff,", sum.diff: ",sum.diff),collapse=""))
		title(sub=paste(c("id: ",x$id,"   pep_seq: ",x$pep_seq),collapse=""),line=3)
		text(x=-0.9,y=100, x$filt0)
	}
	dev.off()
	stopCluster(cl)
}

cat(paste0(Sys.time(),paste(c("   now creating condens.3"),collapse=""),"\n"))

xx			<- create.condens.3(combinedResults3,combineDir,cutoff.intens,timestamp,project,rawNames)

cat(paste0(Sys.time(),paste(c("   done project: ",project),collapse=""),"\n"))
flush.console()

sink()
sink(type = "message")
close(error.logfile.con)
x			<- file(paste(c(logFolder,"\\done_",mainName,"_",project,".txt"),collapse=""),"w")
close(x)

