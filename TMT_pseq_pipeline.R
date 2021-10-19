#Pipeline for SLN, zero imputation, and IRS on TMT data followed by DE analysis using PoissonSeq
#Last updated: Dec 16, 2020

TMT_pseq_pipeline = function(workdir, datafile, metadatafile, exp, SLN, PTM, stat, qval,compsfile){

#make sure the exp name is syntatically valid
exp = make.names(exp)

setwd(workdir)
message("Loading data files...")

#load the data file
#make the peptide ids the row names
data = read.csv(datafile,row.names="id",stringsAsFactors=FALSE)
colnames(data)[1]="Proteins"

#read in the metadata
metadataorg = read.table(metadatafile,sep="\t",header=TRUE,stringsAsFactors=FALSE)
#remove blanks
noblanks = metadataorg[!grepl('blank',metadataorg$name,ignore.case=TRUE),]
metadata=noblanks
plex = length(metadata$sample[metadata$run==1])
reps = max(metadata$run)
nums = metadata$rep
numrefs = sum(grepl("Ref",metadata$name,ignore.case=TRUE))/reps

message("Beginning data processing...")
if (PTM == "P"){
  #get rid of peptides with no phospho sites
  data = data[!(data$Number.of.Phospho..STY.==""),]
  #get the raw intensities for differential abundance
  intensities = data[,grepl('Reporter.intensity',colnames(data))]
  rawintensities = intensities[!grepl('corrected',colnames(intensities))]
  rawintensities = rawintensities[!grepl('count',colnames(rawintensities))]
  #rawintensities = intensities[grepl('corrected',colnames(intensities))]
  if(exp=="None"){
    intensitiesK=rawintensities
  }else{
    intensitiesK = rawintensities[,grepl(exp,colnames(rawintensities))] 
  }
  #remove blanks
  metadata3 = rep(metadataorg$name,each=3)
  intensitiesK = intensitiesK[,!grepl('blank',metadata3,ignore.case=TRUE)]
  #splitnames = strsplit(colnames(intensitiesK),'NLB')
  #splitnames = unlist(splitnames)
  #splitnames2 = strsplit(splitnames[seq(2,924,by=2)],"___")
  #splitnames2 = unlist(splitnames2)
  #intensitiesK = intensitiesK[,order(strtoi(splitnames2[seq(1,923,by=2)]))]
  
  message("Separating multiplicities...")
  #first we need to make a new table where we 1) take only the columns we need and 2) separate peptides by multiplicity
  #this assumes columns are always named the same thing
  currentrow=1
  for (i in 1:length(row.names(data))){
    myprotein = data[i,]
    #check number of phospho sites
    #if there is a semicolon, that means we have multiple sites and need to split into multiple rows
    if (length(grep(";",myprotein$Number.of.Phospho..STY.))>0){
      #get the info for this protein
      mydata = data.frame(myprotein$Proteins,myprotein$Positions.within.proteins,myprotein$Protein,myprotein$Fasta.headers,
                          myprotein$Localization.prob,myprotein$Number.of.Phospho..STY.,
                          myprotein$Sequence.window, myprotein$Modification.window, myprotein$Peptide.window.coverage,
                          myprotein$Phospho..STY..Probabilities, row.names(myprotein), stringsAsFactors=FALSE)
      colnames(mydata)=c("Proteins","Positions.within.proteins","Protein","Fasta.headers",
                         "Localization.prob",
                         "Phospho.Site","Sequence.window",'Modification.window',
                         'Peptide.window.coverage','Phospho.site.probs','original.id')
      #split the phospho sites by semincolon
      sites = strsplit(myprotein$Number.of.Phospho..STY.,";")
      mysites = sites[[1]]
      for (j in 1:length(mysites)){
        #if the site # is >3, we need to make it 3
        if (mysites[j]>3){
          currentsite = 3
        } else{
          currentsite=mysites[j]
        }
        myintensities = intensitiesK[i,grepl(paste('_',currentsite,sep=""),colnames(intensitiesK))]
        colnames(myintensities)=paste(rep("R",length(nums)),nums,sep="_")
        if (i==1){
          newdata = data.frame(mydata,myintensities,stringsAsFactors=FALSE)
          row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",mysites[j],sep="")
          currentrow = currentrow+1
        } else{
          newdata[currentrow,] =  data.frame(mydata,myintensities)
          row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",mysites[j],sep="")
          currentrow = currentrow+1
        }
      }
    }else{
      #get the info for this protein
      mydata = data.frame(myprotein$Proteins,myprotein$Positions.within.proteins,myprotein$Protein,myprotein$Fasta.headers,
                          myprotein$Localization.prob,myprotein$Number.of.Phospho..STY.,
                          myprotein$Sequence.window, myprotein$Modification.window, myprotein$Peptide.window.coverage,
                          myprotein$Phospho..STY..Probabilities, row.names(myprotein), stringsAsFactors=FALSE)
      colnames(mydata)=c("Proteins","Positions.within.proteins","Protein","Fasta.headers",
                         "Localization.prob",
                         "Phospho.Site","Sequence.window",'Modification.window',
                         'Peptide.window.coverage','Phospho.site.probs','original.id')
      #get the correct intensity values for this protein
      if (myprotein$Number.of.Phospho..STY.>3){
        currentsite = 3
      } else{
        currentsite=myprotein$Number.of.Phospho..STY.
      }
      myintensities = intensitiesK[i,grepl(paste('_',currentsite,sep=""),colnames(intensitiesK))]
      colnames(myintensities)=paste(rep("R",length(nums)),nums,sep="_")
      if (i==1){
        newdata = data.frame(mydata,myintensities,stringsAsFactors=FALSE)
        row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",myprotein$Number.of.Phospho..STY.,sep="")
        currentrow = currentrow+1
      } else{
        newdata[currentrow,] =  data.frame(mydata,myintensities)
        row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",myprotein$Number.of.Phospho..STY.,sep="")
        currentrow = currentrow+1
      }
    }
  }
}else if(PTM=="U"){
  #get rid of peptides with no phospho sites
  data = data[!(data$Number.of.GlyGly..K.==""),]
  #get the raw intensities for differential abundance
  intensities = data[,grepl('Reporter.intensity',colnames(data))]
  rawintensities = intensities[!grepl('corrected',colnames(intensities))]
  rawintensities = rawintensities[!grepl('count',colnames(rawintensities))]
  #rawintensities = intensities[grepl('corrected',colnames(intensities))]
  if(exp=="None"){
    intensitiesK=rawintensities
  }else{
    intensitiesK = rawintensities[,grepl(exp,colnames(rawintensities))] 
  }
  #remove blanks
  metadata3 = rep(metadataorg$name,each=3)
  intensitiesK = intensitiesK[,!grepl('blank',metadata3,ignore.case=TRUE)]
  #splitnames = strsplit(colnames(intensitiesK),'NLB')
  #splitnames = unlist(splitnames)
  #splitnames2 = strsplit(splitnames[seq(2,924,by=2)],"___")
  #splitnames2 = unlist(splitnames2)
  #intensitiesK = intensitiesK[,order(strtoi(splitnames2[seq(1,923,by=2)]))]
  
  message("Separating multiplicities...")
  #first we need to make a new table where we 1) take only the columns we need and 2) separate peptides by multiplicity
  #this assumes columns are always named the same thing
  currentrow=1
  for (i in 1:length(row.names(data))){
    myprotein = data[i,]
    #check number of phospho sites
    #if there is a semicolon, that means we have multiple sites and need to split into multiple rows
    if (length(grep(";",myprotein$Number.of.GlyGly..K.))>0){
      #get the info for this protein
      mydata = data.frame(myprotein$Proteins,myprotein$Positions.within.proteins,myprotein$Protein,myprotein$Fasta.headers,
                          myprotein$Localization.prob,myprotein$Number.of.GlyGly..K.,
                          myprotein$Sequence.window, myprotein$Modification.window, myprotein$Peptide.window.coverage,
                          myprotein$GlyGly..K..Probabilities, row.names(myprotein), stringsAsFactors=FALSE)
      colnames(mydata)=c("Proteins","Positions.within.proteins","Protein","Fasta.headers",
                         "Localization.prob",
                         "GlyGly.Site","Sequence.window",'Modification.window',
                         'Peptide.window.coverage','GlyGly.site.probs','original.id')
      #split the phospho sites by semincolon
      sites = strsplit(myprotein$Number.of.GlyGly..K.,";")
      mysites = sites[[1]]
      for (j in 1:length(mysites)){
        #if the site # is >3, we need to make it 3
        if (mysites[j]>3){
          currentsite = 3
        } else{
          currentsite=mysites[j]
        }
        myintensities = intensitiesK[i,grepl(paste('_',currentsite,sep=""),colnames(intensitiesK))]
        colnames(myintensities)=paste(rep("R",length(nums)),nums,sep="_")
        if (i==1){
          newdata = data.frame(mydata,myintensities,stringsAsFactors=FALSE)
          row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",mysites[j],sep="")
          currentrow = currentrow+1
        } else{
          newdata[currentrow,] =  data.frame(mydata,myintensities)
          row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",mysites[j],sep="")
          currentrow = currentrow+1
        }
      }
    }else{
      #get the info for this protein
      mydata = data.frame(myprotein$Proteins,myprotein$Positions.within.proteins,myprotein$Protein,myprotein$Fasta.headers,
                          myprotein$Localization.prob,myprotein$Number.of.GlyGly..K.,
                          myprotein$Sequence.window, myprotein$Modification.window, myprotein$Peptide.window.coverage,
                          myprotein$GlyGly..K..Probabilities, row.names(myprotein), stringsAsFactors=FALSE)
      colnames(mydata)=c("Proteins","Positions.within.proteins","Protein","Fasta.headers",
                         "Localization.prob",
                         "GlyGly.Site","Sequence.window",'Modification.window',
                         'Peptide.window.coverage','GlyGly.site.probs','original.id')
      #get the correct intensity values for this protein
      if (myprotein$Number.of.GlyGly..K.>3){
        currentsite = 3
      } else{
        currentsite=myprotein$Number.of.GlyGly..K.
      }
      myintensities = intensitiesK[i,grepl(paste('_',currentsite,sep=""),colnames(intensitiesK))]
      colnames(myintensities)=paste(rep("R",length(nums)),nums,sep="_")
      if (i==1){
        newdata = data.frame(mydata,myintensities,stringsAsFactors=FALSE)
        row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",myprotein$Number.of.GlyGly..K.,sep="")
        currentrow = currentrow+1
      } else{
        newdata[currentrow,] =  data.frame(mydata,myintensities)
        row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",myprotein$Number.of.GlyGly..K.,sep="")
        currentrow = currentrow+1
      }
    }
  }
}else{
  #get the raw intensities for differential abundance
  intensities = data[,grepl('Reporter.intensity',colnames(data))]
  rawintensities = intensities[!grepl('corrected',colnames(intensities))]
  rawintensities = rawintensities[!grepl('count',colnames(rawintensities))]
  intensitiesK = rawintensities[,grepl(exp,colnames(rawintensities))]
  #remove blanks
  if(exp=="None"){
    intensitiesK=rawintensities
  }else{
    intensitiesK = rawintensities[,grepl(exp,colnames(rawintensities))] 
  }
  #intensitiesK = intensitiesK[,metadata$sample]
  #splitnames = strsplit(colnames(intensitiesK),'NLB')
  #splitnames = unlist(splitnames)
  #intensitiesK = intensitiesK[,order(strtoi(splitnames[seq(2,308,by=2)]))]
  #get info for all proteins
  mydata = data.frame(data$Proteins,data$Majority.protein.IDs,data$Peptide.counts..all,
                      data$Peptide.counts..unique.,data$Fasta.headers,data$Number.of.proteins,data$Peptides,
                      data$Unique.peptides, data$Sequence.coverage....,data$Unique.sequence.coverage....,
                      data$Sequence.lengths, data$Peptide.IDs, data$Evidence.IDs,
                      row.names(data), stringsAsFactors=FALSE)
  colnames(mydata)=c("Proteins","Majority.protein.IDs","Peptide.counts.all","Peptide.counts.unique",
                     "Fasta.headers",
                     "Number.of.proteins","Peptides","Unique.peptides",'Sequence.coverage',
                     'Unique.sequence.coverage','Sequence.lengths','Peptide.IDs','Evidence.IDs','id')
  newdata = data.frame(mydata,intensitiesK)
}

#replace column names with sample names from the metadata
colnames(newdata)[(dim(mydata)[2]+1):dim(newdata)[2]]=paste(metadata$name,metadata$rep,sep="_")

#get rid of all contaminants
if (PTM=="P" | PTM=="U"){
  newdata = newdata[!grepl("CON",newdata$Protein,ignore.case=TRUE),]
  newdata = newdata[!grepl("REV",newdata$Protein, ignore.case=TRUE),]
}else{
  newdata = newdata[!grepl("CON",newdata$Majority.protein.IDs,ignore.case=TRUE),]
  newdata = newdata[!grepl("REV",newdata$Majority.protein.IDs, ignore.case=TRUE),]
}

#finally add log2 transformed values because apparently people like those
intensities = newdata[,(dim(mydata)[2]+1):dim(newdata)[2]]
newdatalog = data.frame(newdata,log2(intensities+1),stringsAsFactors=FALSE)
colnames(newdatalog)=c(colnames(newdata),paste("log2(",metadata$name,"_",metadata$rep,")",sep=""))

#save table
write.csv(newdatalog,'prenormalized_data.csv')

message("Removing proteins in less than 50% of the runs...")
#now we can normalize the intensity data and do DE analysis with PoissonSeq
#first we need to get rid of proteins that are not expressed in at least half of the runs
#to do this we can just pull the first lanes from all the samples as those are the reference lanes
#if there are no references, we need to impute them
if (numrefs==0){
  refs = rowMeans(intensities)
} else{
  refs = newdata[,grepl('Ref',colnames(newdata),ignore.case=TRUE)]
}
if (numrefs>1){
  for (i in 1:reps){
      if (i==1){
        refsums = data.frame(rowSums(refs[,(numrefs*(i-1)+1):(numrefs*i)]))
      }else{
        refsums = data.frame(refsums,rowSums(refs[,(numrefs*(i-1)+1):(numrefs*i)]))
      }
    }
    #count how many times the protein is zero in the reference lane
    zerosums = rowSums(refsums[,]==0)
    #remove the proteins that are zero in >=50% of the runs
    nozeros = newdata[zerosums<=(dim(refsums)[2]/2),]
} else{
  #remove the proteins that are zero in all samples
    nozeros = newdata[refs>0,]
}
#save
write.csv(nozeros,'prenormalized_data_in_at_least_half_of_runs.csv')
finalintensities = nozeros[,(dim(mydata)[2]+1):dim(newdata)[2]]

#plot boxplot before normalization
#QC plots
colors = c("blue","green","orange","yellow","red","purple","white","blue","green","orange","yellow","red","purple","white")
png(filename='boxplot_log2.png',width=5000,height=2000,res=300)
invisible(b <- boxplot(log2(finalintensities+1),col=colors[unlist(lapply(1:reps,function(x) rep(x,plex)))],ylab="log2(Intensity)",cex.axis=0.75,las=2))
print(b)
dev.off()

if(SLN=="YES"){
  message("Sample loading normalization...")
  #peform sample loading normalization (SLN)
  if (reps>1){
    for (i in 1:reps){
      #get intensity values for this run
      myints = finalintensities[,(1+(i-1)*plex):(i*plex)]
      #column normalize
      sums = colSums(myints)
      meansums = mean(sums)
      normfactor = sums/meansums
      myintensitiesnorm = myints
      for (j in 1:dim(myints)[2]){
        myintensitiesnorm[,j] = ceiling(myints[,j]/normfactor[j])
      }
      #make a table of normalized intensities
      if (i==1){
        normintensities = myintensitiesnorm
      }else{
        normintensities = data.frame(normintensities,myintensitiesnorm)
      }
    }
    write.csv(normintensities,"sample_loading_normalization.csv")
  }else{
    myints = finalintensities
    #column normalize
    sums = colSums(myints)
    meansums = mean(sums)
    normfactor = sums/meansums
    myintensitiesnorm = myints
    for (j in 1:dim(myints)[2]){
      myintensitiesnorm[,j] = ceiling(myints[,j]/normfactor[j])
    }
    normintensities=myintensitiesnorm
    write.csv(normintensities,"sample_loading_normalization.csv")
  }
}else{
  normintensities = finalintensities
}

#zero imputation
#first we need to get the sd of the technical replicates so we can simulate missing tech reps
#first we need to look at the SD of the log2 transformed tech reps to see tech rep error
# samples = plex-numrefs
# if (reps>1){
#   message("Imputing missing values...")
#   if (numrefs==0){
#     refs = rowMeans(normintensities)
#   }else{
#     refs = normintensities[,grepl('Ref',colnames(normintensities),ignore.case=TRUE)]
#   }
#   for (i in 1:reps){
#     if (i==1){
#       refsums = data.frame(rowSums(refs[,(numrefs*(i-1)+1):(numrefs*i)]))
#     }else{
#       refsums = data.frame(refsums,rowSums(refs[,(numrefs*(i-1)+1):(numrefs*i)]))
#     }
#   }
#   zerosums = rowSums(refsums[,]==0)
#   refszero = as.matrix(log2(refs[zerosums==0,]))
#   for (i in 1:reps){
#     if (i==1){
#       refstech = data.frame(apply(refszero[,(numrefs*(i-1)+1):(numrefs*i)],1,sd))
#       refscomb = data.frame(rowMeans(refszero[,(numrefs*(i-1)+1):(numrefs*i)]))
#     }else{
#       refstech = data.frame(refstech,apply(refszero[,(numrefs*(i-1)+1):(numrefs*i)],1,sd))
#       refscomb = data.frame(refscomb,rowMeans(refszero[,(numrefs*(i-1)+1):(numrefs*i)]))
#     }
#   }
#   refstechnorm = stack(refstech/refscomb)
#   myvar = mean(na.omit(refstechnorm[,1]))
#   
#   #now we can impute references
#   logrefs = log2(refs+1)
#   logrefsimp = logrefs
#   impvals = t(rep(1,2))
#   set.seed(1)
#   for (i in 1:dim(logrefs)[1]){
#     #check for zero values
#     if (sum(logrefs[i,]==0)==numrefs){
#       mylogrefs = logrefs[i,]
#       mylogrefsnz = mylogrefs[mylogrefs>0]
#       #impute first ref
#       impvals[,1] = rnorm(1,mean(mylogrefsnz),sd(mylogrefsnz))
#       #impute follow technical reps
#       for (j in 2:numrefs){
#         impvals[,j] = impvals[,1]+myvar*impvals[,1]
#       }
#       mylogrefs[mylogrefs==0] = impvals
#       logrefsimp[i,] = mylogrefs
#     }
#   }
#   allrefsimp = floor(2^logrefsimp)
#   normintensitiesimp = cbind(normintensities[,!grepl('Ref',colnames(normintensities),ignore.case=TRUE)],allrefsimp)
#   
#   #now we need to impute the non reference values
#   #to do this we will again log2 transform everything and then use the results of the other 2 replicates to approximate the biological effect
#   normintensitiesimplog = log2(normintensitiesimp+1)
#   normintensitiesimpall = normintensitiesimplog
#   for (i in 1:dim(normintensitiesimp)[1]){
#     mylogexpr = normintensitiesimplog[i,]
#     if(sum(mylogexpr==0)>0){
#       #first get correct reference
#       missingrep = mylogexpr[,mylogexpr==0]
#       notmissing = mylogexpr[,mylogexpr!=0]
#       #if we have missing values within a run, get rid of this protein
#       if (length(missingrep)!=samples){
#         normintensitiesimpall[i,]=NA
#       }else{
#         missingrepind = colnames(missingrep)
#         missingrepnum = strsplit(missingrepind[1],"_")
#         missingrepnum = missingrepnum[[1]][2]
#         myref = mylogexpr[,grepl("Ref",colnames(mylogexpr),ignore.case=TRUE)]
#         #now we need to approximate the biological effect for each sample
#         #for this we will compare the FC of each sample to the reference
#         for (j in 1:reps){
#           if (j==1){
#             avgrefs = data.frame(rowMeans(myref[,(numrefs*(j-1)+1):(numrefs*j)]))
#           }else{
#             avgrefs = data.frame(avgrefs,rowMeans(myref[,(numrefs*(j-1)+1):(numrefs*j)]))
#           }
#         }
#         colnames(avgrefs)=c("ref1","ref2","ref3")
#         notmissingrefs = avgrefs[,!grepl(missingrepnum,colnames(avgrefs))]
#         missingref = avgrefs[,grepl(missingrepnum,colnames(avgrefs))]
#         #first 8 are divided by first ref, and so on
#         numnotmissing = (length(notmissing)-numrefs*reps)/samples
#         for (j in 1:numnotmissing){
#           if (j==1){
#             notmissingFC = data.frame(notmissing[,(samples*(j-1)+1):(samples*j)]/notmissingrefs[[j]])
#           }else{
#             notmissingFC = data.frame(notmissingFC, notmissing[,(samples*(j-1)+1):(samples*j)]/notmissingrefs[[j]])
#           }
#         }
#         #get mean and SD for each
#         #then approximate missing values by drawing from a distribution and multiplying by the reference
#         for (j in 1:samples){
#           mysamples = notmissingFC[,grep(metadata$name[j],colnames(notmissingFC))]
#           mymean = apply(mysamples,1,mean)
#           mysd = apply(mysamples,1,sd)
#           missingrep[,j] = rnorm(1,mymean,mysd)*missingref
#         }
#         #save
#         mylogexpr[,mylogexpr==0]=missingrep
#         normintensitiesimpall[i,]=mylogexpr
#       }
#     }
#   }
#   normintensitiesimpall = na.omit(normintensitiesimpall) #get rid of proteins that had missing vals within a run
#   normintensitiesimpall = 2^normintensitiesimpall #undo log transform
#   write.csv(normintensitiesimpall,'missing_value_imputation.csv')
# }
# else{
#   normintensitiesimpall=normintensities
# }
normintensitiesimpall=normintensities

#perform IRS
if (reps>1){
  message("Internal reference normalization...")
  refs = as.matrix(normintensitiesimpall[,grepl('Ref',colnames(normintensitiesimpall),ignore.case=TRUE)])
  if (numrefs>1){
    for (i in 1:reps){
      if (i==1){
        refscomb = data.frame(rowMeans(refs[,(numrefs*(i-1)+1):(numrefs*i)]))
      }else{
        refscomb = data.frame(refscomb,rowMeans(refs[,(numrefs*(i-1)+1):(numrefs*i)]))
      }
    }
  }else{
    refscomb=refs
  }
  refslog = log(refscomb+1)
  zerosums = rowSums(refscomb[,]==0)
  irsaverage <- rowSums(refslog)/(dim(refscomb)[2]-zerosums)
  irsaverage = exp(irsaverage)
  normfactorref = irsaverage/refscomb
  for (i in 1:dim(normfactorref)[2]){
    normfactorref[normfactorref[,i]=="Inf",i]=NA
  }
  for (i in 1:dim(refscomb)[2]){
    #get intensity values for this run
    myintensities = normintensitiesimpall[,colnames(normintensitiesimpall)%in%paste(metadata$name[metadata$run==i],"_",metadata$rep[metadata$run==i],sep="")]
    mynormfactorref = normfactorref[,i]
    #row normalize each protein against the reference
    for (j in 1:dim(myintensities)[2]){
      myintensities[,j] = ceiling(myintensities[,j]*mynormfactorref)
    }
    #make a table of normalized intensities
    if (i==1){
      normintensitiesIRS = myintensities
    }else{
      normintensitiesIRS = data.frame(normintensitiesIRS,myintensities)
    }
  }
  write.csv(normintensitiesIRS,'IRS_normalized_values.csv')
  finalimpintensitiesIRS=normintensitiesIRS
  # if (numrefs==0){
  #   refs = rowMeans(normintensitiesimpall)
  # } else{
  #   refs = as.matrix(normintensitiesimpall[,grepl('Ref',colnames(normintensitiesimpall),ignore.case=TRUE)])
  # }
  # if (numrefs>1){
  #   for (i in 1:reps){
  #     if (i==1){
  #       refscomb = data.frame(rowMeans(refs[,(numrefs*(i-1)+1):(numrefs*i)]))
  #     }else{
  #       refscomb = data.frame(refscomb,rowMeans(refs[,(numrefs*(i-1)+1):(numrefs*i)]))
  #     }
  #     irsaverage <- apply(refscomb, 1, function(x) exp(mean(log(x))))
  #     normfactorref = irsaverage/refscomb
  #   }
  # }else{
  #   refscomb=refs
  #   
  # }
  # irsaverage <- apply(refscomb, 1, function(x) exp(mean(log(x))))
  # normfactorref = irsaverage/refscomb
  #   for (i in 1:reps){
  #     #get intensity values for this run
  #     myints = normintensitiesimpall[,(1+(i-1)*plex):(i*plex)]
  #     mynormfactorref = normfactorref[,i]
  #     #row normalize each protein against the reference
  #     for (j in 1:dim(myints)[2]){
  #       myints[,j] = ceiling(myints[,j]*mynormfactorref)
  #     }
  #     #make a table of normalized intensities
  #     if (i==1){
  #       finalimpintensitiesIRS = myints
  #     }else{
  #       finalimpintensitiesIRS = data.frame(finalimpintensitiesIRS,myints)
  #     }
  #   }
  # write.csv(finalimpintensitiesIRS,'IRS_normalized_values.csv')
}else{
  finalimpintensitiesIRS=normintensitiesimpall
}

#QC plots
colors = c("blue","green","orange","yellow","red","purple","white","blue","green","orange","yellow","red","purple","white")
png(filename='boxplot_log2_norm.png',width=5000,height=2000,res=300)
invisible(b <- boxplot(log2(finalimpintensitiesIRS+1),col=colors[unlist(lapply(1:reps,function(x) rep(x,plex)))],ylab="log2(Intensity)",cex.axis=0.75,las=2))
print(b)
dev.off()

#clustering
h = hclust(dist(colMeans(na.omit(finalimpintensitiesIRS))))
png(filename='clusters.png',width=5000,height=2000,res=300)
plot(h)
dev.off()

#pca
if (numrefs>0){
  pcaresults = prcomp(t(na.omit(finalimpintensitiesIRS[,!grepl("Ref",colnames(finalimpintensitiesIRS),ignore.case=TRUE)])))
  #plot PCA
  g <- ggbiplot(pcaresults, groups = metadata$name[!grepl("Ref",metadata$name,ignore.case=TRUE)],
                labels=colnames(finalimpintensitiesIRS)[grep("Ref",colnames(finalimpintensitiesIRS),invert=TRUE,ignore.case=TRUE)],
                var.axes=FALSE,labels.size=3,ellipse=TRUE)
  if (length(unique(metadata$name))>2){
    g <- g+scale_color_brewer(palette="Set1")
  }
  g<- g+ theme(text = element_text(size=14))
  png(filename='PCA.png',width=5000,height=2000,res=300)
  print(g)
  dev.off()
}else{
  pcaresults = prcomp(t(finalimpintensitiesIRS))
  #plot PCA
  g <- ggbiplot(pcaresults, groups = metadata$name,
                labels=colnames(finalimpintensitiesIRS),
                var.axes=FALSE,labels.size=3,ellipse=TRUE)
 # if (length(unique(metadata$name))>2){
#    g <- g+scale_color_brewer(palette="Paired")
#  }
  g<- g+ theme(text = element_text(size=14))
  png(filename='PCA.png',width=4000,height=2000,res=300)
  print(g)
  dev.off()
}

#make list of all pairwise comparisons
# mysamples = unique(metadata$name)
# mysamples = mysamples[!grepl("Ref",mysamples,ignore.case=TRUE)]
# comps = rep("",choose(length(mysamples),2))
# currentrow = 1
# for (i in 1:(length(mysamples)-1)){
#   for (j in (i+1):length(mysamples)){
#     comps[currentrow] = paste(mysamples[i],"_vs_",mysamples[j],sep="")
#     currentrow = currentrow+1
#   }
# }

#read in list of comparisons
comps = read.xlsx(compsfile)

#perform PoissonSeq
message("Differential expression analysis...")
pseqdata = finalimpintensitiesIRS
newwb <- createWorkbook()
newwb2 <- createWorkbook()
for (i in 1:dim(comps)[1]){
  #get the intensities for this comparison
  sepcomps = strsplit(comps[i,1],"_vs_")
  intensities1 = pseqdata[,grepl(sepcomps[[1]][1],colnames(pseqdata))]
  intensities2 = pseqdata[,grepl(sepcomps[[1]][2],colnames(pseqdata))]
  #intensities1 = na.omit(intensities1)
  #intensities2 = na.omit(intensities2)
  #make indicator variable y
  y= c(rep(1,dim(intensities1)[2]),rep(2,dim(intensities2)[2]))
  #perform PSeq
  pdata = data.frame(intensities1,intensities2)
  pdata = na.omit(pdata)
  pseq<- PS.Main(dat=list(n=pdata,y=y,type="twoclass",pair=FALSE,gname=row.names(pdata)),para=list(ct.sum=0,ct.mean=0))
  #get the actual fc
  pseq = pseq[order(pseq$gname),]
  pdata = pdata[order(row.names(pdata)),]
  #make sure everything in pdata is in pseq
  pdata = pdata[row.names(pdata)%in%pseq$gname,]
  myFC = data.frame(rowMeans(pdata[,y==2])/rowMeans(pdata[,y==1]),row.names=row.names(pdata))
  pseq[,7]=log2(myFC)
  colnames(pseq)[7]="log2FC"
  pseqdata = pseqdata[order(row.names(pseqdata)),]
  nozeros = nozeros[order(row.names(nozeros)),]
  myresults = data.frame(pseq[,c(1:5,7)],nozeros[row.names(nozeros)%in%pseq$gname,1:dim(mydata)[2]],pseqdata[row.names(pseqdata)%in%pseq$gname,])
  #save
  mycomp = comps[i,1]
  if (nchar(mycomp)>31){
    mysheet = abbreviate(mycomp,minlength=31)
  }else{
    mysheet=comps[i,1]
  }
  addWorksheet(wb = newwb2, sheetName = mysheet, gridLines = TRUE)
  writeDataTable(wb=newwb2, sheet=mysheet,x=myresults,tableStyle="none",
                 rowNames=TRUE,withFilter=FALSE,
                 bandedRows=FALSE,bandedCols=FALSE)
  #make volcano plot
  signum = sum(pseq$pval<qval)
  if (stat=="q"){
    png(filename=paste(paste(comps[i,1],"_volcano_plot_",qval,".png",sep="")),width=2500,height=2000,res=300)
    e <- EnhancedVolcano(pseq,rownames(pseq),'log2FC','fdr',ylim=c(0,3),xlim=c(-3,3),pointSize=1,labSize=0,FCcutoff=log2(1.1),pCutoff=qval,
                         title=paste(comps[i,1],"(",sum(pseq$fdr<qval),")",sep=""),
                         col=c('grey30','grey60','royalblue','red2'),
                         legendLabels=c('FC<1.1, q>0.1','FC>1.1, q>0.1','FC<1.1, q<0.1','FC>1.1, q<0.1'),
                         legendLabSize=10, ylab = bquote(~-Log[10]~italic(q)))
    plot(e)
    dev.off() 
  }else{
    png(filename=paste(paste(comps[i,1],"_volcano_plot_",qval,".png",sep="")),width=2500,height=2000,res=300)
    e <- EnhancedVolcano(pseq,rownames(pseq),'log2FC','pval',ylim=c(0,3),xlim=c(-3,3),pointSize=1,labSize=0,FCcutoff=log2(1.1),pCutoff=qval,
                         title=paste(comps[i,1],"(",signum,")",sep=""),
                         col=c('grey30','grey60','royalblue','red2'),
                         legendLabels=c('FC<1.1, p>0.05','FC>1.1, p>0.05','FC<1.1, p<0.05','FC>1.1, p<0.05'),
                         legendLabSize=10, ylab = bquote(~-Log[10]~italic(p)))
    plot(e)
    dev.off()
  }
  #make qvalue histogram
  if (stat=="q"){
    png(filename=paste(paste(comps[i,1],"_qval_hist.png",sep="")),width=2000,height=2000,res=300)
    h <- hist(pseq$fdr,breaks=100)
    plot(h)
    dev.off()
  }else{
    png(filename=paste(paste(comps[i,1],"_pval_hist.png",sep="")),width=2000,height=2000,res=300)
    h <- hist(pseq$pval,breaks=100)
    plot(h)
    dev.off()
  }
  #get differentially expressed genes and save
  if (stat=="q"){
    mypros = pseq[pseq$fdr<qval,c(1:5,7)]
  }else{
    mypros = pseq[pseq$pval<qval,c(1:5,7)]
  }
  mypros = mypros[order(mypros$gname),]
  myresults = data.frame(mypros,nozeros[row.names(nozeros)%in%mypros$gname,1:dim(mydata)[2]],pseqdata[row.names(pseqdata)%in%mypros$gname,])
  #save
  addWorksheet(wb = newwb, sheetName = mysheet, gridLines = TRUE)
  writeDataTable(wb=newwb, sheet=mysheet,x=myresults,tableStyle="none",
                 rowNames=TRUE,withFilter=FALSE,
                 bandedRows=FALSE,bandedCols=FALSE)
}

#write workbook
if (stat=="q"){
  saveWorkbook(newwb, paste("Pseq_all_comps_q",qval,".xlsx",sep=""),overwrite=TRUE)
}else{
  saveWorkbook(newwb, paste("Pseq_all_comps_p",qval,".xlsx",sep=""),overwrite=TRUE)
}
saveWorkbook(newwb2, "Pseq_all_comps.xlsx",overwrite=TRUE)
message("Finished!")
}
