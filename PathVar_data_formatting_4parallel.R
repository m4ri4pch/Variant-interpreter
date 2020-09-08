# Data formatting for PathVar in Hipathia Web
# Script by: MPC
# Mon May  6 15:55:06 2019
# Updated Fri May 31 10:26:30 2019 for parallel processing by tissue

#My paths
path_to_GTEx<-"path_to_file" #indicate the path where gtex expression file is located

#Input parameters

Tissue<-"Blood" #input the tissue to assess
#geneList<-c("BRCA1", "FANCD2") #input the gene or genes to evaluate in the pathways, genes to KO with PathAct
geneList<-c("672", "2177") #input the gene or genes to be evaluated in entrez

#Data from GTEx

load(file = file.path(path_to_GTEx,"gtexInHipathia.rda"))
cat("Loading gtex data.\n")
load(file= file.path(path_to_GTEx,"tissue_gtex.rda"))
#updated Thu May 30 16:35:11 2019
load(file= file.path(path_to_GTEx,"allgenes.rda")) #Load Hipathia genes

#Select tissues

#Search for custom tissue 
if(Tissue=="custom"){
  stopifnot(length(list.files(pattern="custom"))==1)
  file=list.files(pattern="custom")
  custom_exp<-read.delim(file,stringsAsFactors = F,header=T,dec=".")
  #updated Thu May 30 16:35:11 2019
  custom_exp<-custom_exp[rownames(custom_exp)%in%allg,] #Filter genes in Hipathia
  gtex_Sel<-custom_exp
}else{
  gtex_Sel<-as.data.frame(gtex[,t==Tissue])
}

#Load functions:

#Function to perform in silico Knock outs. It returns a list with each element being a tissue and containing its design and expression matrix

autoKO<-function(x,genes){
  c<-x
  colnames(c)<-paste(colnames(x),"KO",sep="_")
  c[rownames(c)%in%genes,]<-0.0001
  r<-cbind(x,c) 
  d<-data.frame("type"=c(rep("normal",times=ncol(x)),rep("mutated",times=ncol(x))),stringsAsFactors = F,row.names = c(colnames(x),colnames(c)))
  return(list("design"=d,"expMatrix"=r))
  
}

#Function to save expression and design files from a list object

saveMatrix<-function(x){
  if(!dir.exists("./tissues")){
    dir.create("./tissues")
  }
  dir.create(file.path("./tissues",Tissue),showWarnings = F)
  cat("Files are saved here:\n")
    pd<-paste(file.path("./tissues",Tissue),"/design.txt",sep="")
  print(pd)
    write.table(x$design,file=pd,col.names = F,row.names = T,sep="\t", quote = F)
    pm<-paste(file.path("./tissues",Tissue),"/expression_matrix.txt",sep="")
  print(pm)
    write.table(x$expMatrix,file=pm,col.names = T,row.names = T,sep="\t",dec=".", quote=F)
    
}


#In silico Knock outs
gtexHipath<-autoKO(gtex_Sel,genes=geneList)


#Create tissue folder and save expression matrix and design as txt TAB sep files 

saveMatrix(gtexHipath)

#Remove non-necessary objects
#updated Thu May 30 16:35:11 2019
rm(list=setdiff(ls(),c("gtexHipath"))) #Add needed objects' names

####################################