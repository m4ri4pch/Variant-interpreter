# Data formatting for PathVar in Hipathia Web
# Script by: MPC
# Mon May  6 15:55:06 2019

#My paths
path_to_GTEx<-"path_to_file" #indicate the path where gtex expression file is located

#Input parameters. This can be passed as parameters by an external config file

TissueList<-c("Blood","Blood Vessel") #input the tissue or tissues to assess
#geneList<-c("BRCA1", "FANCD2") #input the gene or genes to evaluate in the pathways, genes to KO with PathAct
geneList<-c("672", "2177") #input the gene or genes to be evaluated in entrez
custom<-T #indicates (TRUE or FALSE) whether we are including custom tissue expression

#Data from GTEx

load(file = file.path(path_to_GTEx,"gtex.rda"))
cat("Loading gtex data.\n")
load(file= file.path(path_to_GTEx,"tissue_gtex.rda"))

#Select tissues
load_pathways
gtex_Sel<-lapply(TissueList, function(x) as.data.frame(gtex[,t==x]))
names(gtex_Sel)<-TissueList

#Search for custom tissue 
if(custom){
  stopifnot(length(list.files(pattern="custom"))==1)
  file=list.files(pattern="custom")
  custom_exp<-read.delim(file,stringsAsFactors = F,header=T,dec=".")
  gtex_Sel$custom<-custom_exp
}

#Load functions:

#Function to perform in silico Knock outs. It returns a list with each element being a tissue and containing its design and expression matrix

autoKO<-function(x,genes){
  
  c<-x
  colnames(c)<-paste(colnames(x),"KO",sep="_")
  c[rownames(c)%in%genes,]<-c[rownames(c)%in%genes,]*0.01
  r<-cbind(x,c) 
  d<-data.frame("type"=c(rep("normal",times=ncol(x)),rep("mutated",times=ncol(x))),stringsAsFactors = F,row.names = c(colnames(x),colnames(c)))
  return(list("design"=d,"expMatrix"=r))

}

#Function to save expression and design files from a list object

saveMatrix<-function(x){
  if(!dir.exists("./tissues")){
    dir.create("./tissues")
  }
  tissues<-names(x)
  lapply(tissues,function(y) dir.create(file.path("./tissues",y)))
  cat("Files are saved here:\n")
  for(i in 1:length(names(x))){ x[[i]]$Tissue <- tissues[i] }
  aux<-function(y){
    
    pd<-paste(file.path("./tissues",y$Tissue),"/design.txt",sep="")
    print(pd)
    pm<-paste(file.path("./tissues",y$Tissue),"/expression_matrix.txt",sep="")
    print(pm)
    
    exp<-y$expMatrix
    print("exp")
    print(head(exp))
    write.table(y$design,file=pd,col.names = F,row.names = T,sep="\t", quote = F)
    write.table(y$expMatrix,file=pm,col.names = T,row.names = T,sep="\t",dec=".", quote=F)
    
    
    
  }
  lapply(x,aux)
}

#In silico Knock outs
gtexHipath<-lapply(gtex_Sel,autoKO,genes=geneList)


#Create tissue folder

saveMatrix(gtexHipath)

####################################