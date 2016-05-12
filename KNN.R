library(readr)
library(dplyr)
library(stringr)


###########data wrangling############

data<-read.table("expression_data.csv",sep="\t",header=TRUE)
data[is.na(data)]<-0

rownames(data)<-data[,1]
data<-data%>%select(-PlasmoDBID)

corMatrix<-cor(t(data),method="spearman")

pc<-prcomp(corMatrix)

pc$sdev[1:10]^2/sum(pc$sdev^2)


###we found 4 is good enough (explain 95% of the variance)

pcData<-pc$x[,1:4]
pcData<-cbind(pcData,rownames(pcData))

colnames(pcData)[5]<-"PlasmoDBID"

#write.table(pcData,"pca.csv",sep="\t",quote = FALSE)



goData<-read.table("gene_association.Pfalciparum.1.4.2016",sep="\t",header=FALSE)
goData<-goData%>%select(V2,V3,V5)

colnames(goData)<-c("PlasmoDBID","geneId","GOTermID")


idData<-read_tsv("ids.csv")

idData<-idData%>%select(c(1,2,3))

idData<-tbl_df(idData)


pcData<-tbl_df(as.data.frame(pcData))

npcData<-pcData%>%left_join(idData,by=c("PlasmoDBID"="PlasmoDBID (v6.3)"))
colnames(npcData)[6]<-"PlasmoDBID2"
npcData<-npcData%>%mutate(id=ifelse(is.na(npcData[,6]),PlasmoDBID,PlasmoDBID2))
npcData<-npcData%>%select(c(1,2,3,4,8))
colnames(npcData)[5]<-"PlasmoDBID"

label<-goData[,1]%>%as.character()%>%str_split("\\.",n=3)
labels<-list()
for(i in 1:length(label)){
  labels[i]<-label[[i]][1]
}
labels<-unlist(labels)
goData[,1]<-labels
goData2<-goData%>%left_join(npcData)

goData2<-goData2%>%filter(complete.cases(.))

#write.table(goData2[1:3],"plasmoid2GoId.csv",sep="\t",row.names = FALSE,quote = FALSE)

GOTermIDs<-goData2%>%select(GOTermID)%>%distinct()

geneData<-tbl_df(read.csv("Top100Bot100_All_DrugsV2.csv"))
s<-list()
for(i in 0:30){
  for(j in 101:200){
    s[i*200+j-100]<-i*200+j
  }
}

s<-s%>%unlist

geneData<-geneData[s,]

completeData<-geneData%>%filter(complete.cases(.))%>%filter(!grepl("unknown function",description))

geneData<-geneData%>%filter(!is.na(id))%>%filter(description=="conserved Plasmodium protein, unknown function")

completeData<-

dataWithoutGO<-geneData%>%left_join(goData,by=c("id"="PlasmoDBID"))

dataWithGO<-dataWithoutGO%>%filter(!is.na(GOTermID))
dataWithoutGO<-dataWithoutGO%>%filter(is.na(GOTermID))


smallGOData<-read.table("GO.csv",sep="\t",header = FALSE,stringsAsFactors = FALSE)
smallGOData<-smallGOData%>%as.list()%>%unlist()

newData<-goData2%>%filter(GOTermID %in% smallGOData)
newData<-newData%>%left_join(dataWithoutGO,by=c("PlasmoDBID"="id"))

newData<-newData%>%mutate(known=ifelse(!is.na(description),0,1))%>%select(-c(8,9,10,11))
colnames(newData)[2:3]<-c("geneId","GOTermID")

#write.table(newData,"total_data.csv",quote = FALSE,row.names = FALSE,sep = "\t")





######### build KNN model ###############

newData<-read.csv("total_data.csv",sep="\t")


knownData<-newData%>%filter(known==1)
unknownData<-newData%>%filter(known==0)
knownData[,4:7]<-lapply(knownData[,4:7],as.numeric)
unknownData[,4:7]<-lapply(unknownData[,4:7],as.numeric)


knownData<-knownData%>%group_by(GOTermID)%>%mutate(n=n())%>%ungroup()%>%filter(n>=10)
clusterID<-read.table("cluster.csv",sep="\t",header = TRUE,stringsAsFactors = FALSE)

knownData<-knownData%>%left_join(clusterID)


library(caret)
library(pROC)
#inTrain <- createDataPartition(y = knownData$cluster,p=0.8)$Resample

#knnFit<-knn3(GOTermID~PC1+PC2+PC3+PC4,knownData,k=5)
#knnPred<-predict(knnFit,newData=unknownData)

control<-trainControl(method='cv',number=10)
res<-train(as.factor(cluster)~PC1+PC2+PC3+PC4,
           data=knownData,
           method="knn",
           trControl=control,
           tuneLength=1,
           tuneGrid=data.frame(k=seq(2,100,5)),
           metric="Accuracy")
plot(res)

## select() interface:
## Objects in this package can be accessed using the select() interface
## from the AnnotationDbi package. See ?select for details.
## Bimap interface:
# Convert the object to a list
xx <- as.list(GO.db::GOBPCHILDREN)
# Remove GO IDs that do not have any children
xx <- xx[!is.na(xx)]
if(length(xx) > 0){
  # Get the parent GO IDs for the first elents of xx
  goids <- xx[[1]]
  # Find out the GO terms for the first parent goid
  GOID(GOTERM[[goids[1]]])
  Term(GOTERM[[goids[1]]])
  Synonym(GOTERM[[goids[1]]])
  Secondary(GOTERM[[goids[1]]])
  Definition(GOTERM[[goids[1]]])
  Ontology(GOTERM[[goids[1]]])
}




