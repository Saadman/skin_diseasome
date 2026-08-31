library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(readxl)
library(biomaRt)
library(ggplot2)
library(igraph)
data(gcSample)
library(UpSetR)
library(disgenet2r)
#str(gcSample) 
ck <- compareCluster(gcSample, fun =enrichKEGG)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID") #maps gene ids to gene symbols
head(ck) 
dotplot(ck)

cnetplot(ck)



#For all unique differentially expressed genes for each disease
diseasome<-read_excel("DEG/GEX_datasets/Venn_allDisease_final.xlsx")

PS_3<-unique(diseasome$PS[!is.na(diseasome$PS)])
AD_3<-unique(diseasome$AD[!is.na(diseasome$AD)])
CSU_3<-unique(diseasome$CSU[!is.na(diseasome$CSU)])
SLS_3<-unique(diseasome$SLS[!is.na(diseasome$SLS)])
three_common<-Reduce(intersect, list(PS_3,AD_3,CSU_3))
two_common<-Reduce(intersect, list(PS_3,AD_3))
SLS_out<-Reduce(intersect, list(three_common,SLS_3))
only_three<-three_common[!three_common%in%SLS_out]
only_two<-two_common[!two_common%in%c(SLS_3,CSU_3)]

only_three<-unname(mapIds(org.Hs.eg.db, only_three, 'ENTREZID', 'SYMBOL'))
only_two<-unname(mapIds(org.Hs.eg.db, only_two, 'ENTREZID', 'SYMBOL'))
dis_list<-list(
  Three_common_CSU_AD_PS=only_three
  #Two_common_AD_PS=only_two
  
)

#dk <- compareCluster(dis_list, fun ="enrichGO",pvalueCutoff=0.05,OrgDb='org.Hs.eg.db')
dk <- compareCluster(dis_list, fun ="enrichGO",OrgDb='org.Hs.eg.db',ont="BP") #biological process
write.csv(dk,"disease_GO_enrichment_table_BiologicalProcess_ThreeCommon.csv")
dotplot(dk)

PS_map<-mapIds(org.Hs.eg.db, diseasome$PS[!is.na(diseasome$PS)], 'ENTREZID', 'SYMBOL')
Psoriasis<-unname(PS_map)
Atopic_Dermatitis<-unname(mapIds(org.Hs.eg.db, diseasome$AD[!is.na(diseasome$AD)], 'ENTREZID', 'SYMBOL'))
Chronic_Urticaria<-unname(mapIds(org.Hs.eg.db, diseasome$CSU[!is.na(diseasome$CSU)], 'ENTREZID', 'SYMBOL'))
Vitiligo<-unname(mapIds(org.Hs.eg.db, diseasome$VL[!is.na(diseasome$VL)], 'ENTREZID', 'SYMBOL'))
Hidradenitis<-unname(mapIds(org.Hs.eg.db, diseasome$HS[!is.na(diseasome$HS)], 'ENTREZID', 'SYMBOL'))
SLS<-unname(mapIds(org.Hs.eg.db, diseasome$SLS[!is.na(diseasome$SLS)], 'ENTREZID', 'SYMBOL'))


##WRITE entrez ids to txt ##
write.table(Psoriasis[!is.na(Psoriasis)], file = "PS_final.txt", sep = "",row.names = FALSE,col.names = FALSE) #6066 out of 6523
write.table(Atopic_Dermatitis[!is.na(Atopic_Dermatitis)], file = "AD_final.txt", sep ="",row.names = FALSE,col.names = FALSE)#1456 out of 1515 
write.table(Chronic_Urticaria[!is.na(Chronic_Urticaria)], file = "CSU_final.txt", sep ="",row.names = FALSE,col.names = FALSE)#948 out of 977
write.table(Vitiligo[!is.na(Vitiligo)], file = "VL_final.txt", sep ="",row.names = FALSE,col.names = FALSE) #2918 out of  3041
write.table(Hidradenitis[!is.na(Hidradenitis)], file = "HS_final.txt", sep ="",row.names = FALSE,col.names = FALSE) # 6324 out of 6801
write.table(SLS[!is.na(SLS)], file = "SLS_final.txt", sep ="",row.names = FALSE,col.names = FALSE)  #34 out of 34
####
dis_list<-list(
  Psoriasis=Psoriasis,
  Atopic_Dermatitis=Atopic_Dermatitis,
  Chronic_Urticaria=Chronic_Urticaria,
  Vitiligo=Vitiligo,
  Hidradenitis=Hidradenitis,
  SLS=SLS
  
)



#dk <- compareCluster(dis_list, fun ="enrichGO",pvalueCutoff=0.05,OrgDb='org.Hs.eg.db')
dk <- compareCluster(dis_list, fun ="enrichGO",OrgDb='org.Hs.eg.db',ont="BP") #biological process



dk <- setReadable(dk, OrgDb = org.Hs.eg.db, keyType="ENTREZID") #maps gene ids to gene symbols
write.csv(dk,"disease_GO_enrichment_table_BiologicalProcess_final.csv	")

dotplot(dk)

##Upsetplot##
upset_list<-list(
  
  Psoriasis=diseasome$PS[!is.na(diseasome$PS)],
  Atopic_Dermatitis=diseasome$AD[!is.na(diseasome$AD)],
  Chronic_Urticaria=diseasome$CSU[!is.na(diseasome$CSU)],
  Vitiligo=diseasome$VL[!is.na(diseasome$VL)],
  Hidradenitis=diseasome$HS[!is.na(diseasome$HS)],
  SLS=diseasome$SLS[!is.na(diseasome$SLS)]

  
)

upset(fromList(upset_list),order.by = "freq",nsets=6)
####

cnetplot(dk)
cnetplot(dk,node_label="category", 
         cex_label_category = 1)







##For only genes unique to particular diseases (upstate plot one dot columns)








diseaseome_up<-read_excel("DEG/unique_disease_upsetplotGenes.xlsx")




Psoriasis<-unname(mapIds(org.Hs.eg.db, diseaseome_up$PS[!is.na(diseaseome_up$PS)], 'ENTREZID', 'SYMBOL'))
Atopic_Dermatitis<-unname(mapIds(org.Hs.eg.db, diseaseome_up$AD[!is.na(diseaseome_up$AD)], 'ENTREZID', 'SYMBOL'))
Chronic_Urticaria<-unname(mapIds(org.Hs.eg.db, diseaseome_up$CSU[!is.na(diseaseome_up$CSU)], 'ENTREZID', 'SYMBOL'))
SLS<-unname(mapIds(org.Hs.eg.db, diseaseome_up$SLS[!is.na(diseaseome_up$SLS)], 'ENTREZID', 'SYMBOL'))
Shared_Genes<-unname(mapIds(org.Hs.eg.db, diseaseome_up$Common[!is.na(diseaseome_up$Common)], 'ENTREZID', 'SYMBOL'))

dis_list<-list(
  Psoriasis=Psoriasis,
  Atopic_Dermatitis=Atopic_Dermatitis,
  Chronic_Urticaria=Chronic_Urticaria,
  SLS=SLS,
  Shared_Genes=Shared_Genes
  
  
)



#dk <- compareCluster(dis_list, fun ="enrichGO",pvalueCutoff=0.05,OrgDb='org.Hs.eg.db')
dk <- compareCluster(dis_list, fun ="enrichGO",OrgDb='org.Hs.eg.db',ont="BP") #biological process



dk <- setReadable(dk, OrgDb = org.Hs.eg.db, keyType="ENTREZID") #maps gene ids to gene symbols
write.csv(dk,"disease_GO_enrichment_table_BiologicalProcess_upsetPlot.csv	")

dotplot(dk)

cnetplot(dk)
cnetplot(dk,node_label="category", 
         cex_label_category = 1)


#DOSE for DO similiarity
a <- c("DOID:14095", "DOID:5844", "DOID:2044", "DOID:8432", "DOID:9146",
       "DOID:10588", "DOID:3209", "DOID:848", "DOID:3341", "DOID:252")
b <- c("DOID:9409", "DOID:2491", "DOID:4467", "DOID:3498", "DOID:11256")
doSim(a, a, measure="Wang")


#DOID from malacards: CSU=DOID:0080747, Psoriasis:DOID:8893, irritant contact dermatitis:DOID:2772,Atopic dermatitis:DOID:3310

dis_sem<-c("DOID:0080747","DOID:8893","DOID:2772","DOID:3310")

s<-doSim(dis_sem, dis_sem, measure="Wang") #no csu
doSim(dis_sem, dis_sem, measure="Resnik") #no csu
doSim(dis_sem, dis_sem, measure="Lin") 
colnames(s)<-c("CSU","PS","ICD","AD")
rownames(s)<-c("CSU","PS","ICD","AD")
simplot(s,
        color.low="white", color.high="red",
        labs=TRUE, digits=2, labs.size=5,
        font.size=14, xlab="", ylab="")



#DOSE for gene cluster similarity

PS_entrez<-unname(mapIds(org.Hs.eg.db, PS_3, 'ENTREZID', 'SYMBOL'))
fwrite(list(as.integer(PS_entrez)), file = "PS.txt")
AD_entrez<-unname(mapIds(org.Hs.eg.db, AD_3, 'ENTREZID', 'SYMBOL'))
fwrite(list(as.integer(AD_entrez)), file = "AD.txt")
CSU_entrez<-unname(mapIds(org.Hs.eg.db, CSU_3, 'ENTREZID', 'SYMBOL'))
fwrite(list(as.integer(CSU_entrez)), file = "CSU.txt")
SLS_entrez<-unname(mapIds(org.Hs.eg.db, SLS_3, 'ENTREZID', 'SYMBOL'))
fwrite(list(as.integer(SLS_entrez)), file = "SLS.txt")



clusters <- list(PS=PS_entrez, AD=AD_entrez, CSU=CSU_entrez,SLS=SLS_entrez)
DOSE::mclusterSim(clusters, measure="Wang", combine="BMA")

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
myAttributes_a<-c("entrezgene_id","hgnc_symbol")
#myValues_a<-HGNC_symbols$SYMBOL
#myValues_a<-names(PS_map[is.na(PS_map)])
myValues_a<- h
#myFilter_a<-"hgnc_symbol"
myFilter_a<-"hgnc_symbol"
## assemble and query the mart
geneInfo_eQTL_a<- getBM(attributes =  myAttributes_a, filters =  myFilter_a,
                        values =  myValues_a, mart = ensembl)



###diseasome network viz##

nodes <- read.csv("go_enrichment/nodes.csv", header=T, as.is=T)

links <- read.csv("go_enrichment/links.csv", header=T, as.is=T)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
V(net)$size<-(V(net)$size/max(V(net)$size)*100)
V(net)$size[3]<-80
V(net)$size[6]<-14
E(net)$width <- 1+E(net)$Network_separation_s_AB*5

colrs <- sample(categorical_pal(50),6)
#"#CC79A7" "#0072B2" "#009E73" "#E69F00" "#999999" "#D55E00"
V(net)$color <- colrs[V(net)$color]

l <- layout_in_circle(net)



plot(net, layout=l,vertex.label=NA)
legend(x=-1.5, y=-1.5, c("AD","CSU", "HS","PS","VL","SLS"), pch=21,
       
       col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
