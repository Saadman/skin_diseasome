library(readxl)

library(UpSetR)
library(ComplexUpset)
library(VennDiagram)
library(ggplot2)

# Run with working directory = repo root (e.g. Rscript scripts/Common_DEG_list.R from repo root)
DEG_XLSX <- "DEG/GEX_datasets/all_skin_disease_DEGs.xlsx"
OUT_DIR <- "DEG/GEX_datasets"
IMG_DIR <- "images/high_res"

skin_DEG<-excel_sheets(DEG_XLSX)

HS<-read_xlsx(DEG_XLSX,sheet=1)
PS<-read_xlsx(DEG_XLSX,sheet=2)
VL<-read_xlsx(DEG_XLSX,sheet=3)
AD<-read_xlsx(DEG_XLSX,sheet=4)
SLS<-read_xlsx(DEG_XLSX,sheet=5)
CSU<-read_xlsx(DEG_XLSX,sheet=6)

##Reduce(intersect,list(CSU$Gene.symbol,HS$Gene.symbol,PS$Gene.symbol,VL$Gene.symbol,AD$Gene.symbol))


PS_AD_HS_VL<-merge(PS_AD_HS,VL,by="Gene.symbol")
PS_AD_HS_CSU<-merge(PS_AD_HS,CSU,"Gene.symbol")
PS_AD_HS_VL_CSU<-merge(PS_AD_HS_VL,CSU,by="Gene.symbol")
write.csv(PS_AD_HS_VL,file.path(OUT_DIR,"83_common_DEG.csv"))
write.csv(PS_AD_HS_CSU,file.path(OUT_DIR,"106_common_DEG.csv"))
write.csv(PS_AD_HS_VL_CSU,file.path(OUT_DIR,"14_common_DEG.csv"))
#PS and HS
PS_HS<-merge(PS,HS,by="Gene.symbol")

colnames(PS_HS)<-sub(".y",".HS",gsub(".x",".PS",colnames(PS_HS)),fixed = TRUE)
##write.csv(PS_HS,file.path(OUT_DIR,"PS_HS_common_DEG.csv"))

#PS and AD
PS_AD<-merge(PS,AD,by="Gene.symbol")

colnames(PS_AD)<-gsub(".y",".AD",gsub(".x",".PS",colnames(PS_AD)),fixed = TRUE)
##write.csv(PS_AD,file.path(OUT_DIR,"PS_AD_common_DEG.csv"))

#PS and VL

PS_VL<-merge(PS,VL,by="Gene.symbol")

colnames(PS_VL)<-gsub(".y",".VL",gsub(".x",".PS",colnames(PS_VL)),fixed=TRUE)
##write.csv(PS_VL,file.path(OUT_DIR,"PS_VL_common_DEG.csv"))



##VL and HS
VL_HS<-merge(VL,HS,by="Gene.symbol")
colnames(VL_HS)<-gsub(".y",".HS",gsub(".x",".VL",colnames(VL_HS)),fixed=TRUE)
write.csv(VL_HS,file.path(OUT_DIR,"VL_HS_common_DEG.csv"))

#PS and CSU
PS_CSU<-merge(PS,CSU,by="Gene.symbol")
colnames(PS_CSU)<-gsub(".y",".CSU",gsub(".x",".PS",colnames(PS_CSU)),fixed=TRUE)
write.csv(PS_CSU,file.path(OUT_DIR,"PS_CSU_common_DEG.csv"))


#PS, AD, HS

PS_AD_HS<-merge(PS_AD,HS,by="Gene.symbol")
colnames(PS_AD_HS)[16:22]<-paste0(colnames(PS_AD_HS)[16:22],".HS")
##write.csv(PS_AD_HS,file.path(OUT_DIR,"PS_AD_HS_common_DEG.csv"))

##PS,VL,HS
PS_VL_HS<-merge(PS_VL,HS,by="Gene.symbol")
colnames(PS_VL_HS)[16:22]<-paste0(colnames(PS_VL_HS)[16:22],".HS")
##write.csv(PS_VL_HS,file.path(OUT_DIR,"PS_VL_HS_common_DEG.csv"))

length(Reduce(intersect,list(PS$Gene.symbol,AD$Gene.symbol,HS$Gene.symbol)))

#Upset plot
list<-list(SLS=SLS$SYMBOL,Chronic_Urticaria=CSU$Gene.symbol,Atopic_Dermatitis=AD$Gene.symbol,
           Vitiligo=VL$Gene.symbol,Psoriasis=PS$Gene.symbol,Hidradenitis=HS$Gene.symbol)
upset_list<-UpSetR::upset(fromList(list),order.by ="freq",nsets=6)


upset_list$New_data
#set_size(8, 3)

jpeg(filename=file.path(IMG_DIR,"upsetplot_DEG.jpg"))
ComplexUpset::upset(as.data.frame(upset_list$New_data),colnames(upset_list$New_data),mode="intersect",set_sizes=FALSE
                    #base_annotations = list(  annotate(
                    #  geom='text', x=Inf, y=Inf,
                     # label=paste('Total:', nrow(as.data.frame(upset_list$New_data))),
                      #vjust=, hjust=1
                   # ))
                   #,min_size=10
                    
                    )
dev.off()
venn.diagram(
  x = list(PS$Gene.symbol, HS$Gene.symbol, AD$Gene.symbol),
  category.names = c("PS" , "HS " , "AD"),
  filename = file.path(IMG_DIR,'14_venn_diagramm.tiff'),
  output=TRUE
)

