library(disgenet2r)
library(readxl)
library(tidyverse)
disgenet_api_key <- Sys.getenv("DISGENET_API_KEY")
if (disgenet_api_key == "") {
  stop("Set DISGENET_API_KEY in your environment before running this script (e.g. in ~/.Renviron).")
}

Sys.setenv(DISGENET_API_KEY = disgenet_api_key)
data1 <- gene2evidence( gene = "LEPR",
                        vocabulary = "HGNC",
                        disease ="C3554225",
                        database = "ALL",
                        score =c(0.3,1))
data1
results <- extract(data1)

# Run with working directory = repo root (e.g. Rscript scripts/disgenet.R from repo root)
OUT_DIR <- "disgenet"

diseasome<-read_excel("DEG/GEX_datasets/Venn_allDisease_final.xlsx")

Psoriasis=diseasome$PS[!is.na(diseasome$PS)] #6522 unique genes
Atopic_Dermatitis=diseasome$AD[!is.na(diseasome$AD)] #1527 unique genes
Chronic_Urticaria=diseasome$CSU[!is.na(diseasome$CSU)] #976 unique genes
Vitiligo=diseasome$VL[!is.na(diseasome$VL)] #3040 unique genes
Hidradenitis=diseasome$HS[!is.na(diseasome$HS)]# 6801 unique genes
SLS=diseasome$SLS[!is.na(diseasome$SLS)] #33 unique genes
common_skin_genes<-Reduce(intersect,list(Psoriasis,Atopic_Dermatitis,Chronic_Urticaria,Vitiligo,Hidradenitis)) #SLS has nothing in common

Reduce(intersect,list(Psoriasis,Atopic_Dermatitis,Chronic_Urticaria,Vitiligo,Hidradenitis))

####Get common elements in each overlap of skin diseases##
AD_PS_HS<-Reduce(intersect,list(Vitiligo,Psoriasis))
overlap<-calculate.overlap(
x=list(
  "p"=Psoriasis,
  "a"=Atopic_Dermatitis,
  #"c"=Chronic_Urticaria,
  #"v"=Vitiligo,
   "h"=Hidradenitis
  #"s"=SLS
)
);
upsetP<-upset(fromList(x))

#####
skin_genes_diseases<- gene2disease(
  gene     = common_skin_genes,
 
  score =c(0.2, 1),
  verbose  = FALSE
)

plot( skin_genes_diseases,
      class = "Network",
      prop = 1)

plot( skin_genes_diseases,
      class  ="Heatmap",
      limit  = 20, nchars = 50 )

plot( skin_genes_diseases,
      class="DiseaseClass", nchars=60,limit  = 100)

#OMIM CODE
HS_OM<-142690

data9 <- disease2disease_by_gene(
  ##Hidradenitis,ICD, Vitligo,psoriasis,atopic dermatitis, chronic idiopathic/spontaneous urticaria,
  disease  = c("C0085160","C0162823","C0042900","C0033860", "C1853965", "C0578870"),
  database = "ALL", ndiseases = 5,verbose = FALSE )
d9_diseaseome<-plot( data9, 
      class = "Diseasome",
      layout="layout.lgl" ,
      prop = 5 )

plot( data9, 
      class = "Network",
      #layout="layout.lgl" ,
      prop = 0 )





##gene2disease for each skin disease##


##PS
top200_PS<-c(head(Psoriasis,100),tail(Psoriasis,100))
#top200_PS<-Psoriasis
skin_genes_PS<- gene2disease(
  gene     = top200_PS,
  
  database="CURATED",
  score =c(0.2, 1),
  verbose  = FALSE
)

PS_df<-skin_genes_PS@qresult
PS_df$from<-rep('PS',47)

##HS###
top200_HS<-c(head(Hidradenitis,100),tail(Hidradenitis,100))

skin_genes_HS<- gene2disease(
  gene     = top200_HS,
  
  database="CURATED",
  score =c(0.7, 1),
  verbose  = FALSE
)

HS_df<-skin_genes_HS@qresult
HS_df$from<-rep('HS',44)

##VL###

top200_VL<-c(head(Vitiligo,100),tail(Vitiligo,100))

skin_genes_VL<- gene2disease(
  gene     = top200_VL,
  
  database="CURATED",
  score =c(0.7, 1),
  verbose  = FALSE
)

VL_df<-skin_genes_VL@qresult
VL_df$from<-rep('VL',43)

##AD###

top200_AD<-c(head(Atopic_Dermatitis,100),tail(Atopic_Dermatitis,100))

skin_genes_AD<- gene2disease(
  gene     = top200_AD,
  
  database="CURATED",
  score =c(0.7, 1),
  verbose  = FALSE
)

AD_df<-skin_genes_AD@qresult
AD_df$from<-rep('AD',53)

#CSU##
top200_CSU<-c(head(Chronic_Urticaria,100),tail(Chronic_Urticaria,100))

skin_genes_CSU<- gene2disease(
  gene     = top200_CSU,
  
  database="CURATED",
  score =c(0.7, 1),
  verbose  = FALSE
)

CSU_df<-skin_genes_CSU@qresult
CSU_df$from<-rep('CSU',48)

##SLS##
skin_genes_SLS<- gene2disease(
  gene     = SLS,
  
  database="CURATED",
  score =c(0.7, 1),
  verbose  = FALSE
)
SLS_df<-skin_genes_SLS@qresult
SLS_df$from<-rep('SLS',6)

master_skin_df<-rbindlist(list(SLS_df,CSU_df,AD_df,PS_df,HS_df,VL_df),use.names = TRUE)
write.csv(master_skin_df,file.path(OUT_DIR,"skin_diseaseom_v1.csv"))


##Disease enrichment

#An FDR of 25% indicates that the result is likely to be valid 3 out of 4 times, which is reasonable in the setting of exploratory discovery where one is interested in finding candidate hypothesis to be further validated as a results of future research. Given the lack of coherence in most expression datasets and the relatively small number of gene sets being analyzed, using a more stringent FDR cutoff may lead you to overlook potentially significant results. 



##PS###
PS_enrich <-disease_enrichment( entities =Psoriasis, vocabulary = "HGNC",
                                 database = "CURATED" ) #6522 genes found
PS_df<-PS_enrich @qresult
# adding suffix to column names 
colnames(PS_df) <- paste("PS",colnames(PS_df),sep="_")
colnames(PS_df)[2]<-"Description"

PS_df$PS_from<-rep('Psoriasis',6047 )

PS_df_final<-PS_df[which(PS_df$PS_FDR<0.25),] #33 rows/diseases



##AD##
AD_enrich <-disease_enrichment( entities =Atopic_Dermatitis, vocabulary = "HGNC",
                                       database = "CURATED" ) #1527 genes captured

AD_df<-AD_enrich @qresult
# adding suffix to column names 
colnames(AD_df) <- paste("AD",colnames(AD_df),sep="_")
colnames(AD_df)[2]<-"Description"

AD_df$AD_from<-rep('Atopic Dermatitis',2822)

AD_df_final<-AD_df[which(AD_df$AD_FDR<0.25),] #145 rows/diseases


##HS##

HS_enrich <-disease_enrichment( entities =Hidradenitis, vocabulary = "HGNC",
                                database = "CURATED" ) #6801 genes captured


HS_df<-HS_enrich @qresult
# adding suffix to column names 
colnames(HS_df) <- paste("HS",colnames(HS_df),sep="_")
colnames(HS_df)[2]<-"Description"

HS_df$HS_from<-rep('Hidradenitis Suppurativa',6109)

HS_df_final<-HS_df[which(HS_df$HS_FDR<0.25),] #136 diseases/rows


##VL##

VL_enrich <-disease_enrichment( entities =Vitiligo, vocabulary = "HGNC",
                                database = "CURATED" ) #3040 genes captured

VL_df<-VL_enrich @qresult
# adding suffix to column names 
colnames(VL_df) <- paste("VL",colnames(VL_df),sep="_")
colnames(VL_df)[2]<-"Description"

VL_df$VL_from<-rep('Vitiligo',3705 )

VL_df_final<-VL_df[which(VL_df$VL_FDR<0.25),] #12 diseases/rows



##CSU##


CSU_enrich <-disease_enrichment( entities =Chronic_Urticaria, vocabulary = "HGNC",
                                database = "CURATED" ) #976 genes captured

CSU_df<-CSU_enrich @qresult
# adding suffix to column names 
colnames(CSU_df) <- paste("CSU",colnames(CSU_df),sep="_")
colnames(CSU_df)[2]<-"Description"

CSU_df$CSU_from<-rep('Chronic Urticaria',2190 )

CSU_df_final<-CSU_df[which(CSU_df$CSU_FDR<0.25),] # 749 diseases/rows


##SLS##

SLS_enrich <-disease_enrichment( entities =SLS, vocabulary = "HGNC",
                                 database = "CURATED" ) #33 genes captured

SLS_df<-SLS_enrich @qresult
# adding suffix to column names 
colnames(SLS_df) <- paste("SLS",colnames(SLS_df),sep="_")
colnames(SLS_df)[2]<-"Description"

SLS_df$SLS_from<-rep('SLS',306)

SLS_df_final<-SLS_df[which(SLS_df$SLS_FDR<0.25),] # 259 diseases/rows

##Merge all disease enrichment for FDR < 0.25

#put all data frames into list
enrich_fdr_list <- list(PS_df_final, AD_df_final, HS_df_final,VL_df_final,CSU_df_final,SLS_df_final)      

#merge all data frames together
merged_enrich_fdr<-Reduce(function(x, y) merge(x, y, all=TRUE), enrich_fdr_list)   #1030 diseases total   85 (no disease common across all conditions). all=TRUE is outer join meaning all disease across datasets are captured
write.csv(merged_enrich_fdr,file.path(OUT_DIR,"skin_disease_enrichment_FDR_131122.csv"))



##Merge all disease enrichment for FDR < 0.25

#put all data frames into list
enrich_all_list <- list(PS_df, AD_df, HS_df,VL_df,CSU_df,SLS_df)

#merge all data frames together
merged_enrich_all<-Reduce(function(x, y) merge(x, y, all=FALSE), enrich_all_list)   #242 disease 85 columns #IMPORTANT for publication##

write.csv(merged_enrich_all,file.path(OUT_DIR,"skin_disease_enrichment_all_partialJoin_131122.csv"))



library("SPARQL")
dis2com <- disease2compound(
  disease = c("C0020538"),
  database = "CURATED" )
dis2com
qr <- extract(dis2com)
head(qr, 3)


