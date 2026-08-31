# Run with working directory = repo root (e.g. Rscript scripts/diseasome_cyto_filePrep.R from repo root)
# Input is produced by scripts/disgenet.R (merged_enrich_fdr -> disgenet/skin_disease_enrichment_FDR_131122.csv)
FDR_diseases<-read.csv("disgenet/skin_disease_enrichment_FDR_131122.csv")



all_diseases<-c("Psoriasis","Atopic Dermatitis","Hidradenitis Suppurativa","Vitiligo","Chronic Urticaria","SLS")
abb_names<-c('PS','AD','HS','VL','CSU','SLS')

source<-c()
target<-c()
target_code<-c()
target_class<-c()
ratio<-c()
FDR<-c()

for (diseases in 1:length(all_diseases)){
  #[1] "Description"           "PS_disease_class_name" "PS_FDR"                "PS_Ratio"              "PS_ID"
  abb_disease<-paste0(abb_names[diseases],"_from")
  col_num<-which(colnames(FDR_diseases)==abb_disease)
  single_ds<-FDR_diseases[which(FDR_diseases[,col_num]==all_diseases[diseases]),c("Description",paste0(abb_names[diseases],"_disease_class_name"),paste0(abb_names[diseases],"_FDR"),paste0(abb_names[diseases],"_Ratio"),paste0(abb_names[diseases],"_ID"))]
  
  #print(dim(single_ds))
  for (i in 1:nrow(single_ds)){
        #strsplit(single_ds$PS_disease_class_name)
        
        class_types<-strsplit(single_ds[i,2],";")[[1]]
        
        #for (j in 1:length(class_types)){
          source<-c(source,all_diseases[diseases])
          target<-c(target,single_ds[i,1])
          target_code<-c(target_code,single_ds[i,5])
          #target_class<-c(target_class,class_types[j])
          target_class<-c(target_class,single_ds[i,2])
          ratio<-c(ratio,eval(parse(text = single_ds[i,4])))
          FDR<-c(FDR,single_ds[i,3])
          
       # }
        
       
  }
  

}

cyto_df<-data.frame(source,target,target_code,target_class,ratio,FDR)
cyto_df_sub<-cyto_df[which(cyto_df$FDR<0.20),]

#cyto_df$target<-make.unique(cyto_df$target, sep = ".")
write.table(cyto_df,"networksim/combinedDiseaseClass_cytoscape_disease_details.txt",sep='\t',row.names = FALSE,col.names = TRUE)
write.table(cyto_df_sub,"networksim/combinedDiseaseClass_cytoscape_disease_details_sub.txt",sep='\t',row.names = FALSE,col.names = TRUE)

