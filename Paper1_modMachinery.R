####load files, nanocount and protAtlas
library(tidyverse)
library(ggplot2)
library(vioplot)
library(dplyr)
library(viridis)
library(gplots)
library(pheatmap)
library(Metrics)
library(VariantAnnotation)
library(GenomicRanges)
library(VennDiagram)
library(pheatmap)
library(readxl)
library(ggseqlogo)
library(ComplexHeatmap)
library(Cairo)
library(rtracklayer)

###### NanoCount N1
ActivatedNano1<-read.csv("/Dropbox/NanoCountN1/Activated_DRS.gencodev43.nanocount.N1.tsv", sep = '\t',header = T)
NaiveNano1<-read.csv("/Dropbox/NanoCountN1/Naive_DRS.gencodev43.nanocount.N1.tsv", sep = '\t',header = T)
JurkatNano1<-read.csv("/Dropbox/NanoCountN1/Jurkat_DRS.gencodev43.nanocount.N1.tsv", sep = '\t',header = T)
ImmortalizedNano1<-read.csv("/Dropbox/nanocount_six_cells_N1/immortalized.gencodev43.N1.tsv",sep='\t', header=T) #all six cell lines merged together
#add gene name
if(T){
  JurkatNano1$gene<-""
  NaiveNano1$gene<-"" 
  ActivatedNano1$gene <-""
  ImmortalizedNano1$gene <-""
  for (i in c(1:nrow(JurkatNano1))) {
    JurkatNano1$gene[i] <- strsplit(JurkatNano1$transcript_name[i], "\\|") [[1]][6]
  }
  for (i in c(1:nrow(NaiveNano1))) {
    NaiveNano1$gene[i] <- strsplit(NaiveNano1$transcript_name[i], "\\|") [[1]][6]
  }
  for (i in c(1:nrow(ActivatedNano1))) {
    ActivatedNano1$gene[i] <- strsplit(ActivatedNano1$transcript_name[i], "\\|") [[1]][6]
  }
  for (i in c(1:nrow(ImmortalizedNano1))) {
    ImmortalizedNano1$gene[i] <- strsplit(ImmortalizedNano1$transcript_name[i], "\\|") [[1]][6]
  }

}

#Immortalized unmerged Cell Lines N1
A549Nano1<-read.csv("/Dropbox/six_cell_lines/nanocount_six_cells_N1/A549_Direct_merged.gencodev43.N1.tsv", sep = '\t',header = T)
HepG2Nano1<-read.csv("/Dropbox/six_cell_lines/nanocount_six_cells_N1/HepG2_Direct_merged.gencodev43.N1.tsv", sep = '\t',header = T)
HeLaNano1<-read.csv("/Dropbox/six_cell_lines/nanocount_six_cells_N1/HeLa_Direct_merged.gencodev43.N1.tsv", sep = '\t',header = T)
NTERANano1<-read.csv("/Dropbox/six_cell_lines/nanocount_six_cells_N1/NTERA_Direct_merged.gencodev43.N1.tsv", sep = '\t',header = T)
SHSY5YNano1<-read.csv("/Dropbox/six_cell_lines/nanocount_six_cells_N1/SHSY5Y_Direct_merged.gencodev43.N1.tsv", sep = '\t',header = T)
#add gene names
if(T){
  A549Nano1$gene<-""
  for (i in c(1:nrow(A549Nano1))) {
    A549Nano1$gene[i] <- strsplit(A549Nano1$transcript_name[i], "\\|") [[1]][6]
  }
  A549Nano1 <- A549Nano1 %>%
    group_by(gene) %>%
    summarise(tpm = sum(tpm, na.rm = TRUE))
  
  HepG2Nano1$gene<-""
  for (i in c(1:nrow(HepG2Nano1))) {
    HepG2Nano1$gene[i] <- strsplit(HepG2Nano1$transcript_name[i], "\\|") [[1]][6]
  }
  HepG2Nano1 <- HepG2Nano1 %>%
    group_by(gene) %>%
    summarise(tpm = sum(tpm, na.rm = TRUE))
  
  HeLaNano1$gene<-""
  for (i in c(1:nrow(HeLaNano1))) {
    HeLaNano1$gene[i] <- strsplit(HeLaNano1$transcript_name[i], "\\|") [[1]][6]
  }
  HeLaNano1 <- HeLaNano1 %>%
    group_by(gene) %>%
    summarise(tpm = sum(tpm, na.rm = TRUE))
  
  NTERANano1$gene<-""
  for (i in c(1:nrow(NTERANano1))) {
    NTERANano1$gene[i] <- strsplit(NTERANano1$transcript_name[i], "\\|") [[1]][6]
  }
  NTERANano1 <- NTERANano1 %>%
    group_by(gene) %>%
    summarise(tpm = sum(tpm, na.rm = TRUE))
 
  SHSY5YNano1$gene<-""
  for (i in c(1:nrow(SHSY5YNano1))) {
    SHSY5YNano1$gene[i] <- strsplit(SHSY5YNano1$transcript_name[i], "\\|") [[1]][6]
  }
  SHSY5YNano1 <- SHSY5YNano1 %>%
    group_by(gene) %>%
    summarise(tpm = sum(tpm, na.rm = TRUE)) 
}

#### Protein Atlas
#### TPM Illumina protein atlas primary
protAtlasPrimary<-read.csv("/Dropbox/proteinAtlas_rna_single_cell_type.tsv",sep = '\t',header = T)
### TPM Illumina protein Atlas immortalized
protAtlas<-read.csv("/Dropbox/proteinAtlas_rna_celline.tsv", sep = '\t',header = T)

#mod-enzymes and random enzymes definition
if(T){
  
  immortalized_cellines<-c("SK-OV-3","PANC-1","Raji", "SK-MEL-1","hTERT-HME1","Jurkat E6.1","CACO-2","SuSa","Hep-G2","A-549")
  primary_cellines<-c("Ovarian stromal cells","Ductal cells","B-cells","Melanocytes","Breast glandular cells","T-cells","Proximal enterocytes","Sertoli cells","Hepatocytes","Alveolar cells type 2")
  
  pus_enzymes <- c("PUS1", "PUSL1", "PUS3", "TRUB1", "TRUB2",
                   "DKC1", "PUS7", "PUS7L", "RPUSD1", "RPUSD2",
                   "RPUSD3", "RPUSD4", "PUS10")
  m6a_enzymes_writers <- c("METTL3","MTA70","METTL4","METTL5",
                   "METTL14","METTL16","METT10D",
                   "RBM15","OTT","RBM15B","OTT3",
                   "ZC3H13", "VIRMA","HAKAI","CBLL1","WTAP","PCIF1","CAPAM") 
  m6a_enzymes_erasers <- c('ALKBH5',"FTO")
  m6a_enzymes_readers <- c('EIF3A',"HNRNPA2B1","HNRNPC","YTHDC1","YTHDC2","YTHDC3","ELAVL1","IGF2BP1","IGF2BP2","IGF2BP3")
  m6a_enzymes<-c(m6a_enzymes_writers,m6a_enzymes_erasers)
  #https://www.nature.com/articles/s41392-020-00450-x  KIAA1429 is VIRMA, HAKAI is CBLL1, METTL16 is METT10D
  #METTL3 IS MTA70, PCIF1 is CAPAM, HAKAI is CBLL1, RBM15 is OTT, RBM15B is OTT3
  #motif DRACH (D=G/A/U, R=G/A, H=A/U/C)
  ai_enzymes <- c( "ADAR1","ADAR", "ADAR2", "ADARB1", "ADAR3", "ADARB2","RED2")
  dhu_enzymes <- c("DUS","DUS1","DUS2","DUS3","DUS4")
  
  psi_and_m6a_enzymes<-c(pus_enzymes,m6a_enzymes_writers,m6a_enzymes_erasers)
  all_mod_enzymes<-c(pus_enzymes,m6a_enzymes_writers,m6a_enzymes_erasers,ai_enzymes,dhu_enzymes) 
  
  all_mod_writers<-c(pus_enzymes,m6a_enzymes_writers,ai_enzymes) 
  
  
  housekeeping<-c("GAPDH","RPL37")
  
  all_mod_enzymes_housekeeping<-c(pus_enzymes,m6a_enzymes_writers,m6a_enzymes_erasers,ai_enzymes,housekeeping) 
  
  #random genes
  #find extremes: min and max TPM over all the protAtlas cell lines in the immortalized sample
  TPM_min_RNA_mod<-min(protAtlas$nTPM[which((protAtlas$Gene.name %in% all_mod_enzymes) & (protAtlas$Cell.line %in% immortalized_cellines) & protAtlas$nTPM>0  )] )
  TPM_max_RNA_mod<-max(protAtlas$nTPM[which((protAtlas$Gene.name %in% all_mod_enzymes) & (protAtlas$Cell.line %in% immortalized_cellines))] )
  
  common_genes<-intersect(protAtlas$Gene.name, protAtlasPrimary$Gene.name)
  Ovary<-merge(protAtlas[which(protAtlas$Cell.line == immortalized_cellines[1] & protAtlas$Gene.name %in% common_genes),],protAtlasPrimary[which(protAtlasPrimary$Cell.type == primary_cellines[1] & protAtlasPrimary$Gene.name %in% common_genes),], by="Gene.name") #immortalized-primary
  Ovary<- Ovary[apply(Ovary!=0, 1, all),]
  Ovary$FC_ovary<-log2(Ovary$nTPM.x)-log2(Ovary$nTPM.y) #imm/primary
  names(Ovary)[names(Ovary) == "nTPM.x"] <- "nTPM_OvaryImmortalized"
  names(Ovary)[names(Ovary) == "nTPM.y"] <- "nTPM_OvaryPrimary"
  ovary_25<-quantile(Ovary$diff, probs = 0.5)
  ovary_75<-quantile(Ovary$diff, probs = 0.95)
  Ovary_diff<-Ovary[which(Ovary$diff<ovary_25 | Ovary$diff>ovary_75),]

  Ductal<-merge(protAtlas[which(protAtlas$Cell.line == immortalized_cellines[2] & protAtlas$Gene.name %in% common_genes),],protAtlasPrimary[which(protAtlasPrimary$Cell.type == primary_cellines[2] & protAtlasPrimary$Gene.name %in% common_genes),], by="Gene.name") #immortalized-primary
  Ductal<- Ductal[apply(Ductal!=0, 1, all),]
  Ductal$FC_ductal<-log2(Ductal$nTPM.x)-log2(Ductal$nTPM.y) #imm/primary
  names(Ductal)[names(Ductal) == "nTPM.x"] <- "nTPM_ductalImmortalized"
  names(Ductal)[names(Ductal) == "nTPM.y"] <- "nTPM_ductalPrimary"
  d_25<-quantile(Ductal$diff, probs = 0.5)
  d_75<-quantile(Ductal$diff, probs = 0.95)
  Ductal_diff<-Ductal[which(Ductal$diff<d_25 | Ductal$diff>d_75),]
  
  Bcell<-merge(protAtlas[which(protAtlas$Cell.line == immortalized_cellines[3] & protAtlas$Gene.name %in% common_genes),],protAtlasPrimary[which(protAtlasPrimary$Cell.type == primary_cellines[3] & protAtlasPrimary$Gene.name %in% common_genes),], by="Gene.name") #immortalized-primary
  Bcell<- Bcell[apply(Bcell!=0, 1, all),]
  Bcell$FC_bcell<-log2(Bcell$nTPM.x)-log2(Bcell$nTPM.y) #imm/primary
  names(Bcell)[names(Bcell) == "nTPM.x"] <- "nTPM_BcellImmortalized"
  names(Bcell)[names(Bcell) == "nTPM.y"] <- "nTPM_BcellPrimary"
  b_25<-quantile(Bcell$diff, probs = 0.5)
  b_75<-quantile(Bcell$diff, probs = 0.95)
  Bcell_diff<-Bcell[which(Bcell$diff<b_25 | Bcell$diff>b_75),]
  
  Melan<-merge(protAtlas[which(protAtlas$Cell.line == immortalized_cellines[4] & protAtlas$Gene.name %in% common_genes),],protAtlasPrimary[which(protAtlasPrimary$Cell.type == primary_cellines[4] & protAtlasPrimary$Gene.name %in% common_genes),], by="Gene.name") #immortalized-primary
  Melan<- Melan[apply(Melan!=0, 1, all),]
  Melan$FC_melan<-log2(Melan$nTPM.x)-log2(Melan$nTPM.y) #imm/primary
  names(Melan)[names(Melan) == "nTPM.x"] <- "nTPM_MelanImmortalized"
  names(Melan)[names(Melan) == "nTPM.y"] <- "nTPM_MelanPrimary"
  e_25<-quantile(Melan$diff, probs = 0.5)
  e_75<-quantile(Melan$diff, probs = 0.95)
  Melan_diff<-Melan[which(Melan$diff<e_25 | Melan$diff>e_75),]
  
  Breast<-merge(protAtlas[which(protAtlas$Cell.line == immortalized_cellines[5] & protAtlas$Gene.name %in% common_genes),],protAtlasPrimary[which(protAtlasPrimary$Cell.type == primary_cellines[5] & protAtlasPrimary$Gene.name %in% common_genes),], by="Gene.name") #immortalized-primary
  Breast<- Breast[apply(Breast!=0, 1, all),]
  Breast$diff<-Breast$nTPM.x-Breast$nTPM.y
  Breast$FC_breast<-log2(Breast$nTPM.x)-log2(Breast$nTPM.y) #imm/primary
  names(Breast)[names(Breast) == "nTPM.x"] <- "nTPM_BreastImmortalized"
  names(Breast)[names(Breast) == "nTPM.y"] <- "nTPM_BreastPrimary"
  br_25<-quantile(Breast$diff, probs = 0.5)
  br_75<-quantile(Breast$diff, probs = 0.95)
  Breast_diff<-Breast[which(Breast$diff<br_25 | Breast$diff>br_75),]
 
  common_filt<-0 
  common_filt<-merge(Breast,Melan, by="Gene.name", all=TRUE)
  common_filt<-merge(common_filt, Bcell, by="Gene.name", all=TRUE)
  common_filt<-merge(common_filt,Ductal,by="Gene.name", all=TRUE)
  common_filt<-merge(common_filt,Ovary,by="Gene.name", all=TRUE)
  
  # Set the seed for reproducibility
  set.seed(123)
  # Generate 28 random numbers between 1 and length of the total protAtlas library
  n_genes<- 3*length(all_mod_enzymes)  #length(table(all_mod_enzymes))-1
  random_numbers <- sample(1:nrow(common_filt), n_genes, replace = FALSE)
  random_genes_df<-common_filt[random_numbers,]
  random_genes<-unique(JurkatNano1$gene[which(JurkatNano1$gene %in% random_genes_df$Gene.name)]) #keeping just those random genes we have in a549 cell line nanopore data
  set.seed(12)
  random_numbers<-sample(1:length(random_genes),length(all_mod_enzymes),replace=FALSE)
  random_genes<-random_genes[random_numbers]
  random_genes_df2<-random_genes_df[which(random_genes_df$Gene.name %in% random_genes),]
}

### plots fig1
if(T){
  
############################mod machinery with Illumina data ONLY from protAtlas
  cell_i<-6 #goes from 1 to 6. 1 to 5 are immortalized general cell lines, 6 is jurkat
  if(T){
    
    #plot trend for known cell lines
    primary_cell_line=primary_cellines[cell_i]
    cell_line=immortalized_cellines[cell_i]
    
    #keeping just the genes with tpm!=0 in both primary and immortalzied libraries.
    di<-protAtlas[which(protAtlas$Cell.line==cell_line),]
    dp<-protAtlasPrimary[which(protAtlasPrimary$Cell.type==primary_cell_line),]
    data<-merge(di,dp,by="Gene.name")
    data<- data[apply(data!=0, 1, all),]
    #using immortalized only as the random genes are coming from the immortalized set
    mean_imm_RME<-mean(log10(data$nTPM.x[which(data$Gene.name %in% all_mod_enzymes)]))
    mean_imm_nRME<-mean(log10(data$nTPM.x[which( !(data$Gene.name %in% all_mod_enzymes))]))
    mean_prim_RME<-mean(log10(data$nTPM.y[which(data$Gene.name %in% all_mod_enzymes)]))
    mean_prim_nRME<-mean(log10(data$nTPM.y[which( !(data$Gene.name %in% all_mod_enzymes))]))
    sd_imm_RME<-sd(log10(data$nTPM.x[which(data$Gene.name %in% all_mod_enzymes)]))
    sd_prim_RME<-sd(log10(data$nTPM.y[which(data$Gene.name %in% all_mod_enzymes)]))
    sd_imm_nRME<-sd(log10(data$nTPM.x[which( !(data$Gene.name %in% all_mod_enzymes))]))
    sd_prim_nRME<-sd(log10(data$nTPM.y[which( !(data$Gene.name %in% all_mod_enzymes))]))
        
    labels_genes<-c(all_mod_enzymes)
    pus_df<-data.frame(gene=labels_genes, "pa_immTPM"=0, "pa_primaryTPM"=0)
    for (i in c(1:nrow(pus_df))){
      pus_df$pa_immTPM[i]<-protAtlas$nTPM[which(protAtlas$Cell.line == cell_line & protAtlas$Gene.name==pus_df$gene[i])][1]
      pus_df$pa_primaryTPM[i]<-protAtlasPrimary$nTPM[which(protAtlasPrimary$Cell.type == primary_cell_line & protAtlasPrimary$Gene.name==pus_df$gene[i])][1]
    }
    pus_df<-na.omit(pus_df)
    pus_df_rna_mod<- pus_df[apply(pus_df!=0, 1, all),]
    pus_df_rna_mod$I_minus_P<-abs(log10(pus_df_rna_mod$pa_immTPM)-log10(pus_df_rna_mod$pa_primaryTPM))
 
    
    labels_genes<-c(random_genes,housekeeping)
    pus_df<-data.frame(gene=labels_genes, "pa_immTPM"=0, "pa_primaryTPM"=0)
    for (i in c(1:nrow(pus_df))){
      pus_df$pa_immTPM[i]<-protAtlas$nTPM[which(protAtlas$Cell.line == cell_line & protAtlas$Gene.name==pus_df$gene[i])][1]
      pus_df$pa_primaryTPM[i]<-protAtlasPrimary$nTPM[which(protAtlasPrimary$Cell.type == primary_cell_line & protAtlasPrimary$Gene.name==pus_df$gene[i])][1]
    }
    pus_df<-na.omit(pus_df)
    pus_df<- pus_df[apply(pus_df!=0, 1, all),]
    pus_df$I_minus_P<-abs(log10(pus_df$pa_immTPM)-log10(pus_df$pa_primaryTPM))
    pus_df_random<-pus_df[order(pus_df$I_minus_P,decreasing=FALSE),]
    
    pus_df<-rbind(pus_df_rna_mod,pus_df_random)

    
    #boxplot/violin
    if(T){
      
      file=paste0("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Fig1_violin_",primary_cellines[cell_i],".pdf")
      
      imm<-protAtlas[which(protAtlas$Cell.line==immortalized_cellines[cell_i]),]
      names(imm)[names(imm)=="nTPM"] <- "nTPM_immortalized"
      prim<-protAtlasPrimary[which(protAtlasPrimary$Cell.type == primary_cellines[cell_i]),]
      names(prim)[names(prim)=="nTPM"] <- "nTPM_primary"
      merged<-merge(imm, prim, by = "Gene.name")
      
      merged$nTPM_Immort_min_Primary<-(merged$nTPM_immortalized-merged$nTPM_primary)
      merged$log2FC<-log2((merged$nTPM_immortalized+0.1)/(merged$nTPM_primary+0.1))

      #plot
      if(T){
        # Save the plot as a PDF file
        pdf(file, width = 9, height = 7)
        
        vioplot(merged$nTPM_Immort_min_Primary,col=adjustcolor("grey",alpha=0.5),rectCol = "white",border="black",
                horizontal =  TRUE, lineCol = "white", plotCentre = "line", ylim=c(-50000,50000))
        
        stripchart(merged$nTPM_Immort_min_Primary, method="jitter", col="black", vertical=FALSE, pch=1, add=TRUE,ylim=c(-50000,50000))  
        
        par(new=TRUE)
        boxplot(merged$nTPM_Immort_min_Primary,
                main = "ProtAtlas TPMs",
                xlab = "TPM",
                ylab = "",
                col = "grey",
                border = "black",
                horizontal = TRUE,
                outline = TRUE,
                range=0,
                ylim=c(-50000,50000)
               
        )
        
        stripchart(merged$nTPM_Immort_min_Primary[which(merged$Gene.name %in% all_mod_enzymes)], method="jitter", col=adjustcolor("green3",alpha=0.8), vertical=FALSE, pch=16,ylim=c(-50000,50000),add=TRUE)
        
        dev.off()
        
        
      }
    }
    
    
    #median log2FC calculation
    fold_changes_pa<-data.frame("percentiles"=c("0","25","50","75","100"), "SK-OV-3"=0,
                                                   "PANC-1"=0,
                                                   "Raji"=0,
                                                   "SK-MEL-1"=0,
                                                   "hTERT-HME1"=0,
                                                   "Jurkat E6.1"=0,check.names = FALSE )
    
    for (cell_i in c(1:6)){
    
      primary_cell_line=primary_cellines[cell_i]
      cell_line=immortalized_cellines[cell_i]
      
      imm<-protAtlas[which(protAtlas$Cell.line==immortalized_cellines[cell_i]),]
      names(imm)[names(imm)=="nTPM"] <- "nTPM_immortalized"
      prim<-protAtlasPrimary[which(protAtlasPrimary$Cell.type == primary_cellines[cell_i]),]
      names(prim)[names(prim)=="nTPM"] <- "nTPM_primary"
      merged<-merge(imm, prim, by = "Gene.name")

      merged$nTPM_Immort_min_Primary<-(merged$nTPM_immortalized-merged$nTPM_primary)
      merged$log2FC<-log2((merged$nTPM_immortalized+0.1)/(merged$nTPM_primary+0.1))
        
      p0=quantile(merged$log2FC, probs = 0, na.rm = TRUE)
      p25=quantile(merged$log2FC, probs = 0.25, na.rm = TRUE)
      p50=quantile(merged$log2FC, probs = 0.5, na.rm = TRUE)
      p75=quantile(merged$log2FC, probs = 0.75, na.rm = TRUE)
      p100=quantile(merged$log2FC, probs = 1, na.rm = TRUE)
      
      p<-c(p0,p25,p50,p75,p100)
      fold_changes_pa[cell_line]<-p
      
    }
      
    fold_changes_Nano<-data.frame("percentiles"=c("0","25","50","75","100"),"Jurkat E6.1"=0,check.names = FALSE )
    imm<- JurkatNano1
    names(imm)[names(imm)=="tpm"] <- "TPM_immortalized"
    prim<-NaiveNano1
    names(prim)[names(prim)=="tpm"] <- "TPM_primary"
    merged<-merge(imm, prim, by = "gene")
    
    merged$log2FC<-log2((merged$TPM_immortalized+0.1)/(merged$TPM_primary+0.1))
    
    p0=quantile(merged$log2FC, probs = 0, na.rm = TRUE)
    p25=quantile(merged$log2FC, probs = 0.25, na.rm = TRUE)
    p50=quantile(merged$log2FC, probs = 0.5, na.rm = TRUE)
    p75=quantile(merged$log2FC, probs = 0.75, na.rm = TRUE)
    p100=quantile(merged$log2FC, probs = 1, na.rm = TRUE)
    
    p<-c(p0,p25,p50,p75,p100)
    fold_changes_Nano[cell_line]<-p
    
    
    data_matrix <- cbind(fold_changes_Nano$`Jurkat E6.1`, fold_changes_pa$`Jurkat E6.1`,fold_changes_pa$`PANC-1`,fold_changes_pa$Raji,fold_changes_pa$`SK-MEL-1`,fold_changes_pa$`hTERT-HME1`,fold_changes_pa$`SK-OV-3`)
    colnames(data_matrix) <- c("Jurkat vs T cells Nanopore","Jukat vs T cells ProteinAtlas", "PANC-1 vs Ductal cells", "Rajii vs B-cells", "SK-MEL-1 vs Melanocytes", "hTERT_HME1 vs Breast glandular cells","SK-OV-3 vs Ovary")
    rownames(data_matrix) <- fold_changes_pa$percentiles
    
    
    #plot percentiles
    if(T){
      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/heatmapModEnzymes_percentiles.pdf"
      min.x = 0
      max.x = 100
      min.y = 0
      max.y = 10
      line.thickness = 1
      label.size = 1
      text.size = 0.2
      point.size = 1
      tck.length = 0.01
      tick.thickness = 1
      transparency = 0.2
      pdf(file.name,width = 5,height = 6)
      par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2)
      
      mako_palette<-viridis(20)
      mako_palette[20] <- "grey"
      
      heatmap.2(data_matrix,
                Rowv = FALSE, #false,
                Colv = FALSE ,
                col = mako_palette,
                dendrogram = "none", #none,
                key = TRUE,
                key.title = "Color Key",
                density.info = "none",
                trace ="none",
                margins=c(8,4), #0 starts from the bottom right, the first value controls if it's up/down. The second value left/right
                cexRow=0.3,
                cexCol = 0.4,
                lhei = c(1,4),
                lwid = c(2,3)
      )
      
      dev.off()
    }
      
      
      
     }    
       
    

###########################big heatmap for all 5 cell lines and all enzymes and random  
  if(T){
    
    labels_genes<-c(pus_enzymes, m6a_enzymes_writers,m6a_enzymes_erasers, ai_enzymes, dhu_enzymes)
    pus_df_big<-data.frame(gene=labels_genes,"SK-OV-3"=0,
                                             "PANC-1"=0,
                                             "Raji"=0,
                                             "SK-MEL-1"=0,
                                             "hTERT-HME1"=0,
                                             "Jurkat E6.1"=0,check.names = FALSE )
    #number of cell lines we use from Human Protein Atlas
    for (c in c(1:6)){
     
      primary_cell_line=primary_cellines[c]
      cell_line=immortalized_cellines[c]

      pus_df<-data.frame(gene=labels_genes, "pa_immTPM"=0, "pa_primaryTPM"=0)
      
      for (i in c(1:nrow(pus_df))){
        pus_df$pa_immTPM[i]<-protAtlas$nTPM[which(protAtlas$Cell.line == cell_line & protAtlas$Gene.name==pus_df$gene[i])][1]
        pus_df$pa_primaryTPM[i]<-protAtlasPrimary$nTPM[which(protAtlasPrimary$Cell.type == primary_cell_line & protAtlasPrimary$Gene.name==pus_df$gene[i])][1]
      }
      pus_df$logFC<-log2(pus_df$pa_immTPM)-log2(pus_df$pa_primaryTPM)
      pus_df_big[cell_line]<-pus_df$logFC
    
      }

    pus_df_big$JurkatNano<-0
    for (i in c(1:nrow(pus_df_big))){
      pus_df_big$JurkatNano[i]<- log2(sum(JurkatNano1$tpm[which(JurkatNano1$gene==pus_df$gene[i])])) - log2(sum(NaiveNano1$tpm[which(NaiveNano1$gene==pus_df$gene[i])]))
    }
    pus_df_big_rnamod<-pus_df_big

    pus_df_big_rnamod$dissimilarity<-0
    for (g in c(1:nrow(pus_df_big_rnamod))){
      data<-c(pus_df_big_rnamod$`SK-OV-3`[g],pus_df_big_rnamod$`PANC-1`[g],pus_df_big_rnamod$Raji[g],pus_df_big_rnamod$`SK-MEL-1`[g],pus_df_big_rnamod$`hTERT-HME1`[g])
      pus_df_big_rnamod$dissimilarity[g]<- sd(data)
    }  
    pus_df_big_rnamod<-na.omit(pus_df_big_rnamod)
    
    
    labels_genes<-random_genes
    pus_df_big<-data.frame(gene=labels_genes, 
                           "SK-OV-3"=0,
                           "PANC-1"=0,
                           "Raji"=0,
                           "SK-MEL-1"=0,
                           "hTERT-HME1"=0, 
                           "Jurkat E6.1"=0,
                           check.names = FALSE )
    #number of cell lines we use from Human Protein Atlas
    for (c in c(1:6)){
      
      primary_cell_line=primary_cellines[c]
      cell_line=immortalized_cellines[c]
      
      pus_df<-data.frame(gene=labels_genes, "pa_immTPM"=0, "pa_primaryTPM"=0)
      
      for (i in c(1:nrow(pus_df))){
        pus_df$pa_immTPM[i]<-protAtlas$nTPM[which(protAtlas$Cell.line == cell_line & protAtlas$Gene.name==pus_df$gene[i])][1]
        pus_df$pa_primaryTPM[i]<-protAtlasPrimary$nTPM[which(protAtlasPrimary$Cell.type == primary_cell_line & protAtlasPrimary$Gene.name==pus_df$gene[i])][1]
      }
      
      pus_df$logFC<-log2(pus_df$pa_immTPM)-log2(pus_df$pa_primaryTPM)
      pus_df_big[cell_line]<-pus_df$logFC
      
    }
    
    pus_df_big$JurkatNano<-0
    for (i in c(1:nrow(pus_df_big))){
      pus_df_big$JurkatNano[i]<- log2(sum(JurkatNano1$tpm[which(JurkatNano1$gene==pus_df$gene[i])])) - log2(sum(NaiveNano1$tpm[which(NaiveNano1$gene==pus_df$gene[i])]))
    }

    pus_df_big_rand<-pus_df_big
  
    pus_df_big_rand<- pus_df_big_rand %>%
      mutate(across(everything(), ~ ifelse(is.na(.),Inf, .)))
    
    
    #bigger std means bigger variation from the mean value, so bigger dissimilarity 
    pus_df_big_rand$dissimilarity<-0
    for (g in c(1:nrow(pus_df_big_rand))){
      data<-c(pus_df_big_rand$`SK-OV-3`[g],pus_df_big_rand$`PANC-1`[g],pus_df_big_rand$Raji[g],pus_df_big_rand$`SK-MEL-1`[g],pus_df_big_rand$`hTERT-HME1`[g])
      pus_df_big_rand$dissimilarity[g]<- sd(data,na.rm = TRUE)
    }  
    pus_df_big_rand<-pus_df_big_rand[order(pus_df_big_rand$dissimilarity,decreasing=TRUE),]

    
    pus_df_big<-rbind(pus_df_big_rnamod,pus_df_big_rand)
    
    pus_df_big <- pus_df_big %>%
    mutate(across(everything(), ~ ifelse(is.infinite(.), 100, .)))
    pus_df_big<- pus_df_big %>%
    mutate(across(everything(), ~ ifelse(is.na(.), 100, .)))

    #percentiles 
    ############
    if(T){
      
      #median log2FC calculation
      fold_changes_pa<-data.frame("gene"=c("0","10","20","30","40","50","60","70","80","90","100"), "SK-OV-3"=0,
                                  "PANC-1"=0,
                                  "Raji"=0,
                                  "SK-MEL-1"=0,
                                  "hTERT-HME1"=0,
                                  "Jurkat E6.1"=0,check.names = FALSE )
      
      for (cell_i in c(1:6)){
        
        primary_cell_line=primary_cellines[cell_i]
        cell_line=immortalized_cellines[cell_i]
        
        imm<-protAtlas[which(protAtlas$Cell.line==immortalized_cellines[cell_i]),]
        names(imm)[names(imm)=="nTPM"] <- "nTPM_immortalized"
        prim<-protAtlasPrimary[which(protAtlasPrimary$Cell.type == primary_cellines[cell_i]),]
        names(prim)[names(prim)=="nTPM"] <- "nTPM_primary"
        merged<-merge(imm, prim, by = "Gene.name")
        # Remove rows that contain zeros
        merged <- merged[rowSums(merged == 0) == 0, ]
  
        merged$log2FC<-log2((merged$nTPM_immortalized)/(merged$nTPM_primary))
        
        p0=quantile(merged$log2FC,   probs = 0, na.rm = TRUE)
        p10=quantile(merged$log2FC,  probs = 0.10, na.rm = TRUE)
        p20=quantile(merged$log2FC,  probs = 0.20, na.rm = TRUE)
        p30=quantile(merged$log2FC,  probs = 0.30, na.rm = TRUE)
        p40=quantile(merged$log2FC,  probs = 0.40, na.rm = TRUE)
        p50=quantile(merged$log2FC,  probs = 0.5, na.rm = TRUE)
        p60=quantile(merged$log2FC,  probs = 0.60, na.rm = TRUE)
        p70=quantile(merged$log2FC,  probs = 0.70, na.rm = TRUE)
        p80=quantile(merged$log2FC,  probs = 0.80, na.rm = TRUE)
        p90=quantile(merged$log2FC,  probs = 0.90, na.rm = TRUE)
        p100=quantile(merged$log2FC, probs = 1, na.rm = TRUE)
        
        p<-c(p0,p10,p20,p30,p40,p50,p60,p70,p80,p90,p100)
        fold_changes_pa[cell_line]<-p
        
      }
      
      #median fold change calculation Nanopore
      fold_changes_Nano<-data.frame("gene"=c("0","10","20","30","40","50","60","70","80","90","100"),"JurkatNano"=0,check.names = FALSE )
      imm<- JurkatNano1
      names(imm)[names(imm)=="tpm"] <- "TPM_immortalized"
      prim<-NaiveNano1
      names(prim)[names(prim)=="tpm"] <- "TPM_primary"
      merged<-merge(imm, prim, by = "gene")
      merged <- merged[rowSums(merged == 0) == 0, ]
      merged$log2FC<-log2((merged$TPM_immortalized)/(merged$TPM_primary))
      
      p0=quantile(merged$log2FC, probs = 0, na.rm = TRUE)
      p10=quantile(merged$log2FC, probs = 0.10, na.rm = TRUE)
      p20=quantile(merged$log2FC, probs = 0.20, na.rm = TRUE)
      p30=quantile(merged$log2FC, probs = 0.30, na.rm = TRUE)
      p40=quantile(merged$log2FC, probs = 0.40, na.rm = TRUE)
      p50=quantile(merged$log2FC, probs = 0.5, na.rm = TRUE)
      p60=quantile(merged$log2FC, probs = 0.60, na.rm = TRUE)
      p70=quantile(merged$log2FC, probs = 0.70, na.rm = TRUE)
      p80=quantile(merged$log2FC, probs = 0.80, na.rm = TRUE)
      p90=quantile(merged$log2FC, probs = 0.90, na.rm = TRUE)
      p100=quantile(merged$log2FC, probs = 1, na.rm = TRUE)
      
      p<-c(p0,p10,p20,p30,p40,p50,p60,p70,p80,p90,p100)
      fold_changes_Nano<-p
      
      df <- cbind(fold_changes_pa, JurkatNano = fold_changes_Nano)
      df$dissimilarity<-0
      
      

    }

    pus_df_big_percentiles<-rbind(pus_df_big_rnamod,df)  
    data_matrix <- cbind( pus_df_big_percentiles$`Jurkat E6.1`,pus_df_big_percentiles$`PANC-1`,pus_df_big_percentiles$Raji,pus_df_big_percentiles$`SK-MEL-1`,pus_df_big_percentiles$`hTERT-HME1`,pus_df_big_percentiles$`SK-OV-3`)
    colnames(data_matrix) <- c("Jukat vs T cells ProteinAtlas", "PANC-1 vs Ductal cells", "Rajii vs B-cells", "SK-MEL-1 vs Melanocytes", "hTERT_HME1 vs Breast glandular cells","SK-OV-3 vs Ovary")
    rownames(data_matrix) <- pus_df_big_percentiles$gene


    #plot
    if(T){
      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/heatmapModEnzymes.pdf"
      min.x = 0
      max.x = 100
      min.y = 0
      max.y = 10
      line.thickness = 1
      label.size = 1
      text.size = 0.2
      point.size = 1
      tck.length = 0.01
      tick.thickness = 1
      transparency = 0.2
      pdf(file.name,width = 5,height = 6)
      par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2)
      
      mako_palette<-viridis(120)
      mako_palette[120] <- "grey"
      
      heatmap.2(data_matrix,
                Rowv = FALSE, #false,
                Colv = FALSE ,
                col = mako_palette,
                dendrogram = "none", #none,
                key = TRUE,
                key.title = "Color Key",
                density.info = "none",
                trace ="none",
                margins=c(8,4), #0 starts from the bottom right, the first value controls if it's up/down. The second value left/right
                cexRow=0.3,
                cexCol = 0.4,
                lhei = c(1,4),
                lwid = c(2,3)
      )
      
      dev.off()
    }
    
    
    pus_df_big_percentiles<-rbind(pus_df_big_rnamod,df)  
    data_matrix <- cbind(pus_df_big_percentiles$JurkatNano, pus_df_big_percentiles$`JurkatNano` )
    colnames(data_matrix) <- c("Jukat vs T cells Nanopore", "Jurkat vs T cells Nanopore")
    rownames(data_matrix) <- pus_df_big_percentiles$gene
    #plot
    if(T){
      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/heatmapModEnzymes_Nanopore.pdf"
      min.x = 0
      max.x = 100
      min.y = 0
      max.y = 10
      line.thickness = 1
      label.size = 1
      text.size = 0.2
      point.size = 1
      tck.length = 0.01
      tick.thickness = 1
      transparency = 0.2
      pdf(file.name,width = 5,height = 6)
      par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2)
      
      mako_palette<-mako(120)
      
      heatmap.2(data_matrix,
                Rowv = FALSE, #false,
                Colv = FALSE ,
                col = mako_palette,
                dendrogram = "none", #none,
                key = TRUE,
                key.title = "Color Key",
                density.info = "none",
                trace ="none",
                margins=c(8,4), #0 starts from the bottom right, the first value controls if it's up/down. The second value left/right
                cexRow=0.3,
                cexCol = 0.4,
                lhei = c(1,4),
                lwid = c(2,3)
      )
      
      dev.off()
    }
    
    
    
  }
  
######################compare prot Atlas Tcell-Jurkat vs nanopore Tcell-Jurkat 
  #abs plot mod machinery
  if(T){
    primary_cell_line="T-cells"
    cell_line="Jurkat E6.1"
    pus_df_rna_mod<-0
    pus_df_random<-0
    
    labels_genes<-c(all_mod_enzymes, random_genes, housekeeping)
    pus_df<-data.frame(gene=labels_genes, "paTPM"=0, "pa_primaryTPM"=0, "nanocountTPM"=0,"nanocount_primaryTPM"=0)
    for (i in c(1:nrow(pus_df))){
      pus_df$paTPM[i]<-protAtlas$nTPM[which(protAtlas$Cell.line == cell_line & protAtlas$Gene.name==pus_df$gene[i])][1]
      pus_df$pa_primaryTPM[i]<-protAtlasPrimary$nTPM[which(protAtlasPrimary$Cell.type == primary_cell_line & protAtlasPrimary$Gene.name==pus_df$gene[i])][1]
      pus_df$nanocountTPM[i]<-JurkatNano1$tpm[which(JurkatNano1$gene==pus_df$gene[i])][1]
      pus_df$nanocount_primaryTPM[i]<-sum(NaiveNano1$tpm[which(NaiveNano1$gene==pus_df$gene[i])])[1]
    }
    pus_df<-na.omit(pus_df)
    pus_df<-pus_df[apply(pus_df!=0, 1, all),]
   
    
    #plot absolute values
    if(T){
      min.x = 0
      max.x = length(pus_df$gene)
      min.y = -1
      max.y = 4
      line.thickness = 1
      label.size = 1
      text.size = 0.8
      point.size = 0.5
      tck.length = 0.01
      tick.thickness = 1
      transparency = 0.8
      hk="GAPDH"
      file.name = paste0("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Fig1_trend_ModRand_",cell_line,".pdf")
      pdf(file.name,width = 10,height = 5)
      par(mar = c(4, 4, 1, 1)) 
      plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
           ylab = NA,xlab=NA,ann = F,axes=F)
      
      space=0.2
      y_pa<-log10(pus_df$paTPM) #/(pus_df$paTPM[which(pus_df$gene==hk)])  #protAtlas Immortal
      x_pa<-c(1:length(y_pa))
      y_pa_primary<-log10(pus_df$pa_primaryTPM) #protAtlas Primary
      x_pa_primary<-c(1:length(y_pa_primary))
      for(i in c(1:length(y_our))){
        lines(x=c(x_pa[i],x_pa[i]-space), y=c(-1.1,min(y_pa[i],y_pa_primary[i])), col="grey90",lty=2)
      }
      
      points(x=x_pa-space, y=y_pa, pch=16, col=adjustcolor("blue2",alpha=transparency), cex=point.size)
      points(x=x_pa_primary-space, y=y_pa_primary, pch=16, col=adjustcolor("orange",alpha=transparency), cex=point.size)  
      for(i in c(1:length(y_our))){
        lines(x=c(x_pa_primary[i]-space,x_pa[i]-space),y=c(y_pa[i],y_pa_primary[i]), col="darkgrey")
      }
      
      y_our<- log10(pus_df$nanocountTPM) #/(pus_df$nanocountTPM[which(pus_df$gene==hk)])  #nanocount Immortal
      x_our<-c(1:length(y_our))
      y_our_primary<-log10(pus_df$nanocount_primaryTPM) #nanocount Primary
      x_our_primary<-c(1:length(y_our_primary))
      for(i in c(1:length(y_our))){
        lines(x=c(x_our[i],x_our[i]+space), y=c(-1.1,min(y_our[i],y_our_primary[i])), col="grey90",lty=2)
       }
      points(x=x_our+space, y=y_our, pch=16, col=adjustcolor("blue2",alpha=transparency), cex=point.size)
      points(x=x_our_primary+space, y=y_our_primary, pch=16, col=adjustcolor("orange",alpha=transparency), cex=point.size)
      for(i in c(1:length(y_our))){
        lines(x=c(x_our[i]+space,x_our_primary[i]+space),y=c(y_our[i],y_our_primary[i]), col="darkcyan")
      }

      axis(side = 2,at = c(min.y:max.y) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
           tck=-tck.length,las=1,lwd=line.thickness, line=0)
      axis(side = 2,at = c(min.y:max.y),labels = c(min.y:max.y),lwd.ticks = 0 ,cex.axis=label.size/2 ,tck= 0 ,las=1 ,lwd = 0, line= -0.5 )
      axis(side = 1,at = seq(0,length(pus_df$gene)+1,by=1) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
           tck=-tck.length,las=1,lwd=line.thickness, line=0)
      axis(side = 1,at = seq(1,length(pus_df$gene),by=1),labels = pus_df$gene ,lwd.ticks = 0 ,cex.axis=label.size/2 ,tck= 0 ,las=2 ,lwd = 0, line= -0.5 )
      legend("topleft", c(paste0("Immortalized ",cell_line),paste0("Primary ",primary_cell_line),"Protein Atlas" ,"Nanocount"),
             cex= text.size, pt.cex=point.size, y.intersp=1,x.intersp=0.01,
             pch=c(16,16),lty=c(0,0,1,1),
             lwd=c(2.5,2.5,2.5),col=c("blue2", "orange","darkgrey","darkcyan"), bty = "n", bg="white")
      mtext(text = "Abundance log10(TPMs)",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
      dev.off()
      
      
    }

  }
  
  #KDE/Histogram of delta for all transcripts in transcriptome Nanopore (Jurkat and Naive only)
  if(T){
    
    imm<-JurkatNano1
    imm <- imm %>%
      group_by(gene) %>%
      summarise(tpm = sum(tpm, na.rm = TRUE))
    names(imm)[names(imm)=="tpm"] <- "nTPM_immortalized"
    
    prim<-NaiveNano1
    prim <- prim %>%
      group_by(gene) %>%
      summarise(tpm = sum(tpm, na.rm = TRUE))
    names(prim)[names(prim)=="tpm"] <- "nTPM_primary"
    
    merged<-merge(imm, prim, by = "gene")
    
    
    pus_df<-data.frame(gene=merged$gene, nTPM_Immort_min_Primary=(merged$nTPM_immortalized-merged$nTPM_primary))
    
    
    if(T){
      min.y=0
      max.y=12000
      min.x=-300
      max.x=300
      file.name = paste0("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Fig1_Nanopore_delta_KDE_zoom_Tcell.pdf")
      pdf(file.name,width = 6,height = 5)
      data<-data.frame(x=(pus_df$nTPM_Immort_min_Primary), y=pus_df$gene)

      h<-hist(data$x[which(abs(data$x)<max.x)], main = paste0(immortalized_cellines[cell_i], ' vs ', primary_cellines[cell_i]), breaks=8,
              border="lightgrey",xlab = "TPMs (immortalized-primary) ", ylab = "Frequency", ylim=c(min.y,max.y), xlim=c(min.x,max.x))
      data_mod<-data[which(data$y %in% all_mod_enzymes),]
      
      #histogram countour
      lines(rep(h$breaks, each=2)[-c(1,2*length(h$breaks))], 
            rep(h$counts, each=2), lwd=1,col="darkgrey")
      
      #number of positive and negative mod-enzymes 
      l_pos=max(data_mod$x)
      h_pos=50*length(which(data_mod$x>=0))
      l_neg=min(data_mod$x)-5
      h_neg=50*length(which(data_mod$x<0))
      rect(xleft=l_neg, ybottom=max.y-500, xright=0, ytop=max.y-500+h_neg, density = NULL,
           col = adjustcolor("forestgreen",alpha=0.2), border = NA, lty = par("lty"), lwd = par("lwd"))
      rect(xleft=0, ybottom=max.y-500, xright=l_pos, ytop=max.y-500+h_pos, density = NULL,
           col = adjustcolor("forestgreen",alpha=0.2), border = NA, lty = par("lty"), lwd = par("lwd"))
      
      abline(h=max.y-500)
      
      for (j in c(1:length(data_mod$x))){
        points(x=data_mod$x[j],y=max.y-500, col="forestgreen", lw=1, cex=0.5)
        text(x=data_mod$x[j],y=max.y,labels=data_mod$y[j], srt=90,cex=0.3)
        lines(x=c(data_mod$x[j],data_mod$x[j]),y=c(min.y,max.y-500),lwd=0.5 ,lend=1, lty=2, col = adjustcolor("forestgreen",alpha=1)) 
        
      }
      
      dev.off()
      
    }
    
    if(T){
      
      file.name = paste0("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Fig1_Nanopore_delta_KDE_all_Tcell.pdf")
      pdf(file.name,width = 6,height = 5)
      data<-data.frame(x=(pus_df$nTPM_Immort_min_Primary), y=pus_df$gene )
      breaks_large <- seq(min(data$x), max(data$x), by = 1000) 
      data_range <- range(data$x)
      buffer <- 0.1 * diff(data_range)  # Adjust the buffer space as needed
      breaks_large <- seq(data_range[1] - buffer, data_range[2] + buffer, by = 800) 
      
      h<-hist(data$x, main = paste0(immortalized_cellines[cell_i], ' vs ', primary_cellines[cell_i]), breaks=200,
              xlab = "TPMs (immortalized-primary) ", ylab = "Frequency", ylim=c(min.y,max.y),border = "lightgrey")
      data_mod<-data[which(data$y %in% all_mod_enzymes),]
      
      lines(rep(h$breaks, each=2)[-c(1,2*length(h$breaks))], 
            rep(h$counts, each=2), lwd=1,col="darkgrey")
      dev.off()
      
      
    } 
    
    
  }
  
  #boxplot/violin for nanopore primary and immortalized T cells
  if (T){
   
    #nanocount
    if(T){
      file="/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Fig1_boxplot_Tcell_nanocount.pdf"
      
      imm<-JurkatNano1
      imm <- imm %>%
        group_by(gene) %>%
        summarise(tpm = sum(tpm, na.rm = TRUE))
      names(imm)[names(imm)=="tpm"] <- "TPM_immortalized"
      
      prim<-NaiveNano1
      prim <- prim %>%
        group_by(gene) %>%
        summarise(tpm = sum(tpm, na.rm = TRUE))
      names(prim)[names(prim)=="tpm"] <- "TPM_primary"
      
      merged<-merge(imm, prim, by = "gene")
      merged$TPM_Immort_min_Primary<-merged$TPM_immortalized-merged$TPM_primary
      
      data<-data.frame(merged)
      
      #plot
      if(T){
        
         # Save the plot as a PDF file
         pdf(file, width = 9, height = 7)

         vioplot(merged$TPM_Immort_min_Primary,col=adjustcolor("grey",alpha=0.5),rectCol = "white",border="black",
                horizontal =  TRUE, lineCol = "white", plotCentre = "line", ylim=c(-50000,50000))

         stripchart(merged$TPM_Immort_min_Primary, method="jitter", col="black", vertical=FALSE, pch=1, add=TRUE, ylim=c(-50000,50000))  
      
         axis(1)
         par(new=TRUE)
         boxplot(merged$TPM_Immort_min_Primary,
                      main = "Nanocount TPMs",
                      xlab = "TPM",
                      ylab = "",
                      col = "black",
                      border = "black",
                      horizontal = TRUE,
                      outline = TRUE,
                      range=0,  ylim=c(-50000,50000)
                      
         )
         stripchart(merged$TPM_Immort_min_Primary[which(merged$gene %in% all_mod_enzymes)],method="jitter", col=adjustcolor("green3",alpha=0.8), vertical=FALSE, pch=16, add=TRUE,  ylim=c(-50000,50000))
         
        
        dev.off()
        
        
        }
    }
    
    #illumina protein atlas
    if(T){
      file="/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Fig1_boxplot_Tcell_pA.pdf"
      
      imm<-protAtlas[which(protAtlas$Cell.line=="Jurkat E6.1"),]
      imm <- imm %>%
        group_by(Gene.name) %>%
        summarise(tpm = sum(nTPM, na.rm = TRUE))
      names(imm)[names(imm)=="tpm"] <- "TPM_immortalized"
      
      prim<-protAtlasPrimary[which(protAtlasPrimary$Cell.type=="T-cells"),]
      prim <- prim %>%
        group_by(Gene.name) %>%
        summarise(tpm = sum(nTPM, na.rm = TRUE))
      names(prim)[names(prim)=="tpm"] <- "TPM_primary"
      
      merged<-merge(imm, prim, by = "Gene.name")
      merged$TPM_Immort_min_Primary<-merged$TPM_immortalized-merged$TPM_primary
      

      #plot
      if(T){
        
        # Save the plot as a PDF file
        pdf(file, width = 9, height = 7)
        
        
        vioplot(merged$TPM_Immort_min_Primary,col=adjustcolor("grey",alpha=0.5),rectCol = "white",border="black", ylim=c(-50000,50000),
                horizontal =  TRUE, lineCol = "white", plotCentre = "line")
        
        stripchart(merged$TPM_Immort_min_Primary, method="jitter", col="black", vertical=FALSE, pch=1, add=TRUE, ylim=c(-50000,50000))  
        
        par(new=TRUE)
        boxplot(merged$TPM_Immort_min_Primary,
                main = "Nanocount TPMs",
                xlab = "TPM",
                ylab = "",
                col = "grey",
                border = "black",
                horizontal = TRUE,
                outline = TRUE,
                range=0,
                ylim=c(-50000,50000)
        )
        
        stripchart(merged$TPM_Immort_min_Primary[which(merged$Gene.name %in% all_mod_enzymes)], method="jitter", col=adjustcolor("green3",alpha=0.8), vertical=FALSE, pch=16, add=TRUE, ylim=c(-50000,50000))
        
        
        dev.off()
        
        
      }
    }
    
    #scatter plot RMSE between prot atlas and nanopore Jurkat-Naive data
    if(T){
      
      imm<-JurkatNano1
      imm <- imm %>%
        group_by(gene) %>%
        summarise(tpm = sum(tpm, na.rm = TRUE))
      names(imm)[names(imm)=="tpm"] <- "TPM_immortalized"
      
      prim<-NaiveNano1
      prim <- prim %>%
        group_by(gene) %>%
        summarise(tpm = sum(tpm, na.rm = TRUE))
      names(prim)[names(prim)=="tpm"] <- "TPM_primary"
      
      merged_nanopore<-merge(imm, prim, by = "gene")
      merged_nanopore$TPM_Immort_min_Primary_nano<- abs(merged_nanopore$TPM_immortalized-merged_nanopore$TPM_primary)
      
      
      imm<-protAtlas[which(protAtlas$Cell.line=="Jurkat E6.1"),]
      imm <- imm %>%
        group_by(Gene.name) %>%
        summarise(tpm = sum(nTPM, na.rm = TRUE))
      names(imm)[names(imm)=="tpm"] <- "TPM_immortalized"
      
      prim<-protAtlasPrimary[which(protAtlasPrimary$Cell.type=="T-cells"),]
      prim <- prim %>%
        group_by(Gene.name) %>%
        summarise(tpm = sum(nTPM, na.rm = TRUE))
      names(prim)[names(prim)=="tpm"] <- "TPM_primary"
      
      merged_protAtlas<-merge(imm, prim, by = "Gene.name")
      names(merged_protAtlas)[names(merged_protAtlas) == "Gene.name"]<-"gene"
      merged_protAtlas$TPM_Immort_min_Primary_pa<- abs(merged_protAtlas$TPM_immortalized-merged_protAtlas$TPM_primary)
      
      
      merged_df <- merge(merged_nanopore[, c("gene", "TPM_Immort_min_Primary_nano")], merged_protAtlas[, c("gene", "TPM_Immort_min_Primary_pa")], by = "gene")
    
      
      #scatter plot
      if(T){
        min.x = -6
        max.x = 16
        min.y = -6
        max.y = 16
        line.thickness = 1
        label.size = 1
        text.size = 1
        point.size = 0.5
        tck.length = 0.01
        tick.thickness = 1
        transparency = 1
        
        ##############create Figures folder first###############
        file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/correlationPA_Nano.pdf"
        pdf(file.name,width = 5, height = 6.6)
        par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
        plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
             ylab = NA,xlab=NA,ann = F,axes=F)
        
        mat<-merged_df

        points(x=(log2(mat$TPM_Immort_min_Primary_nano)), y=(log2((mat$TPM_Immort_min_Primary_pa)) ), pch=16, col=adjustcolor("gray13",alpha=transparency), cex=point.size-0.1)

        RMSE<-rmse(mat$TPM_Immort_min_Primary_nano,mat$TPM_Immort_min_Primary_pa)
        max_value <- max(c(mat$TPM_Immort_min_Primary_nano, mat$TPM_Immort_min_Primary_pa))
        min_value <- min(c(mat$TPM_Immort_min_Primary_nano, mat$TPM_Immort_min_Primary_pa))
        range_value <- max_value - min_value
        RMSE <- (RMSE / range_value) * 100
        
        
        correlation_coefficient <- cor((mat$TPM_Immort_min_Primary_nano), (mat$TPM_Immort_min_Primary_pa), method="pearson")
        
        # Print the correlation coefficient and R-squared value
        print(paste("RMSE N-J:", RMSE))
        
        # Add correlation line
        abline(lm( log2(mat$TPM_Immort_min_Primary_nano+1) ~ log2(mat$TPM_Immort_min_Primary_pa+1) ), col = "royalblue1", lty=2, lwd=2 )
        text(paste0("RMSE=",round(RMSE,digits=2)),x=2,y=13, col="black")
        
        
        axis(side = 2,at = seq(min.y,max.y,by=2) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
             tck=-tck.length,las=1,lwd=line.thickness, line=0)
        axis(side = 2,at = seq(min.y,max.y,by=2),labels = seq(min.y,max.y,by=2),lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line= -0.5 )
        axis(side = 1,at = seq(min.x,max.x,by=2)  ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
             tck=-tck.length,las=1,lwd=line.thickness, line=0)
        axis(side = 1,at = seq(min.x,max.x,by=2) ,labels = seq(min.x,max.x,by=2),lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line= -0.5 )
        mtext(text = "log2(TPM Jurkat-Naive) ProteinAtlas",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
        mtext(text = "log2(TPM Jurkat-Naive) Nanopore", side = 1, cex = text.size, line = 2.3, family = "Helvetica")
        
        
        dev.off()
        
        
      }
      
      
      
      
    }
    
  }
  
  
  #DELTA plot
  if(T){
    
    #enz<-c(all_mod_enzymes, random_genes,housekeeping) 
    enz<-c(all_mod_enzymes)
    pus_df<-data.frame("gene"=enz, "TPM_Jurkat_pA"=0, "TPM_Tcell_pA"=0, "TPM_Jurkat_nano"=0, "TPM_Naive_nano"=0)
    
    for (i in c(1:nrow(pus_df))) {
      pus_df$TPM_Tcell_pA[i]   <- protAtlasPrimary$nTPM[which(protAtlasPrimary$Gene.name == pus_df$gene[i] & protAtlasPrimary$Cell.type=="T-cells")][1]
      pus_df$TPM_Jurkat_pA[i]    <- protAtlas$nTPM[which(protAtlas$Gene.name == pus_df$gene[i] & protAtlas$Cell.line=="Jurkat E6.1")][1]
      pus_df$TPM_Jurkat_nano[i] <- sum(JurkatNano1$tpm[which(JurkatNano1$gene == pus_df$gene[i])])[1]
      pus_df$TPM_Naive_nano[i]  <- sum(NaiveNano1$tpm[which(NaiveNano1$gene == pus_df$gene[i])])[1]
    }
    
    pus_df <- na.omit(pus_df)
    
    pus_df$log_fc_nano<-log2(pus_df$TPM_Jurkat_nano)-log2(pus_df$TPM_Naive_nano)
    pus_df$log_fc_pa<-log2(pus_df$TPM_Jurkat_pA)-log2(pus_df$TPM_Tcell_pA)
    
    
    #plot
    if(T){
      min.x = 0
      max.x = length(pus_df$gene)
      min.y = -10
      max.y =  10
      line.thickness = 1
      label.size = 1
      text.size = 0.8
      point.size = 1
      tck.length = 0.01
      tick.thickness = 1
      transparency = 1
      ##############create Figures folder first###############
      file.name = paste0("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Fig1_deltas_JurkatTcell.pdf")
      pdf(file.name,width = 10,height = 6.6)
      plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
           ylab = NA,xlab=NA,ann = F,axes=F)
      labels_genes<-pus_df$gene

      for (i in c(1:nrow(pus_df))) {
        cnano<-0
        if(pus_df$TPM_Jurkat_nano[i]-pus_df$TPM_Naive_nano[i]<0){
          cnano<-15
        }
        cpA<-0
        if(pus_df$TPM_Jurkat_pA[i]-pus_df$TPM_Tcell_pA[i]<0){
          cpA<-15
        }
        lines(x=c(i+.1,i+.1),y=c(0,pus_df$log_fc_nano[i]),lwd=1.5 ,lend=1, lty=1, col = adjustcolor("darkcyan",alpha=0.8)) 
        lines(x=c(i+.1,i+.1),y=c(0,pus_df$log_fc_pa[i]),lwd=1.5 ,lend=1, lty=1, col = adjustcolor("grey30",alpha=0.8))
        points(x=i+0.1, y=(pus_df$log_fc_nano[i]), pch=cnano, col=adjustcolor("darkcyan",alpha=0.8), cex=point.size+0.2)
        points(x=i+0.1, y=(pus_df$log_fc_pa[i]), pch=cpA, col=adjustcolor("grey30",alpha=0.8), cex=point.size) #protAtlas
      }
      
      abline(h=0)
      axis(side = 2,at = seq(min.y,max.y,by=10) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
           tck=-tck.length,las=1,lwd=line.thickness, line=0)
      axis(side = 2,at = seq(min.y,max.y,by=10),labels = seq(min.y,max.y,by=10),lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line= -0.5 )
      axis(side = 1,at = c(0:length(labels_genes)+1) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
           tck=-tck.length,las=1,lwd=line.thickness, line=0)
      axis(side = 1,at = c(1:length(labels_genes)),labels = labels_genes ,lwd.ticks = 0 ,cex.axis=label.size/2 ,tck= 0 ,las=2 ,lwd = 0, line= -0.5 )
      legend("topright", c("protAtlas","nanocount"),
             cex= text.size, pt.cex=point.size, y.intersp=1,x.intersp=0.01,
             pch=c(15,15),lty=c(0,0),
             col=c("grey30","darkcyan"),
             lwd=c(1,1), bty = "n")
      mtext(text = "log2FC (Jurkat/Tcell)  ",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
   
      dev.off()
      
      
      
      
    }
    
    
    
  }
  
  #scatterplot all tpms + t cell markers
  if(T){
    #scatterplot all+tcell markers for tpm>2
    
    # Adding a prefix to TPM columns to distinguish cell lines
    NaiveNano1_copy<-NaiveNano1
    JurkatNano1_copy<-JurkatNano1
    ImmortalizedNano1_copy<-ImmortalizedNano1
    colnames(NaiveNano1_copy)[2:ncol(NaiveNano1_copy)] <- paste0("Naive_", colnames(NaiveNano1_copy)[2:ncol(NaiveNano1_copy)])
    colnames(JurkatNano1_copy)[2:ncol(JurkatNano1_copy)] <- paste0("Jurkat_", colnames(JurkatNano1_copy)[2:ncol(JurkatNano1_copy)])
    colnames(ImmortalizedNano1_copy)[2:ncol(ImmortalizedNano1_copy)] <- paste0("Immortalized_", colnames(ImmortalizedNano1_copy)[2:ncol(ImmortalizedNano1_copy)])
    
    merged_NaiveJurkatNanocountN1<-0
    merged_NaiveJurkatNanocountN1 <- merge(NaiveNano1_copy, JurkatNano1_copy, by = "transcript_name", all = TRUE)
    merged_NaiveJurkatNanocountN1$gene <-""
    for (i in c(1:nrow(merged_NaiveJurkatNanocountN1))) {
      merged_NaiveJurkatNanocountN1$gene[i] <- strsplit(merged_NaiveJurkatNanocountN1$transcript_name[i], "\\|") [[1]][6]
    }
    
    
    merged_NaiveJurkatNanocountN1_noNA<-data.frame(gene=merged_NaiveJurkatNanocountN1$gene,Naive_tpm=merged_NaiveJurkatNanocountN1$Naive_tpm,Jurkat_tpm=merged_NaiveJurkatNanocountN1$Jurkat_tpm)
    # Subset dataframe to exclude rows with NA values in Naive or Jurkat TPM columns
    merged_NaiveJurkatNanocountN1_noNA <- merged_NaiveJurkatNanocountN1_noNA[complete.cases(merged_NaiveJurkatNanocountN1_noNA[, c("Naive_tpm", "Jurkat_tpm")]), ]
    #sum all isoforms tpm together in order to retain gene information only
    merged_NaiveJurkatNanocountN1_noNA <- aggregate(. ~ gene, data = merged_NaiveJurkatNanocountN1_noNA, FUN = sum, na.rm = TRUE)
    
    
    up <- c("CD3","CD4", "CD8","CCR4","CD194","CCR5","CD195","CCR6","CD196","CCR7","CD197","CCR10","CD127","CD27","CD28","CD38","CD45RA","CD45RO","CD45", "CD58","LFA3","CD69",
            "CTLA4","CD152","CXCR3","CD183","FAS","CD95","HLA-DR","IL2RA","IL2RB","ITGAE","ITGAL","KLRB1","NCAM1","PECAM1","PTGDR2","SELL","IFNG","IL10",
            "IL13","IL17A","IL2","IL21","IL22","IL25","IL26","IL4","IL5","IL9","TGFB1","TNF","AHR",'EOMES',"FOXO4",'FOXP1',"FOXP3",
            "GATA3",'IRF4',"LEF1",'PRDM1',"RORC",'STAT4',"TBX21","TCF7","GZMA","CD25","CD122","CD103",'CD11a','CD161',"CD56",'CD31',
            "CD294","CRTH2","CD26L",'TGFB',"TNF",'MUM1',"BLIMP1","RORY")
    merged_NaiveJurkatNanocountN1_noNA$marker<-"black"
    merged_NaiveJurkatNanocountN1_noNA$marker[which(merged_NaiveJurkatNanocountN1_noNA$gene %in% up)]<-"yellow"
    
    mat_filt<-merged_NaiveJurkatNanocountN1_noNA[which(merged_NaiveJurkatNanocountN1_noNA$Naive_tpm>2 & merged_NaiveJurkatNanocountN1_noNA$Jurkat_tpm >2),]
    rownames(mat_filt)<-NULL
    rows_to_remove <- grepl("^RP", mat_filt$gene)
  

    #scatter plot
    if(T){
      min.x = round(log2(2))
      max.x = 15
      min.y = round(log2(2))
      max.y = 15
      line.thickness = 1
      label.size = 1
      text.size = 1
      point.size = 0.5
      tck.length = 0.01
      tick.thickness = 1
      transparency = 1

      ##############create Figures folder first###############
      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Fig1_scatterPlot_NanoCountN1.pdf"
      pdf(file.name,width = 5, height = 6.6)
      par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
      plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
           ylab = NA,xlab=NA,ann = F,axes=F)
      
      mat<-mat_filt
      #annotate abundant transcripts
      not_abundant<-mat[which(log2(mat$Naive_tpm)<10 & log2(mat$Jurkat_tpm)<10),]
  
      mat$marker[which(log2(mat$Naive_tpm)>=10 & log2(mat$Jurkat_tpm)>=10)]<-"grey"
      #text(abundant$gene,x=log2(abundant$Naive_tpm)+0.    if(T){
      
      points(x=log2(mat$Naive_tpm), y=log2(mat$Jurkat_tpm), pch=16, col=adjustcolor("gray13",alpha=transparency), cex=point.size-0.1)
      points(x=log2(mat$Naive_tpm[which(mat$marker=="yellow")]), y=log2(mat$Jurkat_tpm[which(mat$marker=="yellow")]), pch=16, col=adjustcolor("red2",alpha=transparency), cex=point.size*2)
      points(x=log2(mat$Naive_tpm[which(mat$marker=="grey")]), y=log2(mat$Jurkat_tpm[which(mat$marker=="grey")]), pch=16, col=adjustcolor("grey",alpha=transparency), cex=point.size)
      
      mat_noRP<-not_abundant
      # Calculate Pearson correlation coefficient for log-transformed data
      correlation_coefficient <- cor((mat_noRP$Naive_tpm), (mat_noRP$Jurkat_tpm))
      RMSE<-rmse(mat_noRP$Naive_tpm,mat_noRP$Jurkat_tpm)
      max_value <- max(c(mat_noRP$Naive_tpm,mat_noRP$Jurkat_tpm))
      min_value <- min(c(mat_noRP$Naive_tpm,mat_noRP$Jurkat_tpm))
      range_value <- max_value - min_value
      RMSE <- (RMSE / range_value) * 100
      
      
      # Calculate R-squared value
      r_squared <- correlation_coefficient^2
      

      # Print the correlation coefficient and R-squared value
      print(paste("Correlation Coefficient N-J:", correlation_coefficient))
      print(paste("RMSE N-J:", RMSE))
      print(paste("R-squared Value:", r_squared))
      
      # Add correlation line
      abline(lm(log2(mat_noRP$Jurkat_tpm) ~ log2(mat_noRP$Naive_tpm)), col = "royalblue1", lty=2, lwd=2 )
      text(paste0("RMSE=",round(RMSE,digits=2)),x=2,y=13, col="black")
      
      
      axis(side = 2,at = seq(min.y,max.y,by=2) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
           tck=-tck.length,las=1,lwd=line.thickness, line=0)
      axis(side = 2,at = seq(min.y,max.y,by=2),labels = seq(min.y,max.y,by=2),lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line= -0.5 )
      axis(side = 1,at = seq(min.x,max.x,by=2)  ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
           tck=-tck.length,las=1,lwd=line.thickness, line=0)
      axis(side = 1,at = seq(min.x,max.x,by=2) ,labels = seq(min.x,max.x,by=2),lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line= -0.5 )
      legend("topright", c("T-cell markers"),
             cex= text.size, pt.cex=point.size*2, y.intersp=1,x.intersp=0.01,
             pch=c(16,16,16),lty=c(0,0,0),
             lwd=c(2.5,2.5,2.5),col=c("red"), bty = "n")
      mtext(text = "log2(TPMs Jurkat)",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
      mtext(text = "log2(TPMs Naive)", side = 1, cex = text.size, line = 2.3, family = "Helvetica")
      
      
      dev.off()
      
      
  }
  
  }
  
  #heatmap Jurkat Naive PA and Nanopore 
  if(T){
    
    all_enzymes<-c(pus_enzymes, m6a_enzymes_writers, m6a_enzymes_erasers,ai_enzymes)
    syn_enz=all_enzymes
    pus_df<- data.frame("gene" = syn_enz,
                        "NaiveTPM"=0,
                        "JurkatTPM"=0,
                        "NaivePA"=0,
                        "JurkatPA"=0)
    #Add tpms
    pus_df$NaiveTPM<-0
    pus_df$JurkatTPM<-0
    
    for (i in c(1:nrow(pus_df))){
      pus_df$NaiveTPM[i]<-sum(NaiveNano1$tpm[which(NaiveNano1$gene==pus_df$gene[i])])
      pus_df$JurkatTPM[i]<-sum(JurkatNano1$tpm[which(JurkatNano1$gene==pus_df$gene[i])])
      pus_df$NaivePA[i]<-protAtlasPrimary$nTPM[which(protAtlasPrimary$Cell.type=="T-cells" & protAtlasPrimary$Gene.name==pus_df$gene[i])][1]
      pus_df$JurkatPA[i]<-protAtlas$nTPM[which(protAtlas$Cell.line=="Jurkat E6.1" & protAtlas$Gene.name==pus_df$gene[i])][1]
    }
    
    pus_df<-pus_df[apply(pus_df!=0, 1,all),]
    pus_df$similarity_JN_HPA<- pus_df$JurkatPA-pus_df$NaivePA
    pus_df$similarity_JN_Nano<- pus_df$JurkatTPM-pus_df$NaiveTPM

    #heatmap
    if(T){
    
    colormy<-mako(250)
    colormy<-rev(colormy) 
    pdf("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/HeatmapNaiveJurkatImmortalized_machinery_Nanocount.pdf", width = 8, height = 10)
    pheatmap(pus_df[,(2:5)], cluster_rows = F, cluster_cols = F, color = colormy,
             labels_col = colnames(pus_df[,(2:5)]),name = "TPM", show_rownames = TRUE, labels_row = pus_df$gene
    )  
    dev.off()
  }

    #table pnas
    if(T){ 
      pnas_filter=pus_df
      min.x = -50
      max.x = 70
      min.y = -10
      max.y = nrow(pnas_filter)+ 5
      line.thickness = 0.8
      txt.size = 0.2
      point.size = 0.3
      label.size = 0.5
      anot.pos = -20
      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Fig2_heatmapSimilarityNJ.pdf"
      pdf(file=file.name)
      
      par(mfrow=c(1,1),lend=1)
      plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
           ylab = NA,xlab=NA,ann = F,axes=F)
      
      for (i in c(1:nrow(pnas_filter))) {
        ### Add the similarity PA
        text(x = anot.pos ,y = nrow(pnas_filter)-i, round(pnas_filter$similarity_JN_HPA[i],digits = 1) ,cex = 1.2*txt.size, pos = 2, srt = 0)

        ### Add the similarity Nanopore
        text(x = anot.pos+20 ,y = nrow(pnas_filter)-i,
             labels = round(pnas_filter$similarity_JN_Nano[i],digits = 1) ,cex = 1.2*txt.size, pos = 2, srt = 0) 
      }
      
      ### Add the axis
      text(x = c(anot.pos,anot.pos+10), y = (nrow(pnas_filter))*c(1,1), labels = c("HPA","Nanopore"),
           srt = 0,cex = 1.2*txt.size, pos = 4, col = "black",family = "Helvetica", font = 3)
      
      dev.off()
    }
    
    
}
    
}  

### figures with psi sites analysis, fig 2 

#panIVT enrichment, the filtered and enriched panIVT libraries will be written in Jurkat_Naive_filtered_pan_noPvalue.csv
if(T){

  Naive_pairedIVT <- read.csv("/Dropbox/merged_pvals_Naive.csv", header = T)
  Naive_pairedIVT$ACP<-paste0( Naive_pairedIVT$Annotation, Naive_pairedIVT$chr, Naive_pairedIVT$position)
  Jurkat_pairedIVT <- read.csv("/Dropbox/merged_pvals_Jurkat.csv", header = T)
  Jurkat_pairedIVT$ACP<-paste0(Jurkat_pairedIVT$Annotation,Jurkat_pairedIVT$chr,Jurkat_pairedIVT$position)
  
  panIVT_Naive<-read.csv("/Dropbox/merged_pvals_Naive_panIVT.csv", header = T)
  panIVT_Naive$ACP<-paste0(panIVT_Naive$Annotation,panIVT_Naive$chr,panIVT_Naive$position)
  panIVT_Jurkat<-read.csv("/Dropbox/merged_pvals_JurkatIVT.csv", header = T)
  panIVT_Jurkat$ACP<-paste0(panIVT_Jurkat$Annotation,panIVT_Jurkat$chr,panIVT_Jurkat$position)
  
  
  #data preparation short way for Fig1 without p-val filtration
  #for Jurkat and Naive together, with >30 filtration on %mm.DRS in Jurkat OR Naive but no pvalue filtration. panIVT enrichment
  if(T){
    
    JurkatDRS_filter <- Jurkat_pairedIVT[which(Jurkat_pairedIVT$N_reads_JurkatDRS > 10  ),]
    NaiveDRS_filter <- Naive_pairedIVT[which(Naive_pairedIVT$N_reads_NaiveDRS > 10 ),]
    
    JurkatNaiveCommon<-intersect(JurkatDRS_filter$ACP,NaiveDRS_filter$ACP)
    #I need to supplement with panIVT just Naive or Jurkat, as the panIVT is the same.
    CommonSites_Jurkat<-JurkatDRS_filter[which(JurkatDRS_filter$ACP %in% JurkatNaiveCommon),]
    CommonSites_Naive <- NaiveDRS_filter[which(NaiveDRS_filter$ACP %in% JurkatNaiveCommon),]
    
    # Create a combined data frame for easier plotting
    CommonSites <- data.frame(
      Annotation=CommonSites_Jurkat$Annotation,
      chr=CommonSites_Jurkat$chr,
      position=CommonSites_Jurkat$position,
      strand=CommonSites_Jurkat$strand,
      kmer=CommonSites_Jurkat$kmer,
      ACP=CommonSites_Jurkat$ACP,
      T_JurkatDRS=CommonSites_Jurkat$T_JurkatDRS,
      C_JurkatDRS=CommonSites_Jurkat$C_JurkatDRS,
      T_JurkatIVT=CommonSites_Jurkat$T_JurkatIVT,
      C_JurkatIVT=CommonSites_Jurkat$C_JurkatIVT,
      T_NaiveDRS=CommonSites_Naive$T_NaiveDRS,
      C_NaiveDRS=CommonSites_Naive$C_NaiveDRS,
      T_NaiveIVT=CommonSites_Naive$T_NaiveIVT,
      C_NaiveIVT=CommonSites_Naive$C_NaiveIVT,
      N_reads_JurkatDRS=CommonSites_Jurkat$N_reads_JurkatDRS,
      N_reads_JurkatIVT=CommonSites_Jurkat$N_reads_JurkatIVT,
      N_reads_NaiveDRS=CommonSites_Naive$N_reads_NaiveDRS,
      N_reads_NaiveIVT=CommonSites_Naive$N_reads_NaiveIVT,
      mm.JurkatDRS=CommonSites_Jurkat$mm.JurkatDRS,
      mm.JurkatIVT=CommonSites_Jurkat$mm.JurkatIVT,
      mm.NaiveDRS=CommonSites_Naive$mm.NaiveDRS,
      mm.NaiveIVT=CommonSites_Naive$mm.NaiveIVT,
      p.value.JurkatDRS=CommonSites_Jurkat$p.value.JurkatDRS,
      p.value.NaiveDRS=CommonSites_Naive$p.value.NaiveDRS
      
    )
    
    #filtering so both jurkat and Naive have never a %mm=0 in both (means it's not psi in both)
    CommonSitesOR<-CommonSites[which(CommonSites$mm.JurkatDRS>0 | CommonSites$mm.NaiveDRS>0),]
    
    # Filtering for p-value to never be =1 in both (menas it's not psi in both for sure)
    CommonSitesORp <- CommonSites %>% 
      filter(!(p.value.JurkatDRS == 1 & p.value.NaiveDRS == 1))
    CommonSitesORp <-CommonSitesORp %>%
      filter(!(p.value.JurkatDRS == Inf | p.value.NaiveDRS == Inf))
    CommonSitesORp <-CommonSitesORp %>%
      filter(!(is.na(p.value.JurkatDRS) | is.na(p.value.NaiveDRS) ))
    
    #just keeping sites that have >30% mm in Naive or Jurkat 
    CommonSitesORp<-CommonSitesORp[which(CommonSitesORp$mm.NaiveDRS>30 | CommonSitesORp$mm.JurkatDRS >30),]
    
    ### Substitute mmIVT==0 with panIVT
    CommonSitesORp$ACP<-paste0(CommonSitesORp$Annotation,CommonSitesORp$chr,CommonSitesORp$position)
    ACP<-CommonSitesORp$ACP[which(CommonSitesORp$N_reads_JurkatIVT < 10 | CommonSitesORp$N_reads_NaiveIVT <10)]
    #filter the panIVT to have #readsIVT and mm according to our criteria
    panIVT_Jurkat_filt<-panIVT_Jurkat[which(panIVT_Jurkat$N_reads_panIVT>10 & panIVT_Jurkat$mm.panIVT<10),]
    panIVT_Naive_filt<-panIVT_Naive[which(panIVT_Naive$N_reads_panIVT>10 & panIVT_Naive$mm.panIVT<10),]
    #find the intersection of the 2 dataframes from jurkat an naive pan IVT. The positions may differ but the mm and #reads are same in J and N
    panIVT_JurkatNaive_idx<-intersect(panIVT_Jurkat_filt$ACP, panIVT_Naive_filt$ACP)
    panIVT_JurkatNaive_df<-panIVT_Jurkat_filt[which(panIVT_Jurkat_filt$ACP %in% panIVT_JurkatNaive_idx),]
    panIVT_JurkatNaive<-data.frame(
      ACP=  panIVT_JurkatNaive_df$ACP,
      mm.panIVT= panIVT_JurkatNaive_df$mm.panIVT,
      T_panIVT= panIVT_JurkatNaive_df$T_panIVT,
      C_panIVT= panIVT_JurkatNaive_df$C_panIVT,
      N_reads_panIVT= panIVT_JurkatNaive_df$N_reads_panIVT
    )
    #filter to search in the merged panIVt having at least one C or T at the interesting position
    panIVT_JurkatNaive<-panIVT_JurkatNaive[which( (panIVT_JurkatNaive$T_panIVT+panIVT_JurkatNaive$C_panIVT)!=0 ),]
    ACP_pan<-panIVT_JurkatNaive$ACP[which(panIVT_JurkatNaive$ACP %in% ACP)]
    pan<-panIVT_JurkatNaive[which(panIVT_JurkatNaive$ACP %in% ACP),]
    for (i in c(1:length(ACP_pan))) {
      if (ACP_pan[i] %in% CommonSitesORp$ACP){
        idx<-which(CommonSitesORp$ACP==ACP_pan[i])
        CommonSitesORp$N_reads_JurkatIVT[idx]<-panIVT_JurkatNaive$N_reads_panIVT[which(panIVT_JurkatNaive$ACP %in% ACP_pan[i])]
        CommonSitesORp$N_reads_NaiveIVT[idx]<-panIVT_JurkatNaive$N_reads_panIVT[which(panIVT_JurkatNaive$ACP %in% ACP_pan[i])]
        CommonSitesORp$T_JurkatIVT[idx]<-panIVT_JurkatNaive$T_panIVT[which(panIVT_JurkatNaive$ACP %in% ACP_pan[i])]
        CommonSitesORp$C_JurkatIVT[idx]<-panIVT_JurkatNaive$C_panIVT[which(panIVT_JurkatNaive$ACP %in% ACP_pan[i])]
        CommonSitesORp$T_NaiveIVT[idx]<-panIVT_JurkatNaive$T_panIVT[which(panIVT_JurkatNaive$ACP %in% ACP_pan[i])]
        CommonSitesORp$C_NaiveIVT[idx]<-panIVT_JurkatNaive$C_panIVT[which(panIVT_JurkatNaive$ACP %in% ACP_pan[i])]
        CommonSitesORp$mm.JurkatIVT[idx]<-panIVT_JurkatNaive$mm.panIVT[which(panIVT_JurkatNaive$ACP %in% ACP_pan[i])]
        CommonSitesORp$mm.NaiveIVT[idx]<-panIVT_JurkatNaive$mm.panIVT[which(panIVT_JurkatNaive$ACP %in% ACP_pan[i])]
      }
      
      period= round(length(ACP_pan) / 100)
      if (i%%period == 0) { 
        print(paste0(round(i/length(ACP_pan) * 100, 2),"%")) 
      }
    }
    
    CommonSitesORp <- CommonSitesORp[which(CommonSitesORp$N_reads_JurkatIVT > 10),]  #let's keep only sites with total number of reads>10.
    CommonSitesORp <- CommonSitesORp[which(CommonSitesORp$N_reads_NaiveIVT > 10),]
    CommonSitesORp<-CommonSitesORp[which((CommonSitesORp$T_JurkatIVT+CommonSitesORp$C_JurkatIVT)!=0),] #but let's also keep only those sites with at least some T and C in that positions
    CommonSitesORp<-CommonSitesORp[which((CommonSitesORp$T_NaiveIVT+CommonSitesORp$C_NaiveIVT)!=0),]
    
    ### Add the color column based on the kmer
    CommonSitesORp$color <- "black"
    for (i in c(1:nrow(CommonSitesORp))) {
      if (CommonSitesORp$kmer[i]  %in%  c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")) { #UNUAR
        CommonSitesORp$color[i] <- "blue"
      }
      if (CommonSitesORp$kmer[i] %in%  c("GTTCA", "GTTCT", "GTTCC", "GTTCG")) { #GUUCN
        CommonSitesORp$color[i] <- "red"
      }
      else {}
    }

    
  }    
  write.csv(CommonSitesORp,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveJurkat_panIVT/Jurkat_panIVT/init_gene_pileup/finalized_pileup/Jurkat_Naive_filtered_pan_noPvalue.csv",row.names = F)
    
  
  #data preparation for panIVT with pval<0.001. without p-val <0.001 filtration takes too long. Files with p-value filtration are:
  #Jurkat_filtered_pan.csv
  #Naive_filtered_pan.csv
  if (T) {
    ################# For Jurkat
    if (T) {

      ### Filter to total_number_reads_TRUB1_KD & total_number_reads_SCR < 10
      JurkatDRS_filter <- Jurkat_pairedIVT[which(Jurkat_pairedIVT$N_reads_JurkatDRS > 10 & Jurkat_pairedIVT$p.value<0.001 ),]
  
      ### Substitute mmIVT==0 with panIVT
      JurkatDRS_filter$ACP<-paste0(JurkatDRS_filter$Annotation,JurkatDRS_filter$chr,JurkatDRS_filter$position)
      ACP<-JurkatDRS_filter$ACP[which(JurkatDRS_filter$N_reads_JurkatIVT < 10)]
      panIVT_Jurkat_filt<-panIVT_Jurkat[which(panIVT_Jurkat$N_reads_panIVT>10 & panIVT_Jurkat$mm.panIVT<10),]
      ACP_pan<-panIVT_Jurkat_filt$ACP[which(panIVT_Jurkat_filt$ACP %in% ACP)]
      ACP_pan_mm<-panIVT_Jurkat_filt$mm.panIVT[which(panIVT_Jurkat_filt$ACP %in% ACP)]
      pan<-panIVT_Jurkat_filt[which(panIVT_Jurkat_filt$ACP %in% ACP),]
      for (i in c(1:length(ACP_pan))) {
        if (ACP_pan[i] %in% JurkatDRS_filter$ACP){
          idx<-which(JurkatDRS_filter$ACP==ACP_pan[i])
          JurkatDRS_filter$N_reads_JurkatIVT[idx]<-panIVT_Jurkat_filt$N_reads_panIVT[which(panIVT_Jurkat_filt$ACP %in% ACP_pan[i])]
          JurkatDRS_filter$T_JurkatIVT[idx]<-panIVT_Jurkat_filt$T_panIVT[which(panIVT_Jurkat_filt$ACP %in% ACP_pan[i])]
          JurkatDRS_filter$C_JurkatIVT[idx]<-panIVT_Jurkat_filt$C_panIVT[which(panIVT_Jurkat_filt$ACP %in% ACP_pan[i])]
          JurkatDRS_filter$mm.JurkatIVT[idx]<-panIVT_Jurkat_filt$mm.panIVT[which(panIVT_Jurkat_filt$ACP %in% ACP_pan[i])]
        }
        
        period= round(length(ACP_pan) / 100)
        if (i%%period == 0) { 
          print(paste0(round(i/length(ACP_pan) * 100, 2),"%")) 
        }
      }
      
      JurkatDRS_filter <- JurkatDRS_filter[which(JurkatDRS_filter$N_reads_JurkatIVT > 10),]  #let's keep only sites with total number of reads>10
      JurkatDRS_filter<-JurkatDRS_filter[which((JurkatDRS_filter$T_JurkatIVT+JurkatDRS_filter$C_JurkatIVT)!=0),] #but let's also keep only those sites with at least some T and C in that positions
      
      ### Add total_TRUB1_KD & total_SCR
      JurkatDRS_filter$total_DRS <- JurkatDRS_filter$T_JurkatDRS + JurkatDRS_filter$C_JurkatDRS
      JurkatDRS_filter$total_IVT <- JurkatDRS_filter$T_JurkatIVT + JurkatDRS_filter$C_JurkatIVT
      
      ### Add the color column based on the kmer
      JurkatDRS_filter$color <- "black"
      for (i in c(1:nrow(JurkatDRS_filter))) {
        if (JurkatDRS_filter$kmer[i]  %in%  c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")) { #UNUAR
          JurkatDRS_filter$color[i] <- "blue"
        }
        if (JurkatDRS_filter$kmer[i] %in%  c("GTTCA", "GTTCT", "GTTCC", "GTTCG")) { #GUUCN
          JurkatDRS_filter$color[i] <- "red"
        }
        else {}
      }
      
    }
    write.csv(JurkatDRS_filter,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveJurkat_panIVT/Jurkat_panIVT/init_gene_pileup/finalized_pileup/Jurkat_filtered_pan.csv",row.names = F)
    
    ################# For Naive
    if (T) {
      
      ### Filter to Nreads_DRS
      NaiveDRS_filter <- Naive_pairedIVT[which(Naive_pairedIVT$N_reads_NaiveDRS > 10 & Naive_pairedIVT$p.value<0.001),]
    
      ### Substitute mmIVT==0 with panIVT
      NaiveDRS_filter$ACP<-paste0( NaiveDRS_filter$Annotation, NaiveDRS_filter$chr, NaiveDRS_filter$position)
      ACP<- NaiveDRS_filter$ACP[which(NaiveDRS_filter$N_reads_NaiveIVT < 10)]
      ACP_pan<-panIVT_Naive$ACP[which(panIVT_Naive$ACP %in% ACP)]
      ACP_pan_mm<-panIVT_Naive$mm.panIVT[which(panIVT_Naive$ACP %in% ACP)]
      
      for (i in c(1:length(ACP_pan))) {
        if (ACP_pan[i] %in% NaiveDRS_filter$ACP){
          idx<-which(NaiveDRS_filter$ACP==ACP_pan[i])
          NaiveDRS_filter$N_reads_NaiveIVT[idx]<-panIVT_Naive$N_reads_panIVT[which(panIVT_Naive$ACP %in% ACP_pan[i])]
          NaiveDRS_filter$T_IVT[idx]<-panIVT_Naive$T_panIVT[which(panIVT_Naive$ACP %in% ACP_pan[i])]
          NaiveDRS_filter$C_IVT[idx]<-panIVT_Naive$C_panIVT[which(panIVT_Naive$ACP %in% ACP_pan[i])]
          NaiveDRS_filter$mm.NaiveIVT[idx]<-panIVT_Naive$mm.panIVT[which(panIVT_Naive$ACP %in% ACP_pan[i])]
        }
        
        period= round(length(ACP_pan) / 100)
        if (i%%period == 0) { 
          print(paste0(round(i/length(ACP_pan) * 100, 2),"%")) 
        }
        
      }

      NaiveDRS_filter <-  NaiveDRS_filter[which( NaiveDRS_filter$N_reads_NaiveIVT > 10),]
      NaiveDRS_filter<-NaiveDRS_filter[which((NaiveDRS_filter$T_NaiveIVT +NaiveDRS_filter$C_NaiveIVT)!=0),]
      
      ### Add total_PUS7_KD & total_PUS7_KD
      NaiveDRS_filter$total_DRS <-  NaiveDRS_filter$T_NaiveDRS + NaiveDRS_filter$C_NaiveDRS
      NaiveDRS_filter$total_IVT <-  NaiveDRS_filter$T_NaiveIVT + NaiveDRS_filter$C_NaiveIVT

      ### Add the color column based on the kmer
      NaiveDRS_filter$color <- "black"
      for (i in c(1:nrow(NaiveDRS_filter))) {
        if (NaiveDRS_filter$kmer[i]  %in%  c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")) { #UNUAR
          NaiveDRS_filter$color[i] <- "blue"
        }
        if (NaiveDRS_filter$kmer[i] %in%  c("GTTCA", "GTTCT", "GTTCC", "GTTCG")) {
          NaiveDRS_filter$color[i] <- "red"
        }
        else {}
      }
      
      
    }
    write.csv(NaiveDRS_filter,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveJurkat_panIVT/Naive_panIVT/init_gene_pileup/finalized_pileup/Naive_filtered_pan.csv",row.names = F)
    
}

  }

### plots figure co-expressed transcripts in Naive and Jurkat
if(T){
 if(T){
  
  #write Naive and Jurkat sites
  if(T){
  #supplemented with panIVT and p-val<0.001. NOT USING ANYMORE CAUSE WE DON'T WANT TO FILTER FOR P_VALUE
  # Naive_orig <- read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveJurkat_panIVT/Naive_panIVT/init_gene_pileup/finalized_pileup/Naive_filtered_pan.csv", header = T)
  # Jurkat_orig <- read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveJurkat_panIVT/Jurkat_panIVT/init_gene_pileup/finalized_pileup/Jurkat_filtered_pan.csv", header = T)
  
  #NaiveJurkat contains sites common to both datasets, panIVT enriched
  NaiveJurkat<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveJurkat_panIVT/Jurkat_panIVT/init_gene_pileup/finalized_pileup/Jurkat_Naive_filtered_pan_noPvalue.csv",header = T)
  Naive_orig<-NaiveJurkat
  Jurkat_orig<-NaiveJurkat
  
  #ALL THESE HAVE A P_VAL FILTRATION OF 0.001
  Naive_filt<-Naive_orig
  
  Jurkat_filt<-Jurkat_orig

  #common values are all those sites with high coverage in DRS and IVT and a low %mm in IVT (to exclude SNVs)
  common_values <- rbind(Naive_filt, Jurkat_filt)
  common_values <- common_values[which(common_values$mm.JurkatIVT<10 & common_values$mm.NaiveIVT<10),]
  common_values <- common_values[which(common_values$N_reads_JurkatIVT>10 & common_values$N_reads_NaiveIVT>10),]
  common_values <- common_values[which(common_values$N_reads_JurkatDRS>20 & common_values$N_reads_NaiveDRS>20),]
  common_values <- common_values[which(common_values$p.value.JurkatDRS<0.01 | common_values$p.value.NaiveDRS<0.01),]

  
  #vcf
  if(T){
    # Search my sites in the SNV dataset
    rna_modifications <- data.frame(
      chr = common_values$chr,  # Replace with your actual chromosome data
      position = common_values$position,  # Replace with your actual positions
      annotation = common_values$Annotation  # Replace with your actual annotations
    )
    
    rna_mod_gr <- GRanges(
      seqnames = rna_modifications$chr,
      ranges = IRanges(start = rna_modifications$position, end = rna_modifications$position),
      annotation = rna_modifications$annotation
    )
    
    vcf_file <- "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/vcf_GRCh38p10/Jurkat_pileup.vcf"  # Replace with the path to your VCF file
    vcf <- readVcf(vcf_file, "hg18") 
    snv_gr_Jurkat <- rowRanges(vcf)
    vcf_file <- "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/vcf_GRCh38p10/Naive_pileup.vcf"  # Replace with the path to your VCF file
    vcf <- readVcf(vcf_file, "hg18") 
    snv_gr_Naive <- rowRanges(vcf)
    #concatenating both SNV from Jurkat and Naive in one large object
    snv_gr<-c(snv_gr_Jurkat,snv_gr_Naive)
    
    # Get the indices of the overlapping sites between my sites and the snv_gr containing all the SNV indexes
    overlaps <- findOverlaps(rna_mod_gr, snv_gr)
    overlapping_indices <- queryHits(overlaps)
    overlapping_sites <- rna_mod_gr[overlapping_indices]
    
    # Get all the non-overlapping sites that are NOT SNVs
    all_indices <- seq_len(nrow(rna_modifications))
    overlapping_indices <- queryHits(overlaps)
    non_overlapping_indices <- setdiff(all_indices, overlapping_indices)
    non_overlapping_sites <- rna_modifications[non_overlapping_indices, ]
    non_overlapping_ACP<- paste0(non_overlapping_sites$annotation, non_overlapping_sites$chr,non_overlapping_sites$position)
    #removing all SNVs from the list
    common_values<-common_values[which(common_values$ACP %in% non_overlapping_ACP),]
    
  }
  
  #venn diagram sets FILES WRITING, unfiltered for mm<10 to assess absence/presence.  
  if(T){

    #separating cause I want to find common and unique sites to jurkat and naive
    Naive_filt<-data.frame(chr=common_values$chr,position=common_values$position,Annotation=common_values$Annotation,ACP=common_values$ACP,N_reads_NaiveDRS=common_values$N_reads_NaiveDRS,N_reads_NaiveIVT=common_values$N_reads_NaiveIVT,mm.NaiveDRS=common_values$mm.NaiveDRS,mm.NaiveIVT=common_values$mm.NaiveIVT,p.value.NaiveDRS=common_values$p.value.NaiveDRS)
    Jurkat_filt<-data.frame(chr=common_values$chr,position=common_values$position,Annotation=common_values$Annotation,ACP=common_values$ACP,N_reads_JurkatDRS=common_values$N_reads_JurkatDRS,N_reads_JurkatIVT=common_values$N_reads_JurkatIVT,mm.JurkatDRS=common_values$mm.JurkatDRS,mm.JurkatIVT=common_values$mm.JurkatIVT,p.value.JurkatDRS=common_values$p.value.JurkatDRS)
  
    Naive_filt_v<-Naive_filt[which(Naive_filt$mm.NaiveDRS>=30 & Naive_filt$mm.NaiveIVT<10 ),]
    write.csv(Naive_filt_v$ACP,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Naive_set_venn.csv",row.names=FALSE)
  
    Jurkat_filt_v<-Jurkat_filt[which(Jurkat_filt$mm.JurkatDRS>=30 & Jurkat_filt$mm.JurkatIVT<10 ),]
    write.csv(Jurkat_filt_v$ACP,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Jurkat_set_venn.csv",row.names=FALSE)
  
    }
  }
  
  #psi panel, common sites defined as %mm>30 in one cell line and %mm>=10 in the other 
  #unique sites defined as %mm>30 in one cell line and %mm<10 in the other 
  if(T){
    #read venn files with ACP names only 
    Naive_venn<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Naive_set_venn.csv")
    Jurkat_venn<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Jurkat_set_venn.csv")
    #read file with all the info about Naive and Jurkat sites
    NaiveJurkat<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveJurkat_panIVT/Jurkat_panIVT/init_gene_pileup/finalized_pileup/Jurkat_Naive_filtered_pan_noPvalue.csv",header = T)
    
    #read manually excluded sites names
    excluded<-c("DIAPH1chr5141573627","SELLchr1169701554","CEP295chr1193721658","ZNF337chr2025674348",
                "IKZF1chr750404512","SSNA1chr9137189053","ACOT8chr2045857191",
                "JADE2chr5134578694","TRIM14chr998087768","UCHL3chr1375605736","POP7chr7100707137",
                "RNF185chr2231206868","CTDSP1chr2218403055","UBE2Hchr7129831273","TMUB1chr7151082501",
                "DRAM2chr1111139690","ERAL1chr1728858165")
    
    # Find the intersections
    intersect_1_2 <- intersect(Naive_venn$x, Jurkat_venn$x)  # Intersection of Group 1 and Group 2
    
    # Unique elements in each group
    unique_N <- setdiff(Naive_venn$x, Jurkat_venn$x)  # Unique to Group 1
    unique_J <- setdiff(Jurkat_venn$x, Naive_venn$x)  # Unique to Group 2

    #plot
    all_genes<-c(unique_J,intersect_1_2,unique_N)
    df_Naive <- NaiveJurkat[which(NaiveJurkat$ACP %in% all_genes),]
    df_Jurkat<- NaiveJurkat[which(NaiveJurkat$ACP %in% all_genes),]
    
    #This is the cutoff for defining the horizontal lines in the plots
    delta_cutoff=20
    
    # Create a combined data frame for easier plotting
    df_combined <- data.frame(
      Index = 1:nrow(df_Naive),  # Assuming same number of rows
      ACP = df_Naive$ACP,
      chr= df_Naive$chr,
      position=df_Naive$position,
      Annotation <- df_Naive$Annotation,
      kmer <- df_Naive$kmer,
      mm.Naive = df_Naive$mm.NaiveDRS,
      mm.Jurkat = df_Jurkat$mm.JurkatDRS,
      mm.NaiveIVT = df_Naive$mm.NaiveIVT,
      mm.JurkatIVT = df_Jurkat$mm.JurkatIVT,
      delta = df_Naive$mm.NaiveDRS-df_Jurkat$mm.JurkatDRS,
      N_reads_JurkatIVT = df_Jurkat$N_reads_JurkatIVT,
      N_reads_NaiveIVT = df_Naive$N_reads_NaiveIVT,
      N_reads_JurkatDRS = df_Jurkat$N_reads_JurkatDRS,
      N_reads_NaiveDRS = df_Naive$N_reads_NaiveDRS
    )
    
    df_combined<-df_combined[!df_combined$ACP %in% excluded, ] #remove manually selected sites
    
    df_combined$color<-"cyan"
    df_combined$color[which(df_combined$mm.Naive>30 & df_combined$mm.Jurkat<10)]<-"orange"
    df_combined$color[which(df_combined$mm.Naive<10 & df_combined$mm.Jurkat>30)]<-"blue"
    
    df_combined_shared<-df_combined[which(df_combined$color=="cyan"),]
    write.csv(df_combined_shared,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/CommonSitesNaiveJurkat.csv",row.names = F)
    
    df_combined <- df_combined %>%
      arrange(color,Index)
    
    median_cyan <- median(df_combined$delta[df_combined$color == "cyan"])
    median_lightblue <- median(df_combined$delta[df_combined$color == "blue"])
    median_orange <- median(df_combined$delta[df_combined$color =="orange"])
   
    df_combined$Index <- factor(df_combined$Index, levels = df_combined$Index[order(df_combined$color, df_combined$Index)])
    
    #order in Ascending order according to delta values N-J
    if(T){
      
      cyan_ordered <- df_combined %>%
        filter(color == "cyan") %>%
        arrange(delta)
      
      orange_ordered <- df_combined %>%
        filter(color == "orange") %>%
        arrange(delta)
      
      blue_ordered <- df_combined %>%
        filter(color == "blue") %>%
        arrange(delta)
      
      # Combine the ordered cyan points with the rest of the data
      df_combined <- bind_rows(
        blue_ordered,
        cyan_ordered, 
        orange_ordered
      )
      
      # Update the Index factor to reflect the new order
      df_combined$Index <- factor(df_combined$Index, levels = df_combined$Index)
      
      
    }
    
    #for Venn diagram
    common<-c(df_combined$ACP[which((df_combined$mm.Naive>30 & df_combined$mm.Jurkat>=10)|(df_combined$mm.Naive>=10 & df_combined$mm.Jurkat>30))])
    naive_sites <-  c(df_combined$ACP[which(df_combined$mm.Naive>30 & df_combined$mm.Jurkat<10 )], common)
    #same as df_combined$ACP[which(df_combined$color=="orange" | df_combined$color=="cyan")]
    jurkat_sites <- c(df_combined$ACP[which(df_combined$mm.Jurkat>30 & df_combined$mm.Naive<10 )], common)
    #same as df_combined$ACP[which(df_combined$color=="blue" | df_combined$color=="cyan")]
    n<-df_combined[which(df_combined$mm.Naive>30 & df_combined$mm.Jurkat<10 ),]
    j<-df_combined[which(df_combined$mm.Jurkat>30 & df_combined$mm.Naive<10 ),]
    cm<-df_combined[which(df_combined$ACP %in% common),]
    
    # Plot with individual points, delta values Naive-Jurkat
    if(F){
    ggplot(df_combined, aes(x = Index, y = delta)) +
      # Add median lines
      geom_hline(yintercept = median_cyan, color = "cyan", linetype = "dashed", size = 0.25) +
      geom_hline(yintercept = delta_cutoff, color = "red", linetype = "dashed", size = 0.25) +
      geom_hline(yintercept = -delta_cutoff, color = "red", linetype = "dashed", size = 0.25) +
      # Plot points with colors from the 'color' column
      geom_point(aes(color = color), size = 0.5) +
      scale_color_identity() +  # Use the color column for colors
      # Labels and theme
      labs(title = "Delta Values Colored by Group", x = "Index", y = "Delta") +
      theme_minimal() +  # Use minimal theme
      theme(
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank()   # Remove minor grid lines
      )
    }
    
    # Plot delta values, Naive-Jurkat traces
    if(T){
      ggplot(df_combined, aes(x = Annotation)) +
        # Plot points and lines for mm.Naive (blue)
        geom_line(aes(y = delta, group = 1, color = color)) +
        # Add median lines
        geom_hline(yintercept = 10, color = "red", linetype = "dashed", size = 0.25) +
        geom_hline(yintercept = -10, color = "red", linetype = "dashed", size = 0.25) +
        # Add labels and theme
        labs(x = "Index", y = "Values", title = "mm.Naive and mm.Jurkat Values") +
        ylim(-100,100)+
        theme_minimal() +  # Use minimal theme
        theme(
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),   # Remove minor grid lines
        )
      }
    
    # Plot %mm values traces
    if(T){
      ggplot(df_combined, aes(x = Annotation)) +
        # Plot points and lines for mm.Naive (blue)
        geom_line(aes(y = mm.Naive, group = 1), color = "orange") +
        # Add threshold lines
         geom_hline(yintercept = 10, color = "black", linetype = "dashed", size = 1) +
         geom_hline(yintercept = 30, color = "red", linetype = "dashed", size = 1) +
        # geom_hline(yintercept = median_orange, color = "orange", linetype = "dashed", size = 1) +
        geom_line(aes(y = mm.Jurkat, group = 1), color = "royalblue") +
        # Add labels and theme
        labs(x = "Index", y = "Values", title = "mm.Naive and mm.Jurkat Values") +
        theme_minimal() +  # Use minimal theme
        theme(
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),   # Remove minor grid lines
          axis.text.x = element_text(angle = 80, vjust = 0.5, hjust = 1, size=6) 
        )
       }

    #pie chart % variable and stable occupancy for shared sites
    if(T){
        
        data<-data.frame(
          Naive<-df_combined$mm.Naive[which(df_combined$color=="cyan")],
          Jurkat<-df_combined$mm.Jurkat[which(df_combined$color=="cyan")]
        )
        colnames(data)<-c("Naive","Jurkat")
        # Calculate delta
        data$delta <- data$Jurkat - data$Naive
        
        # Categorize as stable or variable
        data$category <- ifelse(abs(data$delta) <= 20, "Stable", "Variable")
        
        # Count stable and variable points
        counts <- table(data$category)
        
        # Plot pie chart
        pie(counts, labels = paste0(names(counts)," ",counts), main = "Stability of %mm_delta",
            col = c("skyblue", "skyblue3"))
        
        
      }
    
    # Create a Venn Diagram for CondA>30 and CondB<10 and CondA>30 and CondB>10 so for unique and shared
    if(T){
      # Define scaling factors based on size (for rough visual approximation)
      naive_size <- length(naive_sites)
      jurkat_size <- length(jurkat_sites)
      
      # Calculate cex and cat.cex based on size difference (adjust scaling_factor as needed)
      scaling_factor <- 0.01
      base_cex <- 1.5
      adjusted_cex <- base_cex + abs(naive_size - jurkat_size) * scaling_factor
      
      # Set higher values for larger dataset
      overall_cex <- if (naive_size > jurkat_size) adjusted_cex else base_cex
      overall_cat_cex <- if (naive_size > jurkat_size) adjusted_cex else base_cex
      
      venn.plot <- venn.diagram(
        x = list(Naive = naive_sites, Jurkat = jurkat_sites),
        category.names = c("Naive", "Jurkat"),
        filename = NULL, # NULL to plot directly to RStudio
        fill = c("orange", "blue"),  # Colors for each group
        alpha = 0.5,                 # Transparency level for fill colors
        cex = overall_cex,  # Use scaled cex values for numbers in each group
        cat.cex = overall_cat_cex,  # Use scaled cat.cex values for labels
        cat.pos = c(-20, 20),         # Positioning of category names
        cat.dist = 0.05,              # Distance of category names from the circles
        cat.col = c( "orange", "blue") # Color of category labels
      )
      # Display the plot
      grid.draw(venn.plot)
   
       }

    #histogram fig1 on Common NaiveJurkat sites %mm>30 in one matched by %mm>10 in the other
    if(T){
      cm$Jurkat_minus_Naive<-cm$mm.Jurkat-cm$mm.Naive
      
      ggplot(cm, aes(x = Jurkat_minus_Naive)) +
        geom_histogram(binwidth = 1, fill="skyblue", ) +
        geom_vline(xintercept = 0, color="black", size=1) + # Black line at 0
        geom_vline(xintercept = mean(cm$Jurkat_minus_Naive), color="red", size=1) + # Red line at mean of mm
        labs(title="Histogram of %mm delta", x="(Jurkat - Naive) mismatch", y="Count")+
        xlim(-100,100)
      ggsave("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/histogramCommon.pdf")
    }
    
    #seqLogos for variable/stable shared sites and unique 
    if (plot_figs==T){
      
      common<-df_combined[which(df_combined$color=="cyan"),]
      common$category <- ifelse(abs(common$delta) <= 20, "Stable", "Variable")
      
      uniqueN<-df_combined[which(df_combined$color=="orange"),]
      uniqueJ<-df_combined[which(df_combined$color=="blue"),]
      
      #for unique
      kmers<-sapply( uniqueN$kmer,function(x) str_replace_all(x,"T","U"))
      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/logo_UniqueNaive.pdf"
      
      kmers<-sapply( uniqueJ$kmer,function(x) str_replace_all(x,"T","U"))
      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/logo_UniqueJurkat.pdf"
      
      #for common
      #delta <=-20, delta >=20 
      # -20<delta<=-10,  -10<delta<10,  10<=delta<20
    
      kmers<-sapply( common$kmer[which(common$delta<=-20)],function(x) str_replace_all(x,"T","U"))
      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/logo_CommonNaive-20.pdf"
      
      kmers<-sapply( common$kmer[which(common$delta>=20)],function(x) str_replace_all(x,"T","U"))
      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/logo_CommonNaive20.pdf"
      
      kmers<-sapply(common$kmer[which(common$delta>-20 & common$delta<=-10 )],function(x) str_replace_all(x,"T","U") )
      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/logo_CommonNaiveJurkat-20-10mm.pdf"
      
      #highly stable
      kmers<-sapply(common$kmer[which(common$delta>-10 & common$delta<10 )],function(x) str_replace_all(x,"T","U") )
      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/logo_CommonNaiveJurkat-1010mm.pdf"
      
      kmers<-sapply(common$kmer[which(common$delta>=10 & common$delta<20 )],function(x) str_replace_all(x,"T","U") )
      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/logo_CommonNaiveJurkat1020mm.pdf"
      
      # kmers<-sapply(common$kmer[which(common$delta>-5 & common$delta<5 )],function(x) str_replace_all(x,"T","U") )
      # file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/logo_CommonNaiveJurkat-55mm.pdf"

      
      my_list<-list("new motifs" = kmers)
      min.x = 0
      max.x = 200
      min.y = 0
      max.y = 20
      line.thickness = 1
      label.size = 1
      text.size = 0.8
      point.size = 1
      tck.length = 0.01
      tick.thickness = 1 
      transparency = 0.2
      Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
      CairoPDF(file.name, width = 6, height = 3, family = "Helvetica")
      par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
      plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
           ylab = NA,xlab=NA,ann = F,axes=F)
      
      ggseqlogo::ggseqlogo(my_list,ncol=1,method="prob")
      dev.off()
      
      #highly stable sites, percentage of highly stable sites containing the GUUCN motif
      if(T){
        
        #highly stable
        kmers<-sapply(common$kmer[which(common$delta>-10 & common$delta<10 )],function(x) str_replace_all(x,"T","U") )
        file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/logo_CommonNaiveJurkat-1010mm.pdf"
        k<-N_GUUCN<-length( which(kmers %in% c("GUUCA", "GUUCU", "GUUCC", "GUUCG") ) )
        n<-length(kmers)
        N_GUUCN/length(kmers) #percentage of GUUCN sites
        
        
        kmers_non_stable<-sapply(c(common$kmer[which(common$delta< -10)], common$kmer[which(common$delta>10)]),function(x) str_replace_all(x,"T","U") )
        N_GUUCN<-length( which(kmers_non_stable %in% c("GUUCA", "GUUCU", "GUUCC", "GUUCG") ) )
        p_GUUCN<-N_GUUCN/length(kmers_non_stable) #percentage of GUUCN sites
        
        binom.test(x = k, n = n, p = p_GUUCN, alternative = "greater")
        
        
      }
      
      #percentage of cytidine at the -1 site, to repeat for each kmer common group
      if(T){
        
        count_C_second <- sum(substr(kmers, 2, 2) == "C")
        percent_C_second <- (count_C_second / length(kmers)) * 100
        percent_C_second
        
      }
      
      
      
    }
  
  }
  
  #show scatterplot with Jurkat and Naive common genes and a box for >20 reads to show how we filter
  if(T){
    
    NaiveJurkat<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveJurkat_panIVT/Jurkat_panIVT/init_gene_pileup/finalized_pileup/Jurkat_Naive_filtered_pan_noPvalue.csv",header = T)
    
    ggplot(NaiveJurkat) +
      # Plot NaiveDRS in orange
      geom_point(aes(x = Annotation , y = log10(N_reads_JurkatDRS)), color = "orange", shape=16, size=0.1) +
      # Plot JurkatDRS in blue
      geom_point(aes(x = Annotation, y = log10(N_reads_NaiveDRS)), color = "blue", alpha=0.5, shape=16, size=0.1) +
      geom_hline(yintercept = log10(20), color = "black", linetype = "dashed", size = 0.25)+
      # Adding labels and titles
      labs(title = "Scatterplot of mm.NaiveDRS vs mm.JurkatDRS",
           x = "Index",
           y = "N_reads") +
      theme_minimal()+
      theme(
        panel.background = element_rect(fill = "white"),  # Set background to white
        axis.line = element_line(color = "black"),        # Black axis lines
        panel.grid.major = element_blank(),               # Remove major grid lines
        panel.grid.minor = element_blank(),               # Remove minor grid lines
        plot.background = element_rect(fill = "white"),    # Ensure the plot background is white
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6)
      )
    ggsave("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/N_reads_cutoff.pdf", width = 4, height=1)
    
    
    ggplot(NaiveJurkat) +
      # Plot NaiveDRS in orange with transparency (alpha = 0.5)
      geom_point(aes(x = log10(N_reads_JurkatDRS), y = mm.JurkatDRS), color = "orange", alpha = 0.5) +
      # Plot JurkatDRS in blue with transparency (alpha = 0.5)
      geom_point(aes(x = log10(N_reads_NaiveDRS), y = mm.NaiveDRS), color = "blue", alpha = 0.5) +
      geom_vline(xintercept = log10(20), color="black", linetype="dashed")+
      geom_hline(yintercept = 10, color="black", linetype="dashed")+
      # Adding labels and titles
      labs(title = "Scatterplot of mm.NaiveDRS vs mm.JurkatDRS",
           x = "mm.NaiveDRS",
           y = "mm.JurkatDRS") +
      theme_minimal()+
      theme(
        panel.background = element_rect(fill = "white"),  # Set background to white
        axis.line = element_line(color = "black"),        # Black axis lines
        panel.grid.major = element_blank(),               # Remove major grid lines
        panel.grid.minor = element_blank(),               # Remove minor grid lines
        plot.background = element_rect(fill = "white")    # Ensure the plot background is white
      )
    
  }
  
  #heatmaps for presence and absence
  if(T){
    #presence
    if(T){
      
      nread=20
      mm=30

      common_values_filt<-common_values[which(common_values$p.value.JurkatDRS<0.001 & common_values$p.value.NaiveDRS<0.001),] #In this case I want both libraries with low pvalue cause I'm looking for hypermod sites
      common_values_filt<-common_values_filt[which(common_values_filt$N_reads_JurkatDRS>nread & common_values_filt$N_reads_NaiveDRS>nread),]
      common_values_filt<-common_values_filt[which(common_values_filt$mm.JurkatDRS>mm | common_values_filt$mm.NaiveDRS>mm),]
      
      ALL_NaiveJurkat <-common_values_filt #maybe select for difference values, so calculate the similarity here and not later on line 943 
      
      ALL_NaiveJurkat$color<-"black"
      for (i in c(1:length( ALL_NaiveJurkat$ACP))) {
        if ( ALL_NaiveJurkat$kmer.x[i]  %in%  c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")) { #UNUAR
          ALL_NaiveJurkat$color[i] <- "blue"
        }
        else if ( ALL_NaiveJurkat$kmer.x[i] %in%  c("GTTCA", "GTTCT", "GTTCC", "GTTCG")) {
          ALL_NaiveJurkat$color[i] <- "red"
        }
      }
      
      # Naive Jurkat
      ALL_NaiveJurkat$similarity<-ALL_NaiveJurkat$mm.JurkatDRS-ALL_NaiveJurkat$mm.NaiveDRS
      ALL_NaiveJurkat <- ALL_NaiveJurkat[order(ALL_NaiveJurkat$similarity ,decreasing=TRUE), ] #order
      data_matrix <- cbind(ALL_NaiveJurkat$mm.NaiveDRS, ALL_NaiveJurkat$mm.JurkatDRS)
      colnames(data_matrix) <- c("Naive T cells", "Jurkat")
      AC<-paste0(ALL_NaiveJurkat$Annotation.x,", ",ALL_NaiveJurkat$chr.x)
      ACS<-paste0(AC,", ",ALL_NaiveJurkat$position.x)
      ACP<-paste(ACS," ",format( round(ALL_NaiveJurkat$similarity, digits=2 ), nslow=2))
      rownames(data_matrix) <- ACP
      
      fakerow<-c(0,0)
      fakerowname<-"VOID, chrS, 11111"
      data_matrix<-rbind(data_matrix,fakerow)
      rownames(data_matrix)[nrow(data_matrix)]<-fakerowname
      fakerow<-c(100,100)
      data_matrix<-rbind(data_matrix,fakerow)
      rownames(data_matrix)[nrow(data_matrix)]<-fakerowname
      
      #plot HEATMAP
      if(T){
        file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/heatmapNaiveJurkat.pdf"
        min.x = 0
        max.x = 100
        min.y = 0
        max.y = 800
        line.thickness = 1
        label.size = 1
        text.size = 0.8
        point.size = 1
        tck.length = 0.01
        tick.thickness = 1
        transparency = 0.2
        pdf(file.name,width = 5,height = 6)
        par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2)
        
        heatmap.2(data_matrix,
                  Rowv = FALSE,
                  Colv = FALSE,
                  #our colors and breaks
                  #col = colorRampPalette(colors)(length(breaks)-1),
                  col = mako(250),
                  #breaks = breaks,
                  dendrogram = "none",
                  key = FALSE,
                  key.title = "Color Key",
                  density.info = "none",
                  trace ="none",
                  margins=c(5,9), #0 starts from the bottom right, the first value controls if it's up/down. The second value left/right
                  cexRow=0.3,
                  cexCol = 0.8,
                  lhei = c(0.5,5),
                  lwid = c(0.5,5))
        dev.off()
      }
      
      #table pnas
      if(T){ 
        pnas_filter=ALL_NaiveJurkat
        
        min.x = -50
        max.x = 70
        min.y = -10
        max.y = nrow(pnas_filter)+ 5
        line.thickness = 0.8
        txt.size = 0.2
        point.size = 0.3
        label.size = 0.5
        anot.pos = -20
        file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Fig2_heatmapPresence_list.pdf"
        pdf(file=file.name)
        
        par(mfrow=c(1,1),lend=1)
        plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
             ylab = NA,xlab=NA,ann = F,axes=F)
        
        for (i in c(1:nrow(pnas_filter))) {
          ### Add the gene names
          text(x = anot.pos ,y = nrow(pnas_filter)-i, labels = pnas_filter$Annotation.x[i],
               srt = 0,cex = 1.2*txt.size, pos = 4, col = pnas_filter$color[i], family = "Helvetica", font = 3)
          ### Add the chromosome
          text(x = anot.pos+20 ,y = nrow(pnas_filter)-i,
               labels = pnas_filter$chr.x[i] ,cex = 1.2*txt.size, pos = 2, srt = 0) 
          ### Add the position
          text(x = anot.pos+30 ,y = nrow(pnas_filter)-i,
               labels = pnas_filter$position.x[i],cex = 1.2*txt.size, pos = 2, srt = 0) 
          ### Add the similarity
          text(x = anot.pos+40 ,y = nrow(pnas_filter)-i,
               labels = round(pnas_filter$similarity[i],digits = 2) ,cex = 1.2*txt.size, pos = 2, srt = 0) 
        }
        
        ### Add the axis
        text(x = c(anot.pos,anot.pos+10,anot.pos+20,anot.pos+30), y = (nrow(pnas_filter))*c(1,1,1), labels = c("Gene","Chr","pos","Jurkat-Naive %mm"),
             srt = 0,cex = 1.2*txt.size, pos = 4, col = "black",family = "Helvetica", font = 3)
        
        dev.off()
      }
      
    }
    
    #absence
    if(T){
      
      #sites definition and utr/cds assignment
      if(T){
    
      #I'm focusing on the extremes, so I need to pick the "blue" and "orange" sites
      #In this case I'm ok with one significant only, so I use data with no pvalue filtration 
      Naive_only_ACP<-df_combined$ACP[which(df_combined$color=="orange")]
      Jurkat_only_ACP<-df_combined$ACP[which(df_combined$color=="blue")]
      all_unique_ACP<-c(Naive_only_ACP,Jurkat_only_ACP)
      
      ALL_NaiveJurkat<-NaiveJurkat[which(NaiveJurkat$ACP %in% all_unique_ACP),]
      
      
      ALL_NaiveJurkat$color<-"black"
      for (i in c(1:length( ALL_NaiveJurkat$ACP))) {
        if ( ALL_NaiveJurkat$kmer[i]  %in%  c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")) { #UNUAR
          ALL_NaiveJurkat$color[i] <- "blue"
        }
        else if ( ALL_NaiveJurkat$kmer[i] %in%  c("GTTCA", "GTTCT", "GTTCC", "GTTCG")) {
          ALL_NaiveJurkat$color[i] <- "red"
        }
      }
      
      #ORTHOGONAL
      if(T){
        orthogonal<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/orthogonal.hg38_noPsiSeq.csv", header=T)
        orthogonal$ACP<-paste0(orthogonal$Annotation,orthogonal$chr,orthogonal$position)
        
        orthogonal_clean <- orthogonal[, c("ACP", "ortho")]
        
        # Pivot ortho values into columns
        orthogonal_new <- orthogonal_clean %>%
          mutate(ortho = replace_na(ortho, "None")) %>% # Replace NA with "None"
          distinct() %>%  # Remove duplicates
          pivot_wider(names_from = ortho, values_from = ortho, 
                      values_fn = ~ 1, values_fill = 0)
        
        
        ALL_NaiveJurkat$BID.seq<-0
        ALL_NaiveJurkat$CeU.seq<-0
        ALL_NaiveJurkat$PRAISE<-0
        ALL_NaiveJurkat$RBS.seq<-0
        
        for (i in c(1:length(ALL_NaiveJurkat))){
          
          ACP<-ALL_NaiveJurkat$ACP[i]
          
          if(ACP %in% orthogonal_new$ACP==TRUE){
            idx<- which(orthogonal_new$ACP %in% ACP )
            ALL_NaiveJurkat$CeU.seq<-orthogonal_new$CeU.seq[idx]
            ALL_NaiveJurkat$BID.seq[i]<-orthogonal_new$BID.seq[idx]
            ALL_NaiveJurkat$PRAISE[i]<-orthogonal_new$PRAISE[idx]
            ALL_NaiveJurkat$RBS.seq[i]<-orthogonal_new$RBS.seq[idx]
          }
          
        }
       
        
      
        
      }
      
      
      write.csv(ALL_NaiveJurkat,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/UniqueSitesNaiveJurkat.csv",row.names = F)
    
      }  
    
      #plots 
      if(T){
      
      #this is after the vcf file filtration  
      ALL_NaiveJurkat<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/UniqueSitesNaiveJurkat.csv",header=T)
      
      # Naive Jurkat
      ALL_NaiveJurkat$similarity<-ALL_NaiveJurkat$mm.JurkatDRS-ALL_NaiveJurkat$mm.NaiveDRS
      ALL_NaiveJurkat <- ALL_NaiveJurkat[order(ALL_NaiveJurkat$similarity ,decreasing=TRUE), ] #order
      data_matrix_high <- cbind(ALL_NaiveJurkat$mm.NaiveDRS[1:10], ALL_NaiveJurkat$mm.JurkatDRS[1:10])
      data_matrix_low <- cbind(ALL_NaiveJurkat$mm.NaiveDRS[(nrow(ALL_NaiveJurkat)-10):nrow(ALL_NaiveJurkat)],ALL_NaiveJurkat$mm.JurkatDRS[(nrow(ALL_NaiveJurkat)-10):nrow(ALL_NaiveJurkat)])
      data_matrix<-rbind(data_matrix_high,data_matrix_low)
      colnames(data_matrix) <- c("Naive T cells", "Jurkat")
      AC<-paste0(ALL_NaiveJurkat$Annotation,", ",ALL_NaiveJurkat$chr)
      ACS<-paste0(AC,", ",ALL_NaiveJurkat$position)
      ACP<-paste(ACS," ",format( round(ALL_NaiveJurkat$similarity, digits=2 ), nslow=2))
      ACP_low<-ACP[1:10]
      ACP_high<-ACP[(length(ACP)-10):length(ACP)]
      ACP<-c(ACP_low,ACP_high)
      rownames(data_matrix) <- ACP
      
      fakerow<-c(0,0)
      fakerowname<-"VOID, chrS, 11111"
      data_matrix<-rbind(data_matrix,fakerow)
      rownames(data_matrix)[nrow(data_matrix)]<-fakerowname
      fakerow<-c(100,100)
      data_matrix<-rbind(data_matrix,fakerow)
      rownames(data_matrix)[nrow(data_matrix)]<-fakerowname
      
      #plot heatmap
      if(T){
        file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/heatmapNaiveJurkat_absence.pdf"
        min.x = 0
        max.x = 100
        min.y = 0
        max.y = 800
        line.thickness = 1
        label.size = 1
        text.size = 0.8
        point.size = 1
        tck.length = 0.01
        tick.thickness = 1
        transparency = 0.2
        pdf(file.name,width = 5,height = 6)
        par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2)

        heatmap.2(data_matrix,
                  Rowv = FALSE,
                  Colv = FALSE,
                  col = mako(100),
                  dendrogram = "none",
                  key = TRUE,
                  key.title = "Color Key",
                  density.info = "none",
                  trace ="none",
                  margins=c(9,9), #0 starts from the bottom right, the first value controls if it's up/down. The second value left/right
                  cexRow=0.4,
                  cexCol = 0.8,
                  lhei = c(2,5),
                  lwid = c(3,5)
                  )
        
        dev.off()
      }
      
      #table pnas for heatmap
      if(T){ 
        pnas_filter=ALL_NaiveJurkat
        
        pnas_filter$kmer<-sapply(pnas_filter$kmer,function(x) str_replace_all(x,"T","U"))
        
        min.x = -50
        max.x = 70
        min.y = -10
        max.y = nrow(pnas_filter)+ 5
        line.thickness = 0.8
        txt.size = 0.2
        point.size = 0.3
        label.size = 0.5
        anot.pos = -20
        file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Fig2_heatmapAbsence_list.pdf"
        pdf(file=file.name)
        
        par(mfrow=c(1,1),lend=1)
        plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
             ylab = NA,xlab=NA,ann = F,axes=F)
        
        for (i in c(1:nrow(pnas_filter))) {
          ### Add the gene names
          text(x = anot.pos ,y = nrow(pnas_filter)-i, labels = pnas_filter$Annotation[i],
               srt = 0,cex = 1.2*txt.size, pos = 4, col = pnas_filter$color[i], family = "Helvetica", font = 3)
          ### Add the chromosome
          text(x = anot.pos+20 ,y = nrow(pnas_filter)-i,
               labels = pnas_filter$chr[i] ,cex = 1.2*txt.size, pos = 2, srt = 0) 
          ### Add the position
          text(x = anot.pos+30 ,y = nrow(pnas_filter)-i,
               labels = pnas_filter$position[i],cex = 1.2*txt.size, pos = 2, srt = 0) 
          ### Add the similarity
          text(x = anot.pos+40 ,y = nrow(pnas_filter)-i,
               labels = round(pnas_filter$similarity[i],digits = 2) ,cex = 1.2*txt.size, pos = 2, srt = 0) 
          ### Add kmer
          text(x = anot.pos+50 ,y = nrow(pnas_filter)-i,
               labels = pnas_filter$kmer[i] ,cex = 1.2*txt.size, pos = 2, srt = 0) 
        }
        
        ### Add the axis
        text(x = c(anot.pos,anot.pos+10,anot.pos+20,anot.pos+30, anot.pos+40), y = (nrow(pnas_filter))*c(1,1,1,1), labels = c("Gene","Chr","pos","Jurkat-Naive %mm","kmer"),
             srt = 0,cex = 1.2*txt.size, pos = 4, col = "black",family = "Helvetica", font = 3)
        
        dev.off()
      }
      
      #cds-utr unique sites
      if(T){
        
        library(rtracklayer)
        g = readGFF("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/gencode/gencode.v27.annotation.gff3")
        gff3=as.data.frame(g)
        gff3 <- gff3[which(gff3$gene_type == "protein_coding" & ((gff3$type == "CDS") | (gff3$type == "three_prime_UTR") | (gff3$type == "five_prime_UTR") | (gff3$type == "start_codon") | (gff3$type == "stop_codon") | (gff3$type == "stop_codon_redefined_as_selenocysteine"))), c("type", "seqid","start", "end", "gene_name")]
        
        #for unique sites
        ALL_NaiveJurkat<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/UniqueSitesNaiveJurkat.csv",header=T)
        ALL_NaiveJurkat_idx<-df_combined$ACP[which(df_combined$color!="cyan")] #orange and blue, unique
        ALL_NaiveJurkat<-ALL_NaiveJurkat[which(ALL_NaiveJurkat$ACP %in% ALL_NaiveJurkat_idx),]
        file.name="/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/UniqueSitesNaiveJurkatCDSUTR.csv"
        #for common sites
        # ALL_NaiveJurkat<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/CommonSitesNaiveJurkat.csv",header=T)
        # ALL_NaiveJurkat_idx<-df_combined$ACP[which(df_combined$color=="cyan")] #cyan, shared
        # ALL_NaiveJurkat<-ALL_NaiveJurkat[which(ALL_NaiveJurkat$ACP %in% ALL_NaiveJurkat_idx),]
        # file.name="/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/CommonSitesNaiveJurkatCDSUTR.csv"
        
        #utr-cds attributions
        if(T){
          library(rtracklayer)
          ALL_NaiveJurkat$CDS = F
          ALL_NaiveJurkat$Three_prime_UTR = F
          ALL_NaiveJurkat$Five_Prime_UTR = F
          ALL_NaiveJurkat$start_codon = F
          ALL_NaiveJurkat$stop_codon = F
          ALL_NaiveJurkat$stop_codon_redefined_as_selenocysteine = F
          
          for (row in 1:nrow(ALL_NaiveJurkat)){
            
            type = gff3$type[which(gff3$start <= ALL_NaiveJurkat$position[row] & 
                                     gff3$end >= ALL_NaiveJurkat$position[row] &
                                     gff3$seqid == ALL_NaiveJurkat$chr[row])]
            
            if (length(type)>0){
              if ("CDS" %in% type) {
                ALL_NaiveJurkat$CDS[row] = T
              }
              if ("three_prime_UTR" %in% type) {
                ALL_NaiveJurkat$Three_prime_UTR[row] = T
              }
              if ("five_prime_UTR" %in% type) {
                ALL_NaiveJurkat$Five_Prime_UTR[row] = T
              }
              if ("start_codon" %in% type) {
                ALL_NaiveJurkat$start_codon[row] = T
              }
              if ("stop_codon" %in% type) {
                ALL_NaiveJurkat$stop_codon[row] = T
              }
              if ("stop_codon_redefined_as_selenocysteine" %in% type) {
                ALL_NaiveJurkat$stop_codon_redefined_as_selenocysteine[row] = T
              }
            }
          }
          
          
          start.codon <- ALL_NaiveJurkat[which(ALL_NaiveJurkat$start_codon == 1),]
          
          all.false <- ALL_NaiveJurkat[which(ALL_NaiveJurkat$CDS == 0 &
                                                 ALL_NaiveJurkat$Three_prime_UTR == 0 &
                                                 ALL_NaiveJurkat$Five_Prime_UTR == 0 &
                                                 ALL_NaiveJurkat$start_codon == 0 &
                                                 ALL_NaiveJurkat$stop_codon == 0),]
          
        }
        write.csv(ALL_NaiveJurkat,file.name,row.names = F)  
        
        #plot unique sites 
        #plot cds utr, exclude positions from plot if found in multiple regions as they are on different isoforms
        if(T){
          
          ALL_NaiveJurkat<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/UniqueSitesNaiveJurkatCDSUTR.csv",header = T)
          
          min.x = 0
          max.x = 6
          min.y = 0
          max.y = 100
          label.size = 1
          text.size = 0.8
          point.size = 1
          tck.length = 0.01
          tick.thickness = 1
          transparency = 0.2
          file.name ="/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Unique_sites_UTR-CDS.pdf"
          pdf(file.name, width = 6, height = 6, family = "Helvetica")
          plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
               ylab = NA,xlab=NA,ann = F,axes=F)
          #box(lwd = line.thickness)
          #### Add the lines
          par(mfrow = c(1, 2)) 
          
          #normalizing for number of reads. NO NORMALIZATIONS IN THIS CASE
          max_reads= 1#max(c(sum(ALL_NaiveJurkat$total_TRUB1_KD),sum(PUS7KD_sites$total_PUS7_KD)))
          factor_norm=1#000
          spacer=0.1
          
          ALL_NaiveJurkat_JurkatOnly<-ALL_NaiveJurkat[which(ALL_NaiveJurkat$mm.JurkatDRS>30),]
          ALL_NaiveJurkat_NaiveOnly<-ALL_NaiveJurkat[which(ALL_NaiveJurkat$mm.NaiveDRS>30),]

          Prop<-c(length(which(ALL_NaiveJurkat_JurkatOnly$stop_codon==T))/max_reads*factor_norm,
                  length(which(ALL_NaiveJurkat_JurkatOnly$start_codon==T))/max_reads*factor_norm,
                  length(which(ALL_NaiveJurkat_JurkatOnly$Three_prime_UTR==T & ALL_NaiveJurkat_JurkatOnly$Five_Prime_UTR==F & ALL_NaiveJurkat_JurkatOnly$CDS==F ))/max_reads*factor_norm,
                  length(which(ALL_NaiveJurkat_JurkatOnly$CDS==T & ALL_NaiveJurkat_JurkatOnly$Three_prime_UTR==F & ALL_NaiveJurkat_JurkatOnly$Five_Prime_UTR==F ))/max_reads*factor_norm,
                  length(which(ALL_NaiveJurkat_JurkatOnly$Five_Prime_UTR==T & ALL_NaiveJurkat_JurkatOnly$CDS==F & ALL_NaiveJurkat_JurkatOnly$Three_prime_UTR==F))/max_reads*factor_norm)

          # Prop<-c(length(which(ALL_NaiveJurkat_JurkatOnly$stop_codon==T))/max_reads*factor_norm,
          #         length(which(ALL_NaiveJurkat_JurkatOnly$start_codon==T))/max_reads*factor_norm,
          #         length(which(ALL_NaiveJurkat_JurkatOnly$Three_prime_UTR==T))/max_reads*factor_norm,
          #         length(which(ALL_NaiveJurkat_JurkatOnly$CDS==T ))/max_reads*factor_norm,
          #         length(which(ALL_NaiveJurkat_JurkatOnly$Five_Prime_UTR==T))/max_reads*factor_norm)
          
          labels<-c("STOP","START","3'UTR","CDS","5'UTR")
          numbers<-Prop
          pie(Prop, labels = paste0(labels,"-",numbers), border="white", col= c("red","green","lightblue", "royalblue4","royalblue"),cex = 0.7 ) 
    
          
          Prop<-c(length(which(ALL_NaiveJurkat_NaiveOnly$stop_codon==T))/max_reads*factor_norm,
                  length(which(ALL_NaiveJurkat_NaiveOnly$start_codon==T))/max_reads*factor_norm,
                  length(which(ALL_NaiveJurkat_NaiveOnly$Three_prime_UTR==T & ALL_NaiveJurkat_NaiveOnly$Five_Prime_UTR==F & ALL_NaiveJurkat_NaiveOnly$CDS==F ))/max_reads*factor_norm,
                  length(which(ALL_NaiveJurkat_NaiveOnly$CDS==T & ALL_NaiveJurkat_NaiveOnly$Three_prime_UTR==F & ALL_NaiveJurkat_NaiveOnly$Five_Prime_UTR==F ))/max_reads*factor_norm,
                  length(which(ALL_NaiveJurkat_NaiveOnly$Five_Prime_UTR==T & ALL_NaiveJurkat_NaiveOnly$CDS==F & ALL_NaiveJurkat_NaiveOnly$Three_prime_UTR==F))/max_reads*factor_norm)

          # Prop<-c(length(which(ALL_NaiveJurkat_NaiveOnly$stop_codon==T))/max_reads*factor_norm,
          #         length(which(ALL_NaiveJurkat_NaiveOnly$start_codon==T))/max_reads*factor_norm,
          #         length(which(ALL_NaiveJurkat_NaiveOnly$Three_prime_UTR==T ))/max_reads*factor_norm,
          #         length(which(ALL_NaiveJurkat_NaiveOnly$CDS==T ))/max_reads*factor_norm,
          #         length(which(ALL_NaiveJurkat_NaiveOnly$Five_Prime_UTR==T))/max_reads*factor_norm)
          
          labels<-c("STOP","START","3'UTR","CDS","5'UTR")
          numbers<-Prop
          pie(Prop, labels = paste0(labels,"-",numbers), border="white",col= c("red","green","burlywood1","darkorange3","orange"),cex = 0.7)
          
          dev.off() 
        }
        
      }
      
      #how many sites are completely absent (%mm==0)? overlayed barplot, supplementary fig
      if(T){
        
        # Define %mm intervals and calculate the number of sites in each
        interval_breaks <- seq(0, 10, by = 1)
        interval_labels <- c("0",paste0(interval_breaks[-length(interval_breaks)], "≤mm<", interval_breaks[-1]))
        
        data <- data.frame(
          CellType = c(rep("Naive", 11), rep("Jurkat", 11)),
          ExpressionLevel = rep(interval_labels, 2),
          Sites = c( length(ALL_NaiveJurkat$ACP[which(ALL_NaiveJurkat$mm.NaiveDRS==0)]),
                     sapply(1:10, function(i) { nrow(ALL_NaiveJurkat %>% filter(mm.NaiveDRS >= interval_breaks[i] & mm.NaiveDRS < interval_breaks[i+1]))}) ,
                     length(ALL_NaiveJurkat$ACP[which(ALL_NaiveJurkat$mm.JurkatDRS==0)]),
                     sapply(1:10, function(i) { nrow(ALL_NaiveJurkat %>% filter(mm.JurkatDRS >= interval_breaks[i] & mm.JurkatDRS < interval_breaks[i+1]))})
                    ),
          TotalSites = c(
            rep(nrow(ALL_NaiveJurkat %>% filter(mm.NaiveDRS <= 10)), 11),
            rep(nrow(ALL_NaiveJurkat %>% filter(mm.JurkatDRS <= 10)), 11)
          )
        )
        
        data$Sites[2]<-abs(data$Sites[1]-data$Sites[2])
        data$Sites[13]<-abs(data$Sites[12]-data$Sites[13])
        
        # Calculate percentage
        data$Percentage <- (data$Sites / data$TotalSites) * 100
        
        # Assign interval index (1 to 10) for color mapping within each CellType
        data$Interval <- rep(1:11, 2)
        
        # Generate a rainbow palette starting with white for the first interval
        rainbow_colors <- c("white", viridis(10, begin = 0.95, end = 0))
        
        # Create the plot using the custom rainbow palette with white as the first color
        p <- ggplot(data, aes(x = CellType, y = Percentage, fill = factor(Interval))) +
          geom_bar(stat = "identity", position = "stack", color = "black") +
          scale_fill_manual(values = rainbow_colors, 
                            labels = interval_labels,
                            guide = guide_legend(title = "%mm Interval")) +
          labs(title = "Percentage of Sites in Different %mm Ranges for Jurkat and Naive Cells",
               y = "Percentage (%)",
               x = "Cell Type") +
          theme_minimal() +
          theme(legend.title = element_blank())
        print(p)
        ggsave("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/expression_levels_unique_sites_v2.pdf", p, width=8, height =6)

        #old version =0 and >0
        if(T){
          data <- data.frame(
            CellType = c("Naive","Naive","Jurkat","Jurkat"),
            ExpressionLevel =  c("Expression=0", "Expression<10", "Expression=0", "Expression<10"),
            Sites = c(nrow(ALL_NaiveJurkat[which(ALL_NaiveJurkat$mm.NaiveDRS==0),]), nrow(ALL_NaiveJurkat[which(ALL_NaiveJurkat$mm.NaiveDRS<=10),])-nrow(ALL_NaiveJurkat[which(ALL_NaiveJurkat$mm.NaiveDRS==0),]), nrow(ALL_NaiveJurkat[which(ALL_NaiveJurkat$mm.JurkatDRS==0),]), nrow(ALL_NaiveJurkat[which(ALL_NaiveJurkat$mm.JurkatDRS<=10),])-nrow(ALL_NaiveJurkat[which(ALL_NaiveJurkat$mm.JurkatDRS==0),]) ),
            TotalSites = c(nrow(ALL_NaiveJurkat[which(ALL_NaiveJurkat$mm.JurkatDRS>30),]),nrow(ALL_NaiveJurkat[which(ALL_NaiveJurkat$mm.JurkatDRS>30),]),nrow(ALL_NaiveJurkat[which(ALL_NaiveJurkat$mm.NaiveDRS>30),]),nrow(ALL_NaiveJurkat[which(ALL_NaiveJurkat$mm.NaiveDRS>30),]))
          )
          
          # Compute percentages
          data$Percentage <- (data$Sites / data$TotalSites) * 100
          
          
          # Assign colors
          colors <- c("Naive.Expression=0" = "white", "Naive.Expression<10" = "orange",
                      "Jurkat.Expression=0" = "white", "Jurkat.Expression<10" = "royalblue")
          # Define border colors
          border_colors <- c("Naive.Expression=0" = "orange", "Naive.Expression<10" = "orange",
                             "Jurkat.Expression=0" = "royalblue", "Jurkat.Expression<10" = "royalblue")
          
          # Create the plot
          p<-ggplot(data, aes(x = CellType, y = Percentage, fill = interaction(CellType, ExpressionLevel) )) +
            geom_bar(stat = "identity", position = "stack", color=border_colors) +
            scale_fill_manual(values = colors) +
            labs(title = "Expression Levels in Jurkat and Naive Cells",
                 y = "Percentage (%)",
                 x = "Cell Type") +
            theme_minimal()
          
          ggsave("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/expression_levels_unique_sites.pdf", p, width=8, height =6)
        }
        
        
        
      }
      
      }
      
    }
  }
  

  #cds-utr for Stuart
  if(T){
 
    #tot number of DRS reads>20
    Naive<-Naive_orig[which(Naive_orig$N_reads_NaiveDRS>20),]
    Jurkat<-Jurkat_orig[which(Jurkat_orig$N_reads_JurkatDRS>20),]
    #tot number of IVT reads<10
    Naive<-Naive[which(Naive$N_reads_NaiveIVT>10),]
    Jurkat<-Jurkat[which(Jurkat$N_reads_JurkatIVT>10),]
    #mm IVT<10
    Naive<-Naive[which(Naive$mm.NaiveIVT<10),]
    Jurkat<-Jurkat[which(Jurkat$mm.JurkatIVT<10),]
    
    #mod type psi
    Naive$mod_type<-"psi"
    Jurkat$mod_type<-"psi"
    #add fake color between yellow and purple
    Naive$color_rgb<-"53,92,140"
    Jurkat$color_rgb<-"71,21,103"
    
    
    Naive$color<-"black"
    for (i in c(1:length( Naive$Annotation))) {
      if ( Naive$kmer[i]  %in%  c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")) { #UNUAR
        Naive$color[i] <- "blue"
      }
      else if ( Naive$kmer[i] %in%  c("GTTCA", "GTTCT", "GTTCC", "GTTCG")) {
        Naive$color[i] <- "red"
      }
    }
    
    Jurkat$color<-"black"
    for (i in c(1:length( Jurkat$Annotation))) {
      if ( Jurkat$kmer[i]  %in%  c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")) { #UNUAR
        Jurkat$color[i] <- "blue"
      }
      else if ( Jurkat$kmer[i] %in%  c("GTTCA", "GTTCT", "GTTCC", "GTTCG")) {
        Jurkat$color[i] <- "red"
      }
    }
    
    NaiveStu<-data.frame(Naive$chr,as.numeric(Naive$position),as.numeric(Naive$position+1), Naive$mod_type, as.numeric(Naive$p.value.NaiveDRS), Naive$strand, as.numeric(Naive$position), as.numeric(Naive$position+1),Naive$color_rgb, Naive$kmer, Naive$Annotation, Naive$mm.NaiveDRS, Naive$N_reads_NaiveDRS)
    write.table(NaiveStu,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Naive_drs_file.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)  
    JurkatStu<-data.frame(Jurkat$chr, as.numeric(Jurkat$position), as.numeric(Jurkat$position+1), Jurkat$mod_type, as.numeric(Jurkat$p.value.JurkatDRS), Jurkat$strand, as.numeric(Jurkat$position), as.numeric(Jurkat$position+1), Jurkat$color_rgb, Jurkat$kmer, Jurkat$Annotation, Jurkat$mm.JurkatDRS, Jurkat$N_reads_JurkatDRS)
    write.table(JurkatStu,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Jurkat_drs_file.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)  
  
    Naive_hc<-Naive[which(Naive$mm.NaiveDRS>=30),]
    Jurkat_hc<-Jurkat[which(Jurkat$mm.JurkatDRS>=30),]
    Naive_hc<-Naive_hc[which(Naive_hc$p.value.NaiveDRS<0.001),]
    Jurkat_hc<-Jurkat_hc[which(Jurkat_hc$p.value.JurkatDRS<0.001),]

    
    #intersect these hc bed files with bedtools intersect -wa -wb -a gft/file/path -b bed/file/path
    NaiveStu<-data.frame(Naive_hc$chr,as.numeric(Naive_hc$position),as.numeric(Naive_hc$position+1), Naive_hc$mod_type, as.numeric(Naive_hc$p.value.NaiveDRS), Naive_hc$strand, as.numeric(Naive_hc$position), as.numeric(Naive_hc$position+1),Naive_hc$color_rgb, Naive_hc$kmer, Naive_hc$Annotation, Naive_hc$mm.NaiveDRS, Naive_hc$N_reads_NaiveDRS)
    write.table(NaiveStu,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Naive_hc_file.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)  
    JurkatStu<-data.frame(Jurkat_hc$chr, as.numeric(Jurkat_hc$position), as.numeric(Jurkat_hc$position+1), Jurkat_hc$mod_type, as.numeric(Jurkat_hc$p.value.JurkatDRS), Jurkat_hc$strand, as.numeric(Jurkat_hc$position), as.numeric(Jurkat_hc$position+1), Jurkat_hc$color_rgb, Jurkat_hc$kmer, Jurkat_hc$Annotation, Jurkat_hc$mm.JurkatDRS, Jurkat_hc$N_reads_JurkatDRS)
    write.table(JurkatStu,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Jurkat_hc_file.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)  
    
   
    }
  
  
  ###
  ##################Hypermodified TYPE2 
  if(T){
    
    #load hypermodified
    JurkatHyp<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/pysam/JurkatIVTless10mm40_hypermod.csv",header = T)  
    JurkatHyp$ACP<-paste0(JurkatHyp$Annotation,JurkatHyp$chr,JurkatHyp$position)
    JurkatHyp<-JurkatHyp[which(JurkatHyp$N_reads_JurkatIVT>10),]
    JurkatHyp<-JurkatHyp[which(JurkatHyp$N_reads_JurkatDRS>20),]
    JurkatHyp<-JurkatHyp[which(JurkatHyp$mm.JurkatIVT<10),]
    JurkatHyp<-JurkatHyp[which(JurkatHyp$mm.JurkatDRS>=30),]
    NaiveHyp<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/pysam/NaiveIVTless10mm40_hypermod.csv",header = T) 
    NaiveHyp$ACP<-paste0(NaiveHyp$Annotation,NaiveHyp$chr,NaiveHyp$position)
    NaiveHyp<-NaiveHyp[which(NaiveHyp$N_reads_NaiveIVT>10),]
    NaiveHyp<-NaiveHyp[which(NaiveHyp$N_reads_NaiveDRS>20),]
    NaiveHyp<-NaiveHyp[which(NaiveHyp$mm.NaiveIVT<10),]
    NaiveHyp<-NaiveHyp[which(NaiveHyp$mm.NaiveDRS>=30),]
    ActivatedHyp<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/pysam/ActivatedIVTless10mm40_hypermod.csv",header = T)
    ActivatedHyp$ACP<-paste0(ActivatedHyp$Annotation,ActivatedHyp$chr,ActivatedHyp$position)
    ActivatedHyp<-ActivatedHyp[which(ActivatedHyp$N_reads_NaiveIVT>10),]
    ActivatedHyp<-ActivatedHyp[which(ActivatedHyp$N_reads_ActivatedDRS>20),]
    ActivatedHyp<-ActivatedHyp[which(ActivatedHyp$mm.NaiveIVT<10),]
    ActivatedHyp<-ActivatedHyp[which(ActivatedHyp$mm.ActivatedDRS>=30),]
    
    #search in non hypermod sites too some hypermod type 2 sites
    NaiveDRS_filter=read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveDRS_filter.csv",header=T)
    NaiveDRS_filter<-NaiveDRS_filter[which(NaiveDRS_filter$mm.NaiveIVT<=10),] #doubleckeck these!!!!
    NaiveDRS_filter<-NaiveDRS_filter[which(NaiveDRS_filter$N_reads_NaiveIVT>10),]
    NaiveDRS_filter$ACP<-paste0(NaiveDRS_filter$Annotation,NaiveDRS_filter$chr,NaiveDRS_filter$position)
    ActivatedDRS_filter=read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/ActivatedDRS_filter.csv",header=T)
    ActivatedDRS_filter<-ActivatedDRS_filter[which(ActivatedDRS_filter$mm.NaiveIVT<=10),]
    ActivatedDRS_filter<-ActivatedDRS_filter[which(ActivatedDRS_filter$N_reads_NaiveIVT>10),]
    ActivatedDRS_filter$ACP<-paste0(ActivatedDRS_filter$Annotation,ActivatedDRS_filter$chr,ActivatedDRS_filter$position)
    JurkatDRS_filter=read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/JurkatDRS_filter.csv",header=T)
    JurkatDRS_filter<-JurkatDRS_filter[which(JurkatDRS_filter$mm.JurkatIVT<=10),]
    JurkatDRS_filter<-JurkatDRS_filter[which(JurkatDRS_filter$N_reads_JurkatIVT>10),]
    JurkatDRS_filter$ACP<-paste0(JurkatDRS_filter$Annotation,JurkatDRS_filter$chr,JurkatDRS_filter$position)
    

    #hypermodified
    #commonsites<-intersect(intersect(ActivatedHyp$ACP,JurkatHyp$ACP),NaiveHyp$ACP)
    commonsites<-intersect(JurkatHyp$ACP,NaiveHyp$ACP)
    h2_common<-as.data.frame(table(NaiveHyp$Annotation[which(NaiveHyp$ACP %in% commonsites)]))
    h2_common<-h2_common[which(h2_common$Freq>1),]
    #h2_common<-h2_common[order(h2_common$Freq,decreasing = FALSE),]
    h2_common$ACPs<-""
    for (i in c(1:nrow(h2_common))){
      
      pattern <- as.character(h2_common$Var1[i])
      idx<-grepl(pattern,commonsites,fixed=TRUE);
      idx_gene<-which(idx==TRUE)
      names<-paste(commonsites[idx_gene],collapse=" / ")
      
      h2_common$ACP[i]<-names
      
      
    }
    h2_common_hyp<-h2_common
    
    #not hypermodified
    #commonsites<-intersect(intersect(ActivatedDRS_filter$ACP,JurkatDRS_filter$ACP),NaiveDRS_filter$ACP)
    commonsites<-intersect(JurkatDRS_filter$ACP,NaiveDRS_filter$ACP)
    h2_common<-as.data.frame(table(NaiveDRS_filter$Annotation[which(NaiveDRS_filter$ACP %in% commonsites)]))
    h2_common<-h2_common[which(h2_common$Freq>1),]
    #h2_common<-h2_common[order(h2_common$Freq,decreasing = FALSE),]
    h2_common$ACPs<-""
    for (i in c(1:nrow(h2_common))){
      
      pattern <- as.character(h2_common$Var1[i])
      idx<-grepl(pattern,commonsites,fixed=TRUE);
      idx_gene<-which(idx==TRUE)
      names<-paste(commonsites[idx_gene],collapse=" / ")
      
      h2_common$ACP[i]<-names
      
      
    }
    
  }
  
}

 ### plots 6cell line GO terms 
 if(T){
  
  GO_RNAmod<-read_xlsx("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/GO_analysis/GO_term_summary_RNAModification.xlsx")
  GO_RNAmod_gene<-toupper(unique(GO_RNAmod$Symbol))
 
  A549   <-  read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/dhu_A549/guppy6.4.2/WT_enriched_pan_filt.csv", header = T)
  NTERA  <-  read_excel("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/GO_analysis/6cellLines_values.xlsx", sheet = "NTERA")
  JurkatNaive <-  read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/CommonSitesNaiveJurkatCDSUTR.csv",header=T)
  Naive  <- read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveJurkat_panIVT/Naive_panIVT/init_gene_pileup/finalized_pileup/Naive_filtered_pan.csv",header = T)
  Jurkat <- read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveJurkat_panIVT/Jurkat_panIVT/init_gene_pileup/finalized_pileup/Jurkat_filtered_pan.csv",header=T)

  Naive<-Naive[which(Naive$N_reads_NaiveDRS>20 & Naive$p.value.NaiveDRS<0.001),]
  Naive<-Naive[which(Naive$N_reads_NaiveIVT>10 & Naive$mm.NaiveIVT<10),]
  Naive<-Naive[which(Naive$mm.NaiveDRS>10),]
  
  Jurkat<-Jurkat[which(Jurkat$N_reads_JurkatDRS>20 & Jurkat$p.value.JurkatDRS<0.001),]
  Jurkat<-Jurkat[which(Jurkat$N_reads_JurkatIVT>10 & Jurkat$mm.JurkatIVT<10),]
  Jurkat<-Jurkat[which(Jurkat$mm.JurkatDRS>10),]
  
  A549<-A549[which(A549$N_reads_WT>20 & A549$p.value.WT<0.001),]
  A549<-A549[which(A549$N_reads_IVT >10 & A549$mm.IVT <10),]
  A549<-A549[which(A549$mm.WT>10),]
  
  NTERA<-NTERA[which(NTERA$N_reads_Direct>20 & NTERA$p.value.Direct<0.001),]
  NTERA<-NTERA[which(NTERA$N_reads_IVT>10 & NTERA$mm.IVT<10),]
  NTERA<-NTERA[which(NTERA$mm.Direct>10),]
  
  RNA_mod_A549<-A549[which(A549$Annotation %in% GO_RNAmod_gene),]
  RNA_mod_NTERA<-NTERA[which(NTERA$Annotation %in% GO_RNAmod_gene),]
  RNA_mod_JurkatNaive<-JurkatNaive[which(JurkatNaive$Annotation %in% GO_RNAmod_gene),]
  RNA_mod_Naive<-Naive[which(Naive$Annotation %in% GO_RNAmod_gene),]
  RNA_mod_Jurkat<-Jurkat[which(Jurkat$Annotation %in% GO_RNAmod_gene),]
  
  #keep only one occurrence of each gene, with the highest %mm 
  RNA_mod_A549_filt <- RNA_mod_A549 %>%
    group_by(Annotation) %>%
    slice_max(order_by = mm.WT, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  RNA_mod_NTERA_filt <- RNA_mod_NTERA %>%
    group_by(Annotation) %>%
    slice_max(order_by = mm.Direct, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  RNA_mod_Jurkat_filt <- RNA_mod_Jurkat %>%
    group_by(Annotation) %>%
    slice_max(order_by = mm.JurkatDRS, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  RNA_mod_Naive_filt <- RNA_mod_Naive %>%
    group_by(Annotation) %>%
    slice_max(order_by = mm.NaiveDRS, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  #create a dataframe with all RNA mod genes and fill it with values from J,N,A549,NTERA
  df<-data.frame("GO_gene"=GO_RNAmod_gene,"Jurkatmm"= -1,"Naivemm"= -1,"A549mm"= -1,"NTERAmm"= -1)
  
  for (i in c(1:length(GO_RNAmod_gene))){
    gene<-df$GO_gene[i]
    
    # row_on_JN<-which(JurkatNaive$Annotation %in% gene)
    # if ( !is_empty(row_on_JN) ){
    #   df$Jurkatmm[i]<- JurkatNaive$mm.JurkatDRS[row_on_JN]
    #   df$Naivemm[i]<-JurkatNaive$mm.NaiveDRS[row_on_JN]
    # } 
    # 
    
    row_on_J<-which(RNA_mod_Jurkat_filt$Annotation == gene)
    if (!is_empty(row_on_J) ) {
      df$Jurkatmm[i]<- RNA_mod_Jurkat_filt$mm.JurkatDRS[row_on_J]
    }  
      
    row_on_N<-which(RNA_mod_Naive_filt$Annotation == gene)  
    if (!is_empty(row_on_N) ) {
      df$Naivemm[i]<- RNA_mod_Naive_filt$mm.NaiveDRS[row_on_N]
    }  
      
      
    row_on_A549<-which(RNA_mod_A549_filt$Annotation == gene)
    if (!is_empty(row_on_A549) ) {
      df$A549mm[i]<- RNA_mod_A549_filt$mm.WT[row_on_A549]
    }
      
    row_on_NTERA<-which(RNA_mod_NTERA_filt$Annotation == gene)
    if (!is_empty(row_on_NTERA) ) {
      df$NTERAmm[i]<- RNA_mod_NTERA_filt$mm.Direct[row_on_NTERA]
    }
    
  }
  
  

  # Order from rows that the transcript existing in all libraries to rows with unexistent transcripts
  df <- df %>%
    # Create a new column for the count of positive values across the last 4 columns
    mutate(PosCount = rowSums(across(-GO_gene) >= 0)) %>%
    # Arrange the dataframe based on the PosCount, with higher values first
    arrange(desc(PosCount), GO_gene)
  
  # Drop the 'PosCount' column if it's not needed
  df$PosCount<-NULL
  
  
  
  # Set the GO_gene column as row names
  rownames(df) <- df$GO_gene
  df_f <- df[, -1]  # Remove the GO_gene column
  
  pdf("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/heatmap_output_Supplementary_AllJurkatNaive.pdf", width = 6, height = 10)  
  # Plot heatmap
  pheatmap(df_f, 
           cluster_rows = F,       # Cluster rows (genes)
           cluster_cols = F,       # Cluster columns (conditions)
           #scale = "row",             # Scale data by rows (z-score scaling)
           color = c("white","black", rainbow(99, start = 0, end = 0.8)),  # Choose your color scale
           show_rownames = TRUE,      # Show row names (genes)
           show_colnames = TRUE,       # Show column names (conditions)
           fontsize_row =5,
           main = "Psi mm % on transcrips from GO:0009451-RNA modification")

    dev.off()
  
  
  
}

 #common sites orthogonal validation and 6cell lines confirmation
 if(T){
  
  #vcf filtered
  ALL_NaiveJurkat<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/CommonSitesNaiveJurkatCDSUTR.csv",header=T)
  
  #vcf check, did once
  if(F){
    
    # Search my sites in the SNV dataset
    rna_modifications <- data.frame(
      chr = ALL_NaiveJurkat$chr,  # Replace with your actual chromosome data
      position = ALL_NaiveJurkat$position,  # Replace with your actual positions
      annotation = ALL_NaiveJurkat$Annotation  # Replace with your actual annotations
    )
    
    rna_mod_gr <- GRanges(
      seqnames = rna_modifications$chr,
      ranges = IRanges(start = rna_modifications$position, end = rna_modifications$position),
      annotation = rna_modifications$annotation
    )
    
    vcf_file <- "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/vcf_GRCh38p10/Jurkat_pileup.vcf"  # Replace with the path to your VCF file
    vcf <- readVcf(vcf_file, "hg18") 
    snv_gr_Jurkat <- rowRanges(vcf)
    vcf_file <- "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/vcf_GRCh38p10/Naive_pileup.vcf"  # Replace with the path to your VCF file
    vcf <- readVcf(vcf_file, "hg18") 
    snv_gr_Naive <- rowRanges(vcf)
    #concatenating both SNV from Jurkat and Naive in one large object
    snv_gr<-c(snv_gr_Jurkat,snv_gr_Naive)
    
    # Get the indices of the overlapping sites between my sites and the snv_gr containing all the SNV indexes
    overlaps <- findOverlaps(rna_mod_gr, snv_gr)
    overlapping_indices <- queryHits(overlaps)
    overlapping_sites <- rna_mod_gr[overlapping_indices]
    
    # Get all the non-overlapping sites that are NOT SNVs
    all_indices <- seq_len(nrow(rna_modifications))
    overlapping_indices <- queryHits(overlaps)
    non_overlapping_indices <- setdiff(all_indices, overlapping_indices)
    non_overlapping_sites <- rna_modifications[non_overlapping_indices, ]
    non_overlapping_ACP<- paste0(non_overlapping_sites$annotation, non_overlapping_sites$chr,non_overlapping_sites$position)
    #removing all SNVs from the list
    ALL_NaiveJurkat<-ALL_NaiveJurkat[which(ALL_NaiveJurkat$ACP %in% non_overlapping_ACP),]
  }
  
  #ORTHOGONAL
  if(T){
    orthogonal<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/orthogonal.hg38_noPsiSeq.csv", header=T)
    orthogonal$ACP<-paste0(orthogonal$Annotation,orthogonal$chr,orthogonal$position)
    
    orthogonal_clean <- orthogonal[, c("ACP", "ortho")]
    
    # Pivot ortho values into columns
    orthogonal_new <- orthogonal_clean %>%
      mutate(ortho = replace_na(ortho, "None")) %>% # Replace NA with "None"
      distinct() %>%  # Remove duplicates
      pivot_wider(names_from = ortho, values_from = ortho, 
                  values_fn = ~ 1, values_fill = 0)
    
    
    ALL_NaiveJurkat$BID.seq<-0
    ALL_NaiveJurkat$CeU.seq<-0
    ALL_NaiveJurkat$PRAISE<-0
    ALL_NaiveJurkat$RBS.seq<-0
    
    for (i in c(1:length(ALL_NaiveJurkat))){
      
      ACP<-ALL_NaiveJurkat$ACP[i]
      
      if(ACP %in% orthogonal_new$ACP==TRUE){
        idx<- which(orthogonal_new$ACP %in% ACP )
        ALL_NaiveJurkat$CeU.seq<-orthogonal_new$CeU.seq[idx]
        ALL_NaiveJurkat$BID.seq[i]<-orthogonal_new$BID.seq[idx]
        ALL_NaiveJurkat$PRAISE[i]<-orthogonal_new$PRAISE[idx]
        ALL_NaiveJurkat$RBS.seq[i]<-orthogonal_new$RBS.seq[idx]
      }
      
    }
    
  }
  
  #6 Cell lines
  if(T){
    A549   <-  read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/dhu_A549/guppy6.4.2/WT_enriched_pan_filt.csv", header = T)
    NTERA  <-  read_excel("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/GO_analysis/6cellLines_values.xlsx", sheet = "NTERA")
    HeLa  <-  read_excel("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/GO_analysis/6cellLines_values.xlsx", sheet = "HeLa")
    HepG2  <-  read_excel("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/GO_analysis/6cellLines_values.xlsx", sheet = "HepG2")
    SHSY5Y  <-  read_excel("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/GO_analysis/6cellLines_values.xlsx", sheet = "SH-SY5Y")
    
    ACP_common_JN<-ALL_NaiveJurkat$ACP
    
    A549_NJ<-A549[which(A549$ACP %in% ACP_common_JN),]
    NTERA_NJ<-NTERA[which(NTERA$ACP %in% ACP_common_JN),]
    HeLa_NJ<-HeLa[which(HeLa$ACP %in% ACP_common_JN),]
    HepG2_NJ<-HepG2[which(HepG2$ACP %in% ACP_common_JN),]
    SHSY5Y_NJ<-SHSY5Y[which(SHSY5Y$ACP %in% ACP_common_JN),]
  }
  
}
}

### plots figure uniquely expressed transcripts in Naive or Jurkat
if(T){
  
  #nanocount to quantify transcript expression and identify singular transcripts
  if(T){

  NaiveNano1<-read.csv("/Dropbox/NanoCountN1/Naive_DRS.gencodev43.nanocount.N1.tsv", sep = '\t',header = T)
  JurkatNano1<-read.csv("/Dropbox/NanoCountN1/Jurkat_DRS.gencodev43.nanocount.N1.tsv", sep = '\t',header = T)
  
  NanoJurkatNaive<-merge(NaiveNano1,JurkatNano1, by="transcript_name", all=TRUE)

  i<-intersect(NaiveNano1$transcript_name,JurkatNano1$transcript_name)
  N<-setdiff(NaiveNano1$transcript_name,JurkatNano1$transcript_name)
  J<-setdiff(JurkatNano1$transcript_name,NaiveNano1$transcript_name)
  
  tot<-length(i)+length(J)+length(N)
  coexpr<-length(i)/tot
  exprN<-length(N)/tot
  exprJ<-length(J)/tot
  
  NanoJurkatNaive[is.na(NanoJurkatNaive)]<-0

  #heatmap with all transcripts
  if(T){
    NanoJurkatNaive$ordering <- with(NanoJurkatNaive, 
                        ifelse(tpm.x > 1 & tpm.y == 0, 1, 
                               ifelse(tpm.x > 1 & tpm.y > 1, 2, 
                                      ifelse(tpm.y > 1 & tpm.x == 0, 4, 3))))
    
    # Order the dataframe based on the custom ordering
    NanoJurkatNaive <- NanoJurkatNaive[order(NanoJurkatNaive$ordering), ]
    
    # Recreate the heatmap matrix with the new order
    heatmap_data <- NanoJurkatNaive[, c("tpm.x", "tpm.y")]
    rownames(heatmap_data) <- gsub("\\|.*", "", NanoJurkatNaive$transcript_name)  # Use first part of transcript_name as row names
    # ordered_indices <- order(-heatmap_data$tpm.x)  # Get the indices for descending order of tpm.x
    # heatmap_data <- heatmap_data[ordered_indices, ]
    heatmap_data<-heatmap_data
    heatmap_data<-heatmap_data[which(heatmap_data$tpm.x<1000 & heatmap_data$tpm.y<1000),]
    heatmap_matrix <- as.matrix(heatmap_data)
    colnames(heatmap_matrix) <- c("Naive","Jurkat")
  
    
    pdf("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/tpm_heatmap.pdf", width = 5, height = 20) 
      
      heatmap_colors <- c("white", colorRampPalette(c("lightblue", "blue"))(1000))
      Heatmap(heatmap_matrix, 
              name = "TPM",           # Name of the color legend
              col = heatmap_colors,          # Apply custom color function
              row_order = NULL,       # Keep rows in the current order
              row_names_gp = gpar(fontsize = 1),      # Customize row labels
              cluster_rows = FALSE,   # Disable row clustering
              cluster_columns = FALSE # Disable column clustering
      )
      
      dev.off()  # Close the PDF device
  }
  
  #heatmap for unique Naive and Jurkat Transcripts without co-expressed
  if(T){
    NaiveOnly<-NaiveNano1[which(NaiveNano1$transcript_name %in% N),]
    colnames(NaiveOnly)<-c("transcript_name", "rawNaive",             "est_countNaive",       "tpmNaive"   )
    NaiveOnly<-NaiveOnly[which(NaiveOnly$tpmNaive>1),]
    JurkatOnly<-JurkatNano1[which(JurkatNano1$transcript_name %in% J),]
    colnames(JurkatOnly)<-c("transcript_name", "rawJurkat",             "est_countJurkat",       "tpmJurkat"   )
    JurkatOnly<-JurkatOnly[which(JurkatOnly$tpmJurkat>1),]
    NaiveJurkatOnly<-merge(NaiveOnly,JurkatOnly, all=TRUE)
    colnames(NaiveJurkatOnly)<-c("transcript_name", "Annotation", "rawNaive",
                                 "est_countNaive", "tpmNaive", "rawJurkat",
                                 "est_countJurkat","tpmJurkat" )
    NaiveJurkatOnly[is.na(NaiveJurkatOnly)]<-0
    #NaiveJurkatOnly <- NaiveJurkatOnly[order(NaiveJurkatOnly$tpmNaive,decreasing = TRUE), ]
    
    NaiveJurkatOnly$ordering <- with(NaiveJurkatOnly, 
                                     ifelse(tpmNaive > 1 & tpmJurkat == 0, 1, 
                                            ifelse(tpmNaive > 1 & tpmJurkat > 1, 2, 
                                                   ifelse(tpmJurkat > 1 & tpmNaive == 0, 4, 3))))
    
    NaiveJurkatOnly <- NaiveJurkatOnly[order(NaiveJurkatOnly$ordering), ]
    
    # Recreate the heatmap matrix with the new order
    heatmap_data <- NaiveJurkatOnly[, c("tpmNaive", "tpmJurkat")]
    rownames(heatmap_data) <- gsub("\\|.*", "", NaiveJurkatOnly$transcript_name)  # Use first part of transcript_name as row names
    heatmap_data<-heatmap_data[which(heatmap_data$tpmNaive<1000 & heatmap_data$tpmJurkat<1000),]
    heatmap_matrix <- as.matrix(heatmap_data)
    colnames(heatmap_matrix) <- c("Naive","Jurkat")
    
    
    pdf("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/tpm_heatmap_UniqueTranscripts.pdf", width = 3, height = 7) 
    
    heatmap_colors <- c("white", colorRampPalette(c("lightblue", "royalblue"))(20),colorRampPalette(c("royalblue2", "blue4"))(400))
    Heatmap(heatmap_matrix, 
            name = "TPM",           # Name of the color legend
            col = heatmap_colors,          # Apply custom color function
            row_order = NULL,       # Keep rows in the current order
            row_names_gp = gpar(fontsize = 1),      # Customize row labels
            cluster_rows = FALSE,   # Disable row clustering
            cluster_columns = FALSE # Disable column clustering
    )
    
    dev.off()  # Close the PDF device
     
  }
      
  }
  
  #focus on unique transcripts to detect psi sites
  if(T){
    
    #vcf check of unique sites, do once, then just load the files
    if(T){
      #supplemented with panIVT, with p<0.001 
      Naive_orig <- read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveJurkat_panIVT/Naive_panIVT/init_gene_pileup/finalized_pileup/Naive_filtered_pan.csv", header = T)
      Jurkat_orig <- read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveJurkat_panIVT/Jurkat_panIVT/init_gene_pileup/finalized_pileup/Jurkat_filtered_pan.csv", header = T)
       
      Naive_filt<-Naive_orig[which(Naive_orig$N_reads_NaiveDRS>20 & Naive_orig$mm.NaiveDRS>10),]
      Naive_filt<-Naive_filt[which(Naive_filt$N_reads_NaiveIVT>10 & Naive_filt$mm.NaiveIVT<10),]
      
      Jurkat_filt<-Jurkat_orig[which(Jurkat_orig$N_reads_JurkatDRS>20 & Jurkat_orig$mm.JurkatDRS>10),]
      Jurkat_filt<-Jurkat_filt[which(Jurkat_filt$N_reads_JurkatIVT>10 & Jurkat_filt$mm.JurkatIVT<10),]
      
      # Find the intersections
      intersect_1_2 <- intersect(Naive_filt$ACP, Jurkat_filt$ACP)  # Intersection of Group 1 and Group 2
      # Unique elements in each group
      unique_N <- setdiff(Naive_filt$Annotation, Jurkat_filt$Annotation)  # Unique to Group 1
      unique_J <- setdiff(Jurkat_filt$Annotation, Naive_filt$Annotation)  # Unique to Group 2
      
      common<-cbind(Naive_filt[which(Naive_filt$ACP %in% intersect_1_2),],Jurkat_filt[which(Jurkat_filt$ACP %in% intersect_1_2),] )
      Naive_unique_psi<-Naive_filt[which(Naive_filt$Annotation %in% unique_N),]
      Jurkat_unique_psi<-Jurkat_filt[which(Jurkat_filt$Annotation %in% unique_J),]
      
      #vcf check, did once
      if(F){
        
        # Search my sites in the SNV dataset
        rna_modifications <- data.frame(
          chr = Naive_unique_psi$chr,  # Replace with your actual chromosome data
          position = Naive_unique_psi$position,  # Replace with your actual positions
          annotation = Naive_unique_psi$Annotation  # Replace with your actual annotations
        )
        
        rna_mod_gr <- GRanges(
          seqnames = rna_modifications$chr,
          ranges = IRanges(start = rna_modifications$position, end = rna_modifications$position),
          annotation = rna_modifications$annotation
        )
        
        vcf_file <- "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/vcf_GRCh38p10/Naive_pileup.vcf"  # Replace with the path to your VCF file
        vcf <- readVcf(vcf_file, "hg18") 
        snv_gr_Naive <- rowRanges(vcf)
        #concatenating both SNV from Jurkat and Naive in one large object
        snv_gr<-snv_gr_Naive
        
        # Get the indices of the overlapping sites between my sites and the snv_gr containing all the SNV indexes
        overlaps <- findOverlaps(rna_mod_gr, snv_gr)
        overlapping_indices <- queryHits(overlaps)
        overlapping_sites <- rna_mod_gr[overlapping_indices]
        
        # Get all the non-overlapping sites that are NOT SNVs
        all_indices <- seq_len(nrow(rna_modifications))
        overlapping_indices <- queryHits(overlaps)
        non_overlapping_indices <- setdiff(all_indices, overlapping_indices)
        non_overlapping_sites <- rna_modifications[non_overlapping_indices, ]
        non_overlapping_ACP<- paste0(non_overlapping_sites$annotation, non_overlapping_sites$chr,non_overlapping_sites$position)
        #removing all SNVs from the list
        Naive_unique_psi<-Naive_unique_psi[which(Naive_unique_psi$ACP %in% non_overlapping_ACP),]
      }
      #vcf check, did once
      if(F){
        
        # Search my sites in the SNV dataset
        rna_modifications <- data.frame(
          chr = Jurkat_unique_psi$chr,  # Replace with your actual chromosome data
          position = Jurkat_unique_psi$position,  # Replace with your actual positions
          annotation = Jurkat_unique_psi$Annotation  # Replace with your actual annotations
        )
        
        rna_mod_gr <- GRanges(
          seqnames = rna_modifications$chr,
          ranges = IRanges(start = rna_modifications$position, end = rna_modifications$position),
          annotation = rna_modifications$annotation
        )
        
        vcf_file <- "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/vcf_GRCh38p10/Jurkat_pileup.vcf"  # Replace with the path to your VCF file
        vcf <- readVcf(vcf_file, "hg18") 
        snv_gr_Jurkat <- rowRanges(vcf)
        
        #concatenating both SNV from Jurkat and Naive in one large object
        snv_gr<-snv_gr_Jurkat
        
        # Get the indices of the overlapping sites between my sites and the snv_gr containing all the SNV indexes
        overlaps <- findOverlaps(rna_mod_gr, snv_gr)
        overlapping_indices <- queryHits(overlaps)
        overlapping_sites <- rna_mod_gr[overlapping_indices]
        
        # Get all the non-overlapping sites that are NOT SNVs
        all_indices <- seq_len(nrow(rna_modifications))
        overlapping_indices <- queryHits(overlaps)
        non_overlapping_indices <- setdiff(all_indices, overlapping_indices)
        non_overlapping_sites <- rna_modifications[non_overlapping_indices, ]
        non_overlapping_ACP<- paste0(non_overlapping_sites$annotation, non_overlapping_sites$chr,non_overlapping_sites$position)
        #removing all SNVs from the list
        Jurkat_unique_psi<-Jurkat_unique_psi[which(Jurkat_unique_psi$ACP %in% non_overlapping_ACP),]
      }
      
      #search for transcript names present on co-expressed transcripts and delete those.
      df_combined<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/CommonSitesNaiveJurkat.csv",header=TRUE)
      a<-Naive_unique_psi$Annotation[which(Naive_unique_psi$Annotation %in% df_combined$Annotation)]
      b<-Jurkat_unique_psi$Annotation[which(Jurkat_unique_psi$Annotation %in% df_combined$Annotation)]
      # Filter out rows with Annotation in the list
      Naive_unique_psi <- Naive_unique_psi %>% 
        filter(!Annotation %in% a)
      Jurkat_unique_psi <- Jurkat_unique_psi %>%
        filter(!Annotation %in% b)
      
      write.csv(Naive_unique_psi,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Naive_Uniquetranscripts_psi_noSNV.csv",row.names = F)
      write.csv(Jurkat_unique_psi,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Jurkat_Uniquetranscripts_psi_noSNV.csv",row.names = F)
    }
    
    #these have mmDRS>10, filter it to have >30 
    Naive_unique_psi<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Naive_Uniquetranscripts_psi_noSNV.csv",header=T)
    Jurkat_unique_psi<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Jurkat_Uniquetranscripts_psi_noSNV.csv",header=T)
    
    N_u_30<-Naive_unique_psi[which(Naive_unique_psi$mm.NaiveDRS>30),]
    J_u_30<-Jurkat_unique_psi[which(Jurkat_unique_psi$mm.JurkatDRS>30),]
    
    
    #plot Naive #reads vs %mm
    if(T){
      
      Naive_pairedIVT <- read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveTCell/init_gene_pileup/finalized_pileup/merged_pvals.csv", header = T)
      Naive_pairedIVT$ACP<-paste0( Naive_pairedIVT$Annotation, Naive_pairedIVT$chr, Naive_pairedIVT$position)
      Jurkat_pairedIVT <- read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Jurkat/init_gene_pileup/finalized_pileup/Merged_with_P_vals.csv", header = T)
      Jurkat_pairedIVT$ACP<-paste0(Jurkat_pairedIVT$Annotation,Jurkat_pairedIVT$chr,Jurkat_pairedIVT$position)
      
      unique_N_paired <- setdiff(Naive_pairedIVT$ACP, Jurkat_pairedIVT$ACP)  # Unique to Group 1
      unique_J_paired <- setdiff(Jurkat_pairedIVT$ACP, Naive_pairedIVT$ACP)  # Unique to Group 2
      
      Naive_unique_psi_p<-Naive_pairedIVT[which(Naive_pairedIVT$ACP %in% unique_N_paired),]
      Naive_unique_psi_p<-Naive_unique_psi_p[which(Naive_unique_psi_p$p.value.NaiveDRS<0.01),]
      Naive_unique_psi_p<-Naive_unique_psi_p[which(Naive_unique_psi_p$mm.NaiveIVT<10 & Naive_unique_psi_p$N_reads_NaiveIVT>10),]
      Jurkat_unique_psi_p<-Jurkat_pairedIVT[which(Jurkat_pairedIVT$ACP %in% unique_J_paired),]
      Jurkat_unique_psi_p<-Jurkat_unique_psi_p[which(Jurkat_unique_psi_p$p.value.JurkatDRS<0.01),]
      Jurkat_unique_psi_p<-Jurkat_unique_psi_p[which(Jurkat_unique_psi_p$mm.JurkatIVT<10 & Jurkat_unique_psi_p$N_reads_JurkatIVT>10),]
       
     
      
      N_plot<-ggplot() +
        # Plot points from the first dataframe
        geom_point(data = Naive_unique_psi_p, 
                   aes(x = mm.NaiveDRS, y = log10(N_reads_NaiveDRS)), 
                   color = "black", shape = 16, size = 1) +
        # Plot points from the second dataframe
        geom_point(data = Naive_unique_psi, 
                   aes(x = mm.NaiveDRS, y = log10(N_reads_NaiveDRS)), 
                   color = "black", shape = 16, size = 1) +
        # Add horizontal and vertical lines
        geom_hline(yintercept = log10(20), color = "orange", linetype = "dashed", size = 0.25) +
        geom_vline(xintercept = 30, color = "orange", linetype = "dashed", size = 0.25) +
        # Add labels and theme adjustments
        labs(title = "Scatterplot of NaiveDRS ",
             x = "positional occupancy U-to-C %mm",
             y = "N_reads DRS") +
        theme_minimal() +
        theme(
          panel.background = element_rect(fill = "white"),  # Set background to white
          axis.line = element_line(color = "black"),        # Black axis lines
          panel.grid.major = element_blank(),               # Remove major grid lines
          
        )
      ggsave("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/readsVSmmNaive.pdf", width = 10, height=5)
      
      N_plot<-ggplot() +
        # Plot points from the first dataframe
        geom_point(data = Naive_unique_psi_p[which(Naive_unique_psi_p$kmer %in% c("GTTCA","GTTCC","GTTCT","GTTCG")),], 
                   aes(x = mm.NaiveDRS, y = log10(N_reads_NaiveDRS)), 
                   color = "red", shape = 16, size = 1) +
        geom_point(data = Naive_unique_psi_p[which(Naive_unique_psi_p$kmer  %in% c("TGTAG", "TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")),], 
                   aes(x = mm.NaiveDRS, y = log10(N_reads_NaiveDRS)), 
                   color = "green", shape = 16, size = 1) +
        # Plot points from the second dataframe
        geom_point(data = Naive_unique_psi[which(Naive_unique_psi$kmer %in% c("GTTCA","GTTCC","GTTCT","GTTCG")),], 
                   aes(x = mm.NaiveDRS, y = log10(N_reads_NaiveDRS)), 
                   color = "red", shape = 16, size = 1) +
        geom_point(data = Naive_unique_psi[which(Naive_unique_psi$kmer  %in% c("TGTAG", "TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")),], 
                   aes(x = mm.NaiveDRS, y = log10(N_reads_NaiveDRS)), 
                   color = "green", shape = 16, size = 1) +
        
        # Add horizontal and vertical lines
        geom_hline(yintercept = log10(20), color = "orange", linetype = "dashed", size = 0.25) +
        geom_vline(xintercept = 30, color = "orange", linetype = "dashed", size = 0.25) +
        # Add labels and theme adjustments
        labs(title = "Scatterplot of NaiveDRS ",
             x = "positional occupancy U-to-C %mm",
             y = "N_reads DRS") +
        theme_minimal() +
        theme(
          panel.background = element_rect(fill = "white"),  # Set background to white
          axis.line = element_line(color = "black"),        # Black axis lines
          panel.grid.major = element_blank(),               # Remove major grid lines
          
        )
      ggsave("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/readsVSmmNaiveMotifs.pdf", width = 10, height=5)
      
      
      
       ggplot() +
        # Plot points from the first dataframe
        geom_point(data = Jurkat_unique_psi_p, 
                   aes(x = mm.JurkatDRS, y = log10(N_reads_JurkatDRS)), 
                   color = "black", shape = 16, size = 1) +
        # Plot points from the second dataframe
        geom_point(data = Jurkat_unique_psi, 
                   aes(x = mm.JurkatDRS, y = log10(N_reads_JurkatDRS)), 
                   color = "black", shape = 16, size = 1) +
        # Add horizontal and vertical lines
        geom_hline(yintercept = log10(20), color = "royalblue", linetype = "dashed", size = 0.25) +
        geom_vline(xintercept = 30, color = "royalblue", linetype = "dashed", size = 0.25) +
        # Add labels and theme adjustments
        labs(title = "Scatterplot JurkatDRS",
             x = "positional occupancy U-to-C %mm",
             y = "N_reads DRS") +
        theme_minimal() +
        theme(
          panel.background = element_rect(fill = "white"),  # Set background to white
          axis.line = element_line(color = "black"),        # Black axis lines
          panel.grid.major = element_blank(),               # Remove major grid lines
          
        )

      ggsave("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/readsVSmmJurkat.pdf", width = 10, height=5)
      

      ggplot() +
        # Plot points from the first dataframe
        geom_point(data = Jurkat_unique_psi_p[which(Jurkat_unique_psi_p$kmer %in% c("GTTCA","GTTCC","GTTCT","GTTCG")),], 
                   aes(x = mm.JurkatDRS, y = log10(N_reads_JurkatDRS)), 
                   color = "red", shape = 16, size = 1) +
        geom_point(data = Jurkat_unique_psi_p[which(Jurkat_unique_psi_p$kmer  %in% c("TGTAG", "TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")),], 
                   aes(x = mm.JurkatDRS, y = log10(N_reads_JurkatDRS)), 
                   color = "green", shape = 16, size = 1) + 
        # Plot points from the second dataframe
        geom_point(data = Jurkat_unique_psi[which(Jurkat_unique_psi$kmer %in% c("GTTCA","GTTCC","GTTCT","GTTCG")),], 
                   aes(x = mm.JurkatDRS, y = log10(N_reads_JurkatDRS)), 
                   color = "red", shape = 16, size = 1) +
        geom_point(data = Jurkat_unique_psi[which(Jurkat_unique_psi$kmer  %in% c("TGTAG", "TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")),], 
                   aes(x = mm.JurkatDRS, y = log10(N_reads_JurkatDRS)), 
                   color = "green", shape = 16, size = 1) + 
        # Add horizontal and vertical lines
        geom_hline(yintercept = log10(20), color = "royalblue", linetype = "dashed", size = 0.25) +
        geom_vline(xintercept = 30, color = "royalblue", linetype = "dashed", size = 0.25) +
        # Add labels and theme adjustments
        labs(title = "Scatterplot JurkatDRS",
             x = "positional occupancy U-to-C %mm",
             y = "N_reads DRS") +
        theme_minimal() +
        theme(
          panel.background = element_rect(fill = "white"),  # Set background to white
          axis.line = element_line(color = "black"),        # Black axis lines
          panel.grid.major = element_blank(),               # Remove major grid lines
          
        )
      
      ggsave("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/readsVSmmJurkatMotifs.pdf", width = 10, height=5)
      
      
      
      
       
    }
    
    #plot KDE %mm of unique sites 
    if(T){
     
      Naive_unique_psi<-Naive_unique_psi[which(Naive_unique_psi$mm.NaiveDRS>30),]
      Jurkat_unique_psi<-Jurkat_unique_psi[which(Jurkat_unique_psi$mm.JurkatDRS>30),]
      
      ggplot(Naive_unique_psi, aes(x = mm.NaiveDRS, color = "orange", fill = "orange")) +
        geom_density(alpha = 0.4, adjust = 1) + # 'adjust' can be tuned for smoothness
        labs(title = "Kernel Density Estimate",
             x = "Value",
             y = "Density") +
        theme_minimal() +
        scale_fill_manual(values = c("orange")) + 
        scale_color_manual(values = c("orange")) 
      
      ggsave("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/KDENaivemm30.pdf", width = 5, height=5)
      
      ggplot(Jurkat_unique_psi, aes(x = mm.JurkatDRS, color = "royalblue", fill = "royalblue")) +
        geom_density(alpha = 0.4, adjust = 1) + # 'adjust' can be tuned for smoothness
        labs(title = "Kernel Density Estimate",
             x = "Value",
             y = "Density") +
        theme_minimal() +
        scale_fill_manual(values = c("royalblue")) + 
        scale_color_manual(values = c("royalblue")) 
      
      ggsave("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/KDEJurkatmm30.pdf", width = 5, height=5)

      
      
    }
    
    #utr-cds
    if(T){
      
      library(rtracklayer)
      g = readGFF("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/gencode/gencode.v27.annotation.gff3")
      gff3=as.data.frame(g)
      gff3 <- gff3[which(gff3$gene_type == "protein_coding" & ((gff3$type == "CDS") | (gff3$type == "three_prime_UTR") | (gff3$type == "five_prime_UTR") | (gff3$type == "start_codon") | (gff3$type == "stop_codon") | (gff3$type == "stop_codon_redefined_as_selenocysteine"))), c("type", "seqid","start", "end", "gene_name")]
      
      #for unique Naive sites
      file.name="/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/UniqueTranscriptsNaivePsiCDSUTR.csv"
      #utr-cds attributions
      if(T){
        Naive_unique_psi$CDS = F
        Naive_unique_psi$Three_prime_UTR = F
        Naive_unique_psi$Five_Prime_UTR = F
        Naive_unique_psi$start_codon = F
        Naive_unique_psi$stop_codon = F
        Naive_unique_psi$stop_codon_redefined_as_selenocysteine = F
        
        for (row in 1:nrow(Naive_unique_psi)){
          
          type = gff3$type[which(gff3$start <= Naive_unique_psi$position[row] & 
                                   gff3$end >= Naive_unique_psi$position[row] &
                                   gff3$seqid == Naive_unique_psi$chr[row])]
          
          if (length(type)>0){
            if ("CDS" %in% type) {
              Naive_unique_psi$CDS[row] = T
            }
            if ("three_prime_UTR" %in% type) {
              Naive_unique_psi$Three_prime_UTR[row] = T
            }
            if ("five_prime_UTR" %in% type) {
              Naive_unique_psi$Five_Prime_UTR[row] = T
            }
            if ("start_codon" %in% type) {
              Naive_unique_psi$start_codon[row] = T
            }
            if ("stop_codon" %in% type) {
              Naive_unique_psi$stop_codon[row] = T
            }
            if ("stop_codon_redefined_as_selenocysteine" %in% type) {
              Naive_unique_psi$stop_codon_redefined_as_selenocysteine[row] = T
            }
          }
        }
        
        
        start.codon <- Naive_unique_psi[which(Naive_unique_psi$start_codon == 1),]
        
        all.false <- Naive_unique_psi[which(Naive_unique_psi$CDS == 0 &
                                              Naive_unique_psi$Three_prime_UTR == 0 &
                                              Naive_unique_psi$Five_Prime_UTR == 0 &
                                              Naive_unique_psi$start_codon == 0 &
                                              Naive_unique_psi$stop_codon == 0),]
        
      }
      write.csv(Naive_unique_psi,file.name,row.names = F)  
      
      #for unique Jurkat sites
      file.name="/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/UniqueTranscriptsJurkatPsiCDSUTR.csv"
      #utr-cds attributions
      if(T){
        
        Jurkat_unique_psi$CDS = F
        Jurkat_unique_psi$Three_prime_UTR = F
        Jurkat_unique_psi$Five_Prime_UTR = F
        Jurkat_unique_psi$start_codon = F
        Jurkat_unique_psi$stop_codon = F
        Jurkat_unique_psi$stop_codon_redefined_as_selenocysteine = F
        
        for (row in 1:nrow(Jurkat_unique_psi)){
          
          type = gff3$type[which(gff3$start <= Jurkat_unique_psi$position[row] & 
                                   gff3$end >= Jurkat_unique_psi$position[row] &
                                   gff3$seqid == Jurkat_unique_psi$chr[row])]
          
          if (length(type)>0){
            if ("CDS" %in% type) {
              Jurkat_unique_psi$CDS[row] = T
            }
            if ("three_prime_UTR" %in% type) {
              Jurkat_unique_psi$Three_prime_UTR[row] = T
            }
            if ("five_prime_UTR" %in% type) {
              Jurkat_unique_psi$Five_Prime_UTR[row] = T
            }
            if ("start_codon" %in% type) {
              Jurkat_unique_psi$start_codon[row] = T
            }
            if ("stop_codon" %in% type) {
              Jurkat_unique_psi$stop_codon[row] = T
            }
            if ("stop_codon_redefined_as_selenocysteine" %in% type) {
              Jurkat_unique_psi$stop_codon_redefined_as_selenocysteine[row] = T
            }
          }
        }
        
        
        start.codon <- Jurkat_unique_psi[which(Jurkat_unique_psi$start_codon == 1),]
        
        all.false <- Jurkat_unique_psi[which(Jurkat_unique_psi$CDS == 0 &
                                               Jurkat_unique_psi$Three_prime_UTR == 0 &
                                               Jurkat_unique_psi$Five_Prime_UTR == 0 &
                                               Jurkat_unique_psi$start_codon == 0 &
                                               Jurkat_unique_psi$stop_codon == 0),]
        
      }
      write.csv(Jurkat_unique_psi,file.name,row.names = F)  
      
      
      #plot unique sites 
      #plot cds utr, exclude positions from plot if found in multiple regions as they are on different isoforms
      if(T){
      
        min.x = 0
        max.x = 6
        min.y = 0
        max.y = 100
        label.size = 1
        text.size = 0.8
        point.size = 1
        tck.length = 0.01
        tick.thickness = 1
        transparency = 0.2
        file.name ="/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Unique_Transcripts_UTR-CDSmm30.pdf"
        pdf(file.name, width = 6, height = 6, family = "Helvetica")
        plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
             ylab = NA,xlab=NA,ann = F,axes=F)
        #box(lwd = line.thickness)
        #### Add the lines
        par(mfrow = c(1, 2)) 
        
        #normalizing for number of reads. NO NORMALIZATIONS IN THIS CASE
        max_reads= 1#max(c(sum(ALL_NaiveJurkat$total_TRUB1_KD),sum(PUS7KD_sites$total_PUS7_KD)))
        factor_norm=1#000
        spacer=0.1
        
        #Jurkat only
        ALL_NaiveJurkat_JurkatOnly<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/UniqueTranscriptsJurkatPsiCDSUTR.csv",header = T)
        ALL_NaiveJurkat_JurkatOnly<-ALL_NaiveJurkat_JurkatOnly[which(ALL_NaiveJurkat_JurkatOnly$mm.JurkatDRS>30),]
        Prop<-c(length(which(ALL_NaiveJurkat_JurkatOnly$stop_codon==T & ALL_NaiveJurkat_JurkatOnly$Three_prime_UTR==F & ALL_NaiveJurkat_JurkatOnly$Five_Prime_UTR==F & ALL_NaiveJurkat_JurkatOnly$CDS==F))/max_reads*factor_norm,
                length(which(ALL_NaiveJurkat_JurkatOnly$start_codon==T & ALL_NaiveJurkat_JurkatOnly$Three_prime_UTR==F & ALL_NaiveJurkat_JurkatOnly$Five_Prime_UTR==F & ALL_NaiveJurkat_JurkatOnly$CDS==F))/max_reads*factor_norm,
                length(which(ALL_NaiveJurkat_JurkatOnly$Three_prime_UTR==T & ALL_NaiveJurkat_JurkatOnly$Five_Prime_UTR==F & ALL_NaiveJurkat_JurkatOnly$CDS==F ))/max_reads*factor_norm,
                length(which(ALL_NaiveJurkat_JurkatOnly$CDS==T & ALL_NaiveJurkat_JurkatOnly$Three_prime_UTR==F & ALL_NaiveJurkat_JurkatOnly$Five_Prime_UTR==F ))/max_reads*factor_norm,
                length(which(ALL_NaiveJurkat_JurkatOnly$Five_Prime_UTR==T & ALL_NaiveJurkat_JurkatOnly$CDS==F & ALL_NaiveJurkat_JurkatOnly$Three_prime_UTR==F))/max_reads*factor_norm)
        
        labels<-c("STOP","START","3'UTR","CDS","5'UTR")
        numbers<-Prop
        pie(Prop, labels = paste0(labels,"-",numbers), border="white", col= c("red","green","lightblue", "royalblue4","royalblue"),cex = 0.7 ) 
        
        #Naive only
        ALL_NaiveJurkat_NaiveOnly<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/UniqueTranscriptsNaivePsiCDSUTR.csv",header = T)
        ALL_NaiveJurkat_NaiveOnly<-ALL_NaiveJurkat_NaiveOnly[which(ALL_NaiveJurkat_NaiveOnly$mm.NaiveDRS>30),]
        Prop<-c(length(which(ALL_NaiveJurkat_NaiveOnly$stop_codon==T & ALL_NaiveJurkat_NaiveOnly$Three_prime_UTR==F & ALL_NaiveJurkat_NaiveOnly$Five_Prime_UTR==F & ALL_NaiveJurkat_NaiveOnly$CDS==F ))/max_reads*factor_norm,
                length(which(ALL_NaiveJurkat_NaiveOnly$start_codon==T & ALL_NaiveJurkat_NaiveOnly$Three_prime_UTR==F & ALL_NaiveJurkat_NaiveOnly$Five_Prime_UTR==F & ALL_NaiveJurkat_NaiveOnly$CDS==F))/max_reads*factor_norm,
                length(which(ALL_NaiveJurkat_NaiveOnly$Three_prime_UTR==T & ALL_NaiveJurkat_NaiveOnly$Five_Prime_UTR==F & ALL_NaiveJurkat_NaiveOnly$CDS==F ))/max_reads*factor_norm,
                length(which(ALL_NaiveJurkat_NaiveOnly$CDS==T & ALL_NaiveJurkat_NaiveOnly$Three_prime_UTR==F & ALL_NaiveJurkat_NaiveOnly$Five_Prime_UTR==F ))/max_reads*factor_norm,
                length(which(ALL_NaiveJurkat_NaiveOnly$Five_Prime_UTR==T & ALL_NaiveJurkat_NaiveOnly$CDS==F & ALL_NaiveJurkat_NaiveOnly$Three_prime_UTR==F))/max_reads*factor_norm)
        
        labels<-c("STOP","START","3'UTR","CDS","5'UTR")
        numbers<-Prop
        pie(Prop, labels = paste0(labels,"-",numbers), border="white",col= c("red","green","burlywood1","darkorange3","orange"),cex = 0.7)
        
        dev.off() 
      }
      
    }
    
    #seqLogos for variable/stable shared sites and unique 
    if (plot_figs==T){
      
      uniqueN<-Naive_unique_psi
      uniqueJ<-Jurkat_unique_psi
      
      N=30
      
      kmers<-sapply( uniqueN$kmer,function(x) str_replace_all(x,"T","U"))
      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/logo_UniqueTranscriptsNaive.pdf"
      kmers<-sapply( uniqueN$kmer[which(uniqueN$mm.NaiveDRS>N)],function(x) str_replace_all(x,"T","U"))
      file.name = paste0("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/logo_UniqueTranscriptsNaive",N,".pdf")
      
      
      kmers<-sapply( uniqueJ$kmer,function(x) str_replace_all(x,"T","U"))
      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/logo_UniqueTranscriptsJurkat.pdf"
      kmers<-sapply( uniqueJ$kmer[which(uniqueJ$mm.JurkatDRS>30)],function(x) str_replace_all(x,"T","U"))
      file.name = paste0("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/logo_UniqueTranscriptsJurkat",N,".pdf")
      
      my_list<-list("new motifs" = kmers)
      min.x = 0
      max.x = 200
      min.y = 0
      max.y = 20
      line.thickness = 1
      label.size = 1
      text.size = 0.8
      point.size = 1
      tck.length = 0.01
      tick.thickness = 1 
      transparency = 0.2
      Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
      CairoPDF(file.name, width = 6, height = 3, family = "Helvetica")
      par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
      plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
           ylab = NA,xlab=NA,ann = F,axes=F)
      
      ggseqlogo::ggseqlogo(my_list,ncol=1,method="prob")
      dev.off()
      
    }
    
    #files to use for GO analysis on EnrichR
    N_u_30<-Naive_unique_psi[which(Naive_unique_psi$mm.NaiveDRS>30),]
    J_u_30<-Jurkat_unique_psi[which(Jurkat_unique_psi$mm.JurkatDRS>30),]
    write.csv(N_u_30,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Naive_Uniquetranscripts_psi_noSNV_mm30.csv",row.names = FALSE)
    write.csv(J_u_30,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Jurkat_Uniquetranscripts_psi_noSNV_mm30.csv", row.names = FALSE)
    
     
  }
  
}



### plots figure HYPERMOD TYPE 2 with %mm>30 nin DRS sample
if(T){
  
  NaiveBootstrap<-read.csv("/Dropbox/merged_pvals_Naive_bootstrap.csv",header=TRUE)
  NaiveBootstrap$ACP<-paste0(NaiveBootstrap$Annotation,NaiveBootstrap$chr,NaiveBootstrap$position)
  JurkatBootstrap<-read.csv("/Dropbox/merged_pvals_Jurkat_bootstrap.csv",header=TRUE)
  JurkatBootstrap$ACP<-paste0(JurkatBootstrap$Annotation,JurkatBootstrap$chr,JurkatBootstrap$position)
  
  #No bootstrap
  NaiveJurkat<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveJurkat_panIVT/Jurkat_panIVT/init_gene_pileup/finalized_pileup/Jurkat_Naive_filtered_pan_noPvalue.csv",header = T)
  
  g = readGFF("gencode.v27.annotation.gff3") 
  gff3=as.data.frame(g)
  gff3_gene <- gff3[which(gff3$type == "gene"),]
  gff3_gene$length<- abs(gff3_gene$end - gff3_gene$start)
  
  gff3_cds_utr <- gff3[which(gff3$gene_type == "protein_coding" & ((gff3$type == "CDS") | (gff3$type == "three_prime_UTR") | (gff3$type == "five_prime_UTR") | (gff3$type == "start_codon") | (gff3$type == "stop_codon") | (gff3$type == "stop_codon_redefined_as_selenocysteine"))), c("type", "seqid","start", "end", "gene_name")]
  
  #Barplot with %mm>30
  if(T){
    #Naive results table construction
    if(T){
      
      # Initialize a list to store results for each N
      results_Naive <- list()
      lengths_h2 <- list("psi2" = data.frame("Annotation"="","length"=c(0)), "psi3" = data.frame("Annotation"="","length"=c(0)),
                         "psi4" = data.frame("Annotation"="","length"=c(0)), "psi5" = data.frame("Annotation"="","length"=c(0)),
                         "psi6" = data.frame("Annotation"="","length"=c(0)), "psi7" = data.frame("Annotation"="","length"=c(0)),
                         "psi8" = data.frame("Annotation"="","length"=c(0)), "psi9" = data.frame("Annotation"="","length"=c(0)),
                         "psi10" = data.frame("Annotation"="","length"=c(0)), "psi11" = data.frame("Annotation"="","length"=c(0)),
                         "psi12" = data.frame("Annotation"="","length"=c(0)), "psi13" = data.frame("Annotation"="","length"=c(0)),
                         "psi14" = data.frame("Annotation"="","length"=c(0)), "psi15" = data.frame("Annotation"="","length"=c(0)), "psi16" = data.frame("Annotation"="","length"=c(0)))
      
      ACP_h2<-list("psi2"  = data.frame("ACP"=""),  "psi3" = data.frame("ACP"=""),
                   "psi4"  = data.frame("ACP"=""),  "psi5" = data.frame("ACP"=""),
                   "psi6"  = data.frame("ACP"=""),  "psi7" = data.frame("ACP"=""),
                   "psi8"  = data.frame("ACP"=""),  "psi9" = data.frame("ACP"=""),
                   "psi10" = data.frame("ACP"=""), "psi11" = data.frame("ACP"=""),
                   "psi12" = data.frame("ACP"=""), "psi13" = data.frame("ACP"=""),
                   "psi14" = data.frame("ACP"=""), "psi15" = data.frame("ACP"=""),
                   "psi16" = data.frame("ACP"=""))
      
      
      
      # Iterate through N0 to N4
      for (N in c("N0", "N1", "N2", "N3", "N4")) {
        # Dynamically construct column names
        N_reads_col <- paste0("N_reads_", N)
        mm_NaiveDRS_col <- paste0("mm.", N)
        
        # Filter data for the current N
        filtered <- NaiveBootstrap[
          which(
            NaiveBootstrap[[N_reads_col]] > 10 &
              NaiveBootstrap[[mm_NaiveDRS_col]] > 30 &
              NaiveBootstrap$mm.NaiveIVT < 10 &
              NaiveBootstrap$N_reads_NaiveIVT > 10
          ),
        ]
        
        # Create hypermodified table
        hypermodified2Naive <- table(filtered$Annotation)
        
        # Initialize Naive_h2 for this N
        Naive_h2 <- data.frame("psi2" = 0, "psi3" = 0, "psi4" = 0, "psi5" = 0, "psi6" = 0, 
                               "psi7" = 0, "psi8" = 0, "psi9" = 0, "psi10" = 0, "psi11" = 0,
                               "psi12" = 0, "psi13" = 0, "psi14" = 0, "psi15" = 0, "psi16" = 0)
        
        # Fill Naive_h2
        for (i in 2:(length(Naive_h2) + 1)) {
          Naive_h2[i - 1] <- length(which(hypermodified2Naive == i))
          if (i > length(Naive_h2)) {
            Naive_h2[i - 1] <- length(which(hypermodified2Naive >= i))
          }
        }
        
        
        #calculate mean distance between each site
        for (i in 2:(length(lengths_h2) + 1)) {
          psi<-paste0("psi",i)
          gene_names<-names(hypermodified2Naive[hypermodified2Naive == i])
          
          if(length(gene_names)==0){break}
          merged<-data.frame()
          mean_genes <- data.frame("Annotation"=gene_names,"length"="")
          merged_acp<-c()
          for (gene in gene_names){
            coordinates<-filtered$position[which(filtered$Annotation == gene)]
            ACP<-paste0(filtered$Annotation[which(filtered$Annotation == gene)], filtered$chr[which(filtered$Annotation == gene)],filtered$position[which(filtered$Annotation == gene)])
            merged_acp<-c(merged_acp, ACP)
            mean_diff <- abs(outer(coordinates, coordinates, "-"))
            mean_diff <- as.vector(mean_diff[upper.tri(mean_diff)])
            mean_diff <- mean(mean_diff)
            # Append the mean_diff to the mean_genes vector
            mean_genes$length[which(mean_genes$Annotation==gene)] <- mean_diff 
            
          }
          merged<-merge(lengths_h2[[psi]],mean_genes, by="Annotation",all=TRUE)
          
          merged$mean<-rowMeans(cbind(as.numeric(merged$length.x), as.numeric(merged$length.y)), na.rm = TRUE)
          lengths_h2[[psi]]<-data.frame("Annotation"=merged$Annotation,"length"=merged$mean)
          
          ACP_h2[[psi]]<-merged_acp
        }
        
        
        #save the annotations
        
        
        # Save results for this N
        results_Naive[[N]] <- list(
          Naive_h2 = Naive_h2,
          hypermodified2Naive = as.data.frame(hypermodified2Naive)
        )
      }
      
      # Access results for N0, N1, N2, N3, or N4
      # Example: results[["N0"]]$Naive_h2
      #          results[["N0"]]$hypermodified2Naive
      
      ACP_h2_Naive <- Filter(function(df) length(df) > 1, ACP_h2)
      # Filter out empty data frames from the list
      mean_dist_Naive <- Filter(function(df) nrow(df) > 1, lengths_h2)
    
    }
    
    #Jurkat results table construction
    if(T){
      
      # Initialize a list to store results for each N
      results_Jurkat <- list()
      lengths_h2 <- list("psi2"  = data.frame("Annotation"="","length"=c(0)),  "psi3" = data.frame("Annotation"="","length"=c(0)),
                         "psi4"  = data.frame("Annotation"="","length"=c(0)),  "psi5" = data.frame("Annotation"="","length"=c(0)),
                         "psi6"  = data.frame("Annotation"="","length"=c(0)),  "psi7" = data.frame("Annotation"="","length"=c(0)),
                         "psi8"  = data.frame("Annotation"="","length"=c(0)),  "psi9" = data.frame("Annotation"="","length"=c(0)),
                         "psi10" = data.frame("Annotation"="","length"=c(0)), "psi11" = data.frame("Annotation"="","length"=c(0)),
                         "psi12" = data.frame("Annotation"="","length"=c(0)), "psi13" = data.frame("Annotation"="","length"=c(0)),
                         "psi14" = data.frame("Annotation"="","length"=c(0)), "psi15" = data.frame("Annotation"="","length"=c(0)),
                         "psi16" = data.frame("Annotation"="","length"=c(0)))
      
      ACP_h2<-list("psi2"  = data.frame("ACP"=""),  "psi3" = data.frame("ACP"=""),
                  "psi4"  = data.frame("ACP"=""),  "psi5" = data.frame("ACP"=""),
                  "psi6"  = data.frame("ACP"=""),  "psi7" = data.frame("ACP"=""),
                  "psi8"  = data.frame("ACP"=""),  "psi9" = data.frame("ACP"=""),
                  "psi10" = data.frame("ACP"=""), "psi11" = data.frame("ACP"=""),
                  "psi12" = data.frame("ACP"=""), "psi13" = data.frame("ACP"=""),
                  "psi14" = data.frame("ACP"=""), "psi15" = data.frame("ACP"=""),
                  "psi16" = data.frame("ACP"=""))
      
      
      # Iterate through N0 to N4
      for (N in c("J0", "J1", "J2", "J3", "J4")) {
        # Dynamically construct column names
        N_reads_col <- paste0("N_reads_", N)
        mm_JurkatDRS_col <- paste0("mm.", N)
        
        # Filter data for the current N
        filtered <- JurkatBootstrap[
          which(
            JurkatBootstrap[[N_reads_col]] > 10 &
              JurkatBootstrap[[mm_JurkatDRS_col]] > 30 &
              JurkatBootstrap$mm.JurkatIVT < 10 &
              JurkatBootstrap$N_reads_JurkatIVT > 10
          ),
        ]
        
        
        # Create hypermodified table
        hypermodified2Jurkat <- table(filtered$Annotation)
        
        # Initialize Jurkat_h2 for this N
        Jurkat_h2 <- data.frame("psi2" = 0, "psi3" = 0, "psi4" = 0, "psi5" = 0, "psi6" = 0, 
                               "psi7" = 0, "psi8" = 0, "psi9" = 0, "psi10" = 0, "psi11" = 0,
                               "psi12" = 0, "psi13" = 0, "psi14" = 0, "psi15" = 0, "psi16" = 0)
        
        # Fill Jurkat_h2
        for (i in 2:(length(Jurkat_h2) + 1)) {
          Jurkat_h2[i - 1] <- length(which(hypermodified2Jurkat == i))
          if (i > length(Jurkat_h2)) {
            Jurkat_h2[i - 1] <- length(which(hypermodified2Jurkat >= i))
          }
        }
        
        
        #calculate mean length
        for (i in 2:(length(lengths_h2) + 1)) {
          psi<-paste0("psi",i)
          gene_names<-names(hypermodified2Jurkat[hypermodified2Jurkat == i])
          
          if(length(gene_names)==0){break}
          merged<-data.frame()
          mean_genes <- data.frame("Annotation"=gene_names,"length"="","region"="")
          merged_acp<-c()
          for (gene in gene_names){
            coordinates<-filtered$position[which(filtered$Annotation == gene)]
            ACP<-paste0(filtered$Annotation[which(filtered$Annotation == gene)], filtered$chr[which(filtered$Annotation == gene)],filtered$position[which(filtered$Annotation == gene)])
            merged_acp<-c(merged_acp, ACP)
            mean_diff <- abs(outer(coordinates, coordinates, "-"))
            mean_diff <- as.vector(mean_diff[upper.tri(mean_diff)])
            mean_diff <- mean(mean_diff)
            # Append the mean_diff to the mean_genes vector
            mean_genes$length[which(mean_genes$Annotation==gene)] <- mean_diff 
            
          }
          merged<-merge(lengths_h2[[psi]],mean_genes, by="Annotation",all=TRUE)
          
          merged$mean<-rowMeans(cbind(as.numeric(merged$length.x), as.numeric(merged$length.y)), na.rm = TRUE)
          lengths_h2[[psi]]<-data.frame("Annotation"=merged$Annotation,"length"=merged$mean)
          
          ACP_h2[[psi]]<-merged_acp
          
        }
        
        
        # Save results for this N
        results_Jurkat[[N]] <- list(
          Jurkat_h2 = Jurkat_h2,
          hypermodified2Jurkat = as.data.frame(hypermodified2Jurkat)
        )
      }
      
      # Access results for J0, ..., J4
      # Example: results[["J0"]]$Jurkat_h2
      #          results[["J0"]]$hypermodified2Jurkat
      
      ACP_h2_Jurkat <- Filter(function(df) length(df) > 1, ACP_h2)
      mean_dist_Jurkat <- Filter(function(df) nrow(df) > 1, lengths_h2)
      
    }
    
    #barplot for psi2, psi3, psi4 ... hypermod II sites
    if(T){
      
      library(ggplot2)
      library(dplyr)
  
      # Extract mean values for Naive
      naive_means <- sapply(2:16, function(psi) {
        replicate_values <- sapply(results_Naive, function(res) res$Naive_h2[[paste0("psi", psi)]])
        mean(as.numeric(replicate_values))
      })
      # norm_N<-length(NaiveBootstrap$Annotation)
      # naive_means<-(naive_means/norm_N)*10000000
      
      
      # Extract mean values for Jurkat
      jurkat_means <- sapply(2:16, function(psi) {
        replicate_values <- sapply(results_Jurkat, function(res) res$Jurkat_h2[[paste0("psi", psi)]])
        mean(as.numeric(replicate_values))
      })
      # norm_J<-length(JurkatBootstrap$Annotation)/10000000
      
      
      # Combine into a data frame for barplot
      plot_data <- data.frame(
        psi = rep(paste0("psi", 2:16), 2),
        group = rep(c("Naive", "Jurkat"), each = 15),
        mean = c(round(naive_means), round(jurkat_means))
      )
      
      # Add replicate-level values for individual points
      replicates <- rbind(
        data.frame(
          psi = rep(paste0("psi", 2:16), each = length(results_Naive)),
          group = "Naive",
          value = unlist(lapply(2:16, function(psi) {
            sapply(results_Naive, function(res) res$Naive_h2[[paste0("psi", psi)]])
          }))
        ),
        data.frame(
          psi = rep(paste0("psi", 2:16), each = length(results_Jurkat)),
          group = "Jurkat",
          value = unlist(lapply(2:16, function(psi) {
            sapply(results_Jurkat, function(res) res$Jurkat_h2[[paste0("psi", psi)]])
          }))
        )
      )
      
      
      # Define the correct order of psi sites
      psi_order <- paste0("psi", 2:16)
      
      # Update the order in plot_data and replicates
      plot_data$psi <- factor(plot_data$psi, levels = psi_order)
      replicates$psi <- factor(replicates$psi, levels = psi_order)

      
      # Re-plot with corrected order
      ggplot() +
        # Barplot for means
        geom_bar(data = plot_data, aes(x = psi, y = mean, fill = group),
                 stat = "identity", position = position_dodge(width = 0.8), color = "white", width = 0.7) +
        # Points for individual replicates
        geom_point(data = replicates, aes(x = psi, y = value, color = group),
                   position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 0.5) +
        # Aesthetics
        labs(x = "Number of Psi Sites on a transcript", y = "Number of Transcripts") +
        scale_fill_manual(values = c("Naive" = "orange", "Jurkat" = "royalblue")) +
        scale_color_manual(values = c("Naive" = "black", "Jurkat" = "black")) +
        theme_minimal() +
        theme(
          axis.line = element_line(color = "black"),  # Black x and y axis lines
          panel.grid = element_blank(),  # Remove gridlines
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_blank()  # Remove background shading
        )
      
      ggsave("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/HypermodType2_bootstrap.pdf",width = 8, height = 5)    

    }
  }
  
  #distribution of distances between the sites
  if(T){
    
    #Function to plot KDE
    plot_kde <- function(df, name, color) {
      if (nrow(df) > 0) { # Check if the dataframe has rows
        plot(
          density(df$length), 
          main = paste("KDE of Lengths for", name),
          xlab = "Lengths",
          ylab = "Density",
          col = color,
          lwd = 2,
          xlim = c(0, max(df$length)),
          ylim = c(0,0.0004)
        )
      } else {
        message(paste("Skipping", name, "as it has no data"))
      }
    }
    
    # Set up a grid for multiple plots
    par(mfrow = c(3, 3)) # Adjust grid size as needed
    for (name in names(mean_dist_Naive)) {
      plot_kde(mean_dist_Naive[[name]], name,"orange")
    }
    
    # Set up a grid for multiple plots
    par(mfrow = c(3, 3)) # Adjust grid size as needed
    for (name in names(mean_dist_Jurkat)) {
      plot_kde(mean_dist_Jurkat[[name]], name,"royalblue")
    }
    
  
    
  }
  
  #lengths of transcripts 
  if(T){
    
    #Naive
    par(mfrow = c(3, 3)) 
    for (i in 2:length(mean_dist_Naive)) {
      psi <- paste0("psi", i)
      
      uniqueAnnotations <- unique(mean_dist_Naive[[psi]]$Annotation)
      lengths_transcripts <- gff3_gene$length[which(gff3_gene$gene_name %in% uniqueAnnotations)]
      
      # Adjust breaks dynamically for small data
      num_points <- length(lengths_transcripts)
      breaks <- if (num_points < 10) num_points else 50
      
      hist(
        lengths_transcripts, 
        main = paste("Histogram of Lengths for", psi),
        xlab = "Lengths of Transcripts",
        ylab = "Frequency",
        col = "orange",
        border = "white",
        breaks =breaks,
        xlim = c(0, 300000),
        
      )
    }
    
    #Jurkat
    par(mfrow = c(3, 3)) 
    for (i in 2:length(mean_dist_Jurkat)) {
      psi <- paste0("psi", i)
      
      uniqueAnnotations <- unique(mean_dist_Jurkat[[psi]]$Annotation)
      lengths_transcripts <- gff3_gene$length[which(gff3_gene$gene_name %in% uniqueAnnotations)]
      
      # Adjust breaks dynamically for small data
      num_points <- length(lengths_transcripts)
      breaks <- if (num_points < 10) num_points else 50
      
      hist(
        lengths_transcripts, 
        main = paste("Histogram of Lengths for", psi),
        xlab = "Lengths of Transcripts",
        ylab = "Frequency",
        col = "royalblue",
        border = "white",
        breaks =breaks,
        xlim = c(0, 300000),
      
      )
    }
    
  }
  
  #correlation lengths vs mean distances
  if(T){
    
    #Naive
    if(T){
      
    
       for (i in c(1:length(mean_dist_Naive)+1)) {
        psi <- paste0("psi", i)
        
        mean_dist_Naive[[psi]]$length_transcripts<-0
        
        uniqueAnnotations <- unique(mean_dist_Naive[[psi]]$Annotation)
        
        #CALCULATING LENGTH OF EACH TRANSCRIPT
        for (j in c(2:length(uniqueAnnotations))){
          gene<-uniqueAnnotations[j]
          mean_dist_Naive[[psi]]$length_transcripts[which(mean_dist_Naive[[psi]]$Annotation == gene)]<-gff3_gene$length[which(gff3_gene$gene_name == gene)]
          
        }
        
      }
        
        mean_dist_naive_all<-rbind(mean_dist_Naive[["psi2"]], mean_dist_Naive[["psi3"]], mean_dist_Naive[["psi4"]])
        mean_dist_naive_all <- mean_dist_naive_all[mean_dist_naive_all$Annotation != "", ]
        
        # Calculate Pearson correlation on log10-transformed values
        cor_value <- cor(log10(mean_dist_naive_all$length), log10(mean_dist_naive_all$length_transcripts))
        # Calculate R²
        r2_value <- (cor_value^2)* 100
        
        
        p<-ggplot(mean_dist_naive_all, aes(x = log10(length), y = log10(length_transcripts))) +
          geom_point() +
          geom_smooth(method = "lm", color = "orange", se = FALSE) +  
          labs(
            x = "log10(Mean distance between hypermodified type II sites)",
            y = "log10(length of transcript)"
          ) +
          annotate("text", 
                   x = min(log10(mean_dist_naive_all$length)), 
                   y = max(log10(mean_dist_naive_all$length_transcripts)), 
                   label = paste("R² =", round(r2_value, 3), "%"), 
                   hjust = 0, size = 5) +
          theme_minimal()
        ggsave("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Correlation_length_distance_Naive.pdf", width = 8, height = 4)
      
    }
    
    #Jurkat
    if(T){
      
      for (i in c(1:length(mean_dist_Jurkat)+1)) {
        psi <- paste0("psi", i)
        
        mean_dist_Jurkat[[psi]]$length_transcripts<-0
        
        uniqueAnnotations <- unique(mean_dist_Jurkat[[psi]]$Annotation)
        
        #CALCULATING LENGTH OF EACH TRANSCRIPT
        for (j in c(2:length(uniqueAnnotations))){
          gene<-uniqueAnnotations[j]
          mean_dist_Jurkat[[psi]]$length_transcripts[which(mean_dist_Jurkat[[psi]]$Annotation == gene)]<-gff3_gene$length[which(gff3_gene$gene_name == gene)]
          
        }
        
      }
      
      mean_dist_jurkat_all<-rbind(mean_dist_Jurkat[["psi2"]], mean_dist_Jurkat[["psi3"]], mean_dist_Jurkat[["psi4"]])
      mean_dist_jurkat_all <- mean_dist_jurkat_all[mean_dist_jurkat_all$Annotation != "", ]
      
      # Calculate Pearson correlation on log10-transformed values
      cor_value <- cor(log10(mean_dist_jurkat_all$length), log10(mean_dist_jurkat_all$length_transcripts))
      # Calculate R²
      r2_value <- (cor_value^2)* 100
      
      
      p<-ggplot(mean_dist_jurkat_all, aes(x = log10(length), y = log10(length_transcripts))) +
        geom_point() +
        geom_smooth(method = "lm", color = "royalblue", se = FALSE) +  
        labs(
          x = "log10(Mean distance between hypermodified type II sites)",
          y = "log10(length of transcript)"
        ) +
        annotate("text", 
                 x = min(log10(mean_dist_jurkat_all$length)), 
                 y = max(log10(mean_dist_jurkat_all$length_transcripts)), 
                 label = paste("R² =", round(r2_value, 3), "%"), 
                 hjust = 0, size = 5) +
        theme_minimal()
      ggsave("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/Correlation_length_distance_Jurkat.pdf", width = 8, height = 4)
      
    }
    
  }
  
  
  #CDS/UTR localization of sites
  if(T){
    
    #Naive
    par(mfrow = c(3, 3)) 
    for (i in c(1:length(ACP_h2_Naive)+1)) {
      psi <- paste0("psi", i)
      
      uniqueACP <- unique(ACP_h2_Naive[[psi]])
      
      filtered_Naive<-NaiveBootstrap[which(NaiveBootstrap$ACP %in% ACP_h2_Naive[[psi]]),]
      
      #cds/utr assessment
      if(T){
        filtered_Naive$CDS = F
        filtered_Naive$Three_prime_UTR = F
        filtered_Naive$Five_Prime_UTR = F
        filtered_Naive$start_codon = F
        filtered_Naive$stop_codon = F
        filtered_Naive$stop_codon_redefined_as_selenocysteine = F
        
        for (row in 1:nrow(filtered_Naive)){
          
          type = gff3_cds_utr$type[which(gff3_cds_utr$start <= filtered_Naive$position[row] & 
                                           gff3_cds_utr$end >= filtered_Naive$position[row] &
                                           gff3_cds_utr$seqid == filtered_Naive$chr[row])]
          
          if (length(type)>0){
            if ("CDS" %in% type) {
              filtered_Naive$CDS[row] = T
            }
            if ("three_prime_UTR" %in% type) {
              filtered_Naive$Three_prime_UTR[row] = T
            }
            if ("five_prime_UTR" %in% type) {
              filtered_Naive$Five_Prime_UTR[row] = T
            }
            if ("start_codon" %in% type) {
              filtered_Naive$start_codon[row] = T
            }
            if ("stop_codon" %in% type) {
              filtered_Naive$stop_codon[row] = T
            }
            if ("stop_codon_redefined_as_selenocysteine" %in% type) {
              filtered_Naive$stop_codon_redefined_as_selenocysteine[row] = T
            }
          }
        }
        
        # Add a new column 'position' based on the last six columns
        filtered_Naive$region <- apply(filtered_Naive[, c("CDS", "Three_prime_UTR", "Five_Prime_UTR", 
                                                            "start_codon", "stop_codon", 
                                                            "stop_codon_redefined_as_selenocysteine")], 
                                        1, function(row) {
                                          if (row["CDS"]) {
                                            "CDS"
                                          } else if (row["Three_prime_UTR"]) {
                                            "3UTR"
                                          } else if (row["Five_Prime_UTR"]) {
                                            "5UTR"
                                          } else if (row["start_codon"]) {
                                            "START"
                                          } else if (row["stop_codon"]) {
                                            "STOP"
                                          } else if (row["stop_codon_redefined_as_selenocysteine"]) {
                                            "STOP_Selenocysteine"
                                          } else {
                                            NA  # If none are TRUE
                                          }
                                        })
        
        
        
      }
      

      cds_utr_vector <- table(na.omit(filtered_Naive$region))
      
      # Create a pie chart
      pie(
        cds_utr_vector, 
        labels = paste(names(cds_utr_vector), " " ,cds_utr_vector), 
        main = psi, 
        col = c("yellow2", "gold", "orange2", "red", "forestgreen", "gray"), 
        #col = c("lightblue", "royalblue", "darkblue", "red", "forestgreen", "gray"), 
        border = "white"
      )
    }
    
    #Jurkat
    par(mfrow = c(3, 3)) 
    for (i in c(1:length(ACP_h2_Jurkat)+1)) {
      psi <- paste0("psi", i)
      
      uniqueACP <- unique(ACP_h2_Jurkat[[psi]])
      
      filtered_Jurkat<-JurkatBootstrap[which(JurkatBootstrap$ACP %in% ACP_h2_Jurkat[[psi]]),]
      
      #cds/utr assessment
      if(T){
        filtered_Jurkat$CDS = F
        filtered_Jurkat$Three_prime_UTR = F
        filtered_Jurkat$Five_Prime_UTR = F
        filtered_Jurkat$start_codon = F
        filtered_Jurkat$stop_codon = F
        filtered_Jurkat$stop_codon_redefined_as_selenocysteine = F
        
        for (row in 1:nrow(filtered_Jurkat)){
          
          type = gff3_cds_utr$type[which(gff3_cds_utr$start <= filtered_Jurkat$position[row] & 
                                           gff3_cds_utr$end >= filtered_Jurkat$position[row] &
                                           gff3_cds_utr$seqid == filtered_Jurkat$chr[row])]
          
          if (length(type)>0){
            if ("CDS" %in% type) {
              filtered_Jurkat$CDS[row] = T
            }
            if ("three_prime_UTR" %in% type) {
              filtered_Jurkat$Three_prime_UTR[row] = T
            }
            if ("five_prime_UTR" %in% type) {
              filtered_Jurkat$Five_Prime_UTR[row] = T
            }
            if ("start_codon" %in% type) {
              filtered_Jurkat$start_codon[row] = T
            }
            if ("stop_codon" %in% type) {
              filtered_Jurkat$stop_codon[row] = T
            }
            if ("stop_codon_redefined_as_selenocysteine" %in% type) {
              filtered_Jurkat$stop_codon_redefined_as_selenocysteine[row] = T
            }
          }
        }
        
        # Add a new column 'position' based on the last six columns
        filtered_Jurkat$region <- apply(filtered_Jurkat[, c("CDS", "Three_prime_UTR", "Five_Prime_UTR", 
                                                  "start_codon", "stop_codon", 
                                                  "stop_codon_redefined_as_selenocysteine")], 
                                    1, function(row) {
                                      if (row["CDS"]) {
                                        "CDS"
                                      } else if (row["Three_prime_UTR"]) {
                                        "3UTR"
                                      } else if (row["Five_Prime_UTR"]) {
                                        "5UTR"
                                      } else if (row["start_codon"]) {
                                        "START"
                                      } else if (row["stop_codon"]) {
                                        "STOP"
                                      } else if (row["stop_codon_redefined_as_selenocysteine"]) {
                                        "STOP_Selenocysteine"
                                      } else {
                                        NA  # If none are TRUE
                                      }
                                    })
        
        
        
      }
      

    }
 
  }
    
  #GO hyp2 psi with >30%mmto use in EnrichR
  if(T){
    more30_Jurkat<-unique(data.frame("Annotation"=c(mean_dist_Jurkat[["psi2"]]$Annotation,mean_dist_Jurkat[["psi3"]]$Annotation,mean_dist_Jurkat[["psi4"]]$Annotation)))
    more30_Naive<-unique(data.frame("Annotation"=c(mean_dist_Naive[["psi2"]]$Annotation,mean_dist_Naive[["psi3"]]$Annotation,mean_dist_Naive[["psi4"]]$Annotation)))
  }
  
  #tables_hypermod2 with mm>30 and ternary plots
  if(T){
   
    #Jurkat cds/utr and ternary plots
    if(T){
      uniqueJurkat_ACP<-unique(unlist(ACP_h2_Jurkat))
      JurkatListHypII<-JurkatBootstrap[which(JurkatBootstrap$ACP %in% uniqueJurkat_ACP),]
      
      #At least one replicate has >30%mm in that position
      JurkatListHypII <- JurkatListHypII %>%
        rowwise() %>%  # Treat each row independently
        filter(sum(c_across(starts_with("mm.J")) > 30) >= 1) %>%
        ungroup()  # Reset grouping
      
      JurkatListHypII <- JurkatListHypII %>%
        group_by(Annotation) %>%                # Group by the Annotation column
        filter(n() > 2) %>%                     # Keep groups with more than one row
        ungroup()                               # Remove grouping
      
    filtered_Jurkat<-JurkatListHypII
    #cds/utr assessment
    if(T){
      filtered_Jurkat$CDS = F
      filtered_Jurkat$Three_prime_UTR = F
      filtered_Jurkat$Five_Prime_UTR = F
      filtered_Jurkat$start_codon = F
      filtered_Jurkat$stop_codon = F
      filtered_Jurkat$stop_codon_redefined_as_selenocysteine = F
      
      for (row in 1:nrow(filtered_Jurkat)){
        
        type = gff3_cds_utr$type[which(gff3_cds_utr$start <= filtered_Jurkat$position[row] & 
                                         gff3_cds_utr$end >= filtered_Jurkat$position[row] &
                                         gff3_cds_utr$seqid == filtered_Jurkat$chr[row])]
        
        if (length(type)>0){
          if ("CDS" %in% type) {
            filtered_Jurkat$CDS[row] = T
          }
          if ("three_prime_UTR" %in% type) {
            filtered_Jurkat$Three_prime_UTR[row] = T
          }
          if ("five_prime_UTR" %in% type) {
            filtered_Jurkat$Five_Prime_UTR[row] = T
          }
          if ("start_codon" %in% type) {
            filtered_Jurkat$start_codon[row] = T
          }
          if ("stop_codon" %in% type) {
            filtered_Jurkat$stop_codon[row] = T
          }
          if ("stop_codon_redefined_as_selenocysteine" %in% type) {
            filtered_Jurkat$stop_codon_redefined_as_selenocysteine[row] = T
          }
        }
      }
      
      # Add a new column 'position' based on the last six columns
      filtered_Jurkat$region <- apply(filtered_Jurkat[, c("CDS", "Three_prime_UTR", "Five_Prime_UTR", 
                                                          "start_codon", "stop_codon", 
                                                          "stop_codon_redefined_as_selenocysteine")], 
                                      1, function(row) {
                                        if (row["CDS"]) {
                                          "CDS"
                                        } else if (row["Three_prime_UTR"]) {
                                          "3UTR"
                                        } else if (row["Five_Prime_UTR"]) {
                                          "5UTR"
                                        } else if (row["start_codon"]) {
                                          "START"
                                        } else if (row["stop_codon"]) {
                                          "STOP"
                                        } else if (row["stop_codon_redefined_as_selenocysteine"]) {
                                          "STOP_Selenocysteine"
                                        } else {
                                          NA  # If none are TRUE
                                        }
                                      })
      
      
      
    }
    JurkatListHypII<-filtered_Jurkat
    
    write.csv(JurkatListHypII,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/JurkatHypermodifiedII.csv", row.names = FALSE)
    
    #ternary plot
    if(T){
      
      # Load the data
      data <- read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/JurkatHypermodifiedII.csv",header=TRUE) #JurkatListHypII
      
      library(ggtern)
      # Preprocess the data
      target_kmers <- c("GTTCA", "GTTCC", "GTTCG", "GTTCT")
      target_kmers_green <- c("TGTAG", "TCTAG", "TATAG", "TTTAG", "TGTAA", "TCTAA", "TATAA", "TTTAA")
      if(T){
        # Create the 'highlight' column first
        data <- data %>%
          mutate(
            highlight = case_when(
              kmer %in% target_kmers_red ~ "red",
              kmer %in% target_kmers_green ~ "green",
              TRUE ~ "gray"
            )
          )
        
        # Now process the data and create 'ternary_data'
        ternary_data <- data %>%
          group_by(Annotation, region) %>%
          summarize(count = n(), .groups = "drop") %>%
          pivot_wider(names_from = region, values_from = count, values_fill = 0) %>%
          mutate(
            total = `3UTR` + CDS + `5UTR`,
            `3UTR_frac` = `3UTR` / total,
            `CDS_frac` = CDS / total,
            `5UTR_frac` = `5UTR` / total
          )
        
        # Now join the 'highlight' column back using the 'Annotation' column
        ternary_data <- merge(
          ternary_data, 
          data[, c("Annotation", "highlight")], 
          by = "Annotation", 
          all.x = TRUE
        )
        
        ternary_data<-ternary_data %>% filter(total != 1)
        # Reorder ternary_data to have rows with "red" in highlight last
        ternary_data <- ternary_data[order(ternary_data$highlight == "green", 
                                           ternary_data$highlight == "red"), ]
      }
     
      
      # Plot the ternary diagram
      if(T){
      ggtern(data = ternary_data, aes(x = `3UTR_frac`, y = `CDS_frac`, z = `5UTR_frac`)) +
        geom_point(aes(size = total, color=highlight), alpha = 0.5, stroke=0) +
        #geom_text(aes(label = Annotation), hjust = 0.5, vjust = 1.5, size=3) +
        geom_text(aes(label = Annotation), vjust = 1.2, hjust = 0.5, check_overlap = TRUE, size=3) +
        labs(
          title = "Ternary Plot of Regions",
          x = "3' UTR",
          y = "CDS",
          z = "5' UTR",
          size = "Frequency"
        ) +
        theme_minimal() +
        scale_size_continuous(range = c(2, 6)) +
        scale_color_manual(values = c("red" = rgb(1, 0, 0, alpha = 0.5),   # Red with alpha = 0.5
                                        "black" = rgb(0, 0, 0, alpha = 0.5),
                                        "green" = rgb(0, 1, 0, alpha = 0.1))) +  # Black with alpha = 0.8
        limit_tern(T = 1.2, L = 1.2, R = 1.2)
      }
      
    }
    
  }
    
    #Naive cds/utr and ternary plots
    if(T){
    uniqueNaive_ACP<-unique(unlist(ACP_h2_Naive))
    NaiveListHypII<-NaiveBootstrap[which(NaiveBootstrap$ACP %in% uniqueNaive_ACP),]
    
    #At least one replicate has >30%mm in that position
    NaiveListHypII <- NaiveListHypII %>%
      rowwise() %>%  # Treat each row independently
      filter(sum(c_across(starts_with("mm.N")) > 30) >= 1) %>%
      ungroup()  # Reset grouping
    
    NaiveListHypII <- NaiveListHypII %>%
      group_by(Annotation) %>%                # Group by the Annotation column
      filter(n() > 1) %>%                     # Keep groups with more than one row
      ungroup()                               # Remove grouping
    
    filtered_Naive<-NaiveListHypII
    #cds/utr
    if(T){
      filtered_Naive$CDS = F
      filtered_Naive$Three_prime_UTR = F
      filtered_Naive$Five_Prime_UTR = F
      filtered_Naive$start_codon = F
      filtered_Naive$stop_codon = F
      filtered_Naive$stop_codon_redefined_as_selenocysteine = F
      
      for (row in 1:nrow(filtered_Naive)){
        
        type = gff3_cds_utr$type[which(gff3_cds_utr$start <= filtered_Naive$position[row] & 
                                         gff3_cds_utr$end >= filtered_Naive$position[row] &
                                         gff3_cds_utr$seqid == filtered_Naive$chr[row])]
        
        if (length(type)>0){
          if ("CDS" %in% type) {
            filtered_Naive$CDS[row] = T
          }
          if ("three_prime_UTR" %in% type) {
            filtered_Naive$Three_prime_UTR[row] = T
          }
          if ("five_prime_UTR" %in% type) {
            filtered_Naive$Five_Prime_UTR[row] = T
          }
          if ("start_codon" %in% type) {
            filtered_Naive$start_codon[row] = T
          }
          if ("stop_codon" %in% type) {
            filtered_Naive$stop_codon[row] = T
          }
          if ("stop_codon_redefined_as_selenocysteine" %in% type) {
            filtered_Naive$stop_codon_redefined_as_selenocysteine[row] = T
          }
        }
      }
      
      # Add a new column 'position' based on the last six columns
      filtered_Naive$region <- apply(filtered_Naive[, c("CDS", "Three_prime_UTR", "Five_Prime_UTR", 
                                                        "start_codon", "stop_codon", 
                                                        "stop_codon_redefined_as_selenocysteine")], 
                                     1, function(row) {
                                       if (row["CDS"]) {
                                         "CDS"
                                       } else if (row["Three_prime_UTR"]) {
                                         "3UTR"
                                       } else if (row["Five_Prime_UTR"]) {
                                         "5UTR"
                                       } else if (row["start_codon"]) {
                                         "START"
                                       } else if (row["stop_codon"]) {
                                         "STOP"
                                       } else if (row["stop_codon_redefined_as_selenocysteine"]) {
                                         "STOP_Selenocysteine"
                                       } else {
                                         NA  # If none are TRUE
                                       }
                                     })
      
      
      
    }
    NaiveListHypII<-filtered_Naive
    
    write.csv(NaiveListHypII,"/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveHypermodifiedII.csv", row.names = FALSE)
    
    #ternary plot
    if(T){
      
      # Load the data
      data <- read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/Tcell/NaiveHypermodifiedII.csv",header=TRUE) #JurkatListHypII
      
      library(ggtern)
      # Preprocess the data
      target_kmers_red <- c("GTTCA", "GTTCC", "GTTCG", "GTTCT")
      target_kmers_green <- c("TGTAG", "TCTAG", "TATAG", "TTTAG", "TGTAA", "TCTAA", "TATAA", "TTTAA")

      if(T){
        # Create the 'highlight' column first
        data <- data %>%
          mutate(
            highlight = case_when(
              kmer %in% target_kmers_red ~ "red",
              kmer %in% target_kmers_green ~ "green",
              TRUE ~ "gray"
            )
          )
        
        # Now process the data and create 'ternary_data'
        ternary_data <- data %>%
          group_by(Annotation, region) %>%
          summarize(count = n(), .groups = "drop") %>%
          pivot_wider(names_from = region, values_from = count, values_fill = 0) %>%
          mutate(
            total = `3UTR` + CDS + `5UTR`,
            `3UTR_frac` = `3UTR` / total,
            `CDS_frac` = CDS / total,
            `5UTR_frac` = `5UTR` / total
          )
        
        # Now join the 'highlight' column back using the 'Annotation' column
        ternary_data <- merge(
          ternary_data, 
          data[, c("Annotation", "highlight")], 
          by = "Annotation", 
          all.x = TRUE
        )
        
        ternary_data<-ternary_data %>% filter(total != 1)
        # Reorder ternary_data to have rows with "red" in highlight last
        # Reorder ternary_data to have rows with "red" and "green" highlights last
        ternary_data <- ternary_data[order(ternary_data$highlight == "red", 
                                           ternary_data$highlight == "green"), ]
      }
      
      
      # Plot the ternary diagram
      if(T){
        ggtern(data = ternary_data, aes(x = `3UTR_frac`, y = `CDS_frac`, z = `5UTR_frac`)) +
          geom_point(aes(size = total, color=highlight), alpha = 0.5, stroke=0) +
          #geom_text(aes(label = Annotation), hjust = 0.5, vjust = 1.5, size=3) +
          geom_text(aes(label = Annotation), vjust = 1.2, hjust = 0.5, check_overlap = TRUE, size=3) +
          labs(
            title = "Ternary Plot of Regions",
            x = "3' UTR",
            y = "CDS",
            z = "5' UTR",
            size = "Frequency"
          ) +
          theme_minimal() +
          scale_size_continuous(range = c(2, 6)) +
          scale_color_manual(values = c("red" = rgb(1, 0, 0, alpha = 0.5),   # Red with alpha = 0.5
                                        "black" = rgb(0, 0, 0, alpha = 0.5),
                                        "green" = "green2")) +  # Black with alpha = 0.8
          limit_tern(T = 1.2, L = 1.2, R = 1.2)
      }
      
    }
    
  }
  
}

}

