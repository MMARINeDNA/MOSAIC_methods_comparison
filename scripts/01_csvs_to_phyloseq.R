#This script takes the results from the MURI bioinformatic pipeline from two sequencing runs using the same datasets and merges them, brings in the metadata, and creates a phyloseq object (without sequences or ASV-level differentiation, just to the BestTaxon level)

library(phyloseq)
library(tidyverse)

`%ni%` = Negate(`%in%`)

#load data
metadata<-read.csv("data/SR2311_SR2408_eDNA_metadata.csv")
GGBC1<-read.csv("data/taxon_table_lineage_wide.csv")
GGBC2Fs<-read.csv("data/taxon_table_wide-SR2408-GGBC2.csv")
GGBC2all<-read.csv(file = "data/taxon_table_wide - taxon_table_wide.csv")

#associate metadata with GGBC1
#GGBC1 has dots and GGBC2 has underlines separating the sampleID and the rep number

metadata<-metadata%>%
  mutate(SequencingRun=gsub("-","",SequencingRun))%>%
  mutate(Longitude.start=gsub("NC",NA,Longitude.start))%>%
  mutate(Longitude.end=gsub("NC",NA,Longitude.end))%>%
  mutate(Latitude.start=gsub("NC",NA,Latitude.start))%>%
  mutate(Latitude.end=gsub("NC",NA,Latitude.end))

#add placeholder OTU names to dataframes so I can split out the taxonomy table for phyloseq
  GGBC1_OTUs <- c(sprintf("GGBC1_OTU_%02d", seq(1,length(GGBC1$BestTaxon))))
rownames(GGBC1) <- GGBC1_OTUs

GGBC1_taxtab<-GGBC1%>%
  dplyr::select(c(Kingdom,Phylum,Class,Order,Family,Genus,Species,BestTaxon))
GGBC1_taxtab<-as.matrix(GGBC1_taxtab)

seqrun1<-as.data.frame(t(GGBC1))
seqrun1<-seqrun1%>%
  rownames_to_column(var = "rowname")%>%
  mutate(repID=gsub(".","-",rowname,fixed=TRUE))%>%
  filter(repID %in% metadata$SequencingID)%>%
  mutate(FilterID=str_split_i(pattern = "-",i = 1,repID))

GGBC1_ready<-column_to_rownames(seqrun1, var = "repID")
GGBC1_ready<-GGBC1_ready%>%
  dplyr::select(-c("rowname","FilterID"))

GGBC1_otutab<-as.matrix(t(GGBC1_ready))
class(GGBC1_otutab) <- "numeric"
  
GGBC1_sampletab<-metadata%>%
  filter(SequencingRun=="GGBC1")%>%
  column_to_rownames(var="SequencingID")

#Make Fall GGBC1 phyloseq object
GGBC1.ps<-phyloseq(otu_table(GGBC1_otutab,taxa_are_rows = T), sample_data(GGBC1_sampletab),tax_table(GGBC1_taxtab))

GGBC1.ps

#GGBC2



#add placeholder OTU names to dataframes so I can split out the taxonomy table for phyloseq
GGBC2_F_OTUs <- c(sprintf("GGBC2_F_OTU_%02d", seq(1,length(GGBC2Fs$BestTaxon))))
rownames(GGBC2Fs) <- GGBC2_F_OTUs

GGBC2Fs_taxtab<-GGBC2Fs%>%
  dplyr::select(c(Class,BestTaxon))

GGBC2Fs_taxtab<-as.matrix(GGBC2Fs_taxtab)
  

seqrun2Fs<-as.data.frame(t(GGBC2Fs))


seqrun2Fs<-seqrun2Fs%>%
  rownames_to_column(var="X")%>%
  mutate(repID=gsub(".","-",X,fixed=TRUE))%>%
  filter(repID %in% metadata$SequencingID)%>%
  mutate(FilterID=str_split_i(pattern = "-",i = 1,repID))

GGBC2_Fs_otutab<-seqrun2Fs%>%
  column_to_rownames(var = "repID")%>%
  dplyr::select(-c(X,FilterID))

GGBC2_Fs_otutab<-as.matrix(GGBC2_Fs_otutab)
class(GGBC2_Fs_otutab)<-"numeric"

GGBC2_Fs_sampletab<-metadata%>%
  filter(FilterID %in% seqrun2Fs$FilterID)%>%
  filter(SequencingRun=="GGBC2")%>%
  column_to_rownames(var = "SequencingID")

#make phyloseq object of GGBC2 F series
GGBC2Fs.ps<-phyloseq(otu_table(GGBC2_Fs_otutab,taxa_are_rows = F), sample_data(GGBC2_Fs_sampletab),tax_table(GGBC2Fs_taxtab))

GGBC2Fs.ps

testmerge.ps<-merge_phyloseq(GGBC1.ps,GGBC2Fs.ps)
sample_data(testmerge.ps)

#GGBC2 C and Fall series

GGBC2_OTUs <- c(sprintf("GGBC2_OTU_%02d", seq(1,length(GGBC2all$BestTaxon))))
rownames(GGBC2all) <- GGBC2_OTUs

GGBC2_taxtab<-GGBC2all%>%
  dplyr::select(c(Class,BestTaxon))
GGBC2_taxtab<-as.matrix(GGBC2_taxtab)


GGBC2_CandFall<-as.data.frame(t(GGBC2all))%>%
  rownames_to_column(var = "X")%>%
  mutate(ID=gsub("_","-",X,fixed=TRUE))%>%
  mutate(FilterID=str_split_i(pattern = "-",i = 1,ID))%>%
  filter(FilterID %ni% seqrun2Fs$FilterID)%>%
  filter(FilterID %in% metadata$FilterID)%>%
  mutate(Tech.rep=str_split_i(pattern = "-",i = 2,ID))%>%
  mutate(AnalysisID=paste0(FilterID,"-",as.numeric(Tech.rep)))


GGBC2_otutab<-column_to_rownames(GGBC2_CandFall, var = "AnalysisID")
GGBC2_otutab<-GGBC2_otutab%>%
  dplyr::select(-c("X","ID","FilterID","Tech.rep"))

GGBC2_otutab<-as.matrix(t(GGBC2_otutab))
class(GGBC2_otutab) <- "numeric"

GGBC2_samples<-metadata%>%
  filter(FilterID%in%GGBC2_CandFall$FilterID)

GGBC2_sampletab<-left_join(GGBC2_CandFall, GGBC2_samples, by = "FilterID",multiple = "any")

GGBC2_sampletab<-GGBC2_sampletab%>%
  dplyr::select(AnalysisID,Experiment,Cruise,Year,Month,Date..YYYYMMDD.,SampleID,FilterID,Bio.rep,Line,Station,Collection.method,Sample.type,Collection.device,CTD.Cast,CTD.Bottle.Number.s.,Depth..m.,Filtration.time.start,Latitude.start,Longitude.start,Station.start,Filtration.time.end,Latitude.end,Longitude.end,Station.end,Volume.filtered..mL.,Marker,Filter.type,Proportion.of.filter.extracted)%>%
  mutate(Tech.rep=str_split_i(pattern = "-",i = 2,AnalysisID))%>%
  mutate(SequencingRun="GGBC2")%>%
  column_to_rownames(var="AnalysisID")

#make phyloseq object of GGBC2 Fall and C series
GGBC2.ps<-phyloseq(otu_table(GGBC2_otutab,taxa_are_rows = T), sample_data(GGBC2_sampletab),tax_table(GGBC2_taxtab))

GGBC2.ps

MOSAICpapermerge.ps<-merge_phyloseq(GGBC2.ps,GGBC1.ps,GGBC2Fs.ps)

MOSAIC.ps<-tax_glom(MOSAICpapermerge.ps,taxrank = "BestTaxon")

MOSAIC.ps

#fill in missing taxonomy information

library(taxonomizr)

prepareDatabase('accessionTaxa.sql')

combined_taxtab<-as.data.frame(tax_table(MOSAIC.ps))
besttaxa<-combined_taxtab$BestTaxon

taxaId<-getId(besttaxa,'accessionTaxa.sql')
print(taxaId)
taxa<-getTaxonomy(taxaId,'accessionTaxa.sql')
print(taxa)

OTU_names<-rownames(combined_taxtab)
rownames(taxa)<-OTU_names
taxa<-as.data.frame(taxa)%>%
  dplyr::select(-superkingdom)
taxa<-as.matrix(taxa)
tax_table(MOSAIC.ps)<-taxa



saveRDS(MOSAIC.ps,file = "outputs/MOSAIC.ps.RDS")



  
  
