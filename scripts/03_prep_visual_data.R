library(tidyverse)
library(phyloseq)
library(readxl)
library(taxonomizr)
`%ni%` = Negate(`%in%`)

# NCBI taxonomy database (large download on first use); set ACCESSION_TAXA_SQL to reuse an existing one
taxdb <- Sys.getenv("ACCESSION_TAXA_SQL", "accessionTaxa.sql")
if (!file.exists(taxdb)) prepareDatabase(taxdb, getAccessions = FALSE)

#read in data

sightings<-read_xlsx("visual_data/CalCOFI_sightings_CC2023-11_CC2024-08.xlsx")
segments<-read_xlsx("visual_data/CalCOFI_segments_CC2023-11_CC2024-08.xlsx")
common_names<-read_csv("visual_data/CalCOFI_cetacean_codes (1).csv")
sciname_codes<-read_xlsx("visual_data/species_codes_with_scinames.xlsx")

#reconcile variables

sightings$SegmentID<-sightings$segID_ON_Effort_MNA
sciname_codes$Species1<-sciname_codes$`CalCOFI/NOAA Species Abbreviation`


#select only variables for this analysis

sightings<-sightings%>%
  select(c(Cruise,Date_Local_,Date_Time_Local_,DecLat,DecLong,SightingNo,UniqueID,Species1,SegmentID,Best))

#combine dataframes

combo1<-left_join(sightings, segments)
visuals<-left_join(combo1,sciname_codes)

#assign taxonomy (skipped when visual_data/visuals_for_MOSAIC.csv exists: it is not read again in this script,
#and the lookup downloads the NCBI taxonomy)
if (!file.exists("visual_data/visuals_for_MOSAIC.csv")) {
  scinames<-visuals$`Scientific name`
  taxaId<-getId(scinames,taxdb)
  print(taxaId)
  taxa<-getTaxonomy(taxaId,taxdb)
  print(taxa)

  taxa.df<-as.data.frame(taxa)

  visuals1<-cbind(visuals,taxa.df)
  row.names(visuals1)<-visuals$UniqueID

  write.csv(visuals1, "visual_data/visuals_for_MOSAIC.csv")
}

#in excel included taxonomic levels to family for "unidentified dolphin" =delphinidae, and "unidentified large whale" =balaenopteridae.

visuals<-read_xlsx("visual_data/species_codes_with_scinames_and_taxonomy.xlsx")



visuals_wide<-pivot_wider(visuals, names_from = "Species1", id_cols = "rownames",values_from = "Best")
visuals_wide<-(as.data.frame(visuals_wide))
visuals_wide<-column_to_rownames(visuals_wide,var = "rownames")

#11 taxa in the dataset (not including "other", 9 ID to genus level, 8 to species level)

#########

############
library(phyloseq)


vis_count_table <- table(visuals$Species1,visuals$rownames)
vis_count_table<-as.matrix(vis_count_table)

visuals_tax<-sciname_codes
visuals_tax<-column_to_rownames(visuals_tax, var = "Species1")
visuals_tax<-visuals_tax[-1]
tax2ID<-visuals_tax$`Scientific name`

taxaIds1<-getId(tax2ID,taxdb)
print(taxaIds1)
taxa1<-getTaxonomy(taxaIds1,taxdb)
print(taxa1)

vis_tax<-visuals[22:31]
vis_tax<- unique(vis_tax)
vis_tax<-column_to_rownames(vis_tax, var = "CalCOFI/NOAA Species Abbreviation")
vis_tax<-vis_tax[3:9]
vis_tax<-as.matrix(vis_tax)

visuals_obs<-visuals[1:21]
visuals_obs<-column_to_rownames(visuals_obs,var = "rownames")
visuals_obs<-as.data.frame(visuals_obs)

# visuals.ps<-phyloseq(otu_table(vis_count_table,taxa_are_rows = T), sample_data(visuals_obs),tax_table(vis_tax))
# 
# OTU = otu_table(vis_count_table, taxa_are_rows = TRUE)
# TAX = tax_table(vis_tax)
# OTU


#set equivalent variable names with eDNA data
Mammals.pos.ps<-readRDS("outputs/Mammals.pos.ps.RDS")
Mammals.pos.ps <- prune_taxa(taxa_sums(Mammals.pos.ps) > 0, Mammals.pos.ps)
Mammals.pos.ps<-prune_samples(sample_sums(Mammals.pos.ps)>0,Mammals.pos.ps)
Mammals.pos.ps
Mammals_sample_data<-cbind(as.data.frame(sample_data(Mammals.pos.ps)))
Mammals_otu_tab<-cbind(as.data.frame(otu_table(Mammals.pos.ps)))

Mam_by_date<-merge_samples(Mammals.pos.ps, "Date..YYYYMMDD.",fun = "sum")
Mam_by_date
Mam_by_date_otu_tab<-cbind(as.data.frame(otu_table(Mam_by_date)))
Mam_by_date_otu_tab1<-rownames_to_column(Mam_by_date_otu_tab,var = "date")
Mam_by_date_otu_tab1$occurrence_date<-ymd(Mam_by_date_otu_tab1$date)

Mam_DNA<-Mam_by_date_otu_tab1
Mam_sp<-colnames(Mam_DNA)
Mam_sp<-gsub(" 1","",Mam_sp)
colnames(Mam_DNA)<-Mam_sp

#marine mammals detected on all 30 days with eDNA (100%)
#mammals detected on 18 days with visuals (60%)


Mam_DNA<-Mam_DNA%>%
  mutate(survey="eDNA")

### filter eDNA for cetaceans only and MOSAIC only for best comparison with visuals ###

Cetaceans.ps<-subset_taxa(Mammals.pos.ps,order=="Artiodactyla")
Cetaceans.ps <- prune_taxa(taxa_sums(Cetaceans.ps) > 0, Cetaceans.ps)
Cetaceans.ps<-prune_samples(sample_sums(Cetaceans.ps)>0,Cetaceans.ps)


####
write.csv(visuals, "outputs/visuals_to_coalesce.csv")
  
visuals2<-read.csv("visual_data/visuals_coalesced.csv", row.names = "X")

visuals$visID<-visuals2$visID
  
  
visuals$occurrence_date <- ymd(visuals$Date_Local_)
vis_by_date<-visuals%>%
  group_by(occurrence_date,visID) %>%
  summarise(count = n()) %>%
  ungroup() # Ungroup to remove the grouping structure for subsequent operations

vis_by_date<-as.data.frame(vis_by_date)

vis_by_date_wide<-pivot_wider(vis_by_date,id_cols = "occurrence_date",names_from = visID,values_from = count)

vis_by_date_wide<-vis_by_date_wide%>%
  mutate(survey="visual")

####
visveDNA<-full_join(vis_by_date_wide,Mam_DNA)

write.csv(visveDNA,"outputs/visveDNA.csv")

#5 species were found in both (`Grampus griseus`, `Balaenoptera physalus`, `Tursiops truncatus`, `Balaenoptera musculus`, `Megaptera novaeangliae`)

#24 taxa for eDNA (16 to species level)
#11 for visuals (8 to species level)

#on the scale of a day:
#1 of 4 grampus griseus eDNA detections had a visual.  no visual GG without eDNA.

#B.physalus: 5 times in visuals without eDNA corroboration, 3 times in eDNA without visual corroboration, and 3 times in both sample types on the same day

#T.truncatus: 8 eDNA only, 2 visual only, 2 both

#B. musculus: 2 visual only, 1 eDNA only, 1 both

#M. novaeangliae: 3 visual only, 5 eDNA only, 2 both


#### Now do the same, but keep keep MOSAIC method D separate

#set equivalent variable names with eDNA data
Mammals.pos.ps<-readRDS("outputs/Mammals.pos.ps.RDS")
Mammals.pos.ps <- prune_taxa(taxa_sums(Mammals.pos.ps) > 0, Mammals.pos.ps)
Mammals.pos.ps<-prune_samples(sample_sums(Mammals.pos.ps)>0,Mammals.pos.ps)
Mammals.pos.ps

Mammals.pos.MOSAIC.ps<-prune_samples(Mammals.pos.ps@sam_data$Collection.device=="Underway system",Mammals.pos.ps)
Mammals.pos.CTD.ps<-prune_samples(Mammals.pos.ps@sam_data$Collection.device=="CTD",Mammals.pos.ps)


Cetaceans.pos.MOSAIC.ps<-prune_samples(Cetaceans.ps@sam_data$Collection.device=="Underway system",Cetaceans.ps)
Cetaceans.pos.CTD.ps<-prune_samples(Cetaceans.ps@sam_data$Collection.device=="CTD",Cetaceans.ps)



Mammals_CTD_sample_data<-cbind(as.data.frame(sample_data(Mammals.pos.CTD.ps)))
Mammals_CTD_otu_tab<-cbind(as.data.frame(otu_table(Mammals.pos.CTD.ps)))
Mammals_MOSAIC_sample_data<-cbind(as.data.frame(sample_data(Mammals.pos.MOSAIC.ps)))
Mammals_MOSAIC_otu_tab<-cbind(as.data.frame(otu_table(Mammals.pos.MOSAIC.ps)))
Cetaceans_MOSAIC_sample_data<-cbind(as.data.frame(sample_data(Cetaceans.pos.MOSAIC.ps)))
Cetaceans_MOSAIC_otu_tab<-cbind(as.data.frame(otu_table(Cetaceans.pos.MOSAIC.ps)))

Mam_by_date_CTD<-merge_samples(Mammals.pos.CTD.ps, "Date..YYYYMMDD.",fun = "sum")
Mam_by_date_CTD
Mam_by_date_CTD_otu_tab<-cbind(as.data.frame(otu_table(Mam_by_date_CTD)))
Mam_by_date_CTD_otu_tab1<-rownames_to_column(Mam_by_date_CTD_otu_tab,var = "date")
Mam_by_date_CTD_otu_tab1$occurrence_date<-ymd(Mam_by_date_CTD_otu_tab1$date)

Mam_CTD_DNA<-Mam_by_date_CTD_otu_tab1
Mam_CTD_sp<-colnames(Mam_CTD_DNA)
Mam_CTD_sp<-gsub("1","",Mam_CTD_sp)
colnames(Mam_CTD_DNA)<-Mam_CTD_sp


Mam_by_date_MOSAIC<-merge_samples(Mammals.pos.MOSAIC.ps, "Date..YYYYMMDD.",fun = "sum")
Mam_by_date_MOSAIC
Mam_by_date_MOSAIC_otu_tab<-cbind(as.data.frame(otu_table(Mam_by_date_MOSAIC)))
Mam_by_date_MOSAIC_otu_tab1<-rownames_to_column(Mam_by_date_MOSAIC_otu_tab,var = "date")
Mam_by_date_MOSAIC_otu_tab1$occurrence_date<-ymd(Mam_by_date_MOSAIC_otu_tab1$date)

Mam_MOSAIC_DNA<-Mam_by_date_MOSAIC_otu_tab1
Mam_MOSAIC_sp<-colnames(Mam_MOSAIC_DNA)
Mam_MOSAIC_sp<-gsub("1","",Mam_MOSAIC_sp)
colnames(Mam_MOSAIC_DNA)<-Mam_MOSAIC_sp

Mam_MOSAIC_DNA<-Mam_MOSAIC_DNA%>%
  mutate(survey="MOSAIC")%>%
  select(-date)

Mam_CTD_DNA<-Mam_CTD_DNA%>%
  mutate(survey="CTD")%>%
  select(-date)

eDNA_by_method<-full_join(Mam_MOSAIC_DNA,Mam_CTD_DNA)

#marine mammals detected on all 30 days with MOSAIC eDNA (100%)
#mammals detected on 18 days with visuals (60%) of 30 days on the boat
#mammals detected on 7 days with CTD out of 15 days of CTD sampling

#Cetaceans MOSAIC only:
Cet_by_date_MOSAIC<-merge_samples(Cetaceans.pos.MOSAIC.ps, "Date..YYYYMMDD.",fun = "sum")
Cet_by_date_MOSAIC
Cet_by_date_MOSAIC_otu_tab<-cbind(as.data.frame(otu_table(Cet_by_date_MOSAIC)))
Cet_by_date_MOSAIC_otu_tab1<-rownames_to_column(Cet_by_date_MOSAIC_otu_tab,var = "date")
Cet_by_date_MOSAIC_otu_tab1$occurrence_date<-ymd(Cet_by_date_MOSAIC_otu_tab1$date)

Cet_MOSAIC_DNA<-Cet_by_date_MOSAIC_otu_tab1
Cet_MOSAIC_sp<-colnames(Cet_MOSAIC_DNA)
Cet_MOSAIC_sp<-gsub("1","",Cet_MOSAIC_sp)
colnames(Cet_MOSAIC_DNA)<-Cet_MOSAIC_sp

Cet_MOSAIC_DNA<-Cet_MOSAIC_DNA%>%
  mutate(survey="MOSAIC")%>%
  select(-date)

######
Mam_DNA<-Mam_DNA%>%
  mutate(survey="eDNA")

  
  
Mam_by_date_melt<-psmelt(Mam_by_date)


####

visuals$occurrence_date <- ymd(visuals$Date_Local_)
vis_by_date<-visuals%>%
  group_by(occurrence_date,`Scientific name`) %>%
  summarise(count = n()) %>%
  ungroup() # Ungroup to remove the grouping structure for subsequent operations

vis_by_date<-as.data.frame(vis_by_date)

vis_sci_name<-vis_by_date$`Scientific name`

vis_by_date$species<-vis_sci_name

vis_by_date_wide<-pivot_wider(vis_by_date,id_cols = "occurrence_date",names_from = species,values_from = count)

vis_by_date_wide<-vis_by_date_wide%>%
  mutate(survey="visual")

####
visvMOSAIC_eDNA<-full_join(vis_by_date_wide,Mam_MOSAIC_DNA)

write.csv(visvMOSAIC_eDNA,"outputs/visvMOSAIC_eDNA.csv")


visveDNA_method<-full_join(vis_by_date_wide,eDNA_by_method)

visveDNA_cet_MOSAIC<-full_join(vis_by_date_wide,Cet_MOSAIC_DNA)


#export csv 
write.csv(visveDNA_cet_MOSAIC,"outputs/Cetacean_MOSAIC_vs_visuals.csv")

##### Now do the same, but keep keep all eDNA methods separate

#set equivalent variable names with eDNA data
Mammals.pos.ps<-readRDS("outputs/Mammals.pos.ps.RDS")
Mammals.pos.ps <- prune_taxa(taxa_sums(Mammals.pos.ps) > 0, Mammals.pos.ps)
Mammals.pos.ps<-prune_samples(sample_sums(Mammals.pos.ps)>0,Mammals.pos.ps)
Mammals.pos.ps

Mammals.pos.MOSAICA.ps<-prune_samples(Mammals.pos.ps@sam_data$Collection.method=="MOSAIC-A",Mammals.pos.ps)
Mammals.pos.MOSAICB.ps<-prune_samples(Mammals.pos.ps@sam_data$Collection.method=="MOSAIC-B",Mammals.pos.ps)
Mammals.pos.MOSAICC.ps<-prune_samples(Mammals.pos.ps@sam_data$Collection.method=="MOSAIC-C",Mammals.pos.ps)
Mammals.pos.CTD.ps<-prune_samples(Mammals.pos.ps@sam_data$Collection.device=="CTD",Mammals.pos.ps)

Mam_by_date_MOSAICA<-merge_samples(Mammals.pos.MOSAICA.ps, "Date..YYYYMMDD.",fun = "sum")
Mam_by_date_MOSAICA
Mam_by_date_MOSAICA_otu_tab<-cbind(as.data.frame(otu_table(Mam_by_date_MOSAICA)))
Mam_by_date_MOSAICA_otu_tab1<-rownames_to_column(Mam_by_date_MOSAICA_otu_tab,var = "date")
Mam_by_date_MOSAICA_otu_tab1$occurrence_date<-ymd(Mam_by_date_MOSAICA_otu_tab1$date)

Mam_MOSAICA_DNA<-Mam_by_date_MOSAICA_otu_tab1
Mam_MOSAICA_sp<-colnames(Mam_MOSAICA_DNA)
Mam_MOSAICA_sp<-gsub("1","",Mam_MOSAICA_sp)
colnames(Mam_MOSAIC_DNA)<-Mam_MOSAICA_sp

Mam_MOSAICA_DNA<-Mam_MOSAICA_DNA%>%
  mutate(method="MOSAIC-A")%>%
  select(-date)

visvMOSAICA_eDNA<-full_join(vis_by_date_wide,Mam_MOSAICA_DNA)

Mam_by_date_MOSAICB<-merge_samples(Mammals.pos.MOSAICB.ps, "Date..YYYYMMDD.",fun = "sum")
Mam_by_date_MOSAICB
Mam_by_date_MOSAICB_otu_tab<-cbind(as.data.frame(otu_table(Mam_by_date_MOSAICB)))
Mam_by_date_MOSAICB_otu_tab1<-rownames_to_column(Mam_by_date_MOSAICB_otu_tab,var = "date")
Mam_by_date_MOSAICB_otu_tab1$occurrence_date<-ymd(Mam_by_date_MOSAICB_otu_tab1$date)

Mam_MOSAICB_DNA<-Mam_by_date_MOSAICB_otu_tab1
Mam_MOSAICB_sp<-colnames(Mam_MOSAICB_DNA)
Mam_MOSAICB_sp<-gsub("1","",Mam_MOSAICB_sp)
colnames(Mam_MOSAIC_DNA)<-Mam_MOSAICB_sp

Mam_MOSAICB_DNA<-Mam_MOSAICB_DNA%>%
  mutate(method="MOSAIC-B")%>%
  select(-date)

visvMOSAICAB_eDNA<-full_join(visvMOSAICA_eDNA,Mam_MOSAICB_DNA)


Mam_by_date_MOSAICC<-merge_samples(Mammals.pos.MOSAICC.ps, "Date..YYYYMMDD.",fun = "sum")
Mam_by_date_MOSAICC
Mam_by_date_MOSAICC_otu_tab<-cbind(as.data.frame(otu_table(Mam_by_date_MOSAICC)))
Mam_by_date_MOSAICC_otu_tab1<-rownames_to_column(Mam_by_date_MOSAICC_otu_tab,var = "date")
Mam_by_date_MOSAICC_otu_tab1$occurrence_date<-ymd(Mam_by_date_MOSAICC_otu_tab1$date)

Mam_MOSAICC_DNA<-Mam_by_date_MOSAICC_otu_tab1
Mam_MOSAICC_sp<-colnames(Mam_MOSAICC_DNA)
Mam_MOSAICC_sp<-gsub("1","",Mam_MOSAICC_sp)
colnames(Mam_MOSAIC_DNA)<-Mam_MOSAICC_sp

Mam_MOSAICC_DNA<-Mam_MOSAICC_DNA%>%
  mutate(method="MOSAIC-C")%>%
  select(-date)

visvMOSAICABC_eDNA<-full_join(visvMOSAICAB_eDNA,Mam_MOSAICC_DNA)

Mam_CTDmethod_DNA<-Mam_CTD_DNA%>%
  rename(method=survey)

visv4methods_eDNA<-full_join(visvMOSAICABC_eDNA,Mam_CTDmethod_DNA)

write.csv(visv4methods_eDNA,"outputs/visv4methods_eDNA.csv")

                  
