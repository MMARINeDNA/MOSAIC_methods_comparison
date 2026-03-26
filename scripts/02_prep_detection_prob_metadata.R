#This script takes an initial phyloseq object and calculates important metadata like distance of transect, total time filtered, and then outputs csvs per PCR replicate, and per sample in two versions: with ASV counts per taxon, and as just detection 1/0.  

library(geosphere)
library(geodist)
library(tidyverse)
library(phyloseq)
library(metagMisc)
library(microViz)
library(lubridate)
library(stringr)


`%ni%` = Negate(`%in%`)

#load initial phyloseq object
MOSAIC.ps<-readRDS("outputs/MOSAIC.ps.RDS")
#MOSAIC.ps<-readRDS("outputs/MetabaR_cleaned_MOSAIC.ps")

#remove contam taxa
contams<-c("Homo","Bos","Sus","Gallus","Canis","Felis","Meleagris","Mus")
blanks<-subset_samples(MOSAIC.ps,is.na(Collection.device))
blanks <- prune_taxa(taxa_sums(blanks) > 0, blanks)

nonblanks<-c("CTD", "Underway system")
carboys<-("MOSAIC-C")

MOSAIC.ps<-subset_taxa(MOSAIC.ps,genus%ni%contams)
MOSAIC.ps<-subset_taxa(MOSAIC.ps,!is.na(phylum))
MOSAIC.ps<-subset_samples(MOSAIC.ps,Collection.device%in%nonblanks)

#Investigate blanks to see if MetabaR is necessary
blank_taxa<-cbind(as.data.frame(tax_table(blanks)))
blank_samples<-cbind(as.data.frame(sample_data(blanks)))
blank_otu_tab<-cbind(as.data.frame(otu_table(blanks)))

saveRDS(blanks, "outputs/MOSAIC_blanks.ps.RDS")

MOSAIC_sample_data<-cbind(data.frame(sample_data(MOSAIC.ps)))
#add Inshore/Offshore column manually, based on station start and end numbers (which were not entered correctly in the field.)

write.csv(MOSAIC_sample_data,"outputs/MOSAIC_sample_data_to_add_inshore-offshore.csv")
#edit in Excel or similar based on station 60 being the cutoff
MOSAIC_sample_data1<-read.csv("data/MOSAIC_sample_data_inshore.offshore_added.csv")

MOSAIC_sample_data<-column_to_rownames(MOSAIC_sample_data1,var = "X")

sample_data(MOSAIC.ps)<-MOSAIC_sample_data



MOSAIC_sample_data$Volume.filtered..mL.<-as.numeric(MOSAIC_sample_data$Volume.filtered..mL.)
MOSAIC_sample_data<-MOSAIC_sample_data%>%
  mutate(depth.rank=as.factor(
    if_else(Depth..m.==">250",
            3,
            if_else(Depth..m.=="<250",
                    2,
                    if_else(Depth..m.=="733",
                            3,
                            if_else(Depth..m.=="380-441",
                                    3,
                                    if_else(Depth..m.=="322-441",
                                            3,
                                            if_else(Depth..m.=="320-441",
                                                    3,
                                                    1))))))))

#Fix missing end-point data for CTD samples

MOSAIC_sample_data<-MOSAIC_sample_data%>%
  mutate(Longitude.end = if_else(Collection.device == "CTD", Longitude.start, Longitude.end))%>%
  mutate(Latitude.end = if_else(Collection.device == "CTD", Latitude.start, Latitude.end))
  



start<-cbind(MOSAIC_sample_data$Longitude.start,MOSAIC_sample_data$Latitude.start)
end<-cbind(MOSAIC_sample_data$Longitude.end,MOSAIC_sample_data$Latitude.end)
colnames (start) <- colnames (end) <- c ("x", "y")
distances_in_m<-geodist(start,end,paired = TRUE,measure = "geodesic")

MOSAIC_sample_data$dist_km<-distances_in_m/1000

MOSAIC_sample_data<-MOSAIC_sample_data%>%
  mutate(dist_km= 
           if_else(Collection.device=="CTD",
                   0,
                   dist_km))

sample_data(MOSAIC.ps)<-MOSAIC_sample_data

#calculate time filtered

MOSAIC_sample_data1<-MOSAIC_sample_data%>%
  mutate(
    # Fix known typo: "6.05" → "6:05"
    Filtration.time.start = gsub("^(\\d+)\\.(\\d+)$", "\\1:\\2", Filtration.time.start),
    
    # Build full datetimes
    datetime_start = parse_date_time(
      paste(Date..YYYYMMDD., Filtration.time.start),
      orders = "Ymd H:M",
      quiet = TRUE
    ),
    datetime_end = parse_date_time(
      paste(Date..YYYYMMDD., Filtration.time.end),
      orders = "Ymd H:M",
      quiet = TRUE
    ),
    
    # If end time is earlier than start time, the filter ran past midnight —
    # add one day to the end datetime
    datetime_end = if_else(
      !is.na(datetime_start) & !is.na(datetime_end) & datetime_end < datetime_start,
      datetime_end + days(1),
      datetime_end
    ),
    
    # Calculate elapsed time in minutes
    filtration_time_min = as.numeric(
      difftime(datetime_end, datetime_start, units = "mins")
    )
  )

sample_data(MOSAIC.ps)<-MOSAIC_sample_data1

saveRDS(MOSAIC.ps,"outputs/MOSAIC_clean.ps.RDS")

#tax_fix_interactive(MOSAIC.ps)
MOSAIC_taxfix.ps<-MOSAIC.ps %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )
MOSAIC_taxamerge<-tax_glom(MOSAIC_taxfix.ps,taxrank = "species")
MOSAIC_taxnamed.ps<-tax_rename(MOSAIC_taxfix.ps,rank = "species")
taxnames<-taxa_names(MOSAIC_taxnamed.ps)
taxnames<-gsub(pattern = " 1", replacement = "",taxnames)
taxa_names(MOSAIC_taxnamed.ps)<-taxnames

#sanity check
MOSAIC_taxnamed_sample_data<-cbind(as.data.frame(sample_data(MOSAIC_taxnamed.ps)))
MOSAIC_taxnamed_otu_table<-cbind(as.data.frame(otu_table(MOSAIC_taxnamed.ps)))
MOSAIC_taxnamed_taxa_table<-cbind(as.data.frame(tax_table(MOSAIC_taxnamed.ps)))

write.csv(MOSAIC_taxnamed_sample_data, "outputs/MOSAIC_taxnamed_sample_data.csv")
write.csv(MOSAIC_taxnamed_otu_table,"outputs/MOSAIC_taxnamed_otu_table.csv")
write.csv(MOSAIC_taxnamed_taxa_table,"outputs/MOSAIC_taxnamed_taxa_table.csv")

saveRDS(MOSAIC_taxnamed.ps, "outputs/MOSAIC_taxnamed.ps.RDS")


#Target taxa for MarVer =Vertebrates
Vertebrates<-subset_taxa(MOSAIC_taxnamed.ps,phylum=="Chordata" & class!="Appendicularia")
Vertebrates <- prune_taxa(taxa_sums(Vertebrates) > 0, Vertebrates)
#don't prune samples, just leave 0s
Vertebrate_MOSAIC_tax_table<-cbind(as.data.frame(tax_table(Vertebrates)))


saveRDS(Vertebrates, "outputs/Vertebrates.ps.RDS")

#Marine mammals
Mammals<-subset_taxa(Vertebrates, class=="Mammalia")
Mammals <- prune_taxa(taxa_sums(Mammals) > 0, Mammals)
Mammals
#don't prune samples, just leave 0s


#make another version including only samples with positive mammal detections
Mammals.pos.ps<-prune_samples(sample_sums(Mammals)>0, Mammals)
Mammals.pos.ps
saveRDS(Mammals.pos.ps,"outputs/Mammals.pos.ps.RDS")

#export to 3 csvs
F23_S24_MOSAIC_Mammal_sample_data_by_rep<-cbind(as.data.frame(sample_data(Mammals)))
F23_S24_MOSAIC_Mammal_otu_table_by_rep<-cbind(as.data.frame(otu_table(Mammals)))
F23_S24_MOSAIC_Mammal_taxa_table_by_rep<-cbind(as.data.frame(tax_table(Mammals)))




#Merge by sample

#make a new column for Cruise_FilterID_biorep
F23_S24_MOSAIC_Mammal_sample_data_by_rep<-F23_S24_MOSAIC_Mammal_sample_data_by_rep%>%
  mutate(ConsolidatedID=paste0(Cruise,"_",FilterID,"_",Bio.rep))

F23_S24_MOSAIC_Mammal_sample_data_by_rep<-rownames_to_column(F23_S24_MOSAIC_Mammal_sample_data_by_rep,var="UniqueID")

#count reps per sample
rep_counts<-F23_S24_MOSAIC_Mammal_sample_data_by_rep%>%
  group_by(ConsolidatedID)%>%
  count(ConsolidatedID)%>%
  ungroup%>%
  rename(n_tech.reps=n)

F23_S24_MOSAIC_Mammal_sample_data_by_rep<-left_join(F23_S24_MOSAIC_Mammal_sample_data_by_rep,rep_counts)
F23_S24_MOSAIC_Mammal_sample_data_by_rep<-column_to_rownames(F23_S24_MOSAIC_Mammal_sample_data_by_rep,var = "UniqueID")
sample_data(Mammals)<-F23_S24_MOSAIC_Mammal_sample_data_by_rep
saveRDS(Mammals, "outputs/Mammals_by_rep.ps.RDS")

library(speedyseq)
Mammals_by_sample.ps <- merge_samples2(Mammals, "ConsolidatedID") #sums otus across reps


F23_S24_MOSAIC_Mammal_sample_data_by_sample<-cbind(as.data.frame(sample_data(Mammals_by_sample.ps)))
F23_S24_MOSAIC_Mammal_otu_table_by_sample<-cbind(as.data.frame(otu_table(Mammals_by_sample.ps)))
F23_S24_MOSAIC_Mammal_taxa_table_by_sample<-cbind(as.data.frame(tax_table(Mammals_by_sample.ps)))

F23_S24_MOSAIC_Mammal_sample_data_by_sample<-F23_S24_MOSAIC_Mammal_sample_data_by_sample%>%
  mutate(uL.of.template.per.PCR=.2,
         uL.of.template.per.sample=.2*n_tech.reps)

sample_data(Mammals_by_sample.ps)<-F23_S24_MOSAIC_Mammal_sample_data_by_sample
saveRDS(Mammals_by_sample.ps, "outputs/Mammals_by_sample.ps.RDS")

#103 MOSAIC-A samples (325 reps)


#20 MOSAIC-B samples (80 reps)


#6 MOSAIC-C samples (24 reps)


#36 CTD samples (70 reps)


#Make a presence/absence version:
Mammals_by_sample.pa.ps <- transform_sample_counts(Mammals_by_sample.ps, function(abund) 1*(abund > 0))
saveRDS(Mammals_by_sample.pa.ps, "outputs/Mammals_by_sample.pa.ps.RDS")
F23_S24_MOSAIC_Mammal_otu_table_by_sample_pa<-cbind(as.data.frame(otu_table(Mammals_by_sample.ps)))



#Export csvs for presence/absence analyses:
write.csv(F23_S24_MOSAIC_Mammal_sample_data_by_rep,"outputs/F23_S24_MOSAIC_Mammal_sample_data_by_rep.csv")

write.csv(F23_S24_MOSAIC_Mammal_otu_table_by_rep,"outputs/F23_S24_MOSAIC_Mammal_otu_table_by_rep.csv")

write.csv(F23_S24_MOSAIC_Mammal_taxa_table_by_rep,"outputs/F23_S24_MOSAIC_Mammal_taxa_table_by_rep.csv")

write.csv(F23_S24_MOSAIC_Mammal_taxa_table_by_sample,"outputs/F23_S24_MOSAIC_Mammal_taxa_table_by_sample.csv")

write.csv(F23_S24_MOSAIC_Mammal_sample_data_by_sample,"outputs/F23_S24_MOSAIC_Mammal_sample_data_by_sample.csv")

write.csv(F23_S24_MOSAIC_Mammal_otu_table_by_sample,"outputs/F23_S24_MOSAIC_Mammal_otu_table_by_sample.csv")

write.csv(F23_S24_MOSAIC_Mammal_otu_table_by_sample_pa,"outputs/F23_S24_MOSAIC_Mammal_otu_table_by_sample_pa.csv")

