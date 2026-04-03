---
# MIOP terms
methodology_category: sample extraction and purification
project: CalCOFI MOSAIC eDNA Sampling # can be general
purpose: biodiversity assessment objective [OBI:0001969]
analyses: DNA extraction [OBI:0000257]
geographic_location:  North Pacific Ocean [GAZ:00002418]
broad_scale_environmental_context: marine biome [ENVO:00000447], marine photic zone [ENVO:00000209]
local_environmental_context: marine pelagic zone [ENVO:00000208], marine pelagic biome [ENVO:01000023]
target: environmental DNA [NCIT:C169106]
creator: # names separated by commas
materials_required: vortexer [OBI:0400118], centrifuge [OBI:0400106], micropipette [OBI:0001128]
skills_required: sterile technique, pipetting skills, standard molecular technique
time_required: 345
personnel_required: 1
language: en
issued: # YYYY-MM-DD
audience: scientists
publisher: California Cooperative Oceanic Fisheries Investigations # institution
hasVersion: 1
license: CC0 1.0 Universal
maturity level: mature

# FAIRe terms
samp_vol_we_dna_ext: not provided
samp_vol_we_dna_ext_unit: mL
nucl_acid_ext_lysis: enzymatic
nucl_acid_ext_sep: centrifugation
nucl_acid_ext: # https://files.zymoresearch.com/protocols/_d4302_d4306_d4308_zymobiomics-96_magbead_dna_kit.pdf
nucl_acid_ext_kit: #
nucl_acid_ext_modify: # 
dna_cleanup_0_1: 0.0
dna_cleanup_method: not applicable
concentration: not applicable
concentration_method: not applicable
ratioOfAbsorbance260_280: not applicable
pool_dna_num: not applicable
nucl_acid_ext_method_additional: not applicable
---

# Protocol Template -- DNA Extraction

## PROTOCOL INFORMATION

### Minimum Information about an Omics Protocol (MIOP)

- MIOP terms are listed in the YAML frontmatter of this page.
- See [MIOP_definition.md](https://github.com/BeBOP-OBON/0_protocol_collection_template/blob/main/MIOP_definition.md) for list and definitions.

### Making eDNA FAIR (FAIRe)

- FAIRe terms are listed in the YAML frontmatter of this page.
- See <https://fair-edna.github.io/download.html> for the FAIRe checklist and more information.
- See <https://fair-edna.github.io/guidelines.html#missing-values> for guidelines on missing values that can be used for missing FAIRe or MIOP terms.

### Authors

- All authors known to have contributed to the preparation of this protocol, including those who filled in the template.
- Visit https://orcid.org/ to register for an ORCID.
- Date is the date the author first worked on the protocol.

| PREPARED BY  | AFFILIATION  | ORCID        | DATE       |
| ------------ | ------------ | ------------ | ---------- |
| Content Cell | Content Cell | Content Cell | yyyy-mm-dd |
| Content Cell | Content Cell | Content Cell | yyyy-mm-dd |

### Related Protocols

| PROTOCOL NAME | LINK         | VERSION      | RELEASE DATE | INTERNAL/EXTERNAL |
| ------------- | ------------ | ------------ | ------------ | ----------------- |
| MOSAIC pre-2024 PCR protocol  | Content Cell | Content Cell | yyyy-mm-dd   | Internal      |

### Protocol Revision Record

- Version numbers start at 1.0.0 when the protocol is first completed and will increase when changes that impact the outcome of the procedure are made (patches: 1.0.1; minor changes: 1.1.0; major changes: 2.0.0).
- Release date is the date when a given protocol version was finalised.
- Description of revisions includes a brief description of what was changed relative to the previous version.

| VERSION | RELEASE DATE | DESCRIPTION OF REVISIONS |
| ------------- | ------------- | ------------- |
| 1.0.0 | yyyy-mm-dd | Initial release |

### Acronyms and Abbreviations

| ACRONYM / ABBREVIATION | DEFINITION |
| ------------- | ------------- |
| eDNA  | environmental DNA  |
| CalCOFI | California Cooperative Oceanic Fisheries Investigations |

## BACKGROUND

This document describes the required protocol to conduct DNA extraction on filters used to collect environmental DNA. It has been used for CalCOFI eDNA samples. 

### Summary

This protocol 

### Method Description and Rationale

Insert a short description of the functioning principal of the methodology used in the protocol (i.e. how does the method work?). Please note that this is different from the step-by-step description of the protocol procedure.
Insert a short statement explaining why the specific methodology used in the protocol has been selected (e.g. it is highly reproducible, highly accurate, procedures are easy to execute etc….).

### Spatial Coverage and Environment(s) of Relevance

This protocol is applicable to any environmental DNA sample stored on a flat membrane filter.

## PERSONNEL REQUIRED

1 technician with basic molecular biology lab training. 

### Safety

Use a fume hood when directed and also wear proper PPE. Phenol:Chloroform:Isoamyl alcohol (25:24:1) and chloroform:Isoamyl alcohol (24:1) are hazardous chemicals and should not be inhaled or touched.

### Training Requirements

Basic molecular biology training is sufficient for this protocol including sterile technique, pipetting small volumes and using spin columns.

### Time Needed to Execute the Procedure

Day 1 of the procedure will take approximately 1 hour 45 minutes, and Day 2 will take approximately 4 hours. On Day 1, perform  DNA lysis no later than 3pm if you want to finish by 5pm.

## EQUIPMENT


| DESCRIPTION | PRODUCT NAME AND MODEL | MANUFACTURER | QUANTITY | REMARK |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| **Durable equipment** |
| Centrifuge | Centrifuge | Any brand | 1 | Must fit the standard-size 2 mL spin columns used in the protocol. |
| Mixer/incubator | Hula mixer/incubator | Generic brand | 1 | |
| Vortexer | Vortexer | Generic brand | 1 | With adapter for  5 mL tubes |
| Freezer | -80 ºC Freezer | Generic brand | 1 | |
| Freezer | -20 ºC Freezer | Generic brand | 1 | |
| Tweezers | Metal tweezers | Generic brand | 1 | |
| Scissors | Small metal scissors | Generic lab brand | 1 | |
| Micropipette | 1000 uL | Generic brand | 1 | |
| Micropipette | 100 uL | Generic brand | 1 | |
| **Consumable equipment** |
| 1000 uL pipette tips | generic | any manufacturer | 10 | 10 per sample |
| 100 uL pipette tips | generic | any manufacturer | 3 | 3 per sample |
| Petri dish | generic | any manufacturer | 1 | 1 per sample |
| 1.5 mL tubes | generic | any manufacturer | 10 | 10 per sample |
| 2 mL tubes | generic | any manufacturer | 2 | 2 per sample |
| Parafilm | PARAFILM M | Sigma Aldritch | 1 roll | Generic substitute acceptable |
| **Chemicals** |
| RNaseZap | RNaseZap | Generic brand | 1 bottle |  |
| ATL buffer | Buffer ATL | Qiagen | 1 bottle |  |
| TE buffer | 1XTE buffer | Generic brand | 1 bottle |  |
| Proteinase K | 20 mg mL-2 Proteinase K | Generic brand | 1 bottle | 20 mg mL-1 |
| SDS | 10% SDS | Generic brand | 1 bottle |  |
| Sodium Acetate | 3.0 M pH 5.2 Sodium Acetate | Generic brand | 1 bottle |  |
| Isopropanol | Ice-cold 2-propanol | Generic brand | 1 bottle |  |
| Ethanol | Ice-cold 70% ethanol | Generic brand | 1 bottle |  |
| Ethanol | 100% ethanol | Generic brand | 1 bottle |  |
| Phenol:Chloroform:Isoamyl alcohol | Phenol:Chloroform:Isoamyl alcohol (25:24:1) | Generic brand | 500 uL | 500 uL per sample |
| Chloroform:Isoamyl alcohol | Chloroform:Isoamyl alcohol (24:1) | Generic brand | 500 uL | 500 uL per sample |

## STANDARD OPERATING PROCEDURE

### Preparation

1. Set incubation oven at 56°C
2. Label new 1.5 mL tubes 
3. Disinfect with DNAzap-clean with ethanol (do not clean the PCR cabinet walls with ethanol!) and leave material 15 min under UV (Sterilize tweezers, blade, add all under UV with Petri Dish)
4. Remove samples from -80°C freezer and keep on ice 

### Extraction

#### Day 1 - DNA lysis no later than 3pm if you want to finish by 5pm. This day's work takes about 1 hour 45 minutes for 12 samples.

1. Using clean tweezers and blade, cut ½ filters in small pieces into a petri dish (add remaining filter to original tube on ice immediately and back into -80°C) - tweezers and blades are disinfected with RNA zap and ethanol between each samples
2. . Add 525 µL of ATL Buffer 
7. Add 10 µL of proteinase K and 60 µL of SDS
   a. for 12 samples (X12.1) make a MM of 121uL of Pro. K and 726uL of SDS
   b. for 6 samples (X6.1) make a MM of 61uL of Pro. K and 366uL of SDS
8. Pipette mix 
9. Vortex
10. Parafilm tubes to avoid leakage. 
11. Incubate at 56°C overnight 
12. Label 5 sets of 1.5 mL tubes per samples and 1 set of 2mL tubes (the combination of pooled sample + NaAcetate + isopropanol is 1.5-1.8mL some times)

#### Day 2 - DNA extraction/separation and purification the next morning, taking approximately 4 hours.

1. Rinse filters with 500 µL 1X TE buffer
2. Pipette mix and vortex 
3. Divide extracts equally between 2x1.5mL tubes (~500 µL each)
4. Move the samples and materials under a fume hood
5. Extract DNA with equal volume of P:C:I (~500 µL each)
6. Vortex
7. Centrifuge for 5min at max speed (20°C - 13x1000)
8. Carefully remove tubes from centrifuge
9. Carefully transfer the upper phase into a new tube (use pipette at about 300µL do not pipette the interface- change gloves if any liquid drop was spilled)
10. Add equal volume of C:I (~500 µL each)
11. Vortex
12. Centrifuge for 5min at max speed
13. Carefully remove tubes from centrifuge
14. Carefully transfer the upper phase into a new tube – the 2 upper phase from the same sample are pooled into one tube
15. Add 1/10 vol. of NaAcetate and precipitate with 0.6 vol. of ice-cold isopropanol (in the freezer door, refill if volume gets low for next extraction)
   a. if around 1mL add 100uL NaAcetate and 600uL ice-cold isopropanol
   b. for 800uL add 80uL NaAcetate and 480uL ice-cold isopropanol
   c. for 500uL add 50uL NaAcetate and 300uL ice-cold isopropanol
16. Vortex
17. Incubate at -20°C for at least 1h (DNA precipitation) (If working alone with 12 samples and starting at 9AM, you will get to this step around 11:30a-12p - perfect for lunch!)
18. Set centrifuge at 4°C
19. Centrifuge 20 min at 15 000 rpm at 4°C
20. Remove supernatant (slowly invert the tube in waste bottle, do not invert twice)
21. Wash pellet with 0.5 mL of 70% ice-cold EtOH (in the freezer door, refill if volume gets low for next extraction)
22. Vortex
23. Centrifuge 20 min at 15 000 rpm at 4°C
24. Remove supernatant using pipette tip
25. Dry pellet under the fume hood for 20 min 
26. Dissolve sample into 40 µL of 1XTE buffer 
27. Vortex gently
28. Incubate at 4°C for a few hours (or overnight) and store at -20°C
29. Quantify DNA with Qubit Broad Range kit
30. Run PCR (see related protocols)

#### Clean-up (Some steps lab-specific)
1. Empty all C/I and P/C/I tubes into a waste bottle. 
2. If solid waste bins are full put into ziplock bag
3. If waste bottle is full, ask for new bottle
4. Turn off centrifuge. 
5. Gloves go into the white lab waste bins.
