---
# MIOP terms
methodology_category: Omics Analysis
project: MarVer1 PCR on CalCOFI MOSAIC Samples
purpose: PCR [OBI:0000415]
analyses: PCR [OBI:0000415]
geographic_location: [GAZ:00007899]
broad_scale_environmental_context: marine biome [ENVO:00000447], marine photic zone [ENVO:00000209]
local_environmental_context: marine pelagic zone [ENVO:00000208], marine pelagic biome [ENVO:01000023], marine biome [ENVO:00000447]
environmental_medium: sea water [ENVO:00002149], polymerase chain reaction [OBI:0000415]
target: 12S [NCIT:C128263]
creator: Nastassia Patin, Eldridge Wisely, Ella Crotty # names separated by commas
materials_required: PCR instrument [OBI:0000989], agarose gel electrophoresis system [OBI:0001134]
skills_required: sterile technique, pipetting skills, standard molecular technique
time_required: 180 # minutes (integer)
personnel_required: 1
language: en
issued: 2025-01-21
audience: scientists
publisher: California Cooperative Oceanic Fisheries Investigations
hasVersion: 1
license: CC0 1.0 Universal
maturity level: mature

# FAIRe terms
pcr_0_1: 1 # 1
thermocycler: Applied Biosystems SimpliAmp # name of thermocycler
amplificationReactionVolume: 20
assay_name: MarVer1 
assay_validation: not provided
targetTaxonomicAssay: 12S rRNA gene targeting vertebrates 
targetTaxonomicScope: Vertebrates
target_gene: 12S rRNA 
target_subfragment: not provided
ampliconSize: 197-212 
pcr_primer_forward: CGTGCCAGCCACCGCG 
pcr_primer_reverse: GGGTATCTAATCCYAGTTTG 
pcr_primer_name_forward: MarVer1F
pcr_primer_name_reverse: MarVer1R
pcr_primer_reference_forward: 10.1002/edn3.72
pcr_primer_reference_reverse: 10.1002/edn3.72 
pcr_primer_vol_forward: 1.0
pcr_primer_vol_reverse: 1.0
pcr_primer_conc_forward: 10
pcr_primer_conc_reverse: 10
probeReporter: not applicable
probeQuencher: not applicable
probe_seq: not applicable
probe_ref: not applicable
probe_conc: not applicable
commercial_mm: 2X Phusion Master Mix
custom_mm: not applicable
pcr_dna_vol: 2.0
pcr_rep: 4 # per Eldridge
nucl_acid_amp: # holding for protocols.io link
pcr_cond: initial denaturation:98_0.5;denaturation:98_0.17;annealing:60_0.5;elongation:72_0.75;final elongation:72_10;35
annealingTemp: 60 # 50
pcr_cycles: 35 # 30
pcr_analysis_software: not provided
pcr_method_additional: not provided
---

# Protocol Template -- PCR

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

| PREPARED BY | AFFILIATION | ORCID | DATE |
| ------------ | ------------ | ------------ | ---------- |
| Nastassia Patin | Scripps Institution of Oceanography, University of California San Diego | 0000-0001-8522-7682 | 2026-06-22 |
| Eldridge Wisely | Scripps Institution of Oceanography, University of California San Diego | 0009-0001-2784-4778 | 2026-06-22 |
| Ella Crotty | Scripps Institution of Oceanography, University of California San Diego | 0009-0004-4340-2007 | 2026-06-22 |
### Related Protocols

| PROTOCOL NAME | LINK         | VERSION      | RELEASE DATE | INTERNAL/EXTERNAL |
| ------------- | ------------ | ------------ | ------------ | ----------------- |
| NOAA/AOML PCR Protocol 12S rRNA V5-V6 (MiFish)  | [NOAA/AOML PICR Protocol 12S](https://github.com/aomlomics/protocols/blob/main/markdown/protocol_pcr_ssu12sv5v6_mifish.md) | V5-V6 | 2024-11-16 |  |

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
| SOP  | Standard Operating Procedure  |
| eDNA  | Environmental DNA  |
| PCR  | Polymerase Chain Reaction  |

## BACKGROUND

This protocol describes steps for performing PCR for the 12S MarVer1 marker gene region using eDNA extracted from seawater.

### Summary

The purpose of this protocol is to conduct PCR to optimize marine vertebrate detection using the MarVer1 primer set.

### Method Description and Rationale

This protocol is used for PCR amplifying the 12S MarVer1 marker gene regions of environmental DNA. It is highly reproducible and can easily be adapted for any number of samples (i.e. a full 96-well plate or a few samples).

### Spatial Coverage and Environment(s) of Relevance

This protocol can be used to amplify the 12S marker gene region of any eDNA sample.

## PERSONNEL REQUIRED

One person with molecular laboratory and biology experience.

### Safety

There are no hazardous chemicals or materials involved in this protocol. Standard lab safety techniques should still be used such as wearing PPE to avoid skin or eye contact.

### Training Requirements

Basic molecular biology training is sufficient for this protocol including sterile technique, pipetting small volumes and programming/running thermal cyclers.

### Time Needed to Execute the Procedure

Protocol takes about 3 hours including set-up and thermal cycler run time.

## EQUIPMENT

| DESCRIPTION | PRODUCT NAME AND MODEL | MANUFACTURER | QUANTITY | REMARK |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| **Durable equipment** |
| 100-1000 ul Pipette  | Eppendorf Research Plus Adjustable-Volume Pipette  | Eppendorf  | 1  | Can be substituted with any accurate pipette  |
| 10-100 ul Pipette  | Eppendorf Research Plus Adjustable-Volume Pipette  | Eppendorf  | 1  | Can be substituted with any accurate pipette  |
| 1-10 ul Pipette  | Eppendorf Research Plus Adjustable-Volume Pipette  | Eppendorf  | 1  | Can be substituted with any accurate pipette  |
| Thermal Cycler  | Applied Biosystems SimpliAmp  | Thermo Fisher Scientific  | 1  | Can be substituted for generic  |
| Microwave  | Generic Microwave  | Generic Brand  | 1  |   |
| Flask  | 500 mL Flask  | Generic Brand  | 1  | Used for mixing agarose gel solution  |
| 1 L Glass Container  | 1 L Glass Container  | Generic Brand  | 1  | Used for storing 1x TBE buffer  |
| Gel Tray & Box  | Gel Electrophoresis Box and Tray  | Generic Brand  | 1  |   |
| Gel Combs  | Gel Electrophoresis Combs  | Generic Brand  | 2  | Use the correct number of wells for your samples  |
| **Consumable equipment** |
| Gloves | Nitrile Gloves, Exam Grade, Powder-free  | ULINE  | 1 box  | Can be substituted with generic  |
| Kim Wipes | KimWipe Delicate Task Wipers  | KimTech  | 1 box  | Can be substituted with generic  |
| 8-Count Strip Tubes and Caps | Microamp 8-Tube Strips  | Thermo Fisher Scientific  |   | As many as samples  |
| 1000 uL Filter Tips | TipOne RPT Filter Tips, 1000 uL  | US Scientific  | 1 box  | Can be substituted with generic  |
| 200 uL Filter Tips | 200 uL Filter Tips  | US Scientific  | 2 boxes  | Can be substituted with generic  |
| 10 uL Filter Tips | TipOne RPT Filter Tips, 10 uL  | US Scientific  | 2 boxes  | Can be substituted with generic  |
| Milli-Q quality water | Milli-Q Type 1 Ultrapure Water Systems  | EMD Millipore  | 4 uL per sample well  |   |
| 2X Phusion Mastermix | Phusion® High-Fidelity PCR Master Mix with HF Buffer, M0531S | New England Biolabs  | 10uL per sample well  |   |
| Forward Primer | MarVer1F, CGTGCCAGCCACCGCG  | Any primer manufacturer  |   | Primer must be diluted from 100uM stocks to 10uM |
| Reverse Primer  |  MarVer1R, GGGTATCTAATCCYAGTTTG  | Any primer manufacturer  |   | Primer must be diluted from 100uM stocks to 10uM  |
| TBE Buffer (10x)  | TBE Buffer 10X Solution, Molecular Biology Grade, UltraPure  | Thermo Scientific  | 100 uL  |   |
| Agarose  | Agarose LE, Molecular Biology Grade, UltraPure  | Thermo Scientific  | 4g  | Primer must be diluted from 100uM stocks to 10uM  |
| SYBR Safe  | SYBR Safe DNA Gel Stain  | Invitrogen  | 20uL  | Light sensitive - do not expose to light  |
| Gel stain loading dye  | DNA Gel Loading Dye (6x)  | Thermo Scientific  | 480 ul (per plate)  |   |
| 100bp DNA Ladder  | Generuler 100 bp DNA Ladder  | Thermo Scientific  | 6 ul per lane on gel  | Primer must be diluted from 100uM stocks to 10uM  |
| Parafilm  | Parafilm M Lab Film  | Generic  | 1 roll  | Can substitute with generic brand  |
| **Chemicals** |
| EtOH  | Ethanol  | Generic Brand  | 1 wash bottle  | Must be molecular grade ethanol  |
| **Qubit**  | 
| Qubit Reagents  | Qubit dsDNA Quantification Assay Kit  | Invitrogen  | 1 kit  |   |
| Clear Qubit Assay tubes  | 0.5 mL thin-walled, polypropylene tubes  | Invitrogen  | 98  | Must be correct tubes to allow for fluorometer to read concentration  |

## STANDARD OPERATING PROCEDURE

In the following SOP, please use the exact names of equipment as noted in the table above.

Provide a step-by-step description of the protocol. The identification of difficult steps in the protocol and the provision of recommendations for the execution of those steps are encouraged.

### Preparation

1. Dilute primers from 100 uM to 10 uM if not already at 10uM.
2. Set up PCR under hood by wiping off all surfaces, pipettes, and racks with Bleach and Ethanol and UV sterilization for 5-10 mins.
3. Create the PCR mastermix consisting of the Taq, both forward and reverse primers, MilliQ water, and the eDNA sample. Add 18 uL to each strip tube well.

### PCR

**Primers**: PCR primer sequences

| PCR Primer Name | Direction | Sequence (5’ -> 3’)|
| ----- | ----- | ----- |
| MarVer1F | forward | CGTGCCAGCCACCGCG |
| MarVer1R | reverse | GGGTATCTAATCCYAGTTTG |

**Reaction Mixture**: PCR reagents, volumes, initial and final concentrations

| Reagent | Volume | Initial Concentration | final concentration|
| ----- | ----- | ----- | ----- |
| 2X Phusion Master Mix | 10uL | 2 | 1 |
| 10 uM Forward Primer | 1uL | 10 | 0.5 |
| 10 uM Reverse Primer | 1uL | 10 | 0.5 |
| Nuclease-Free Water | 5uL |  | |
| rAlbumin | 1uL | 10 | 0.5 |
| Template DNA | 2uL |  | |

**PCR Cycling Program**: 

| PCR Step | Temperature | Duration | Repetition |
| ----- | ----- | ----- | ----- |
| Initial Denaturation | 98ºC | 30s | 1x |
| Denaturation | 98ºC | 10s | 35x |
| Annealing | 60ºC | 30s | 35x |
| Extension | 72ºC | 45s | 35x |
| Final Extension | 72ºC | 10min | 1x |
| Hold | 4ºC | Infinity |  |


### Quality Control and PCR cleanup

2% Agarose Gel Following PCR amplification

1. Make stock solution of TBE buffer (1x) in a 1-L glass container by adding 100 ml of stock TBE buffer (10x) to 900 ml DI water.
2. Depending on size of gel tray, mix a proportion of 100 ml of TBE buffer (1x) to 2 g of agarose in a flask. Use scale to weigh agarose.
3. Microwave mixture for 1 minute, followed by 15-30 second intervals. Watch carefully after 1 minute - mixture can bubble out of flask! The agarose should be fully dissolved so that the solution is mostly clear.
4. Wear a protective glove when handling flask as the mixture will be hot. Allow for gel mixture to cool in flask for 5-10 min.
5. While cooling, set up gel tray in gel box. Make sure the tray is oriented properly and sealed tight for gel pouring. Add two gel combs for 20 wells each lane - total of 40 wells (or however many is needed for number of samples).
6. Once cooled to 50°C, add 10 ul of SYBR safe to the mixture and swirl the flask gently to mix (be careful not to create bubbles). SYBR safe is light sensitive so after adding to mixture, immediately close the vial and store in the dark.
7. Pour the gel mixture slowly into the tray and remove any bubbles using a pipette tip.
8. Allow gel to set for 30-45 min.
9. Cut large strips of parafilm or use 8-strip tubes and label each sample as a position on the parafilm/tubes.
10. Pipette 1 ul of blue loading dye onto each sample position or into each tube.
11. Pipette 5 ul of DNA sample into loading dye and pipette to mix 2-3 times.
12. Once the gel is set, load it with the samples (6uL of each post-PCR sample). Also add 1uL of ladder combined with 1uL loading dye and 4uL MilliQ water.
13. Fill the gel box with enough TBE buffer (1x) to fully submerge the gel beneath ~1cm of buffer.
14. Run gel at ~100 V for 40-50 min then visualize on gel reader machine.
15. 
(Optional) To quantify the DNA in the samples post-PCR, the Qubit can be utilized to determine the DNA concentrations. Refer to the Manufacturer protocol for High Senstivity: [Qubit Protocol](https://tools.thermofisher.com/content/sfs/manuals/Qubit_dsDNA_HS_Assay_UG.pdf)

#### Positive Control

No positive controls were used.

#### Negative Control

A negative control should be performed using the same volume of water as eDNA template in the samples (so 2uL of H2O as the negative control). 

### Basic Troubleshooting Guide

- If contamination bands are present on the gel, ensure the lab bench is sterilized and re-run the PCR with new aliquots of reagents and MilliQ water.

## REFERENCES

- [Qubit protocol, Qubit® dsDNA HS Assay Kits For use with the Qubit® Fluorometer (all models) Catalog nos. Q32851, Q32854](https://tools.thermofisher.com/content/sfs/manuals/Qubit_dsDNA_HS_Assay_UG.pdf)

## APPENDIX A: DATASHEETS

## 
