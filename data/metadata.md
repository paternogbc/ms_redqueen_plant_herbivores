# Folder `data/`

**This folder contains 7 files**

| N | file     |      description      |
|:----|----------|:-------------:|
| 1 | full_data.csv | A complete comparatived dataset used in the main analysis |
| 2 | data_exclude_repro.csv  | A comparative dataset generated after excluding records on reproductive organs |
| 3 | data_full_records.csv  | A comparative dataset generated including non-insects |
| 4 | diver_processed_1spe.csv  | Diversity of insects feeding on a single plant species |
| 5 | diver_processed_1gen.csv  | Diversity of insects feeding on a single plant genus |
| 6 | diver_processed_1fam.csv  | Diversity of insects feeding on a single plant family |
| 7 | diver_processed_mfam.csv  | Diversity of insects feeding on more than one plant family |
| 8 | diver_processed_exclude_pollinators.csv  | Diversity of insects after excluding potential pollinators |
| 9 | data_maleness_sd.csv  | Intraspecific variance (SD) for flower maleness  |
| 10 | data_outcrossing_rate.csv  | A comparative dataset with outcrossing rate and flower maleness  |
| 11 | data_lifespan.xlsx  | A comparative dataset of plant life span |
| 12 | phylogenetyic_tree_s1.tre | The study phylogeny generated with scenario S1 |
| 13 | phylogentic_tree_s3_tre | The study phylogeny generated with scenario S3 |
| 14 | phylogenetyic_trees_300_s2.tre | A multiphylo file with 300 trees generated with scenario S2 |
| 15 | study_species_list_taxonomy.csv | The list of of species  and thier taxonomic classification |
| 16 | metadata.md | Metadata file documenting all datasets |


## 1. full_data.csv 

Variables description 

| Column    |      description      | unit |
|-----------|:-------------:|-------:|
| order     | The taxonomic order of the species | |
| family    | The taxonomic family of the species | |
| tip_name  | Binomial species name (genus_epithet) matching the phylogenetic tree | |
| maleness  | Flower maleness. Calculated by dividing the dry biomass of androecium by the dry biomass of androeciem + gynoecium ||
| nspe      | The number of insect species associated to each plant species ||
| nfam      | The number of insect families associated to each plant species ||
| ngui      | The number of insect feeding guilds associated to each plant species ||
| shan      | Shannon diversity index of feeding guilds associated to each plant species ||
| height    | The maximum height of each plant species | (m) |
| sla_imp   | The specific leaf area of each plant species | (mm . mg ^-1) |
| nectar    | The categorical level of nectar offer for each plant species |  |
| color     | The category of flower color for each plant species |  |
| shape     | The category of flower shape for each plant species |  |
| polli     | The category of pollinator group for each plant species |  |
| range     | The area of occupancy for each plant species. Calculated as the number of TK 25 grids out of 3000 |  
