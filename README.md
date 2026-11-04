# LC-FC-Resilience


##########################################################################################################################
### Code for Lawn et al., Locus Coeruleus Connectivity Confers Cognitive Resilience in Preclinical Alzheimer's Disease ###
##########################################################################################################################

## Contents

- `masks` - 2mm LC, shifted LC and VTA masks used in the analyses
- `01_masks_and_atlas` - generating the LC, shifted LC and VTA masks, and network labels for the Schaefer-Tian parcels
- `02_connectivity` - connectivity matrices, ComBat harmonisation, generating FC metrics, weighted degree, cortical thickness, and tSNR
- `03_main_analyses` - mixed-effects models of PACC over time
- `04_sensitivity_specificity` - tissue, shifted LC, and VTA control analyses
- `R` - shared functions for the R scripts 
- `python` - surface plotting

## Data 

The A4 and LEARN data can be requested from the LONI IDA (https://ida.loni.usc.edu). 

## Atlases and templates

- Keren et al. (2015) LC template, 2SD version
- Brainstem Navigator v1.0, VTA-PBP (https://doi.org/10.25790/BML0CM.96)
- Schaefer 2018 parcellation (https://github.com/ThomasYeoLab/CBIG)
- Tian subcortex atlas, S2 (https://github.com/yetianmed/subcortex)

## Citation

Lawn T, Diez I, Bueichekú E, Coughlan G, Rentz DM, Johnson KA, Sperling RA, Sepulcre J, Jacobs HIL. Locus Coeruleus Connectivity Confers Cognitive Resilience in Preclinical Alzheimer's Disease.
