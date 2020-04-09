# edit these lines as needed
### BEGIN SECTION 1 ###
trait_filename <- "katie_data.csv"
pvalue_cutoff <- 1 # 1 = Bonferroni = 1; 2 = "suggested" 0.001 threshold
### END SECTION 1 ###

##############################
### RUN BUT DO NOT EDIT ######
##### THIS SECTION ###########
##############################
##### BEGIN SECTION ##########
##############################
source("Scripts/functions.R")
process_data(trait_filename = trait_filename)
set_threshold(method = 2)
##############################
####### END SECTION ##########
##############################


##############################
### UNCOMMENT AND RUN ########
##### EACH SCRIPT  ###########
#### ONE AT A TIME ###########
##############################
#source("Scripts/1 - Phenotype to GEMMA.R")
#source("Scripts/2 - Make manhattan plots.R")
#source("Scripts/2b (optional) - Make single environment manhattan plots.R")
#source("Scripts/3 - Blocks and heatmaps.R")
#source("Scripts/4 - Trait GWAS colocalization figure.R")
#source("Scripts/4c - (optional) Per chromosome - Trait GWAS colocalization figure.R")
#source("Scripts/5 - List genes in regions.R")
#source("Scripts/6 - Manhattan with blocks.R")
#source("Scripts/7 - Export PVE per trait to table.R")
#source("Scripts/8 - Show blocks on haplotype map.R")


### BEGIN APPENDIX 1 ###
if(FALSE) # change to TRUE to install all required packages
{
  install.packages(c("data.table", "qqman", "tidyverse", "RColorBrewer", "ggpubr", 
                     "grid", "ggrepel", "gridExtra", "cowplot", "wesanderson", "corrr", 
                     "Hmisc", "ggdendro", "urltools", "scales"))
}
### END APPENDIX 1 ###