# GWAS_pipeline
## First use (NOTE: requires Linux (Ubuntu) partition on your computer)
### Downloading the pipeline
#### If synching with GitHub (recommended, but requires you know how to use GitHub)
* Clone https://github.com/ericgoolsby/GWAS_pipeline.git from your Desktop GitHub client (or from command line).
#### If just downloading the repository
* Download the zipped repository here (https://github.com/ericgoolsby/GWAS_pipeline/archive/master.zip). 
### Downloading required software
* Download required software here (https://drive.google.com/open?id=1tQ7ve20fS53lFqpet_-mCxvGBsRZULrt) and unzip it to the Software/ directory. IMPORTANT: MAKE SURE ALL FILES (e.g. XRQv1_412_261_filtered.fam) ARE IN Software/ DIRECTORY (i.e., not in a subfolder).

## Mason/Goolsby Lab modified steps (assuming sunflower GWAS performed in a SINGLE environment)
### Initial data file prep:
1. Put a .csv data file in the data/ directory. The first row should have column names, and each subsequent row corresponds to a SAM line.
2. The first column should have the SAM line number as either integers or in EXACTLY the following format: e.g. SAM002, SAM073, SAM241, etc.
3. Remove all special characters in column (trait) names (but don't worry about underscores -- they will be automatically removed). The first character of a trait name should be a letter, the rest of the trait name should be letters and numbers ONLY.
4. You should leave missing data blank (don't put NA).
5. Make sure your data are numeric (common errors that break the pipelnie: a space instead of a blank cell; putting N/A instead of a blank cell; a comma instead of a decimal place, etc).
### Running the pipeline:
6. Set the working directory to the root (top-level) directory of the pipeline. OR open the R project (Sunflower-GWAS-2.0.Rproj).
7. Open the R script "Scripts/0 - running.R".
8. If necessary, go to the bottom of script 0 and install required packages (APPENDIX 1).
9. Edit the required information in SECTION 1, then run. (note: I have to set pvalue_cutoff <- 2 for the pipline to work, as of 4/9/2020).
10. Run (but do not edit) SECTION 2.
11. Run the scripts in SECTION 3 one-by-one, checking for warnings and errors along the way. Troubleshoot as needed.
12. Results will be found in Tables/ and Plots/.
13. IMPORTANT: The pipeline assumes you grew the plants up in "Wet" and "Dry" environments, and it also analyzes the logdiff of the two. The only valid results are labeled Wet. THE RESULTS FOR DRY AND LOGDIFF ARE MEANINGLESS (BASED ON SIMULATD PLACEHOLDER DATA).


# Sunflower-GWAS-2.0

This is a R only rewrite and extention of the sunflower GWAS pipeline initiated by Rishi Masalia. (Masalia et al., Plos 2018)

This Sunflower-GWAS-2.0 pipeline adds:
- Streamlined folder structure of inputs/outputs/code/software
- GWAS calcuation using GEMMA
- Colocalization visualization using haplotype blocks
- Gene list per haplotype block
- Drawing of manhattan plots with haplotype blocks overlay

UPDATE 2019/05/10: Works like a charm. Needs software data from separate location though

![](Overview.jpg)


Ongoing improvements:

- List R libraries needed
- Make step 0 script that had bits of code to make kinship file, blocks map, etc
- Write guide
- Flag genes that have significant SNPs for traits in gene list output
- PVE per region following Masalia et al 2018
- Heritability per region
- Epistatis graph (R2 between regions on different chromosomes (within cromosomes already captured in LD plot))


Partial guide language
- The blocks are generated from the 1.5 million SNPs using plink. This divides up the genome in chunks that are co-inherited (according to some set thresholds). We now have condensed the genome from 1.5 million independent SNPs to ~20K independent chunks/regions/blocks. Now from the GEMMA output we get a list of SNPs and their p-values. Since not all SNPs are independent (as evidenced from the blocking procedure) the cutoff for significance is based on the 20K independent blocks we have.
- Another place where the block map we've constructed really shines is finding traits that colocate to the same region. If two or more traits have SNPs that lie in the same block we can say they colocate (possibly pleiotropic but most likely close linkage). So in determining if traits have significant SNPs in the same blocks (even if the SNPs are not the same) all we have to do is figure out if they have common blocks.
- A problem with the blocking procedure is that it's sensitive to missplaced SNPs. The algorithm grows a region by walking along the chromosome and stops if a certain fraction of SNPs don't fit the region it is growing (based on D'). What could happen is that two blocks that are adjacent in reality are broken up by a miss placed cluster of SNPs. If there are traits that hit to the first part of the "true" block and traits that hit to the second part of the "true" block we would be wrong in saying these traits are independent.
- To partially protect against this we can look at the SNPs that are significant for a trait at least once, line up those SNPs and re-do the blocking. If the SNPs from two or more blocks now fully go together in a single block we know that the intervening snps could have been missplaced and for the purposes of determining colocalization we should lump the blocks together.
In the triangle plots I'm showing the LD of the SNPs. The blocks the SNPs belong to in the "genome" row, and the new, significant only SNPs, blocks in the "significant" row. 
