# Direct RNA sequencing of primary human T cells reveals the impact of immortalization on mRNA pseudouridine modifications

Immortalized cell lines are commonly used as proxies for primary cells in human biology research. For example, Jurkat leukemic T cells fundamentally contributed to uncovering T cell signaling, activation, and immune responses. However, the immortalization process can alter key cellular properties, and it is widely believed among researchers that RNA modification machinery and sites of modification could be significantly altered by this process. Focusing on pseudouridine (ψ) modifications, which are among the most abundant mRNA modifications, we study here in detail ψ profiles in mRNA from primary and immortalized T cells using direct RNA sequencing (DRS) with mod-p ID analysis. Surprisingly, we find that 87% of ψ-sites were shared between the two cell types on transcripts encoding proteins performing essential cellular processes that include RNA-modification regulation. Furthermore, analysis of the 13% of sites that are unique to each cell type reveals that Jurkat cells contained transcripts linked to immune activation and oncogenesis, while primary T cells contained transcripts associated with calcium signaling and intracellular trafficking. We provide a list of these genes, which should be considered when using immortalized cells to study RNA modifications in immunology contexts. Most differences were driven by whether the mRNA was present or absent in the immortalized or primary cell type. Interestingly, RNA-modification enzyme expression levels were highly conserved in both cell types, suggesting that site-specific differences in ψ levels arise from regulatory processes rather than differences in modification enzyme levels.

In this GitHub repository, you will find all the code used for the analysis/generation of figures. 
You can find all the files required to run the R script on [Dropbox](https://www.dropbox.com/scl/fo/CREATE). We recommend downloading the files from Dropbox and substituting the paths to those on your local machine. 
Due to the size of the data we're hosting fastq raw data for the Direct libraries on NIH NCBI SRA under the accession number [PRJNA1136513](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1136513).
We are happy to share the raw fast5 files: contact the Corresponding author listed in the paper for fast5 files, as these are not accepted by SRA anymore due to size.  

## Summary

The Figures.r R script can be executed in RStudio.

The CSV files containing the p-values analysis for psi detection has been generated using [p-Mod ID](https://github.com/RouhanifardLab/PsiNanopore).


