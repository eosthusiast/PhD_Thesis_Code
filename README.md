# Evaluating the importance of intraspecific variation in fungal phenotypes and genotypes, using Aotearoa New Zealand _Pleurotus_ as an exemplar
This reposity holds various code snippets and scripts used for data analysis and bioinformatics in my 2024 PhD thesis. Within each directory is a brief overview of chapter contents and details about each of the tools used. 

For a brief outline of chapter contents:
## Chapter 2: Multigene phylogeny of _Pleurotus_ in Aotearoa New Zealand supports an indigenous clade of _P. pulmonarius_
Here I've included scripts used for _in silico_ extraction of ITS, LSU, RPB1, RPB2 and Tef gene regions from whole genome assemblies of NZ _Pleurotus_ specimens sequenced in this study, as well as international _Pleurotus_ genomes from NCBI. First, I utilised a BLASTn subject query (Camacho et al., 2009) to extract matches of representative regions of ITS, LSU, RPB1, RPB2, and Tef within each genome assembly using closely related reference sequences from GenBank for each species, where available. Second, I used ThermonucleotideBLAST v2.61 ("tntblast", Gans & Wolinsky, 2008) to extract predicted amplicons from each genome assembly using forward and reverse primers for each target region. Specifically, I used the primer pairs ITS5/ITS4 for ITS (White et al., 1990), LR05/LR5 for LSU (Zervakis et al., 2019), RPB1-Af/RPB1-Cr for RPB1 (Matheny, 2005), bRPB2-3.1F/gRPB2-6R and RPB2-5F/bRPB2-11R1 for RPB2 (Liu et al., 1999), and TEF-983F/TEF-1567R for Tef (Li et al., 2017). The BLASTn approach had a higher success rate but resulted in shorter amplicons than tntblast (Table 2).

## Chapter 3: Phenotypic traits variation is similar within and between five _Pleurotus_ species
Here I've included R scripts used for data handling of phenotypic trait measurements and statistical analyses of the results.

## Chapter 4: Low genome-wide intraspecific variation and high mating type variation reveal local panmixia of _Pleurotus purpureo-olivaceus_ in South Island of Aotearoa New Zealand 
I employed various bioinformatics approaches to assemble and analyse whole genome data of _P. purpureo-olivaceus_.

## Publications
Work in progres...
