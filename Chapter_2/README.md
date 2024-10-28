## Chapter 2: Multigene phylogeny of _Pleurotus_ in Aotearoa New Zealand supports an indigenous clade of _P. pulmonarius_
Here I've included scripts used for _in silico_ extraction of ITS, LSU, RPB1, RPB2 and Tef gene regions from whole genome assemblies of NZ _Pleurotus_ specimens sequenced in this study, as well as international _Pleurotus_ genomes from NCBI. 
First, I utilised a BLASTn subject query (Camacho et al., 2009) to extract matches of representative regions of ITS, LSU, RPB1, RPB2, and Tef within each genome assembly using closely related reference sequences from GenBank for each species, where available. Second, I used ThermonucleotideBLAST v2.61 ("tntblast", Gans & Wolinsky, 2008) to extract predicted amplicons from each genome assembly using forward and reverse primers for each target region. Specifically, I used the primer pairs ITS5/ITS4 for ITS (White et al., 1990), LR05/LR5 for LSU (Zervakis et al., 2019), RPB1-Af/RPB1-Cr for RPB1 (Matheny, 2005), bRPB2-3.1F/gRPB2-6R and RPB2-5F/bRPB2-11R1 for RPB2 (Liu et al., 1999), and TEF-983F/TEF-1567R for Tef (Li et al., 2017). The BLASTn approach had a higher success rate but resulted in shorter amplicons than tntblast (Table 2).

### 1 TNTblast
Run SLURM job that does the following: Runs tntblast over all -d files in the folder Pul, using each of the -i target gene regions within the folder Target_gene_regions, and then extracts the first sequence in the output fasta file (first two lines of the fasta file) to a new fasta file named based on the first two strings before the second _ in the -d filename, with an added _ followed by the first string of the -i target gene. Make sure that not just the fasta filename gets renamed, but also the fasta description (first line after ">"). E.g., the target fasta filename and fasta desription (after ">") for the example code below would be DH1392_Pul_ITS

```
wget https://github.com/jgans/thermonucleotideBLAST/archive/refs/tags/v2.61.tar.gz
sbatch Pul_tntblast.sl
```

### 2 BLASTn subject query
