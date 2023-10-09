# Multiple-Sequence-Alignment-and-Phylogenetic-Tree
Use R to perform multiple sequence alignment using a fasta file and leverage the resulting alignment to create a phylogenetic tree.

This is conducted in R using the following packages:

- Biostrings
- DECIPHER
- msa
- seqinr
- ggseqlogo
- ggmsa
- ape

This R code is used and meant as an example. Example fasta format input data is provided as whole genome shotgun DNA sequence data from _Sola senegalensis_ in July 2021.

There is also simple Python code to pull NCBI taxonomy tree for this project and save it in png format.

1. Paradis E, Schliep K (2019). “ape 5.0: an environment for modern phylogenetics and
  evolutionary analyses in R.” _Bioinformatics_, *35*, 526-528.
  doi:10.1093/bioinformatics/bty633 <https://doi.org/10.1093/bioinformatics/bty633>.
2. Pages H, Aboyoun P, Gentleman R, DebRoy S (2023). _Biostrings: Efficient
  manipulation of biological strings_. doi:10.18129/B9.bioc.Biostrings
  <https://doi.org/10.18129/B9.bioc.Biostrings>, R package version 2.68.1,
  <https://bioconductor.org/packages/Biostrings>.
3. Wright ES (2016). “Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in
  R.” _The R Journal_, *8*(1), 352-359.
4. U. Bodenhofer, E. Bonatesta, C. Horejs-Kainrath, and S. Hochreiter (2015) msa: an R
  package for multiple sequence alignment. Bioinformatics 31(24):3997-9999. DOI:
  10.1093/bioinformatics/btv176.
5. Charif, D. and Lobry, J.R. (2007). seqinr
6. Wagih O (2017). _ggseqlogo: A 'ggplot2' Extension for Drawing Publication-Ready
  Sequence Logos_. R package version 0.1,
  <https://CRAN.R-project.org/package=ggseqlogo>.
7. Guangchuang Yu. (2022). Data Integration, Manipulation and Visualization of
  Phylogenetic Trees (1st edition). Chapman and Hall/CRC.
