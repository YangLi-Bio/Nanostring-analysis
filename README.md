# Nanostring-analysis

Challenges:
1. The popular analysis tool, ROSALIND, cannot perform differentially expressed genes (DEGs) using the free version
2. Nanostring nsolver can identify DEGs via RLF files, which may not exist in most cases

Pipeline:
1. Convert expression matrix file into RCC files using customized R codes
2. Quanlity control (QC) and normalization using ROSALIND or Nanostring nsolver
