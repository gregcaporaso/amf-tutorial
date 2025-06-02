
# Introduction
:::{tip} Are you in the right place?

This is *AMF analysis tutorial*, a dataset-specific tutorial.
If you're looking to get started with QIIME 2, or for general information about the project, there are better resources.

Are you looking for:
- the QIIME 2 homepage? That's [https://qiime2.org](https://qiime2.org).
- learning resources for microbiome marker gene (i.e., amplicon) analysis? See the [QIIME 2 *amplicon distribution* documentation](https://amplicon-docs.readthedocs.io).
- learning resources for microbiome metagenome analysis? See the [MOSHPIT documentation](https://moshpit.readthedocs.io).
- installation instructions, plugins, books, videos, workshops, or resources? See the [QIIME 2 Library](https://library.qiime2.org).
- general help? See the [QIIME 2 Forum](https://forum.qiime2.org)
- information about the author? See [ https://cap-lab.bio](https://cap-lab.bio).

Otherwise, if you're specifically looking for *AMF analysis tutorial*, you're in the right place.
[Read on... 📖](#tutorial)
:::

## Background
In this tutorial you’ll use the amplicon distribution of QIIME 2 to perform an analysis of arbuscular mycorrhizal fungal (AMF) diversity in the rhizosphere microbiome of traditional and modern varieties of rice soil samples from Bangladesh. This tutorial was developed as a part of seminar course BIO 698-063 and BIO498-008 4249 Seminar: Microbial Community Analysis (#4932) Spring 2025 for classroom teaching at School of Earth and Sustainability, Northern Arizona University by Manju M. Gupta, Nancy Jonnson and Catherine Gehring. Links to slides and recording wherever available is provided. This also includes a guest lecture by Greg Caporaso.

Arbuscular mycorrhizal (AM) fungi are ecologically and economically vital symbionts that enhance plant nutrient uptake and contribute to soil health (Smith & Read, 2008; Martin & van der Heijden, 2024). As with other fungi, molecular methods represent the state-of-the-art approach to characterise AM fungal diversity (Tedersoo et al. 2022). Therefore, there is a pressing need for a standardized, stepwise methodology—or a practical tutorial—that guides researchers in the accurate identification of AM fungal taxa from metagenomic datasets.

Initiatives like, GlobalAMF (Větrovský et al., 2023), SPUN, Austarlia (Frew et al. 2024) aim to map the global underground fungal network, underscoring the importance of robust and reproducible workflows for community-level AMF analysis and trained personnel in this direction. 

The dataset used in this tutorial originates from the study by Parvin et al. (2021) and represents a small, curated subset of the original data available on NCBI SRA, chosen to allow rapid processing on a personal computer. This dataset comprises 18S rRNA amplicon sequences targeting AMF in rice root samples, generated using the Illumina MiSeq platform (12,641 reads; 6.3 million bases; ~3.1 MB download). A 350 bp fragment of the small subunit (SSU) rRNA gene was PCR-amplified using Glomeromycota-specific primers AMV4.5NF and AMDGR, each tagged with sample-specific indices to enable dual-index sequencing. 


![image](https://github.com/user-attachments/assets/4eb8f6dd-540d-4cc9-8321-bbac41633e2f)

## License

 *AMF analysis tutorial* (©2025) by [Caporaso Lab](https://cap-lab.bio) is licensed under [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/deed.en).
