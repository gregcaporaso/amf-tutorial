
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
Arbuscular mycorrhizal fungi (AMF) are ecologically and economically vital symbionts that enhance plant nutrient uptake and play a key role in maintaining soil health (Smith & Read, 2010; Martin & van der Heijden, 2024). As with other fungal groups, molecular techniques have become the standard for characterizing AMF diversity (Tedersoo et al., 2022). Global initiatives focused on mapping the distribution and diversity of AMF—such as GlobalAMF (Větrovský et al., 2023), SPUN, and the Australian AMF Mapping Project (Frew et al., 2024)—emphasize the need for accessible, robust, and reproducible workflows. These efforts also highlight the growing demand for trained personnel capable of conducting community-level AMF analyses in this rapidly advancing field.

In this tutorial you’ll use the amplicon distribution of QIIME 2 to perform an analysis of AMF diversity in the rhizosphere microbiome. This tutorial was developed as a part of seminar course BIO 698-063 and BIO498-008 4249 in Spring 2025, Northern Arizona University by Manju M. Gupta, Nancy Johnson and Catherine Gehring. Where available, links to instructional slides and recordings are provided. The tutorial also includes a guest lecture by Greg Caporaso.

The dataset used in this tutorial originates from the study by Parvin et al. (2021) and represents a small, curated subset of the original data available on NCBI SRA, chosen to allow rapid processing on a personal computer. This dataset comprises 18S rRNA amplicon sequences targeting AMF in rice root samples, generated using the Illumina MiSeq platform ( ~3.1 MB download). A 350 bp fragment of the small subunit (SSU) rRNA gene was PCR-amplified using Glomeromycota-specific primers AMV4.5NF and AMDGR, each tagged with sample-specific indices to enable dual-index sequencing. 

References
Smith, S.E. and Read, D.J., 2010. Mycorrhizal symbiosis. Academic press.
Martin, F.M. and van Der Heijden, M.G., 2024. The mycorrhizal symbiosis: research frontiers in genomics, ecology, and agricultural application. New Phytologist, 242(4), pp.1486-1506
Tedersoo et al., 2022. Global patterns in endemicity and vulnerability of soil fungi. Global change biology, 28(22), pp.6696-6710.
Větrovský et al., 2023. GlobalAMFungi: a global database of arbuscular mycorrhizal fungal occurrences from high‐throughput sequencing metabarcoding studies. New Phytologist, 240(5), pp.2151-2163.
Frew et al., 2024. AusAMF: database of arbuscular mycorrhizal fungal communities in Australia. bioRxiv, pp.2024-09.
Parvin et al. (2021) Abundance, characteristics and variation of microplastics in different freshwater fish species from Bangladesh. Science of the Total Environment, 784, p.147137.
![image](https://github.com/user-attachments/assets/a6c7282c-6e94-425b-a64b-f0e2f63089ca)



![image](https://github.com/user-attachments/assets/4eb8f6dd-540d-4cc9-8321-bbac41633e2f)

## License

 *AMF analysis tutorial* (©2025) by [Caporaso Lab](https://cap-lab.bio) is licensed under [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/deed.en).
