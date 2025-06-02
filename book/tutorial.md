(tutorial)=
# Tutorial

Sample metadata
This study investigates differences in AMF communities between modern high-yielding rice varieties (BR28, BR29, BR58) and local traditional varieties (Shampakatar, Ushapari), sampled from eight rice fields with five soil samples each. A total of 200 samples were collected, of which subsample of  143 AMF-relevant samples (79 from modern and 64 from traditional varieties), to test the hypothesis recorded in the metadata file. 
Before starting the analysis, explore the sample metadata to familiarize yourself with the samples used in this study. The following command will download the sample metadata as tab-separated text and save it in the file sample-metadata.tsv. This sample-metadata.tsv file is used throughout the rest of the tutorial. 

Obtaining the data
Explain NCBI SRA if needed 
Importing data

mkdir Rice_metagenomics

Download in directory Rice_metagenomics. The Rice_metagenomics directory contains a fastq subdirectory with 289 files, representing forward (R1) and reverse (R2) reads for 143 samples. The sequences are paired end sequences. Each sample has its own pair of FASTQ files, indicating that the sequences are already demultiplexed. It makes sense to study your fastq files before going ahead.
To import demultiplexed paired-end sequences into QIIME 2, the FASTQ files must follow the Cassava 1.8 format, which uses a specific naming convention. In this format, each file name includes the sample ID, sequencing lane, read direction, and other metadata.
For example, the forward and reverse reads for a single sample might be named:
•	L2S357_15_L001_R1_001.fastq.gz → Forward read
•	L2S357_15_L001_R2_001.fastq.gz → Reverse read
This format allows QIIME 2 to recognize and organize the files correctly during the import step.
The underscore-separated fields in this file name are: 
the sample identifier, L2S357
the barcode sequence or a barcode identifier, 15
the lane number, L001
the direction of the read i.e. R1 (forward) or R2 (reverse), the set number. 

Our data is not in Cassava 1.8 format, so we need to manually import it into QIIME 2 using a manifest file. Our sequences are named as SRR13445888_1.fastq etc.
Converting all 286 files to this format requires a small code included in the resource material creating _manifest_file.txt can be downloaded from here.

A manifest is a tab-separated .txt file that maps each sample name to the absolute file paths of its corresponding FASTQ files (either .fastq or .fastq.gz), including both forward and reverse reads if working with paired-end data.
A portion of manifest file will look like this.

sample-id	forward-absolute-filepath	reverse-absolute-filepath
SRR13445888	/Users/manjum.gupta/Desktop/Rice_metagenomics_data_files/fastq/SRR13445888_1.fastq	/Users/manjum.gupta/Desktop/Rice_metagenomics_data_files/fastq/SRR13445888_2.fastq
SRR13445889	/Users/manjum.gupta/Desktop/Rice_metagenomics_data_files/fastq/SRR13445889_1.fastq	/Users/manjum.gupta/Desktop/Rice_metagenomics_data_files/fastq/SRR13445889_2.fastq
SRR13445890	/Users/manjum.gupta/Desktop/Rice_metagenomics_data_files/fastq/SRR13445890_1.fastq	/Users/manjum.gupta/Desktop/Rice_metagenomics_data_files/fastq/SRR13445890_2.fastq
SRR13445891	/Users/manjum.gupta/Desktop/Rice_metagenomics_data_files/fastq/SRR13445891_1.fastq	/Users/manjum.gupta/Desktop/Rice_metagenomics_data_files/fastq/SRR13445891_2.fastq

We will import the sequences as SampleData[SequencesWithQuality], which is the demultiplexed single-end sequence format. If we wanted to import paired-sequences, we would specify the semantic type SampleData[PairedEndSequencesWithQuality].

  manifest file needs to be created and saved in the fastq files of your system because the path in this file will work only on my system (/users/manjum.gupta/Desktop/Rice_metagenomics_data_files/fastq/SRR13445888_1.fastq) . Run the commands provided to you in creating _manifest_file.txt to create manifest file for path fastq files in your laptop.


a.	Creating a manifest file 

cd fastq
mkdir manifest

echo "# paired-end PHRED 33 fastq manifest file for forward and reverse reads" > manifest1.txt

#Create a new file called manifest1.txt.Insert the text # paired-end PHRED 33 fastq manifest file for forward and reverse reads into that file as a comment (indicated by the # at the start).

echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" >> manifest1.txt

#will append the following line to the manifest1.txt file:
sample-id  forward-absolute-filepath  reverse-absolute-filepath for your system

ls *.fastq | cut -d "_" -f 1 | sort | uniq | while read prefix; do     echo -e "$prefix\t$PWD/${prefix}_1.fastq\t$PWD/${prefix}_2.fastq"; done > manifest2.txt 
 
#The command you've provided generates a manifest file for paired-end sequencing reads.  
cat manifest1.txt manifest2.txt  > manifest/manifest.tsv
#Delete text files
rm *.txt
cd .. 

 tip : # it’s a good idea to view your manifest file by the command cat  command 

cat manifest/manifest.tsv

It will show path of your files.


b.	Importing the data

cd  Rice_metgenomics 

# Activate the qiime environment 

conda activate qiime2-amplicon-2024.5

#Change the path to your manifest file

qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path fastq/manifest/manifest.tsv \
        --output-path demux.qza \
        --input-format PairedEndFastqManifestPhred33V2




Sequence data visualization

In QIIME 2, all data is structured as an Artifact of a specific semantic type. Artifacts contain the data as well as information about the data, including a record of the original data and the tools used to process it. This allows for better tracking of how you actually got to where you are in your analysis. You can learn more about common QIIME 2 Artifacts and semantic types here.

qiime demux summarize \
    --i-data demux.qza \
    --o-visualization visualization.qzv



![image](https://github.com/user-attachments/assets/e2e126e6-429b-4b21-a5e1-c713f62d0c97)

:::{note}
This document was built with its own conda environment.
You can download the environment file that was used from the download link on the top-right of this article.
:::

:::{describe-usage}
:scope: amf-tutorial
project_accession = use.init_artifact_from_url('project-accession',
                                              'https://www.dropbox.com/scl/fi/h470qev6vqrzir3f7yksk/project_id.qza?rlkey=ra7jno6nzs0m2ddkynvgt3w14&st=rmpx5ndn&dl=1')
:::

:::{describe-usage}

project_accession_md = use.view_as_metadata('project_accession_md', project_accession)

use.action(
    use.UsageAction(plugin_id='metadata',
                    action_id='tabulate'),
    use.UsageInputs(input=project_accession_md),
    use.UsageOutputNames(visualization='project-accession'))
:::

:::{describe-usage}
demux = use.init_artifact_from_url('demux',
                                   'https://www.dropbox.com/scl/fi/oklnxxh6wruw9i84m07fr/paired_reads.qza?rlkey=cwf85atfi6j5dyjzas51xorq0&st=9cxj7l5t&dl=1')
:::

:::{describe-usage}
use.action(
    use.UsageAction(plugin_id='demux',
                    action_id='summarize'),
    use.UsageInputs(data=demux),
    use.UsageOutputNames(visualization='demux'))
:::


:::{describe-usage}
metadata = use.init_metadata_from_url('metadata',
                                   'https://www.dropbox.com/scl/fi/zgcnuetdxochydkb0o3bw/metadata_file_rice.tsv?rlkey=fo5ywq8549fn5optv1u1nh4m4&st=hzljxvx3&dl=1')
:::

# Denoising
TODO: TIME RUNNING 11 minutes
:::{describe-usage}
table, denoising_stats, representative_sequences = use.action(
    use.UsageAction(plugin_id='dada2', action_id='denoise_paired'),
    use.UsageInputs(
        demultiplexed_seqs=demux,
        trunc_len_f=240, #different from prov replay [200, 200],[240,240],[220,200], changed to match pptx
        trunc_len_r=220,
        ),
    use.UsageOutputNames(
        table='table',
        denoising_stats='denoising_stats',
        representative_sequences='representative_sequences'
        )
)
:::

:::{describe-usage}
denoising_stats_md = use.view_as_metadata('denoising_stats_md', denoising_stats)

dada2_stats_viz = use.action(
    use.UsageAction(plugin_id='metadata', action_id='tabulate'),
    use.UsageInputs(
        input=denoising_stats_md),
    use.UsageOutputNames(visualization='stats-dada2'))
:::

:::{describe-usage}
rep_seqs_viz = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='tabulate_seqs'),
    use.UsageInputs(
        data=representative_sequences),
    use.UsageOutputNames(visualization='rep-seqs'))
:::

:::{describe-usage}
table_summary = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='summarize'),
    use.UsageInputs(
        table=table),
    use.UsageOutputNames(visualization='table'))
:::

# TODO Reference databases and Taxonomic classification, feature table filtering and building the phylogenetic tree

:::{describe-usage}
def maarjam_refseq_factory():
    from urllib import request
    from qiime2 import Artifact
    # Download the FASTA file
    fp, _ = request.urlretrieve(
        'https://www.dropbox.com/scl/fi/kedvuma4xdl3blf3mb8ig/maarjam_database_SSU_TYPE.qiime.fasta?rlkey=zoj00nmnz97igynnykj0hyqxj&st=xtj7ueg1&dl=1')
    # Import as a QIIME 2 artifact
    return Artifact.import_data(
        'FeatureData[Sequence]', fp)

otus_85_maarjam = use.init_artifact('otus_85_maarjam', maarjam_refseq_factory)
:::


:::{describe-usage}
def maarjam_taxonomy_factory():
    from urllib import request
    from qiime2 import Artifact
    # Download the taxonomy TSV file
    fp, _ = request.urlretrieve(
        'https://www.dropbox.com/scl/fi/qe2z05fh9pe45eytl6f96/maarjam_database_SSU_TYPE.qiime.txt?rlkey=4ma85ti0l94378ctwxf4iwj5n&st=a2oc7ppj&dl=1')
    # Import as a QIIME 2 artifact
    return Artifact.import_data(
        'FeatureData[Taxonomy]', fp, view_type='HeaderlessTSVTaxonomyFormat')

ref_taxonomy_maarjam = use.init_artifact('ref_taxonomy_maarjam', maarjam_taxonomy_factory)

:::

:::{describe-usage}
ref_seqs_maarjam, = use.action(
    use.UsageAction(plugin_id='feature_classifier', action_id='extract_reads'),
    use.UsageInputs(
        sequences=otus_85_maarjam,
        f_primer='AAGCTCGTAGTTGAATTTCG',
        r_primer='CCCAACTATCCCTATTAATCAT',
        trunc_len=250, # mismatched from 120 in provenance replay
        min_length=100,
        max_length=400,
        ),
    use.UsageOutputNames(reads='ref_seqs_maarjam'))
:::
<!-- :::{describe-usage}
classifier_maarjam, = use.action(
    use.UsageAction(plugin_id='feature_classifier', action_id='fit_classifier_naive_bayes'),
    use.UsageInputs(
        reference_reads=ref_seqs_maarjam,
        reference_taxonomy=ref_taxonomy_maarjam,
        ),
    use.UsageOutputNames(classifier='classifier_maarjam'))
::: -->
```
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads ref-seqs-maarjam.qza \
    --i-reference-taxonomy ref-taxonomy-maarjam.qza \
    --o-classifier classifier-maarjam.qza \
```

:::{describe-usage}
classifier_maarjam = use.init_artifact_from_url('classifier_maarjam',
                                   'https://www.dropbox.com/scl/fi/08aeme2zspkrt3q5bjmtu/classifier-maarjam.qza?rlkey=acl2ch723nf13spq25sp48ro4&st=itnbg5bu&dl=1')
:::
<!-- :::{describe-usage}
taxonomy_maarjam, = use.action(
    use.UsageAction(plugin_id='feature_classifier', action_id='classify_sklearn'),
    use.UsageInputs(
        reads=representative_sequences,
        classifier=classifier_maarjam,
        confidence=0.7,
        n_jobs=1,
        read_orientation='auto'),
    use.UsageOutputNames(classification='taxonomy_maarjam'))
::: -->

```
qiime feature-classifier classify-sklearn \
    --i-reads representative-sequences.qza \
    --i-classifier classifier-maarjam.qza \
    --p-confidence 0.7 \
    --o-classification taxonomy-maarjam.qza
```

:::{describe-usage}
taxonomy_maarjam = use.init_artifact_from_url('taxonomy_maarjam',
                                   'https://www.dropbox.com/scl/fi/my2oc0qt0t67u1kgwcrlt/taxonomy-maarjam.qza?rlkey=f3x38pyyhhb2aq7yrdsbabjw9&st=wz36byr7&dl=1')
:::

:::{describe-usage}
taxonomy_maarjam_md = use.view_as_metadata('taxonomy_maarjam_md', taxonomy_maarjam)
use.action(
    use.UsageAction(plugin_id='metadata', action_id='tabulate'),
    use.UsageInputs(
        input=taxonomy_maarjam_md),
    use.UsageOutputNames(
        visualization='taxonomy_maarjam_md'))
:::

# day 8 

```
 qiime feature-classifier classify-consensus-vsearch \
    --i-query representative-sequences.qza \
    --i-reference-reads otus-85-maarjam.qza \
    --i-reference-taxonomy ref-taxonomy-maarjam.qza \
    --p-perc-identity 0.99 \
    --p-top-hits-only \
    --p-maxaccepts 1 \
    --p-strand 'both' \
    --p-unassignable-label 'Unassigned' \
    --o-classification rice-taxonomy-maarjAM.qza \ --o-search-results rice-search-results-maarjAM.qza 
```

:::{describe-usage}
rice_taxonomy_maarjAM, rice_search_results_maarjAM = use.action(
    use.UsageAction(plugin_id='feature_classifier', action_id='classify_consensus_vsearch'),
    use.UsageInputs(
        query=representative_sequences,
        reference_reads=otus_85_maarjam,
        reference_taxonomy=ref_taxonomy_maarjam,
        maxaccepts=1,
        perc_identity=0.99,
        strand='both',
        top_hits_only=True,
        unassignable_label='Unassigned'),
    use.UsageOutputNames(
        classification='rice_taxonomy_maarjAM',
        search_results='rice_search_results_maarjAM'))
:::

:::{describe-usage}
table_rice_unassigned_maarjAM, = use.action(
    use.UsageAction(plugin_id='taxa', action_id='filter_table'),
    use.UsageInputs(
        table=table, #table_rice?
        taxonomy=rice_taxonomy_maarjAM,
        include='Unassigned',
        query_delimiter=',',
        mode='contains'),
    use.UsageOutputNames(filtered_table='table_rice_unassigned_maarjAM'))
:::

:::{describe-usage}
rep_seq_rice_unassigned_maarjAM, = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='filter_seqs'),
    use.UsageInputs(
        data=representative_sequences,
        table=table_rice_unassigned_maarjAM),
    use.UsageOutputNames(filtered_data='rep_seq_rice_unassigned_maarjAM'))
:::

:::{describe-usage}
filtered_unassigned_table, = use.action(
    use.UsageAction(plugin_id='taxa', action_id='filter_table'),
    use.UsageInputs(
        table=table, #table_rice?
        taxonomy=rice_taxonomy_maarjAM,
        include='c__Glomeromycetes',
        query_delimiter=',',
        mode='contains'),
    use.UsageOutputNames(filtered_table='maarjAM_Glomeromycetes_seqs'))
:::

# silva
:::{describe-usage}
def silva_seqs_factory():
    from urllib import request
    from qiime2 import Artifact
    url = 'https://data.qiime2.org/2024.2/common/silva-138-99-seqs.qza'
    fp, _ = request.urlretrieve(url)
    return Artifact.load(fp)

silva_138_99_seqs = use.init_artifact('silva_138_99_seqs', silva_seqs_factory)
:::

:::{describe-usage}
def silva_tax_factory():
    from urllib import request
    from qiime2 import Artifact
    url = 'https://data.qiime2.org/2024.2/common/silva-138-99-tax.qza'
    fp, _ = request.urlretrieve(url)
    return Artifact.load(fp)

silva_138_99_tax = use.init_artifact('silva138_99_tax', silva_tax_factory)
:::

## This is what is taking so long
```
 qiime feature-classifier classify-consensus-vsearch \
    --i-query rep-seq-rice-unassigned-maarjAM.qza \
    --i-reference-reads silva-138-99-seqs.qza \
    --i-reference-taxonomy silva138-99-tax.qza \
    --p-perc-identity 0.99 \
    --p-top-hits-only \
    --p-maxaccepts 1 \
    --p-strand 'both' \
    --p-unassignable-label 'Unassigned' \
    --o-classification rice-taxonomy-assigned-silva.qza --o-search-results rice-search-results-assigned-silva.qza 
```

<!-- :::{describe-usage}
rice_taxonomy_assigned_silva, rice_search_results_assigned_silva = use.action(
    use.UsageAction(plugin_id='feature_classifier', action_id='classify_consensus_vsearch'),
    use.UsageInputs(
        query=rep_seq_rice_unassigned_maarjAM,
        reference_reads=silva_138_99_seqs,
        reference_taxonomy=silva_138_99_tax,
        perc_identity=0.99,
        top_hits_only=True,
        maxaccepts=1,
        strand='both',
        unassignable_label='Unassigned'),
    use.UsageOutputNames(
        classification='rice_taxonomy_assigned_silva',
        search_results='rice_search_results_assigned_silva'))
:::-->
:::{describe-usage}
rice_taxonomy_assigned_silva = use.init_artifact_from_url(
                                'rice_taxonomy_assigned_silva', 'https://www.dropbox.com/scl/fi/q84hfc1ghtrthael2ij77/rice-taxonomy-assigned-silva.qza?rlkey=7tz6odze2rpglj8ubhq5m9t6h&st=qguiyovg&dl=1')
:::

:::{describe-usage}
merged_taxonomy, = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='merge_taxa'),
    use.UsageInputs(
        data=[rice_taxonomy_maarjAM, rice_taxonomy_assigned_silva]),
    use.UsageOutputNames(merged_data='merged_taxonomy'))
:::

:::{describe-usage}
taxa_bar_plots = use.action(
    use.UsageAction(plugin_id='taxa', action_id='barplot'),
    use.UsageInputs(
        table=table,
        taxonomy=merged_taxonomy,
        metadata=metadata),
    use.UsageOutputNames(visualization='taxa_bar_plots'))
:::

# Filtering taxa

:::{describe-usage}
traditional_rice_table, = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='filter_samples'),
    use.UsageInputs(
        table=table,
        metadata=metadata,
        where="[env_broad_scale]='Traditional rice root'"),
    use.UsageOutputNames(filtered_table='traditional_rice_table'))
:::

:::{describe-usage}
traditional_rice_table_viz = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='summarize'),
    use.UsageInputs(
        table=traditional_rice_table),
    use.UsageOutputNames(visualization='traditional_rice_table'))
:::

:::{describe-usage}
modern_rice_table, = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='filter_samples'),
    use.UsageInputs(
        table=table,
        metadata=metadata,
        where="[env_broad_scale]='Modern rice root'"),
    use.UsageOutputNames(filtered_table='modern_rice_table'))
:::

:::{describe-usage}
modern_rice_table_viz = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='summarize'),
    use.UsageInputs(
        table=modern_rice_table),
    use.UsageOutputNames(visualization='modern_rice_table'))
:::

# Differential abundance (ANCOM-BC)

:::{describe-usage}
glomeromycetes_table, = use.action(
    use.UsageAction(plugin_id='taxa', action_id='filter_table'),
    use.UsageInputs(
        table=table,
        taxonomy=merged_taxonomy,
        include='p__Mucoromycota'),
    use.UsageOutputNames(filtered_table='glomeromycetes_table'))
:::

:::{describe-usage}
collapsed_table_level_6, = use.action(
    use.UsageAction(plugin_id='taxa', action_id='collapse'),
    use.UsageInputs(
        table=glomeromycetes_table,
        taxonomy=merged_taxonomy,
        level=6),
    use.UsageOutputNames(collapsed_table='collapsed_table_level_6'))
:::

:::{describe-usage}
l6_ancombc_diffs, = use.action(
    use.UsageAction(plugin_id='composition', action_id='ancombc'),
    use.UsageInputs(
        table=collapsed_table_level_6,
        metadata=metadata,
        formula='env_broad_scale'),
    use.UsageOutputNames(differentials='l6_ancombc_differentials'))
:::

:::{describe-usage}
l6_da_barplot = use.action(
    use.UsageAction(plugin_id='composition', action_id='da-barplot'),
    use.UsageInputs(
        data=l6_ancombc_differentials,
        significance_threshold=0.01,
        level_delimiter=';'),
    use.UsageOutputNames(visualization='l6_da_barplot'))
:::

# alpha rarefaction and core metrics

:::{describe-usage}
rooted_tree, unrooted_tree, aligned_rep_seqs, masked_aligned_rep_seqs = use.action(
    use.UsageAction(plugin_id='phylogeny', action_id='align_to_tree_mafft_fasttree'),
    use.UsageInputs(
        sequences=representative_sequences,
        ),
    use.UsageOutputNames(
        rooted_tree='rooted_tree',
        tree='unrooted_tree',
        alignment='aligned_rep_seqs',
        masked_alignment='masked_aligned_rep_seqs'))
:::

:::{describe-usage}
alpha_rarefaction_plot = use.action(
    use.UsageAction(plugin_id='diversity', action_id='alpha_rarefaction'),
    use.UsageInputs(
        table=table,
        phylogeny=rooted_tree,
        max_depth=13299,
        metadata=metadata,
        ),
    use.UsageOutputNames(visualization='alpha_rarefaction'))
:::

:::{describe-usage}
core_metrics_outputs = use.action(
    use.UsageAction(plugin_id='diversity', action_id='core_metrics_phylogenetic'),
    use.UsageInputs(
        table=table,
        phylogeny=rooted_tree,
        sampling_depth=1602,
        metadata=metadata,
        ),
    use.UsageOutputNames(
        evenness_vector='evenness_vector',
        faith_pd_vector='faith_pd_vector',
        unweighted_unifrac_distance_matrix='unweighted_unifrac_distance_matrix',
        bray_curtis_pcoa_results='bray-curtis-pcoa-results',
        shannon_vector='shannon_vector',
        rarefied_table='rarefied-table',
        weighted_unifrac_distance_matrix='weighted-unifrac-distance-matrix',
        jaccard_pcoa_results='jaccard-pcoa-results',
        unweighted_unifrac_emperor='unweighted-unifrac-emperor',
        weighted_unifrac_pcoa_results='weighted-unifrac-pcoa-results',
        observed_features_vector='observed-features-vector',
        jaccard_distance_matrix='jaccard-distance-matrix',
        jaccard_emperor='jaccard-emperor',
        bray_curtis_emperor='bray-curtis-emperor',
        weighted_unifrac_emperor='weighted-unifrac-emperor',
        bray_curtis_distance_matrix='bray-curtis-distance-matrix',
        unweighted_unifrac_pcoa_results='unweighted-unifrac-pcoa-results'))
:::

:::{describe-usage}
shannon_group_significance = use.action(
    use.UsageAction(plugin_id='diversity', action_id='alpha_group_significance'),
    use.UsageInputs(
        alpha_diversity=core_metrics_outputs.shannon_vector,
        metadata=metadata),
    use.UsageOutputNames(visualization='shannon_group_significance'))
:::

:::{describe-usage}
evenness_group_significance = use.action(
    use.UsageAction(plugin_id='diversity', action_id='alpha_group_significance'),
    use.UsageInputs(
        alpha_diversity=core_metrics_outputs.evenness_vector,
        metadata=metadata),
    use.UsageOutputNames(visualization='evenness_group_significance'))
:::

:::{describe-usage}
faith_pd_group_significance = use.action(
    use.UsageAction(plugin_id='diversity', action_id='alpha_group_significance'),
    use.UsageInputs(
        alpha_diversity=core_metrics_outputs.faith_pd_vector,
        metadata=metadata),
    use.UsageOutputNames(visualization='faith_group'))
:::

:::{describe-usage}

env_broad_col = use.get_metadata_column('env_broad_scale', 'env_broad_scale', metadata)

beta_group_significance_simple = use.action(
    use.UsageAction(plugin_id='diversity', action_id='beta_group_significance'),
    use.UsageInputs(
        distance_matrix=core_metrics_outputs.unweighted_unifrac_distance_matrix,
        metadata=env_broad_col,
        method='permanova',
        pairwise=True,
        permutations=999),
    use.UsageOutputNames(visualization='unweighted_unifrac_subject_group_significance'))
:::

:::{describe-usage}

sample_name_col = use.get_metadata_column('Sample_Name', 'Sample_Name', metadata)

beta_group_significance = use.action(
    use.UsageAction(plugin_id='diversity', action_id='beta_group_significance'),
    use.UsageInputs(
        distance_matrix=core_metrics_outputs.unweighted_unifrac_distance_matrix,
        metadata=sample_name_col,
        method='permanova',
        pairwise=True,
        permutations=999),
    use.UsageOutputNames(visualization='unweighted_unifrac_body_site_significance'))
:::

