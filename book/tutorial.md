(tutorial)=
# Tutorial
:::{note}
This document was built with its own conda environment.
You can download the environment file that was used from the download link on the top-right of this article.
:::

## Sample metadata

This study investigates differences in AMF communities between modern high-yielding rice varieties (BR28, BR29, BR58) and local traditional varieties (Shampakatar, Ushapari), sampled from eight rice fields with five soil samples each. A total of 200 samples were collected, of which subsample of  143 AMF-relevant samples (79 from modern and 64 from traditional varieties), to test the hypothesis recorded in the metadata file. 
Before starting the analysis, explore the sample metadata to familiarize yourself with the samples used in this study. The following command will download the sample metadata as tab-separated text and save it in the file sample-metadata.tsv. This sample-metadata.tsv file is used throughout the rest of the tutorial.

:::{describe-usage}
:scope: amf-tutorial
metadata = use.init_metadata_from_url('metadata',
                                   'https://www.dropbox.com/scl/fi/zgcnuetdxochydkb0o3bw/metadata_file_rice.tsv?rlkey=fo5ywq8549fn5optv1u1nh4m4&st=hzljxvx3&dl=1')
:::


## Obtaining the data
In this tutorial, we will use (q2-fondue)[https://library.qiime2.org/plugins/bokulich-lab/q2-fondue] to download our publically available data and import them into QIIME 2. 

First, Lets download our metadata that has the SRA project accession number.


:::{describe-usage}
project_accession = use.init_artifact_from_url('project-accession',
                                              'https://www.dropbox.com/scl/fi/h470qev6vqrzir3f7yksk/project_id.qza?rlkey=ra7jno6nzs0m2ddkynvgt3w14&st=rmpx5ndn&dl=1')
:::
Then we can visualize it using qiime metadata tabulate. 

:::{describe-usage}

project_accession_md = use.view_as_metadata('project_accession_md', project_accession)

use.action(
    use.UsageAction(plugin_id='metadata',
                    action_id='tabulate'),
    use.UsageInputs(input=project_accession_md),
    use.UsageOutputNames(visualization='project-accession'))
:::

We can use q2-fondue 🫕 to easily import the data using this command: 

```
qiime fondue get_sequences \
    --p-accession-ids project_id.qza \
    --p-email [Insert Your Email] \
    --o-single-reads single-reads-demux.qza \
    --o-paired-reads demux.qza \
    --o-failed-runs failed-runs.qza
```
 or you can download the demux artifact here:

:::{describe-usage}
demux = use.init_artifact_from_url('demux',
                                   'https://www.dropbox.com/scl/fi/oklnxxh6wruw9i84m07fr/paired_reads.qza?rlkey=cwf85atfi6j5dyjzas51xorq0&st=9cxj7l5t&dl=1')
:::
Now lets vizulaze our data using demux summarize

:::note
In QIIME 2, all data is structured as an Artifact of a specific semantic type. Artifacts contain the data as well as information about the data, including a record of the original data and the tools used to process it. This allows for better tracking of how you actually got to where you are in your analysis. You can learn more about common QIIME 2 Artifacts and semantic types here.
:::

:::{describe-usage}
use.action(
    use.UsageAction(plugin_id='demux',
                    action_id='summarize'),
    use.UsageInputs(data=demux),
    use.UsageOutputNames(visualization='demux'))
:::

## #Understanding the Data

If you downloaded your data using q2-fondue, or if you're working with data provided by your sequencing facility, you will typically need three types of files:

1. Metadata file – This mapping file provides contextual information about the experiment, such as the hypothesis, treatments, and sample details.

2. Classifier file – This file is used for taxonomic classification.

3. Sequence files – These are in FASTQ format or imported to qiime2 as .qza files.


To explore your sequencing data and assess quality, open the demux.qzv file in explorer. This visualization will help you examine the sequence quality across reads, which is essential for deciding how much to trim in the next steps.

## Revision Questions

1. What is the minimum and maximum number of reads in our samples?

2. Do any of the samples have fewer than 1,000 sequences?

3. At which position does the median quality score drop below 30?

Note: If any samples have very few sequences (e.g., fewer than 1,000), you may want to omit them from downstream analysis, as they could negatively affect data interpretation.





![image](https://github.com/user-attachments/assets/737d4d72-a40d-4e2b-8815-7c75981a112a)



  
![image](https://github.com/user-attachments/assets/dd99f519-040e-40e7-9776-345b5ffe8a57)


# Denoising Using DADA2
Denoising is the process of correcting errors in the sequencing data and delimitating ASVs( amplicon sequence variants). The Non-biological sequences (e.g., adapters, primers, linker pads, etc.) and errors created by sequencing machines, such as incorrect base calls or random noise which can lead to inaccurate results if not corrected.


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

## Revision Questions
1.  How did we decide    --p-trunc-len-f 240 \ --p-trunc-len-r 220 \ ?
Hint: At what base pair does the median quality drop below 30? 

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



## Revision Questions

1. Do BLAST searches of the representative sequences make sense? Are the features what you would expect like AMF or not?
2. How many features (ASVs) were generated? Are the communities high or low diversity?

Note- If you have a large number (e.g. >50%) of sequences been lost during denoising/filtering? If so, the settings might be too stringent.


# Taxonomic assignment 
## Classifications with fit-classifier-naive-bayes 
### Classification with maarjAM

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
:::{note}
Due to time constraints, this command is not generated in this notebook. The resulting artifacted can be generated by running this command or running wget the artifact below. 👇
:::

```code
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads ref-seqs-maarjam.qza \
    --i-reference-taxonomy ref-taxonomy-maarjam.qza \
    --o-classifier classifier-maarjam.qza 
```

:::{describe-usage}
classifier_maarjam = use.init_artifact_from_url('classifier_maarjam',
                                   'https://www.dropbox.com/scl/fi/08aeme2zspkrt3q5bjmtu/classifier-maarjam.qza?rlkey=acl2ch723nf13spq25sp48ro4&st=itnbg5bu&dl=1')
:::
:::{note}
Due to time constraints, this command is not generated in this notebook. The resulting artifacted can be generated by running this command or running wget the artifact below. 👇
:::

```code
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

## Revision Questions

1.  How many of your ASVs were taxonomically assigned?
   note :view your taxonomy.qzv file
   
## Classifications with classify-consensus-vsearch
### Classification with maarjAM
:::{describe-usage}
rice_taxonomy_maarjAM, rice_search_results_maarjAM = use.action(
    use.UsageAction(plugin_id='feature_classifier', action_id='classify_consensus_vsearch'),
    use.UsageInputs(
        query=representative_sequences,
        reference_reads=otus_85_maarjam,
        reference_taxonomy=ref_taxonomy_maarjam,
        maxaccepts=1,
        perc_identity=0.7,
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

## Revision Questions

1.  How many of your ASVs were taxonomically assigned? 



```



```code
 qiime feature-classifier classify-consensus-vsearch \
    --i-query rep-seq-rice-unassigned-maarjAM.qza \
    --i-reference-reads silva-138-99-seqs.qza \
    --i-reference-taxonomy silva138-99-tax.qza \
    --p-perc-identity 0.99 \
    --p-top-hits-only \
    --p-maxaccepts 1 \
    --p-strand 'both' \
    --p-unassignable-label 'Unassigned' \
    --o-classification rice-taxonomy-assigned-silva.qza \
    --o-search-results rice-search-results-assigned-silva.qza 
```

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

# Note 
Why Vsearch ?

It decides on the basis of global alignment of sequences.
It allows us to pass stricter parameter for identity.
I can get unassignable as output  to run through the other database.



# Filtering Tables

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
l6_ancombc_differentials, = use.action(
    use.UsageAction(plugin_id='composition', action_id='ancombc'),
    use.UsageInputs(
        table=collapsed_table_level_6,
        metadata=metadata,
        formula='env_broad_scale'),
    use.UsageOutputNames(differentials='l6_ancombc_differentials'))
:::

:::{describe-usage}
l6_da_barplot = use.action(
    use.UsageAction(plugin_id='composition', action_id='da_barplot'),
    use.UsageInputs(
        data=l6_ancombc_differentials,
        significance_threshold=0.01,
        level_delimiter=';'),
    use.UsageOutputNames(visualization='l6_da_barplot'))
:::

# Diversity Analysis

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
## Alpha Diversity and Community Richness
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

## Beta Diversity 

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

