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
In this tutorial, we will use [q2-fondue](https://library.qiime2.org/plugins/bokulich-lab/q2-fondue) to download our publically available data and import them into QIIME 2.

First, let’s download the metadata containing the NCBI SRA project, accession number SRR13445888.

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
                                   'https://www.dropbox.com/scl/fi/81fwmwgymuq1ae4o18xpb/subsampled_demux.qza?rlkey=yyzr6inkmt85xv6i54d7qitml&st=3pjj0kd5&dl=1')
:::
Now, let’s visualize our data using demux summarize to assess sequencing quality.

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

### Understanding the Data

If you downloaded your data using q2-fondue, or received it from a sequencing facility, you will typically need three key files:

Metadata file – This mapping file provides contextual information about your experiment, including hypotheses, treatments, and sample details.

Classifier file – Used for assigning taxonomy during analysis.

Sequence files – These are usually in FASTQ format or already imported into QIIME 2 as .qza files.

To begin exploring your data and assessing its quality, open the demux.qzv file in QIIME 2 View. This visualization helps you evaluate read quality, which is critical for determining appropriate trimming parameters in downstream steps.

## Revision Questions

1. What is the minimum and maximum number of reads in our samples?

2. Do any of the samples have fewer than 1,000 sequences?

3. At which position does the median quality score drop below 30?

Note: If any of the samples have very few sequences (e.g., fewer than 1,000), you may want to omit them from downstream analysis, as they could negatively affect data interpretation.

# Denoising Using DADA2
Denoising is the process of correcting errors in the sequencing data and delimitating ASVs (amplicon sequence variants). The Non-biological sequences (e.g., adapters, primers, linker pads, etc.) and errors created by sequencing machines, such as incorrect base calls or random noise which can lead to inaccurate results if not corrected.

:::{describe-usage}
table, denoising_stats, representative_sequences = use.action(
    use.UsageAction(plugin_id='dada2', action_id='denoise_paired'),
    use.UsageInputs(
        demultiplexed_seqs=demux,
        trunc_len_f=240,
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

Note- If a large number (e.g. >50%) of sequences are lost during denoising/filtering the settings may be too stringent.

:::{describe-usage}
rep_seqs_viz = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='tabulate_seqs'),
    use.UsageInputs(
        data=representative_sequences),
    use.UsageOutputNames(visualization='rep-seqs'))
:::
## Revision Questions

1. Do BLAST searches of the representative sequences make sense? Are the features what you would expect like AMF or not?
2. How many features (ASVs) were generated? Are the communities high or low diversity?

:::{describe-usage}
table_summary = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='summarize'),
    use.UsageInputs(
        table=table),
    use.UsageOutputNames(visualization='table'))
:::

# Taxonomic assignment 
## Classifications with fit-classifier-naive-bayes 

To construct a taxonomic classifier using the MaarjAM database, two primary input files are required:

1. A FASTA file containing the reference sequences
2. A taxonomy file mapping those sequences to their taxonomic lineages

These resources can be downloaded from the official MaarjAM database website. Alternatively, pre-imported .qza files are also be downloaded from here:

:::{describe-usage}
def maarjam_refseq_factory():
    from urllib import request
    from qiime2 import Artifact
    # Download the FASTA file
    fp, _ = request.urlretrieve(
        'https://www.dropbox.com/scl/fi/degicmkcyidve4qc3z3jh/maarjam_ref_seq.qza?rlkey=klzt244f0jylqeloxwqp61tb2&st=v6vpxfr0&dl=1')
    # Import as a QIIME 2 artifact
    return Artifact.load(fp)

maarjam_ref_seq = use.init_artifact('maarjam_ref_seq', maarjam_refseq_factory)
:::

:::{describe-usage}
def maarjam_taxonomy_factory():
    from urllib import request
    from qiime2 import Artifact
    # Download the taxonomy TSV file
    fp, _ = request.urlretrieve(
        'https://www.dropbox.com/scl/fi/93gm5lfsn8kli2elemigs/ref-taxonomy.qza?rlkey=afvwcqdjgotn21ermz9rx015c&st=1fhecg2i&dl=1')
    # Import as a QIIME 2 artifact
    return Artifact.load(fp)

ref_taxonomy_maarjam = use.init_artifact('ref_taxonomy_maarjam', maarjam_taxonomy_factory)

:::

:::{note}
Due to time constraints, this command is not generated in this notebook. The resulting artifacts can be generated by running this command or running wget the artifact below. 👇
:::

Note : Building an accurate classifier file is crucial for reliable taxonomic analysis, as it can significantly influence your results. An incorrect or poorly trained classifier may lead to unassigned or misclassified sequences. It's recommended to run your data using both the vsearch and sklearn methods for comparison. Details on both approaches are provided [here].
# comment for chloe please add a link to more details about sklearn and Vsearch methods in qiime docs 


# comment for chloe please note the followings:

```code
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads ref-seqs-maarjam.qza \
    --i-reference-taxonomy ref-taxonomy-maarjam.qza \
    --o-classifier classifier-maarjam.qza 
```

:::{describe-usage}
classifier_maarjam = use.init_artifact_from_url('classifier_maarjam',
                                   'https://www.dropbox.com/scl/fi/bh65ab79wmo9calgiwsxr/classifier-maarjam-1.qza?rlkey=rxegaem82jpclb8e40fp019ru&st=ludlvgxt&dl=1')
:::

:::{note}
Due to time constraints, this command is not generated in this notebook. The resulting artifacts can be generated by running this command or running wget the artifact below. 👇
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
                                   'https://www.dropbox.com/scl/fi/9zpv1fg8nsrcaj0v0oarb/taxonomy-maarjam-subsampled.qza?rlkey=rhzuw7mgz94hz85jx8see82q8&st=3t5cjlc3&dl=1')
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
rice_taxonomy_vsearch, rice_search_results_vsearch = use.action(
    use.UsageAction(plugin_id='feature_classifier', action_id='classify_consensus_vsearch'),
    use.UsageInputs(
        query=representative_sequences,
        reference_reads=maarjam_ref_seq,
        reference_taxonomy=ref_taxonomy_maarjam,
        maxaccepts=1,
        perc_identity=0.7,
        strand='both',
        top_hits_only=True,
        unassignable_label='Unassigned'),
    use.UsageOutputNames(
        classification='rice_taxonomy_vsearch',
        search_results='rice_search_results_vsearch'))
:::

:::{describe-usage}
rice_taxonomy_vsearch = use.view_as_metadata('rice_taxonomy_vsearch', rice_taxonomy_vsearch)
use.action(
    use.UsageAction(plugin_id='metadata', action_id='tabulate'),
    use.UsageInputs(
        input=rice_taxonomy_vsearch),
    use.UsageOutputNames(
        visualization='rice_taxonomy_vsearch'))
:::

# Revision Questions

1.  How many of your ASVs were taxonomically assigned by this method. Is there a difference in outpu between fit-classifier-naive-bayes   and classify-consensus-vsearch?

# Taxonomy Bar Plot 
:::{describe-usage}
taxa_bar_plots = use.action(
    use.UsageAction(plugin_id='taxa', action_id='barplot'),
    use.UsageInputs(
        table=table,
        taxonomy=taxonomy_maarjam,
        metadata=metadata),
    use.UsageOutputNames(visualization='taxa_bar_plots'))
:::

# Revision Questions
 1. What are the dominant phyla in each in each group ?

Note : Visualize the samples at Level 6 (which corresponds to the genus of AMF in this analysis), and then sort the samples by  env_broad_level,  You can add as many taxonomic levels levels you want.
If it's hard to visualize the dominant phyla in each in each group, download the csv file on left hand side and use this to plot a relative abundance chart.

# Filtering Tables

You can filter table if you want to work with specific group of taxa. You can create separate files for traditionnal and modern varieties if you want to.

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
It identifies taxa that are differentially abundant in modern or traditional rice variety.

Accurately identifying features that are differentially abundant across sample types in microbiome data is a challenging problem and an open area of research. Analysis of Compositions of Microbiomes with Bias Correction (ANCOM-BC) is a methodology for differential abundance (DA) testing that corrects bias in microbiome data.

Filter table to selectively pick the Glomeromycota phylum features.

:::{describe-usage}
glomeromycetes_table, = use.action(
    use.UsageAction(plugin_id='taxa', action_id='filter_table'),
    use.UsageInputs(
        table=table,
        taxonomy=taxonomy_maarjam,
        include='p__Mucoromycota'),
    use.UsageOutputNames(filtered_table='glomeromycetes_table'))
:::

:::{describe-usage}
collapsed_table_level_6, = use.action(
    use.UsageAction(plugin_id='taxa', action_id='collapse'),
    use.UsageInputs(
        table=glomeromycetes_table,
        taxonomy=taxonomy_maarjam,
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
# Revision Questions

 1. Which taxa are enriched in traditional varieties?
 Note: try changing sample_name to column 'env_broad_scale'; try with significance threshold of 0.05 also.

# Diversity Analysis

## Phylogenetic tree

Both rooted and unrooted trees can be generated by this command.
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

# Alpha rarefaction Plot

:::{describe-usage}
alpha_rarefaction_plot = use.action(
    use.UsageAction(plugin_id='diversity', action_id='alpha_rarefaction'),
    use.UsageInputs(
        table=table,
        phylogeny=rooted_tree,
        max_depth=5400,
        metadata=metadata,
        steps=4,
        iterations=4
        ),
    use.UsageOutputNames(visualization='alpha_rarefaction'))
:::
# Revision Questions

1. How did we decide as maximum depth 13299 ?

Note : The max depth setting will depend on the number of sequences in your samples. The value that you provide for –p-max-depth should be determined by reviewing the OTU Frequency per sample information presented in the table.qzv file that was created above. Also observe that Do the refraction curves for each sample plateau? If they don’t, the samples haven’t been sequenced deeply enough to capture the full diversity of the AM fungal communities, which is shown on the y-axis. At what sequencing depth (x-axis) do your curves plateau? This value will be important for downstream analyses, particularly for alpha diversity analyses. Considering both features and samples retained is very important : lets go back to view our table

:::{describe-usage}
review_table = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='summarize'),
    use.UsageInputs(
        table=table),
    use.UsageOutputNames(visualization='review_table'))
:::

You may want to increase that value if the lines in the resulting rarefaction plot don’t appear to be leveling out, or decrease that value if you seem to be losing many of your samples due to low total frequencies closer to the minimum sampling depth than the maximum sampling depth.
1602
Retained 206,658 (8.75%) features in 129 (90.21%) samples at the specified sampling depth. So we will use 1602 for further analysis.

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

# Revision Questions
1. Is species richness same as number of ASVs ?

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
    use.UsageOutputNames(visualization='unweighted_unifrac_env_broad_scale_significance'))
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
    use.UsageOutputNames(visualization='unweighted_unifrac_sample_name_significance'))
:::
