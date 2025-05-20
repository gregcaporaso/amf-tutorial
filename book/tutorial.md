(tutorial)=
# Tutorial

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

#denoising 
# Check for primers using cutadapt, ?
:::{describe-usage}
trimmed_sequences = use.action(
    use.UsageAction(plugin_id='cutadapt', action_id='trim_paired'),
    use.UsageInputs(
        demultiplexed_sequences=demux,
        cores=1,
        front_f='TCCTACGGGAGGCAGCAGT',
        front_r='GAGTTTCCCGCAGGTTCAC',
        error_rate=0.2,
        ),
    use.UsageOutputNames(trimmed_sequences='trimmed_sequences'))
:::

:::{describe-usage}
visualization_5 = use.action(
    use.UsageAction(plugin_id='demux', action_id='summarize'),
    use.UsageInputs(
        data=trimmed_sequences,
        n=10000),
    use.UsageOutputNames(visualization='visualization_5'))
:::

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
ref_seqs, = use.action(
    use.UsageAction(plugin_id='feature_classifier', action_id='extract_reads'),
    use.UsageInputs(
        sequences=otus_85_maarjam,
        f_primer='AAGCTCGTAGTTGAATTTCG',
        r_primer='CCCAACTATCCCTATTAATCAT',
        trunc_len=250, # mismatched from 120 in provenance replay
        min_length=100,
        max_length=400,
        ),
    use.UsageOutputNames(reads='ref_seqs'))
:::

qiime feature-classifier extract-reads \
  --i-sequences 85_otus.qza \
  --p-f-primer AAGCTCGTAGTTGAATTTCG \
  --p-r-primer CCCAACTATCCCTATTAATCAT \
  --p-trunc-len 250 \ # mismatched
  --p-min-length 100 \
  --p-max-length 400 \
  --o-reads ref-seqs.qza



:::{describe-usage}
classifier, = use.action(
    use.UsageAction(plugin_id='feature_classifier', action_id='fit_classifier_naive_bayes'),
    use.UsageInputs(
        reference_reads=ref_seqs,
        reference_taxonomy=ref_taxonomy,
        ),
    use.UsageOutputNames(classifier='classifier'))
:::

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier.qza

:::{describe-usage}
taxonomy, = use.action(
    use.UsageAction(plugin_id='feature_classifier', action_id='classify_sklearn'),
    use.UsageInputs(
        reads=representative_sequences,
        classifier=classifier,
        confidence=0.7,
        n_jobs=1,
        read_orientation='auto'),
    use.UsageOutputNames(classification='taxonomy'))
:::

qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads representative-sequences.qza \
  --o-classification taxonomy.qza

:::{describe-usage}
taxonomy_viz = use.action(
    use.UsageAction(plugin_id='metadata', action_id='tabulate'),
    use.UsageInputs(
        input=taxonomy),
    use.UsageOutputNames(
        visualization='taxonomy'))
:::

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime tools view taxonomy.qza

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

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences representative-sequences.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# day 8 
:::{describe-usage}
rice_taxonomy_maarjAM, rice_search_results_maarjAM = use.action(
    use.UsageAction(plugin_id='feature_classifier', action_id='classify_consensus_vsearch'),
    use.UsageInputs(
        query=representative_sequences,
        reference_reads=otus_85,
        reference_taxonomy=ref_taxonomy,
        maxaccepts=1,
        perc_identity=0.99,
        strand='both',
        top_hits_only=True,
        unassignable_label='Unassigned'),
    use.UsageOutputNames(
        classification='rice_taxonomy_maarjAM',
        search_results='rice_search_results_maarjAM'))
:::

qiime feature-classifier classify-consensus-vsearch \
  --i-query rep-seq_rice.qza  \
  --i-reference-reads otus_85._maarjaas.qza \
  --i-reference-taxonomy ref-taxonomy_maarjaas.qza \
  --p-perc-identity 0.99 \
  --p-top-hits-only \
  --p-maxaccepts 1 \
  --p-strand 'both' \
  --p-unassignable-label 'Unassigned' \
  --o-classification rice_taxonomy_maarjAM_aligned_rep_seqs.qza \
  --o-search-results rice_search-results_maarjAM_0.99.qza

:::{describe-usage}
table_rice_unassigned_maarjAM = use.action(
    use.UsageAction(plugin_id='taxa', action_id='filter_table'),
    use.UsageInputs(
        table=table, #table_rice?
        taxonomy=rice_taxonomy_maarjAM,
        include='Unassigned',
        query_delimiter=',',
        mode='contains'),
    use.UsageOutputNames(filtered_table='table_rice_unassigned_maarjAM'))
:::

qiime taxa filter-table \
  --i-table table_rice.qza \
  --i-taxonomy rice_taxonomy_maarjaAM_0.99.qza \
  --p-include "Unassigned" \
  --o-filtered-table table-rice-unassigned_maarjAM.qza

:::{describe-usage}
filtered_unassigned_table = use.action(
    use.UsageAction(plugin_id='taxa', action_id='filter_table'),
    use.UsageInputs(
        table=table, #table_rice?
        taxonomy=rice_taxonomy_maarjAM,
        include='c__Glomeromycetes',
        query_delimiter=',',
        mode='contains'),
    use.UsageOutputNames(filtered_table='maarjAM_Glomeromycetes_seqs'))
:::

qiime taxa filter-table \
  --i-table table_rice.qza \
  --i-taxonomy rice_taxonomy_maarjaAM_0.99.qza \
  --p-include c__Glomeromycetes \
  --o-filtered-table maarjAM_Glomeromycetes-seqs.qza

:::{describe-usage}
rep_seq_rice_unassigned_maarjAM = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='filter_seqs'),
    use.UsageInputs(
        data=representative_sequences,
        table=table_rice_unassigned_maarjAM),
    use.UsageOutputNames(filtered_data='rep_seq_rice_unassigned_maarjAM'))
:::

qiime feature-table filter-seqs \
  --i-data rep-seq_rice.qza \
  --i-table table-rice-unassigned_maarjAM.qza \
  --o-filtered-data rep-seq-rice-unassigned_maarjAM.qza

# silva
:::{describe-usage}
def silva_seqs_factory():
    from urllib import request
    from qiime2 import Artifact
    url = 'https://data.qiime2.org/2024.2/common/silva-138-99-seqs.qza'
    fp, _ = request.urlretrieve(url)
    return Artifact.load(fp)

def silva_tax_factory():
    from urllib import request
    from qiime2 import Artifact
    url = 'https://data.qiime2.org/2024.2/common/silva-138-99-tax.qza'
    fp, _ = request.urlretrieve(url)
    return Artifact.load(fp)

silva_138_99_seqs = use.init_artifact('silva_138_99_seqs', silva_seqs_factory)
silva_138_99_tax = use.init_artifact('silva138_99_tax', silva_tax_factory)
:::

:::{describe-usage}
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
:::

qiime feature-classifier classify-consensus-vsearch \
  --i-query rep-seq-rice-unassigned_maarjAM.qza  \
  --i-reference-reads silva-138-99-seqs.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --p-perc-identity 0.99 \
  --p-top-hits-only \
  --p-maxaccepts 1 \
  --p-strand 'both' \
  --p-unassignable-label 'Unassigned' \
  --o-classification rice_taxonomy_assigned_silva_0.99.qza \
  --o-search-results rice_search-results_assigned_silva_0.99.qza

:::{describe-usage}
merged_taxonomy = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='merge_taxa'),
    use.UsageInputs(
        data=[rice_taxonomy_maarjAM, rice_taxonomy_assigned_silva]),
    use.UsageOutputNames(merged_data='merged_taxonomy'))
:::

qiime feature-table merge-taxa
--i-data rice_taxonomy_maarjAM_0.99.qza rice_taxonomy_assigned_silva_0.99.qza \
--o-merged-data merged-taxonomy.qza

:::{describe-usage}
taxa_bar_plots = use.action(
    use.UsageAction(plugin_id='taxa', action_id='barplot'),
    use.UsageInputs(
        table=table,
        taxonomy=merged_taxonomy,
        metadata=metadata),
    use.UsageOutputNames(visualization='taxa_bar_plots'))
:::

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy merged-taxonomy.qza \
  --m-metadata-file metadata_file_rice.tsv \
  --o-visualization taxa-bar-plots.qzv

# Filtering taxa 

#Traditional rice

:::{describe-usage}
traditional_rice_table, = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='filter_samples'),
    use.UsageInputs(
        table=table,
        metadata=metadata,
        where="[env_broad_scale]='Traditional rice root'"),
    use.UsageOutputNames(filtered_table='traditional_rice_table'))
:::

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metadata_file_rice.tsv \
  --p-where "[env_broad_scale]='Traditional rice root'" \
  --o-filtered-table traditional_rice-table.qza


:::{describe-usage}
traditional_rice_table_viz = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='summarize'),
    use.UsageInputs(
        table=traditional_rice_table),
    use.UsageOutputNames(visualization='traditional_rice_table'))
:::

qiime feature-table summarize \
--i-table traditional_rice-table.qza \
--o-visualization traditional_rice-table.qzv 

# TODO qiime tools view traditional_rice-table.qzv

#modern rice

:::{describe-usage}
modern_rice_table, = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='filter_samples'),
    use.UsageInputs(
        table=table,
        metadata=metadata,
        where="[env_broad_scale]='Modern rice root'"),
    use.UsageOutputNames(filtered_table='modern_rice_table'))
:::

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metadata_file_rice.tsv \
  --p-where "[env_broad_scale]='Modern rice root'" \
  --o-filtered-table modern_rice-table.qza

:::{describe-usage}
modern_rice_table_viz = use.action(
    use.UsageAction(plugin_id='feature_table', action_id='summarize'),
    use.UsageInputs(
        table=modern_rice_table),
    use.UsageOutputNames(visualization='modern_rice_table'))
:::

qiime feature-table summarize \
--i-table Modern_rice-table.qza \
--o-visualization Modern_rice-table.qzv


# TODO qiime tools view Modern_rice-table.qzv

ASV relative abundance bar chart

:::{describe-usage}
taxa_bar_plot = use.action(
    use.UsageAction(plugin_id='taxa', action_id='barplot'),
    use.UsageInputs(
        table=table,
        taxonomy=taxonomy,
        metadata=metadata),
    use.UsageOutputNames(visualization='taxa_bar_plots'))
:::

qiime taxa barplot \
    --i-table table.qza \
    --i-taxonomy taxonomy.qza \
    --m-metadata-file metadata_file_rice.tsv \ 
    --o-visualization taxa-bar-plots.qzv



# TODO qiime tools view taxa-bar-plots.qzv

Differential abundance (ANCOM-BC)

:::{describe-usage}
glomeromycetes_table = use.action(
    use.UsageAction(plugin_id='taxa', action_id='filter_table'),
    use.UsageInputs(
        table=table,
        taxonomy=merged_taxonomy,
        include='p__Mucoromycota'),
    use.UsageOutputNames(filtered_table='Glomeromycetes'))
:::

qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy merged-taxonomy.qza \
  --p-include p__Mucoromycota \
  --o-filtered-table Glomeromycetes.qza

#we observed some contamination

:::{describe-usage}
collapsed_table_lvl6 = use.action(
    use.UsageAction(plugin_id='taxa', action_id='collapse'),
    use.UsageInputs(
        table=Glomeromycetes,
        taxonomy=merged_taxonomy,
        level=6),
    use.UsageOutputNames(collapsed_table='collapsed_table_level_6'))
:::

qiime taxa collapse \
  --i-table Glomeromycetes.qza \
  --i-taxonomy merged-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table collapsed-table-level-6.qza


:::{describe-usage}
l6_ancombc_diffs = use.action(
    use.UsageAction(plugin_id='composition', action_id='ancombc'),
    use.UsageInputs(
        table=collapsed_table_level_6,
        metadata=metadata,
        formula='env_broad_scale'),
    use.UsageOutputNames(differentials='l6_ancombc_differentials'))
:::

qiime composition ancombc \
  --i-table collapsed-table-level-6.qza \
  --m-metadata-file metadata_file_rice.tsv \
  --p-formula 'env_broad_scale' \
  --o-differentials l6-ancombc-differentials.qza

:::{describe-usage}
l6_da_barplot = use.action(
    use.UsageAction(plugin_id='composition', action_id='da-barplot'),
    use.UsageInputs(
        data=l6_ancombc_differentials,
        significance_threshold=0.01,
        level_delimiter=';'),
    use.UsageOutputNames(visualization='l6_da_barplot'))
:::

qiime composition da-barplot \
  --i-data l6-ancombc-differentials.qza \
  --p-significance-threshold 0.01 \
  --p-level-delimiter ';' \
  --o-visualization l6-da-barplot.qzv

# TODO Day 9 alpha rarefaction and core metrics

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

qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 13,299 \
  --m-metadata-file metadata_file_rice.tsv \
  --o-visualization alpha-rarefaction.qzv

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

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 1602\
  --m-metadata-file metadata_file_rice.tsv \
  --output-dir core-metrics-results

:::{describe-usage}
shannon_group_significance = use.action(
    use.UsageAction(plugin_id='diversity', action_id='alpha_group_significance'),
    use.UsageInputs(
        alpha_diversity=shannon_vector,
        metadata=metadata),
    use.UsageOutputNames(visualization='shannon_group_significance'))
:::

qiime diversity alpha-group-significance \
  --i-alpha-diversity shannon_vector.qza \
  --m-metadata-file metadata_file_rice.tsv \
  --o-visualization shannon_group_significance

:::{describe-usage}
evenness_group_significance = use.action(
    use.UsageAction(plugin_id='diversity', action_id='alpha_group_significance'),
    use.UsageInputs(
        alpha_diversity=evenness_vector,
        metadata=metadata),
    use.UsageOutputNames(visualization='evenness_group_significance'))
:::

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file metadata_file_rice.tsv \
  --o-visualization evenness-group-significance.qzv

:::{describe-usage}
faith_pd_group_significance = use.action(
    use.UsageAction(plugin_id='diversity', action_id='alpha_group_significance'),
    use.UsageInputs(
        alpha_diversity=faith_pd_vector,
        metadata=metadata),
    use.UsageOutputNames(visualization='visualization_alphasignificance_faith'))
:::

qiime diversity alpha-group-significance \
    --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
    --m-metadata-file metadata_file_rice.tsv \
    --o-visualization core-metrics-results/ visualization_alphasignificance_faith.qzv

:::{describe-usage}
beta_group_significance_simple = use.action(
    use.UsageAction(plugin_id='diversity', action_id='beta_group_significance'),
    use.UsageInputs(
        distance_matrix=unweighted_unifrac_distance_matrix,
        metadata=metadata,
        metadata_column='env_broad_scale',
        method='permanova',
        no_pairwise=True,
        permutations=999),
    use.UsageOutputNames(visualization='unweighted_unifrac_subject_group_significance'))
:::

qiime diversity beta-group-significance \
  --i-distance-matrix unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata_file_rice.tsv \
  --m-metadata-column env_broad_scale \
  --o-visualization unweighted-unifrac-subject-group-significance.qzv

:::{describe-usage}
beta_group_significance = use.action(
    use.UsageAction(plugin_id='diversity', action_id='beta_group_significance'),
    use.UsageInputs(
        distance_matrix=unweighted_unifrac_distance_matrix,
        metadata=metadata,
        metadata_column='Sample_Name',
        method='permanova',
        pairwise=True,
        permutations=999),
    use.UsageOutputNames(visualization='unweighted_unifrac_body_site_significance'))
:::

qiime diversity beta-group-significance \
  --i-distance-matrix unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata_file_rice.tsv \
  --m-metadata-column Sample_Name \
  --o-visualization unweighted-unifrac-body-site-significance.qzv \
  —p-pairwise
