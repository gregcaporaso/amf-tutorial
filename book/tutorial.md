(tutorial)=
# Tutorial

:::{note}
This document was built with its own conda environment.
You can download the environment file that was used from the download link on the top-right of this article.
:::
# TODO day1-3 her tutorial includes importing instructions, are we including that?
# TODO day4 importing data, do we need to include how to create a manifest file?
# TODO day4 add demux summarize and view
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

# TODO import sample-data-paired-end-sequences-with-quality-2.qza

:::{describe-usage}
trimmed_sequences_0 = use.action(
    use.UsageAction(plugin_id='cutadapt', action_id='trim-paired'),
    use.UsageInputs(
        demultiplexed_sequences='sample-data-paired-end-sequences-with-quality-2.qza',
        cores=1,
        front_f='AAGCTCGTAGTTGAATTTCG', # these are from the prov-replay
        front_r='CCCAACTATCCCTATTAATCAT',
        error_rate=0.2,
        indels=True,
        times=1,
        overlap=3,
        no_match_read_wildcards=True,
        match_adapter_wildcards=True,
        minimum_length=1,
        no_discard_untrimmed=True,
        quality_cutoff_5end=0,
        quality_cutoff_3end=0,
        quality_base=33),
    use.UsageOutputNames(trimmed_sequences='trimmed-sequences.qza'))
:::

:::{describe-usage}
visualization_5 = use.action(
    use.UsageAction(plugin_id='demux', action_id='summarize'),
    use.UsageInputs(
        data='trimmed-sequences.qza',
        n=10000),
    use.UsageOutputNames(visualization='visualization-5.qzv'))
:::

# TODO cutadapt which adapters should be used?

qiime cutadapt trim-paired \
--i-demultiplexed-sequences analysis/seqs/combined.qza \
--p-front-f TCCTACGGGAGGCAGCAGT \ # these are from the powerpoint
--p-front-r GAGTTTCCCGCAGGTTCAC \ # see note above
--p-error-rate 0.20 \
--output-dir analysis/seqs_trimmed_glom \
--verbose
# TODO metadata tabulate?

:::{describe-usage}
table, denoising_stats, representative-sequences.qza = use.action(
    use.UsageAction(plugin_id='dada2', action_id='denoise-paired'),
    use.UsageInputs(
        demultiplexed_seqs='sample-data-paired-end-sequences-with-quality.qza',
        trunc_len_f=220, #different from tutorial see below command
        trunc_len_r=200,
        trim_left_f=0,
        trim_left_r=0,
        max_ee_f=2.0,
        max_ee_r=2.0,
        trunc_q=2,
        min_overlap=12,
        pooling_method='independent',
        chimera_method='consensus',
        min_fold_parent_over_abundance=1.0,
        no_allow_one_off=True,
        n_threads=1,
        n_reads_learn=1000000,
        hashed_feature_ids=True,
        retain_all_samples=True),
    use.UsageOutputNames(
        table='table.qza',
        denoising_stats='denoising-stats.qza',
        representative_sequences='representative-sequences.qza'))
:::

# TODO Denoising with dada2, should I remove all the default parameters above?
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs demux.qza \
    --p-trunc-len-f 240 \ #different from provenance
    --p-trunc-len-r 220 \
    --o-representative-sequences representative-sequences.qza \
    --o-table table.qza \
    --o-denoising-stats denoising-stats.qza

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization stats-dada2.qzv

qiime tools view stats-dada2.qzv

qiime feature-table tabulate-seqs \
--i-data representative-sequences.qza \
--o-visualization rep-seqs.qzv

Saved Visualization to: rep-seqs.qzv
Qiime tools view rep-seqs.qzv

qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv 

Saved Visualization to: table.qzv

qiime tools view table.qzv

# TODO Reference databases and Taxonomic classification, feature table filtering and building the phylogenetic tree
# TODO qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path maarjam_database_SSU_TYPE.qiime.fasta \
  --output-path 85_otus.qza
# TODO qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path maarjam_database_SSU_TYPE.qiime.txt \
  --output-path ref-taxonomy.qza

:::{describe-usage}
ref-seqs.qza = use.action(
    use.UsageAction(plugin_id='feature-classifier', action_id='extract-reads'),
    use.UsageInputs(
        sequences='85_otus.qza',
        f_primer='AAGCTCGTAGTTGAATTTCG',
        r_primer='CCCAACTATCCCTATTAATCAT',
        trim_right=0,
        trunc_len=120, # mismatched from below
        trim_left=0,
        identity=0.8,
        min_length=100,
        max_length=400,
        n_jobs=1,
        batch_size='auto',
        read_orientation='both'),
    use.UsageOutputNames(extracted_reads='reads.qza'))
:::

qiime feature-classifier extract-reads \
  --i-sequences 85_otus.qza \
  --p-f-primer AAGCTCGTAGTTGAATTTCG \
  --p-r-primer CCCAACTATCCCTATTAATCAT \
  --p-trunc-len 250 \ # mismatched
  --p-min-length 100 \
  --p-max-length 400 \
  --o-reads ref-seqs.qza
# TODO import ref-taxonomy.qza

:::{describe-usage}
classifier = use.action(
    use.UsageAction(plugin_id='feature-classifier', action_id='fit-classifier-naive-bayes'),
    use.UsageInputs(
        reference_reads='ref-seqs.qza',
        reference_taxonomy='ref-taxonomy.qza.qza',
        classify__alpha=0.001,
        classify__chunk_size=20000,
        classify__class_prior='null',
        no_classify__fit_prior=True,
        no_feat_ext__alternate_sign=True,
        feat_ext__analyzer='char_wb',
        no_feat_ext__binary=True,
        feat_ext__decode_error='strict',
        feat_ext__encoding='utf-8',
        feat_ext__input='content',
        feat_ext__lowercase=True,
        feat_ext__n_features=8192,
        feat_ext__ngram_range='[7, 7]',
        feat_ext__norm='l2',
        feat_ext__preprocessor='null',
        feat_ext__stop_words='null',
        feat_ext__strip_accents='null',
        feat_ext__token_pattern='(?u)\\b\\w\\w+\\b',
        feat_ext__tokenizer='null',
        no_verbose=True),
    use.UsageOutputNames(classifier='classifier.qza'))
:::

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier.qza

:::{describe-usage}
taxonomy = use.action(
    use.UsageAction(plugin_id='feature-classifier', action_id='classify-sklearn'),
    use.UsageInputs(
        reads='representative-sequences.qza',
        classifier='classifier.qza',
        confidence=0.7,
        n_jobs=1,
        read_orientation='auto'),
    use.UsageOutputNames(classification='taxonomy.qza'))
:::

qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads representative-sequences.qza \
  --o-classification taxonomy.qza

:::{describe-usage}
taxonomy_viz = use.action(
    use.UsageAction(plugin_id='qiime2', action_id='metadata tabulate'),
    use.UsageInputs(
        m_input_file='taxonomy.qza'),
    use.UsageOutputNames(
        visualization='taxonomy.qzv'))
:::

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime tools view taxonomy.qza

:::{describe-usage}
rooted_tree, unrooted-tree, aligned-rep-seqs.qza, masked-aligned-rep-seqs.qza = use.action(
    use.UsageAction(plugin_id='phylogeny', action_id='align-to-tree-mafft-fasttree'),
    use.UsageInputs(
        sequences='representative-sequences.qza',
        n_threads=1,
        mask_max_gap_frequency=1.0,
        mask_min_conservation=0.4,
        no_parttree=True),
    use.UsageOutputNames(
        rooted_tree='rooted-tree.qza',
        tree='unrooted-tree.qza',
        alignment='aligned-rep-seqs.qza',
        masked_alignment='masked-aligned-rep-seqs.qza'))
:::

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences representative-sequences.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# day 8 
:::{describe-usage}
rice_taxonomy_maarjAM_0.99, rice_search-results_maarjAM_0.99 = use.action(
    use.UsageAction(plugin_id='feature-classifier', action_id='classify-consensus-vsearch'),
    use.UsageInputs(
        query='rep-seq_rice.qza',
        reference_reads='85_otus._maarjaas.qza',
        reference_taxonomy='ref-taxonomy_maarjaas.qza',
        maxaccepts=1,
        perc_identity=0.99,
        query_cov=0.8,
        strand='both',
        no_search_exact=True,
        top_hits_only=True,
        maxhits='all',
        maxrejects='all',
        output_no_hits=True,
        weak_id=0.0,
        threads=1,
        min_consensus=0.51,
        unassignable_label='Unassigned'),
    use.UsageOutputNames(
        classification='rice_taxonomy_maarjAM_0.99.qza',
        search_results='rice_search-results_maarjAM_0.99.qza'))
:::

qiime feature-classifier classify-consensus-vsearch \
  --i-query rep-seq_rice.qza  \
  --i-reference-reads 85_otus._maarjaas.qza \
  --i-reference-taxonomy ref-taxonomy_maarjaas.qza \
  --p-perc-identity 0.99 \
  --p-top-hits-only \
  --p-maxaccepts 1 \
  --p-strand 'both' \
  --p-unassignable-label 'Unassigned' \
  --o-classification rice_taxonomy_maarjAM_0.99.qza \
  --o-search-results rice_search-results_maarjAM_0.99.qza

::{describe-usage}
table-rice-unassigned_maarjAM = use.action(
    use.UsageAction(plugin_id='taxa', action_id='filter-table'),
    use.UsageInputs(
        table='table_rice.qza',
        taxonomy='rice_taxonomy_maarjaAM_0.99.qza',
        include='Unassigned',
        query_delimiter=',',
        mode='contains'),
    use.UsageOutputNames(filtered_table='table-rice-unassigned_maarjAM.qza '))
:::

qiime taxa filter-table \
  --i-table table_rice.qza \
  --i-taxonomy rice_taxonomy_maarjaAM_0.99.qza \
  --p-include "Unassigned" \
  --o-filtered-table table-rice-unassigned_maarjAM.qza

::{describe-usage}
filtered_unassigned_table = use.action(
    use.UsageAction(plugin_id='taxa', action_id='filter-table'),
    use.UsageInputs(
        table='table_rice.qza',
        taxonomy='rice_taxonomy_maarjaAM_0.99.qza',
        include='c__Glomeromycetes',
        query_delimiter=',',
        mode='contains'),
    use.UsageOutputNames(filtered_table='maarjAM_Glomeromycetes-seqs.qza'))
:::

qiime taxa filter-table \
  --i-table table_rice.qza \
  --i-taxonomy rice_taxonomy_maarjaAM_0.99.qza \
  --p-include c__Glomeromycetes \
  --o-filtered-table maarjAM_Glomeromycetes-seqs.qza

:::{describe-usage}
rep-seq-rice-unassigned_maarjAM = use.action(
    use.UsageAction(plugin_id='feature-table', action_id='filter-seqs'),
    use.UsageInputs(
        data='rep-seq_rice.qza',
        table='table-rice-unassigned_maarjAM.qza'),
    use.UsageOutputNames(filtered_data='rep-seq-rice-unassigned_maarjAM.qza'))
:::

qiime feature-table filter-seqs \
  --i-data rep-seq_rice.qza \
  --i-table table-rice-unassigned_maarjAM.qza \
  --o-filtered-data rep-seq-rice-unassigned_maarjAM.qza

:::{describe-usage}
rice_taxonomy_assigned_silva_0.99, rice_search-results_assigned_silva_0.99 = use.action(
    use.UsageAction(plugin_id='feature-classifier', action_id='classify-consensus-vsearch'),
    use.UsageInputs(
        query='rep-seq-rice-unassigned_maarjAM.qza',
        reference_reads='silva-138-99-seqs.qza',
        reference_taxonomy='silva-138-99-tax.qza',
        perc_identity=0.99,
        top_hits_only=True,
        maxaccepts=1,
        strand='both',
        unassignable_label='Unassigned'),
    use.UsageOutputNames(
        classification='rice_taxonomy_assigned_silva_0.99.qza',
        search_results='rice_search-results_assigned_silva_0.99.qza'))
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
    use.UsageAction(plugin_id='feature-table', action_id='merge-taxa'),
    use.UsageInputs(
        data=['rice_taxonomy_maarjAM_0.99.qza', 'rice_taxonomy_assigned_silva_0.99.qza']),
    use.UsageOutputNames(merged_data='merged-taxonomy.qza'))
:::

qiime feature-table merge-taxa
--i-data rice_taxonomy_maarjAM_0.99.qza rice_taxonomy_assigned_silva_0.99.qza \
--o-merged-data merged-taxonomy.qza

:::{describe-usage}
axa-bar-plots = use.action(
    use.UsageAction(plugin_id='taxa', action_id='barplot'),
    use.UsageInputs(
        table='table.qza',
        taxonomy='merged-taxonomy.qza',
        metadata_file='metadata_file_rice.tsv'),
    use.UsageOutputNames(visualization='taxa-bar-plots.qzv'))
:::

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy merged-taxonomy.qza \
  --m-metadata-file metadata_file_rice.tsv \
  --o-visualization taxa-bar-plots.qzv

# Filtering taxa 

#Traditional rice

:::{describe-usage}
traditional_rice_table = use.action(
    use.UsageAction(plugin_id='feature-table', action_id='filter-samples'),
    use.UsageInputs(
        table='table.qza',
        metadata_file='metadata_file_rice.tsv',
        where="[env_broad_scale]='Traditional rice root'"),
    use.UsageOutputNames(filtered_table='traditional_rice-table.qza'))
:::

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metadata_file_rice.tsv \
  --p-where "[env_broad_scale]='Traditional rice root'" \
  --o-filtered-table traditional_rice-table.qza


:::{describe-usage}
traditional_rice_table_viz = use.action(
    use.UsageAction(plugin_id='feature-table', action_id='summarize'),
    use.UsageInputs(
        table='traditional_rice-table.qza'),
    use.UsageOutputNames(visualization='traditional_rice-table.qzv'))
:::

qiime feature-table summarize \
--i-table traditional_rice-table.qza \
--o-visualization traditional_rice-table.qzv 

# TODO qiime tools view traditional_rice-table.qzv

#modern rice

:::{describe-usage}
filtered_table_1 = use.action(
    use.UsageAction(plugin_id='feature-table', action_id='filter-samples'),
    use.UsageInputs(
        table='table.qza',
        metadata_file='metadata_file_rice.tsv',
        where="[env_broad_scale]='Modern rice root'",
    use.UsageOutputNames(filtered_table='Modern_rice-table.qza'))
:::

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metadata_file_rice.tsv \
  --p-where "[env_broad_scale]='Modern rice root'" \
  --o-filtered-table Modern_rice-table.qza

:::{describe-usage}
modern_rice_table_viz = use.action(
    use.UsageAction(plugin_id='feature-table', action_id='summarize'),
    use.UsageInputs(
        table='Modern_rice-table.qza'),
    use.UsageOutputNames(visualization='Modern_rice-table.qzv.qzv'))
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
        table='table.qza',
        taxonomy='taxonomy.qza',
        metadata_file='metadata_file_rice.tsv'),
    use.UsageOutputNames(visualization='taxa-bar-plots.qzv'))
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
    use.UsageAction(plugin_id='taxa', action_id='filter-table'),
    use.UsageInputs(
        table='table.qza',
        taxonomy='merged-taxonomy.qza',
        include='p__Mucoromycota'),
    use.UsageOutputNames(filtered_table='Glomeromycetes.qza'))
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
        table='Glomeromycetes.qza',
        taxonomy='merged-taxonomy.qza',
        level=6),
    use.UsageOutputNames(collapsed_table='collapsed-table-level-6.qza'))
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
        table='collapsed-table-level-6.qza',
        metadata_file='metadata_file_rice.tsv',
        formula='env_broad_scale'),
    use.UsageOutputNames(differentials='l6-ancombc-differentials.qza'))
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
        data='l6-ancombc-differentials.qza',
        significance_threshold=0.01,
        level_delimiter=';'),
    use.UsageOutputNames(visualization='l6-da-barplot.qzv'))
:::

qiime composition da-barplot \
  --i-data l6-ancombc-differentials.qza \
  --p-significance-threshold 0.01 \
  --p-level-delimiter ';' \
  --o-visualization l6-da-barplot.qzv

# TODO Day 9 alpha rarefaction and core metrics

:::{describe-usage}
alpha_rarefaction_plot = use.action(
    use.UsageAction(plugin_id='diversity', action_id='alpha-rarefaction'),
    use.UsageInputs(
        table='table.qza',
        phylogeny='rooted-tree.qza',
        max_depth=13299,
        metadata_file='metadata_file_rice.tsv',
        min_depth=1,
        steps=10,
        iterations=10),
    use.UsageOutputNames(visualization='alpha-rarefaction.qzv'))
:::

qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 13,299 \
  --m-metadata-file metadata_file_rice.tsv \
  --o-visualization alpha-rarefaction.qzv

:::{describe-usage}
core_metrics_outputs = use.action(
    use.UsageAction(plugin_id='diversity', action_id='core-metrics-phylogenetic'),
    use.UsageInputs(
        table='table.qza',
        phylogeny='rooted-tree.qza',
        sampling_depth=1602,
        metadata_file='metadata_file_rice.tsv',
        no_with_replacement=True,
        n_jobs_or_threads=1,
        no_ignore_missing_samples=True),
    use.UsageOutputNames(
        evenness_vector='evenness-vector.qza',
        faith_pd_vector='faith-pd-vector.qza',
        unweighted_unifrac_distance_matrix='unweighted-unifrac-distance-matrix.qza',
        bray_curtis_pcoa_results='bray-curtis-pcoa-results.qza',
        shannon_vector='shannon-vector.qza',
        rarefied_table='rarefied-table.qza',
        weighted_unifrac_distance_matrix='weighted-unifrac-distance-matrix.qza',
        jaccard_pcoa_results='jaccard-pcoa-results.qza',
        unweighted_unifrac_emperor='unweighted-unifrac-emperor.qzv',
        weighted_unifrac_pcoa_results='weighted-unifrac-pcoa-results.qza',
        observed_features_vector='observed-features-vector.qza',
        jaccard_distance_matrix='jaccard-distance-matrix.qza',
        jaccard_emperor='jaccard-emperor.qzv',
        bray_curtis_emperor='bray-curtis-emperor.qzv',
        weighted_unifrac_emperor='weighted-unifrac-emperor.qzv',
        bray_curtis_distance_matrix='bray-curtis-distance-matrix.qza',
        unweighted_unifrac_pcoa_results='unweighted-unifrac-pcoa-results.qza'))
:::

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 1602\
  --m-metadata-file metadata_file_rice.tsv \
  --output-dir core-metrics-results

:::{describe-usage}
shannon_group_significance = use.action(
    use.UsageAction(plugin_id='diversity', action_id='alpha-group-significance'),
    use.UsageInputs(
        alpha_diversity='shannon_vector.qza',
        metadata_file='metadata_file_rice.tsv'),
    use.UsageOutputNames(visualization='shannon_group_significance.qzv'))
:::

qiime diversity alpha-group-significance \
  --i-alpha-diversity shannon_vector.qza \
  --m-metadata-file metadata_file_rice.tsv \
  --o-visualization shannon_group_significance.qzv

:::{describe-usage}
evenness_group_significance = use.action(
    use.UsageAction(plugin_id='diversity', action_id='alpha-group-significance'),
    use.UsageInputs(
        alpha_diversity='evenness-vector.qza',
        metadata_file='metadata_file_rice.tsv'),
    use.UsageOutputNames(visualization='evenness-group-significance.qzv'))
:::

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file metadata_file_rice.tsv \
  --o-visualization evenness-group-significance.qzv

:::{describe-usage}
faith_pd_group_significance = use.action(
    use.UsageAction(plugin_id='diversity', action_id='alpha-group-significance'),
    use.UsageInputs(
        alpha_diversity='faith-pd-vector.qza',
        metadata_file='metadata_file_rice.tsv'),
    use.UsageOutputNames(visualization='visualization_alphasignificance_faith.qzv'))
:::

qiime diversity alpha-group-significance \
    --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
    --m-metadata-file metadata_file_rice.tsv \
    --o-visualization core-metrics-results/ visualization_alphasignificance_faith.qzv

:::{describe-usage}
beta_group_significance_simple = use.action(
    use.UsageAction(plugin_id='diversity', action_id='beta-group-significance'),
    use.UsageInputs(
        distance_matrix='unweighted-unifrac-distance-matrix.qza',
        metadata_file='metadata_file_rice.tsv',
        metadata_column='env_broad_scale',
        method='permanova',
        no_pairwise=True,
        permutations=999),
    use.UsageOutputNames(visualization='unweighted-unifrac-subject-group-significance.qzv'))
:::

qiime diversity beta-group-significance \
  --i-distance-matrix unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata_file_rice.tsv \
  --m-metadata-column env_broad_scale \
  --o-visualization unweighted-unifrac-subject-group-significance.qzv

:::{describe-usage}
beta_group_significance = use.action(
    use.UsageAction(plugin_id='diversity', action_id='beta-group-significance'),
    use.UsageInputs(
        distance_matrix='unweighted-unifrac-distance-matrix.qza',
        metadata_file='metadata_file_rice.tsv',
        metadata_column='Sample_Name',
        method='permanova',
        pairwise=True,
        permutations=999),
    use.UsageOutputNames(visualization='unweighted-unifrac-body-site-significance.qzv'))
:::

qiime diversity beta-group-significance \
  --i-distance-matrix unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata_file_rice.tsv \
  --m-metadata-column Sample_Name \
  --o-visualization unweighted-unifrac-body-site-significance.qzv \
  —p-pairwise
