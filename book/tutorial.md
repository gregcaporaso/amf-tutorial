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
