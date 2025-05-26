# Neubauer_et_al-metataxonomic-cow-milk-sara-lps

Repository with the code for reproducing the analysis of the manuscript: Viktoria Neubauer, Siska Aditya, Narciso M. Quijada, Stefanie Urimare Wetzels, Monika Dzieciol, Poulad Pourazad, Qendrim Zebeli and Evelyne Selberherr, *Dynamics of the milk microbial community during subacute ruminal acidosis with or without intramammary lipopolysaccharide challenge in dairy cows* (submitted)

## Index<a name="id0"></a> 
1. [Set up enviornment and variables](#id1)
2. [Quality control, PE-merging chimera removal and ASV table resolution](#id2)
3. [Taxonomic assignment of ASVs](#id3)
4. [Removal of potentially contaminant ASVs](#id4)
5. [Normalization to gene copy number](#id5)
6. [Phylogenetic tree](#id6)
7. [Alpha-rarefaction analysis](#id7)
8. [Alpha-diversity analysis](#id8)
9. [Beta-diversity analysis](#id9)

<br>

## 1. Set up enviornment and variables<a name="id1"></a>

```bash
conda activate qiime2-2023.9
ADAPTERS= #path to the seq-adapters.fasta file contained in "files" directory from this repository
METADATA= #path to the metadata file provided in this repository
MYWD= #define working directory that contains the directory 00-raw-reads with the raw sequencing data
SILVAD= #path to the SILVA database, for the paper: silva-138-99-nb-classifier.qza
```

<br>

[Back to Index](#id0)

<br>

## 2. Quality control, PE-merging chimera removal, ASV table resolution and decontamination<a name="id2"></a>

### 2.1 Remove adapters and inspect quality
```bash
mkdir -p ${MYWD}/00-raw-reads-no-adapt/fastqc
mkdir ${MYWD}/00-raw-reads/fastqc
for SAMPLE in $(ls ${MYWD}/00-raw-reads/*R1_001* | sed "s#${MYWD}/00-raw-reads##" | sed "s/_R1_001.fastq.gz//"); do
  trimmomatic PE -threads 32 -phred33 ${MYWD}/00-raw-reads/${SAMPLE}_R1_001.fastq.gz ${MYWD}/00-raw-reads/${SAMPLE}_R2_001.fastq.gz ${MYWD}/00-raw-reads-no-adapt/${SAMPLE}_L001_R1_001.fastq.gz /dev/null ${MYWD}/00-raw-reads-no-adapt/${SAMPLE}_L001_R2_001.fastq.gz /dev/null ILLUMINACLIP:${ADAPTERS}:1:30:11
  fastqc ${MYWD}/00-raw-reads/${SAMPLE}* -o ${MYWD}/00-raw-reads/fastqc
  fastqc ${MYWD}/${MYWD}/00-raw-reads-no-adapt/${SAMPLE}* -o ${MYWD}/00-raw-reads-no-adapt/fastqc
done
```

### 2.2 Import data to QIIME2

```bash
mkdir ${MYWD}/01-qiime-raw
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ${MYWD}/00-raw-reads-no-adapt --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path ${MYWD}/01-qiime-raw/PE-reads.qza
qiime demux summarize --i-data ${MYWD}/01-qiime-raw/PE-reads.qza --o-visualization ${MYWD}/01-qiime-raw/PE-reads.qzv
```

### 2.3 Run dada2 for QC, PE-merging, chimera removal and ASV resolution
```bash
mkdir ${MYWD}/02-ASV-table
# tried different combinations with different parameters and reported the best
qiime dada2 denoise-paired --i-demultiplexed-seqs ${MYWD}/01-qiime-raw/PE-reads.qza --p-trim-left-f 13 --p-trim-left-r 13 --p-trunc-len-f 247 --p-trunc-len-r 240 --p-max-ee-f 2 --p-max-ee-r 2 --o-representative-sequences ${MYWD}/02-ASV-table/rep-seqs.qza --o-table ${MYWD}/02-ASV-table/table.raw.qza --o-denoising-stats ${MYWD}/02-ASV-table/stats-dada2.qza --p-n-threads 24

# inspect visualization and/or exported files
qiime metadata tabulate --m-input-file ${MYWD}/02-ASV-table/stats-dada2.qza --o-visualization ${MYWD}/02-ASV-table/stats-dada2.qzv
qiime tools export --input-path ${MYWD}/02-ASV-table/stats-dada2.qza --output-path ${MYWD}/02-ASV-table/stats-dada2
qiime feature-table summarize --i-table ${MYWD}/02-ASV-table/table.raw.qza --o-visualization ${MYWD}/02-ASV-table/table.raw.qzv --m-sample-metadata-file ${METDATA}
qiime feature-table tabulate-seqs --i-data ${MYWD}/02-ASV-table/rep-seqs.qza --o-visualization ${MYWD}/02-ASV-table/rep-seqs.qzv
```

### 2.4 Export ASV-table `.qza` to `.txt` for decontamination
```bash

qiime tools export --input-path ${MYWD}/02-ASV-table/table.raw.qza --output-path ${MYWD}/02-ASV-table/exported-table
biom convert -i ${MYWD}/02-ASV-table/exported-table/feature-table.biom -o ${MYWD}/02-ASV-table/exported-table/feature-table.txt --to-tsv
```

<br>

[Back to Index](#id0)

<br>

## 3. Taxonomic assignment of ASVs<a name="id3"></a>

```bash
${MYWD}/03-taxonomy
qiime feature-classifier classify-sklearn --i-classifier ${SILVADB} --i-reads ${MYWD}/02-ASV-table/rep-seqs.qza --o-classification ${MYWD}/03-taxonomy/taxonomy.raw.qza --p-n-jobs 64
```

Removal of taxonomic "unassigned" ASVs or those assigned to:
- Chloroplasts
- Eukaryota
- Mitochondria

```bash
qiime taxa filter-table --i-table ${MYWD}/02-ASV-table/table.raw.qza --i-taxonomy${MYWD}/ 03-taxonomy/taxonomy.raw.qza --p-exclude mitochondria,chloroplast,eukaryota,unassigned --o-filtered-table ${MYWD}/02-ASV-table/table.taxfilt.qza
```

<br>

[Back to Index](#id0)

<br>

## 4. Removal of potentially contaminant ASVs<a name="id4"></a>

[Decontam](https://github.com/benjjneb/decontam) R pipeline was followed, while the main code is in the directory XXX of this repo

Once decontamination is done, import decontaminated ASV-table to QIIME2 
```bash
biom convert -i ${MYWD}/02-ASV-table/exported-table-taxfilt-decontam/feature-table.txt -o ${MYWD}/02-ASV-table/exported-table-taxfilt-decontam/feature-table.biom --to-hdf5
qiime tools import --input-path ${MYWD}/02-ASV-table/exported-table-taxfilt-decontam/feature-table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path ${MYWD}/02-ASV-table/table.taxfilt.decontam.qza
```

Removal of samples belonging to controls

```bash
qiime feature-table filter-samples --i-table ${MYWD}/02-ASV-table/table.taxfilt.decontam.qza --m-metadata-file ${METDATA} --p-where "[FeedingPhase]='NK'" --o-filtered-table ${MYWD}/02-ASV-table/table.taxfilt.decontam.noControls.qza
qiime feature-table summarize --i-table ${MYWD}/02-ASV-table/table.taxfilt.noControls.qza --o-visualization ${MYWD}/02-ASV-table/table.taxfilt.noControls.qzv
```

<br>

[Back to Index](#id0)

<br>

## 5. Normalization to gene copy number<a name="id5"></a>

```bash
qiime gcn-norm copy-num-normalize --i-table ${MYWD}/02-ASV-table/table.taxfilt.noControls.qza --i-taxonomy ${MYWD}/03-taxonomy/taxonomy.raw.qza --o-gcn-norm-table ${MYWD}/02-ASV-table/table.taxfilt.noControls.gcn_norm.qza
```

Remove ASVs with less than 10 counts overall and with a prevalence < 2.

```bash
qiime feature-table filter-features --i-table ${MYWD}/02-ASV-table/table.taxfilt.noControls.gcn_norm.qza --p-min-frequency 10 --o-filtered-table ${MYWD}/02-ASV-table/table.taxfilt.noControls.gcn_norm.ASV_min_10count.qza
#- Remove ASVs only present in one sample
qiime feature-table filter-features --i-table ${MYWD}/02-ASV-table/table.taxfilt.noControls.gcn_norm.ASV_min_10count.qza --p-min-samples 2 --o-filtered-table ${MYWD}/02-ASV-table/table.taxfilt.noControls.gcn_norm.ASV_min_10count.min_prev2.qza
```

This would be final ASV table for donwstream analysis, so you can create a symbolic link to a more "definitive" name for simplify coding, or create a variable for it:

```
ln -s ${MYWD}/02-ASV-table/table.taxfilt.noControls.gcn_norm.ASV_min_10count.min_prev2.qza
Q2TABLE=${MYWD}/02-ASV-table/table.taxfilt.noControls.gcn_norm.ASV_min_10count.min_prev2.qza
```

<br>

[Back to Index](#id0)

<br>

## 6. Phylogenetic tree<a name="id6"></a>

```bash
mkdir ${MYWD}/04-phylogenetic-tree-mafft
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences ${MYWD}/02-ASV-table/rep-seqs.raw.qza --o-alignment ${MYWD}/04-phylogenetic-tree-mafft/aligned-rep-seqs.qza --o-masked-alignment ${MYWD}/04-phylogenetic-tree-mafft/masked-aligned-rep-seqs.qza --o-tree ${MYWD}/04-phylogenetic-tree-mafft/unrooted-tree.qza --o-rooted-tree ${MYWD}/04-phylogenetic-tree-mafft/rooted-tree.qza --p-n-threads 32
```

<br>

[Back to Index](#id0)

<br>

## 7. Alpha-rarefaction analysis<a name="id7"></a>

```bash
mkdir ${MYWD}/05-alpha-rarefaction/
qiime diversity alpha-rarefaction --i-table ${Q2TABLE} --i-phylogeny ${MYWD}/04-phylogenetic-tree-mafft/rooted-tree.qza --p-max-depth 25000 --p-steps 30 --m-metadata-file ${METADATA} --o-visualization ${MYWD}/05-alpha-rarefaction/alpha-rarefaction.qzv --p-metrics chao1 --p-metrics goods_coverage --p-metrics observed_features --p-metrics shannon --p-metrics simpson --p-metrics simpson_e
```
> Adjust `--p-max-depth` to properly see the slope and not too much of the plateau
> The higher `--p-steps` the more resolution and the longer the running time

<br>

Based on the information of the alpha rarefaction curves and the sequencing depth per sample table (`${MYWD}/02-ASV-table/table.raw.qzv`), we define the sequencing depth for rarefaction and beta-diversity analysis:

```
depth=4296
```

<br>

Make rarefaction ASV table:

```bash
qiime feature-table rarefy --i-table ${Q2TABLE} --p-sampling-depth ${depth} --o-rarefied-table ${MYWD}/02-ASV-table/table.taxfilt.noControls.gcn_norm.ASV_min_10count.min_prev2.raref_${depth}.qza
Q2RAREF=${MYWD}/02-ASV-table/table.taxfilt.noControls.gcn_norm.ASV_min_10count.min_prev2.raref_${depth}.qza
```

<br>

[Back to Index](#id0)

<br>

## 8. Alpha-diversity analysis<a name="id8"></a>

### 8.1. Alpha-diversity metrics

```bash
mkdir ${MYWD}/06-alpha-metrics/ 06-alpha-metrics-raref-${depth}/
for alpha in chao1 shannon simpson observed_features goods_coverage simpson_e dominance singles; do
  # without rarefaction
	qiime diversity alpha --i-table ${Q2TABLE} --p-metric ${alpha} --o-alpha-diversity ${MYWD}/06-alpha-metrics/${alpha}-metric.qza
	qiime metadata tabulate --m-input-file ${MYWD}/06-alpha-metrics/${alpha}-metric.qza --o-visualization ${MYWD}/06-alpha-metrics/${alpha}-metric.qzv
  qiime tools export --input-path ${MYWD}/06-alpha-metrics/${alpha}-metric.qza --output-path ${MYWD}/06-alpha-metrics/exported-${alpha}
  # with rarefaction
  qiime diversity alpha --i-table ${Q2RAREF} --p-metric ${alpha} --o-alpha-diversity ${MYWD}/06-alpha-metrics-raref-${depth}/${alpha}-metric.qza
	qiime metadata tabulate --m-input-file ${MYWD}/06-alpha-metrics-raref-${depth}/${alpha}-metric.qza --o-visualization ${MYWD}/06-alpha-metrics-raref-${depth}/${alpha}-metric.qzv
  qiime tools export --input-path ${MYWD}/06-alpha-metrics-raref-${depth}/${alpha}-metric.qza --output-path ${MYWD}/06-alpha-metrics-raref-${depth}/exported-${alpha}
done
```
<br>

### 8.2 Interactive alpha-diversity taxonomic barplots

```bash
mkdir ${MYWD}/07-taxa-barplots/
qiime taxa barplot --i-table ${Q2TABLE} --i-taxonomy ${MYWD}/07-taxonomy/taxonomy.raw.qza --m-metadata-file ${METADATA} --o-visualization ${MYWD}/07-taxa-barplots/taxa-bar-plots.gcn_norm.ASV_min_10count_min_prev2.taxfilt.qzv
qiime taxa barplot --i-table ${Q2RAREF} --i-taxonomy ${MYWD}/07-taxonomy/taxonomy.raw.qza --m-metadata-file ${METADATA} --o-visualization ${MYWD}/07-taxa-barplots/taxa-bar-plots.gcn_norm.ASV_min_10count_min_prev2.taxfilt.raref-${depth}.qzv
```

<br>

[Back to Index](#id0)

<br>

## 9. Beta-diversity analysis<a name="id9"></a>

```bash
qiime diversity core-metrics-phylogenetic --i-phylogeny ${MYWD}/03-phylogenetic-tree-mafft/rooted-tree.qza --i-table ${Q2TABLE} --p-sampling-depth ${depth} --m-metadata-file ${METADATA} --output-dir 08-beta-diversity-GCN-depth-${depth}
for i in $(ls 08-beta-diversity-GCN-depth-${depth}/*qza); do
  mv ${i} $(echo ${i} | sed "s/.qza/_d${depth}.qza/")
done
for i in $(ls 08-beta-diversity-GCN-depth-${depth}/*qzv); do
  mv ${i} $(echo ${i} | sed "s/.qzv/_d${depth}.qzv/")
done
```
<br>

[Back to Index](#id0)

<br>
