# nanoloop: Identify R-loops with Nanopore Reads

- [Introduction](#introduction)
- [Installation](#installation)
- [Parameters](#parameters)
- [Core Functions](#core-functions)
  - [bam_to_tsv](#bam_to_tsv)
  - [tsv_to_plot](#tsv_to_plot)
  - [tsv_to_bed](#tsv_to_bed)
  - [tsv_to_peak](#tsv_to_peak)
  - [bam_to_json](#bam_to_json)
  - [filter_json](#filter_json)
  - [json_to_hotspot](#json_to_hotspot)
  - [stat_hotspot](#stat_hotspot)
  - [cluster_hotspot](#cluster_hotspot)
- [Examples](#examples)
  - [Without PCR](#without-pcr-reference-centric)
    - [Using base quality](#using-base-quality---type-nt_qual)
    - [Using base count](#using-base-count---type-nt_count)
  - [With PCR](#with-pcr-read-centric)
- [Bug reports](#bug-reports)

## Introduction
The `nanoloop` library enables the identification of R-loop regions using nanopore sequencing data.

**R-loops** are three-stranded nucleic acid structures formed when nascent RNA hybridizes with its DNA template, displacing the non-template DNA strand. They play critical roles in transcription regulation, DNA replication, and repair. Dysregulated R-loops are associated with genomic instability and disease. Accurately identifying the genomic locations of R-loop regions provides a critical foundation for subsequent mechanistic investigations. The combination of long-read Nanopore sequencing and A3A treatment enables precise mapping of R-loops. In A3A-treated samples, unprotected single-stranded DNA in R-loops undergoes cytosine (C) deamination to uracil (dU). There are at least two strategies for identifying R-loops, depending on whether the protocol includes an additional PCR step.

### Without PCR (reference-centric)
In the absence of PCR amplification, uracil (dU) residues remain unaltered. When basecalling is performed using standard Dorado models, these dU sites are frequently miscalled as C or T with reduced basecalling quality scores. Additionally, the basecalling qualities of bases neighboring dU sites are often reduced as well. We leverage this signature and implement a reference-centric approach in `nanoloop`to:

+ Quantify and visualize basecalling quality across specified genomic regions

+ Quantify and visualize C-to-T conversion frequencies within those regions

+ Generate MACS3-compatible BED files for downstream peak calling

+ Identify R-loop peaks using a rolling average-based detection algorithm

The reference-centric workflow operates as follows: ![Without PCR](docs/images/without_pcr.svg)

This strategy is termed reference-centric because `nanoloop` scans each position along the reference genome using a sliding window approach. It identifies candidate R-loop regions (peaks) where the window contains an abnormally high fraction of bases with reduced basecalling quality (`--type nt_qual`) or elevated C-to-T conversion rates (`--type nt_count`).

In addition, `nanoloop` can simulate a BED file in which reference positions with a higher fraction of low-quality basecalls or more C-to-T miscalls receive greater coverage. This BED file can then be used as input to popular peak callers like `macs3` for peak calling.

For detailed usage, see the [Example](#without_pcr) section.

### With PCR (read-centric)
Adding a PCR step converts dU into T, effectively creating C-to-T mutations in Nanopore sequencing reads. At the same time, additional mutation patterns can be detected, enabling co-identification of other chromatin features alongside R-loops. 

In addition to the reference-centric strategy, `nanoloop` implements a read-centric approach that locates candidate R-loop regions by identifying mutation hotspots within individual reads. This strategy also uses a sliding-window scan: regions with unusually high mutation density in each read are flagged as hotspots and can be passed directly to `macs3` for peak calling. Main functionalities include:

+ Convert BAM files into JSON format containing relevant mutation information for downstream analyses

+ Filter JSON data based on various metrics

+ Detect mutation hotspots and perform quality control using the JSON data

The strategy operates as follows: ![With PCR](docs/images/with_pcr.svg)

Both the reference-centric and read-centric workflows include built-in plotting functions for quality control. For detailed usage, see the Examples section.

## Installation
```
pip install nanoloop
```

## Parameters
```bash
nanoloop -h
```
```
usage: nanoloop [-h] {bam_to_tsv,tsv_to_plot,tsv_to_bed,tsv_to_peak} ...

nanoloop

positional arguments:
  {bam_to_tsv,tsv_to_plot,tsv_to_bed,tsv_to_peak}
    bam_to_tsv          Parse BAM file and output TSV file
    tsv_to_plot         Parse TSV file and output plot
    tsv_to_bed          Convert TSV file to BED format for MACS3 peak calling. For "--type nt_qual": regions
                        with lower quality scores will generate more simulated read tags, creating peaks in
                        those regions. For "--type nt_count": the number of simulated read tags is proportional
                        to the fraction of C converted to T at each reference C position.
    tsv_to_peak         Call peaks using a sliding window approach

options:
  -h, --help            show this help message and exit
```

You can also run subcommand-specific help, e.g.`nanoloop bam_to_tsv -h`, etc.

## Core Functions

For a complete list of supported options, run `nanoloop <function_name> -h`, replacing <function_name> with the name of the specific subcommand you're interested in.

### `bam_to_tsv`
Converts a BAM file to a bgzipped TSV file. The output reports per-position statistics: base count distribution (if `--type nt_count`) or base quality distribution (if `--type nt_qual`).

### `tsv_to_plot`
Generates plots from TSV files for a given genomic region:

+ `--type nt_qual`: uses TSVs from `nanoloop bam_to_tsv --type nt_qual` and generates plot depicting base quality distribution
+ `--type nt_count`: uses TSVs from `nanoloop bam_to_tsv --type nt_count` and generates plot depicting base count distribution

Note that when `--type nt_count`, the TSV input must be from `bam_to_tsv` with `--type nt_count`; and when `--type nt_qual`, the TSV input must be from `bam_to_tsv` with `--type nt_qual`.

### `tsv_to_bed`
Simulates a MACS3-compatible BED file:

+ `--type nt_qual`: uses TSVs from `nanoloop bam_to_tsv --type nt_qual`, lower average quality yields more tags
+ `--type nt_count`: uses TSVs from `nanoloop bam_to_tsv --type nt_count`, more tags where C-to-T conversion is higher

Note that when `--type nt_count`, the TSV input must be from `bam_to_tsv` with `--type nt_count`; and when `--type nt_qual`, the TSV input must be from `bam_to_tsv` with `--type nt_qual`.

### `tsv_to_peak`
Calls R-loop peaks using a rolling average approach. 

+ Supports both `--type nt_qual` and `--type nt_count`

Note that when `--type nt_count`, the TSV input must be from `bam_to_tsv` with `--type nt_count`; and when `--type nt_qual`, the TSV input must be from `bam_to_tsv` with `--type nt_qual`.

### `bam_to_json`
Parses a BAM file and outputs a JSON file containing per-read information, including read_id, ref_chr, ref_start, ref_end, ref_seq, and detailed mutation data (type, position, and base quality, etc.).

### `filter_json`
Filters a JSON file based on either the raw mutation count or the fraction of mutated bases per read.

### `json_to_hotspot`
Scans each read in the JSON file to detect mutation hotspots using a sliding window approach. Use the `--mutation_type` flag to specify which mutation types to consider.

### `stat_hotspot`
Generates quality control (QC) plots at the hotspot level, including:

+ Hotspot length distribution

+ Number of hotspots per read

+ Number of distinct mutation types per hotspot

### `cluster_hotspot`
Creates two types of clustered heatmaps for QC:

+ Mutation counts per hotspot (clustered by hotspot)

+ Hotspot regions per read (clustered by read)

## Examples

### Without PCR (reference-centric)

The following examples use a downsampled BAM `examples/bam/p1214_no_pcr.bam`, derived from an A3A-treated plasmid. This sample is expected to contain R-loops around its middle region. The filename suffix "no_pcr" indicates that sequencing was performed directly after A3A treatment (without PCR), preserving dU signatures.

`nanoloop` supports both `--type nt_qual` and `--type nt_count` options. Their usage is demosntrated below.

  + `--type nt_qual`: uses the average read quality score at each nucleotide position for peak calling. Regions exhibiting decreased quality scores are identified as potential R-loop candidates.
  + `--type nt_count`: uses the C-to-T conversion frequency at each position for peak calling. Regions with elevated conversion rates are identified as potential R-loop candidates.

#### Using base quality: `--type nt_qual`
This approach leverages the decreased basecalling quality observed at R-loop regions, an artifact caused by the presence of the unconventional deoxyuridien (dU) base.

##### Step 1: generate TSV
```bash
nanoloop bam_to_tsv \
  --bam examples/bam/p1214_no_pcr.bam \
  --ref examples/ref/p1214.fa \
  --type nt_qual \
  --output examples/res/p1214_no_pcr_nt_qual.tsv.gz

zcat < examples/res/p1214_no_pcr_nt_qual.tsv.gz | head
```
```
#chr    start   end     ref_nt  qual_0_10       qual_10_20      qual_20_30      qual_30_40      qual_40_above   qual_avg
p1214   0       1       A       7       275     287     234     166     27.856553147574818
p1214   1       2       A       11      337     307     269     185     27.488728584310188
p1214   2       3       A       13      392     332     283     178     26.808848080133554
p1214   3       4       A       16      408     411     331     173     26.53846153846154
p1214   4       5       A       19      441     523     336     166     26.004040404040403
p1214   5       6       A       41      450     642     373     118     25.34975369458128
p1214   6       7       A       58      441     711     414     89      24.907180385288967
p1214   7       8       A       49      387     702     492     94      25.59570765661253
p1214   8       9       A       43      389     583     546     197     27.357792946530147
```

##### Step 2: visualize quality distribution
```bash
nanoloop tsv_to_plot \
  --tsv examples/res/p1214_no_pcr_nt_qual.tsv.gz \
  --type nt_qual \
  --range p1214:0-10000 \
  --add_qual_avg true \
  --output examples/res/p1214_no_pcr_nt_qual.jpg
```

Plot shows a base quality drop near 5000–6200, suggesting an R-loop

![Quality Distribution Plot](examples/res/p1214_no_pcr_nt_qual.jpg)

##### Step 3: call peaks
Use `nanoloop tsv_to_peak` to call peaks (potential R-loops) using a rolling average approach:
```bash
nanoloop tsv_to_peak \
  --tsv examples/res/p1214_no_pcr_nt_qual.tsv.gz \
  --type nt_qual \
  --output examples/res/p1214_no_pcr_nt_qual_peak.bed.gz
```

The resulting BED file successfully captures the potential R-loop regions:
```bash
zcat < examples/res/p1214_no_pcr_nt_qual_peak.bed.gz | head
```
```
p1214	0	24
p1214	5469	5830
p1214	6018	6035
```

Alternatively, you can convert the TSV file into a MACS3-compatible BED file for peak calling. In this approach:
```bash
# Convert TSV to BED
nanoloop tsv_to_bed \
  --tsv examples/res/p1214_no_pcr_nt_qual.tsv.gz \
  --type nt_qual \
  --output examples/res/p1214_no_pcr_nt_qual.bed.gz

# Call peaks with MACS3
macs3 callpeak -f BED \
  -t examples/res/p1214_no_pcr_nt_qual.bed.gz \
  -n p1214_no_pcr_nt_qual_macs3 \
  -g 8779 \
  --keep-dup all \
  --nomodel \
  --extsize 200 \
  --outdir examples/res/p1214_no_pcr_nt_qual_macs3_peaks
```

Note that in the simulated BED file, each reference position is represented by multiple single-nucleotide read tags. The number of tags at each position reflects:
- The inverse of the average quality score (for `--type nt_qual`)
- The C-to-T conversion frequency (for `--type nt_count`)

To accurately capture these signals, we use `--keep-dup all` to retain all duplicate tags and `--nomodel` to skip fragment length estimation, since our tags are already single-nucleotide in length. The results produced are consistent with the rolling average approach described above.
```bash
cat examples/res/p1214_no_pcr_nt_qual_macs3_peaks/p1214_no_pcr_nt_qual_macs3_peaks.narrowPeak
```
```
p1214	4722	6364	p1214_no_pcr_nt_qual_macs3_peak_1	2414	.	1.57326	245.385	241.431	1044
```

#### Using base count: `--type nt_count`
This approach focuses on the frequency of C-to-T conversions, another signature of R-loops in A3A-treated samples.

##### Step 1: generate TSV
```bash
nanoloop bam_to_tsv \
  --bam examples/bam/p1214_no_pcr.bam \
  --ref examples/ref/p1214.fa \
  --type nt_count \
  --output examples/res/p1214_no_pcr_nt_count.tsv.gz

zcat < examples/res/p1214_no_pcr_nt_count.tsv.gz | head
```
```
#chr    start   end     ref_nt  A       T       C       G       N
p1214   0       1       A       969     0       0       0       0
p1214   1       2       A       1109    0       0       0       0
p1214   2       3       A       1198    0       0       0       0
p1214   3       4       A       1339    0       0       0       0
p1214   4       5       A       1483    0       0       2       0
p1214   5       6       A       1620    0       0       4       0
p1214   6       7       A       1713    0       0       0       0
p1214   7       8       A       1722    0       1       1       0
p1214   8       9       A       1754    0       0       4       0
```

##### Step 2: visualize conversion frequencies
```bash
nanoloop tsv_to_plot \
  --tsv examples/res/p1214_no_pcr_nt_count.tsv.gz \
  --type nt_count \
  --range p1214:0-10000 \
  --add_gc true \
  --output examples/res/p1214_no_pcr_nt_count.jpg
```

The plot shows elevated C-to-T conversion frequencies around 5000-6200, indicating potential R-loop regions:

![Count Distribution Plot](examples/res/p1214_no_pcr_nt_count.jpg)

##### Step 3: call peaks
Use the rolling average approach to identify regions with significant C-to-T conversion:
```bash
nanoloop tsv_to_peak \
  --tsv examples/res/p1214_no_pcr_nt_count.tsv.gz \
  --type nt_count \
  --conversion_cutoff 0.03 \
  --output examples/res/p1214_no_pcr_nt_count_peak.bed.gz
```

The `--conversion_cutoff 0.03` threshold was determined from the visualization in Step 2. The resulting peaks capture the potential R-loop regions:
```bash
zcat < examples/res/p1214_no_pcr_nt_count_peak.bed.gz | head
```
```
p1214   0       24
p1214   5469    5830
p1214   6018    6035
```

Alternatively, you can use MACS3 for peak calling:
```bash
# Convert TSV to BED
nanoloop tsv_to_bed \
  --tsv examples/res/p1214_no_pcr_nt_count.tsv.gz \
  --type nt_count \
  --output examples/res/p1214_no_pcr_nt_count.bed.gz

# Call peaks with MACS3
macs3 callpeak -f BED \
  -t examples/res/p1214_no_pcr_nt_count.bed.gz \
  -n p1214_no_pcr_nt_count_macs3 \
  -g 8779 \
  --keep-dup all \
  --nomodel \
  --extsize 200 \
  --outdir examples/res/p1214_no_pcr_nt_count_macs3_peaks
```

The MACS3 results are consistent with the rolling average approach as well:
```bash
cat examples/res/p1214_no_pcr_nt_count_macs3_peaks/p1214_no_pcr_nt_count_macs3_peaks.narrowPeak
```
```
p1214   5093    6061    p1214_no_pcr_nt_count_macs3_peak_1      2187    .       11.9622 222.658 218.707 671
```

Below is an IGV snapshot comparing peaks called with different approaches:

![IGV snapshot](examples/res/IGV_snapshot.jpg)

### With PCR (read-centric)
The following example uses a downsampled BAM `examples/bam/p1214_with_pcr.bam`, derived from an A3A-treated plasmid sample. This sample is expected to contain R-loops around its middle region. the filename suffix "with_pcr" indicates sequencing was performed after A3A treatment and PCR ammplification, which converts dU to T and introduces other characteristic mutational signatures. 

This is a read-centric approach where `nanoloop` scans each individual read using a sliding windown method. Regions with significantly elevated mutation density (referred to as "hotspots") are identified within each read. The corresponding reference-mapped coordinates of these hotspots are then exported as a BED file, which can be used for peak calling with MACS3.

#### Step 1: generate JSON
First, parse BAM file with `nanoloop bam_to_json` into a JSON file containing information such as read_id, ref_chr, ref_start, ref_end, ref_seq, mutation info (location, type, and base quality of each mutation) on a per-read basis.

```bash
nanoloop bam_to_json \
  --bam examples/bam/p1214_with_pcr.bam \
  --ref examples/ref/p1214.fa \
  --output examples/res/p1214_with_pcr.ndjson.gz
```

#### Step 2: filter JSON 
The resulting JSON file can be further filtered to exclude reads that are unlikely to contain meaningful R-loop signals. For example, reads with too few mutations. This can be done using the `nanoloop filter_json` command.

You can use the `--by count` or `--by frac` options to specify filtering criteria based on the minimum number or fraction of mutations per read, respectively. Additionally, potential false-positive mutations with low basecalling quality can be removed using the `--base_quality_cutoff` option.

```bash
nanoloop filter_json \
  --json examples/res/p1214_with_pcr.ndjson.gz \
  --by count \
  --count_cutoff 10 \
  --base_quality_cutoff 30 \
  --output examples/res/p1214_with_pcr_filtered.ndjson.gz
```

6464 out of 9775 reads passed the filter. They will be used to extract mutation hotspots.

#### Step 3: extract mutation hotspots
You can use `nanoloop json_to_hotspot` to identify mutation hotspots from the filtered JSON file using a sliding window approach. The output is a BED-like file containing the mutation hotspots in each read, that can be passed to MACS3 for peak calling. 

Use `--mutation_type` to specify the types of mutations to consider when calculating the mutation fraction per window. For example:

+ `"all"` includes all mutation types
+ `"CtoT"` includes only C-to-T mutations
+ `"CtoT|CtoG"` includes both C-to-T and C-to-G mutations

Use `--mutation_frac_cutoff` to set the minimum mutation fraction required for a window to be classified as a hotspot. Use `--window_size` to define the size of the sliding window used for hotspot detection. By default, adjacent hotspot windows within `--window_size` distance are merged
```bash
nanoloop json_to_hotspot \
  --json examples/res/p1214_with_pcr_filtered.ndjson.gz \
  --mutation_type CtoT \
  --output examples/res/p1214_with_pcr_filtered_hotspot.bed.gz
```

Below is an example of the output BED file. Note that the first three columns represent the reference coordinates of the hotspot in each read, and the fourth column is the read id.
```bash
zcat < examples/res/p1214_with_pcr_filtered_hotspot.bed.gz | head
```
```
p1214   2647    2672    c2b303ae-93f5-48d3-a76f-5cba776042b5
p1214   5154    5199    b317febd-9ace-4b5e-bafe-ad79ec06669a
p1214   5101    5151    f0a884f9-a99f-481f-8d34-3e8c95b34a01
p1214   5256    5286    5fc1d6fe-bde3-4292-a833-9baa673b58d7
p1214   5628    5673    e930d765-6a87-447a-a730-ecc28338d993
p1214   5276    5321    58ec2b60-da32-4570-9c12-685097a9dc48
p1214   4685    4710    d3c6e465-fb91-449e-9ed1-22043f1cc624
p1214   5116    5141    054dfcee-3310-4e94-b6db-de84e347413f
p1214   5420    5460    ff9eb453-6578-4320-85cf-e5b64a5e2002
p1214   5326    5366    474d6edd-40c3-49d3-b52c-b99a2d62a55a
```

#### Step 4: peak calling with MACS3
A total of 3362 CtoT mutation hotspots were detected in the 6464 reads passed the filter. The BED file can be passed to MACS3 for peak calling.
```bash
macs3 callpeak -f BED \
  -t examples/res/p1214_with_pcr_filtered_hotspot.bed.gz \
  -n p1214_with_pcr_filtered_hotspot_macs3 \
  -g 8779 \
  --keep-dup all \
  --nomodel \
  --extsize 200 \
  --outdir examples/res/p1214_with_pcr_filtered_hotspot_macs3_peaks
```

A single peak is called by MACS3, which is consistent with the reference-based approach above.
```bash
cat examples/res/p1214_with_pcr_filtered_hotspot_macs3_peaks/p1214_with_pcr_filtered_hotspot_macs3_peaks.narrowPeak | head
```
```
p1214   4862    5942    p1214_with_pcr_filtered_hotspot_macs3_peak_1    8661    .       14.6149 869.938 866.131 458
```

#### Step 5: QC plots
You can generate various quality control (QC) plots using the `nanoloop stat_hotspot` and `nanoloop cluster_hotspot` functions. These require additional information to be included in the hotspot BED file, which can be generated with the following command (note the use of `--mutation_type all`):
```bash
nanoloop json_to_hotspot \
  --json examples/res/p1214_with_pcr_filtered.ndjson.gz \
  --output examples/res/p1214_with_pcr_filtered_hotspot_detailed.bed.gz \
  --mutation_type all \
  --ncpus 2 \
  --include_read_id \
  --include_mutation_details \
  --include_ref_seq
```

Then, use `nanoloop stat_hotspot` to produce a series of QC plots related to the detected mutation hotspots:
```bash
nanoloop stat_hotspot \
  --hotspot examples/res/p1214_with_pcr_filtered_hotspot_detailed.bed.gz \
  --output examples/res/stat_p1214_with_pcr_filtered_hotspot
```

This command generates three plots:
+ Hotspot length distribution ![hotspot length distribution](examples/res/stat_p1214_with_pcr_filtered_hotspot/hotspot_length_distribution.svg)

+ Hotspot number per read ![Hotspot number per read](examples/res/stat_p1214_with_pcr_filtered_hotspot/hotspot_number_per_read.svg)

+ Distinct mutation types per hotspot ![Distinct mutation types per hotspot](examples/res/stat_p1214_with_pcr_filtered_hotspot/distinct_mutation_type_per_hotspot.svg)

To generate additional QC plots, use `nanoloop cluster_hotspot`:
```bash
nanoloop cluster_hotspot \
  --hotspot examples/res/p1214_with_pcr_filtered_hotspot_detailed.bed.gz \
  --range p1214:5000-6000 \
  --output examples/res/cluster_p1214_with_pcr_filtered_hotspot
```

This generates two types of heatmaps: 
+ Mutation counts per hotspot (clustered by hotspot)

  + Mutation counts in each hotspot 
  ![Mutation counts in each hotspot](examples/res/cluster_p1214_with_pcr_filtered_hotspot/heatmap_per_hotspot_global_count.svg)

  + Mutation fraction (mutation counts divided by hotspot length) in each hotspot 
  ![Mutation fraction](examples/res/cluster_p1214_with_pcr_filtered_hotspot/heatmap_per_hotspot_global_fraction.svg)

+ Hotspot regions per read (clustered by read)
![Hotspot regions per read](examples/res/cluster_p1214_with_pcr_filtered_hotspot/heatmap_per_read_p1214_5000_6000.svg)

## Bug reports
If you encounter a bug or wish to request a new feature, please open a [new issue](https://github.com/hukai916/nanoloop/issues/new) on the GitHub repository.