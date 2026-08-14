# RESCUE Benchmark Data

This document describes all data used in the RESCUE paper experiments:
sources and versions (§1), how the benchmark query sets and ground truth
were constructed (§2), which experiments use which data (§3), and checksums
for verification (§4).

Everything is derived from **Ensembl release 113** (downloaded 2026-04-03 to
2026-04-06). All processing steps are deterministic; no random sampling or
random seeds are involved at any stage. The full data tree can be regenerated
with the scripts in `scripts/prepare_data/`.

## 1. Data sources

All primary data come from the Ensembl release 113 FTP archive
(`https://ftp.ensembl.org/pub/release-113/`):

| Data | Species | Assembly | Ensembl FTP path (under `pub/release-113/`) |
| --- | --- | --- | --- |
| Genome FASTA (toplevel) | human | GRCh38 | `fasta/homo_sapiens/dna/` |
| Genome FASTA (toplevel) | mouse | GRCm39 | `fasta/mus_musculus/dna/` |
| Genome FASTA (toplevel) | dog | ROS_Cfam_1.0 | `fasta/canis_lupus_familiaris/dna/` |
| Genome FASTA (toplevel) | chicken | bGalGal1.mat.broiler.GRCg7b | `fasta/gallus_gallus/dna/` |
| Genome FASTA (toplevel) | zebrafish | GRCz11 | `fasta/danio_rerio/dna/` |
| Protein FASTA (`pep.all`) | all five species | as above | `fasta/<species>/pep/` |
| Gene annotation GTF | all five species | as above | `gtf/<species>/<Species>.<assembly>.113.gtf.gz` |
| Orthology/paralogy tables | human vs each target | — | `tsv/ensembl-compara/homologies/` |

Derived files produced from these sources:

| Derived file | Built from | Transformation |
| --- | --- | --- |
| `proteins/<species>_canonical.fa` | `pep.all` FASTA + GTF | one protein per gene: the translation of the transcript with the longest total CDS length (Ensembl stable IDs, version suffixes stripped) |
| `gene_structures/<species>_CDS.bed` | GTF | `CDS` features of `protein_coding` genes; 1-based inclusive GTF coordinates converted to 0-based half-open; columns: `chrom start end gene_id transcript_id strand phase` |
| `gene_structures/<species>_exons.bed` | GTF | same, for `exon` features |
| `orthologs/human_vs_<species>_one2one.tsv` | Compara homology dump | rows for the human–target species pair with `homology_type = ortholog_one2one` |
| `orthologs/human_vs_<species>_paralog_rich.tsv` | Compara homology dump | rows with `homology_type ∈ {ortholog_one2many, ortholog_many2many}`, plus an added `n_target_loci` column = number of distinct target genes per human gene |

## 2. Benchmark query sets and ground truth

The benchmark query sets live in `data/v3/<species>/` (one-to-one set) and
`data/v3/<species>_paralog/` (paralog set). Construction is fully
deterministic.

### One-to-one set (`data/v3/<species>/`)

1. Start from all human↔target `ortholog_one2one` pairs (Compara).
2. Keep pairs where the human protein is present in the human canonical
   proteome and the target gene has CDS annotations.
3. Sort the valid pairs by protein sequence identity and select **2,000**
   pairs by evenly spaced sampling over the sorted list (step = N/2000), so
   the sample covers the full identity distribution from distant to easy.
4. Ground truth per query: the CDS exon coordinates of the target gene's
   canonical transcript (longest total CDS), 0-based half-open.
5. Stratification labels (`strata.json`): identity tier
   (`1:1_easy` ≥80%, `1:1_moderate` 60–80%, `1:1_distant` <60%),
   `micro_exon` (any CDS exon < 30 bp), `long_intron` (any intron > 100 kb).

### Paralog set (`data/v3/<species>_paralog/`)

1. Start from `ortholog_one2many`/`ortholog_many2many` pairs; keep human
   genes with `n_target_loci ≥ 2`.
2. The **primary** ground-truth locus is the highest-identity target
   ortholog; all other target loci are recorded as acceptable alternative
   loci (`paralog_info.json`, `all_targets`).
3. Sort candidates by number of paralog targets (descending, ties broken by
   identity) and take the top **500** per species (chicken has only 488
   qualifying genes, so all are kept).

### Ground-truth tables

Per species, three GT tables are derived (schema: `query_id chrom start end
strand gene_id`, gene body span = min/max of CDS exons):

- `<species>.one2one.tsv` — primary locus per one-to-one query.
- `<species>.paralog_strict.tsv` — primary (highest-identity) locus per
  paralog query.
- `<species>.paralog_broad.tsv` — primary locus **plus every acceptable
  co-ortholog locus** from `paralog_info.json` (coordinates from the CDS
  BED), used for 1:many evaluation.

Query counts per species:

| Species | one2one | paralog (strict) | paralog broad rows |
| --- | --- | --- | --- |
| mouse | 2000 | 500 | 5897 |
| dog | 2000 | 500 | 3013 |
| chicken | 2000 | 488 | 4625 |
| zebrafish | 2000 | 500 | 3535 |

## 3. Data usage per experiment

| Experiment | Data used |
| --- | --- |
| End-to-end benchmark (`scripts/exp1/`) | `data/v3/` query proteins + ground truth; genome FASTA; CDS BED |
| Forced-locus diagnosis (`scripts/exp5/`) | `data/v3/` ground-truth loci; genome FASTA |
| Best fixed recipe / oracle analysis (`scripts/exp6/`) | no primary data — analyzes benchmark outputs only |
| RESCUE pipeline (`scripts/tprs/`) | all of the above, plus canonical proteomes and Compara ortholog tables |

## 4. Checksums (MD5)

For verification of a regenerated data tree against the one used in the
paper:

```
32cfe205a8fc667094334f74085dac80  genomes/mouse.fa
74625670c5e784cd11b5ae57d27b84cc  genomes/dog.fa
eeceff90b8c4200216bd5f5b942c65ae  genomes/chicken.fa
3648f1774b9d1c253bce7acf96e67ff7  genomes/zebrafish.fa
a69f594f87056cfa203e38321443a744  proteins/human_canonical.fa
eed4ee7474d21a9305ddf6dec998d793  proteins/mouse_canonical.fa
098c3d56533035afb29069f17615871a  proteins/dog_canonical.fa
d9746e8e9d361b8e5806371f7388fcb4  proteins/chicken_canonical.fa
a077127f24a807cc00cc148b2c17da5c  proteins/zebrafish_canonical.fa
4cce1972265394d3aecad0fd000eeae2  orthologs/human_vs_canis_lupus_familiaris_one2one.tsv
89cdbb94acc6e1798393d742a4557ae7  orthologs/human_vs_canis_lupus_familiaris_paralog_rich.tsv
16a9d8a506db48d9fe852d6c6ff6b8b6  orthologs/human_vs_danio_rerio_one2one.tsv
d886a567a01f9d7cc6e547dbf4a81d71  orthologs/human_vs_danio_rerio_paralog_rich.tsv
7b133cfc72e700512b5f4ed659c445e0  orthologs/human_vs_gallus_gallus_one2one.tsv
59972072e7b4ddb9a484063192b6a7e9  orthologs/human_vs_gallus_gallus_paralog_rich.tsv
5cb1d56889bc178bfa14315a4c49385b  orthologs/human_vs_mus_musculus_one2one.tsv
7950e7eb7f77a81faf235fc17f18f296  orthologs/human_vs_mus_musculus_paralog_rich.tsv
f35e97d53aef70847b7101e653a207ca  gene_structures/mouse_CDS.bed
c980d250acfe4aa83040787e4ebd8318  gene_structures/dog_CDS.bed
7f151040107162b5dca0eae9ec46840a  gene_structures/chicken_CDS.bed
9bbe256955d18c83c31e496b2e6e2dd0  gene_structures/zebrafish_CDS.bed
```

Query-set checksums (`data/v3/`):

```
2027ac6b47ad8ed28f9fc883dcaa83ee  mouse/proteins.fa
f6982e1621acaa300e98eea462e55c30  dog/proteins.fa
bc87cb8e70607b9a87c8f8b443abeaf6  chicken/proteins.fa
7eae53ce18567e4e0f27263c4fd57895  zebrafish/proteins.fa
abeccea768dc8de4a6fd7296256ab655  mouse_paralog/proteins.fa
dfc2d2baa4c446e4e6117dc9add9420a  dog_paralog/proteins.fa
4dd82e009ab2dbd7a5c384d429bc6266  chicken_paralog/proteins.fa
4d17949f89e041b4e5402349513cba8b  zebrafish_paralog/proteins.fa
5edfa40cf803432ad91c256288166efc  mouse/ground_truth.json
2d47e05c41f9a0b034c5a3360bd28414  dog/ground_truth.json
ce40b67494006de901363e9e83c31b87  chicken/ground_truth.json
3a42707b6596907a43553de81581b392  zebrafish/ground_truth.json
c2aaf5588126d0724f9420382387bd16  mouse_paralog/ground_truth.json
f09793ee37f738b11524289bd1eb92b9  dog_paralog/ground_truth.json
e83b6f3a1f381b66a86905f5b0c9a4ff  chicken_paralog/ground_truth.json
18dbee97980af1bd610e06cbe9738741  zebrafish_paralog/ground_truth.json
ed99b2d989ec3ee2ad9083c6261c26c0  mouse/pairs.tsv
a969389df151bc5a5b069be451756d6c  dog/pairs.tsv
49f502aa9a7c835e5d2a42f4f8be0be2  chicken/pairs.tsv
218ceb6e0914cd7261efc56c1b40a8a6  zebrafish/pairs.tsv
633ec11f4fb7595b6f1de39c9ee1122a  mouse_paralog/pairs.tsv
19b4944cacec4ed9359454518f79a147  dog_paralog/pairs.tsv
ab1584859883bde2b58edc153607b22b  chicken_paralog/pairs.tsv
087b6d8c02cc53051764d4a6f6ef9365  zebrafish_paralog/pairs.tsv
```
