# The-Primer-Designer
**SNP-anchored PCR primer designer.**

Given a TSV of variant sites (**RefID, Genotype/Lineage/Variant/SNPname, Position, RefBase, AltBase**) and a reference FASTA, this script designs one or more PCR primer pairs that amplify a region *containing* the SNP. You can fix the PCR **product size range**, **primer length**, **primer Tm**, the **maximum ΔTm** within a pair, and the **number of pairs** to return.

> **Note:** This version does **not** perform off-target or primer-dimer screening. It focuses on robust `primer3` design around the SNP.

---

## Features

* Makes sure that the SNP **inside** the amplicon (`SEQUENCE_TARGET`).
* Adjustable **product length** (min–max), **primer length**, **Tm** bounds, and **ΔTm** ≤ X °C within a pair.
* Can run on multi-threads and run across many variants.

---

## Requirements

* **Python** ≥ 3.8
* Packages:

  * `primer3-py`
  * `biopython`

### Install (pip)

```bash
python -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install primer3-py biopython
```

### Install (conda / mamba)

```bash
mamba create -n primers python=3.10 primer3-py biopython -c conda-forge
mamba activate primers
```

---

## Inputs & Outputs

### Variants TSV (tab-delimited)

Each line (header optional):

```
RefID    Genotype/Lineage/Variant/SNPname    Position    RefBase    AltBase
```

* `RefID` must match a sequence ID in the reference FASTA.
* `Position` is an **1-based INTEGER**

**Example**

```
chr1     2.3.1                               2811222     C          T
Ref_Acc  lineageX                            105433      G          A
```

### Reference FASTA

One or more sequences; IDs must match your `RefID` values.

**Example**

```fasta
>chr1
ACGTTG...
>Ref_Acc
TTAGCA...
```

### Output TSV

```
chrom, genotype, pos, left_seq, right_seq, left_tm, right_tm, product_size
```

---

## Usage

```bash
python ThePrimerDesigner.py \
  -v variants.tsv \
  -r reference.fasta \
  -O primers_out.tsv \
  --product-range 380 420 \
  --primer-min 20 --primer-opt 25 --primer-max 35 \
  --tm-min 59 --tm-opt 62 --tm-max 65 \
  --max-tm-diff 2.0 \
  --num-sets 5 \
  -t 8
```

> You can try using `--product-range` (hyphenated) or `--product_range` (underscored) depending on the version—use the one your CLI help shows.

### Command-line options

```
-v, --variants            TSV of variants (RefID,Genotype,Position,RefBase,AltBase) [required]
-r, --reference           Reference FASTA (IDs must match RefID)                    [required]
-O, --output              Output TSV (default: primers_out.tsv)

--product-range MIN MAX   PCR product size range in bp (default: 380 420)

--primer-min INT          Min primer length (default: 20)
--primer-opt INT          Opt primer length (default: 25)
--primer-max INT          Max primer length (default: 35)

--tm-min FLOAT            Min primer Tm, °C (default: 59)
--tm-opt FLOAT            Opt primer Tm, °C (default: 62)
--tm-max FLOAT            Max primer Tm, °C (default: 65)

--max-tm-diff FLOAT       Max allowed ΔTm within a pair, °C (default: 2.0)
--num-sets INT            Max primer pairs per SNP (default: 5)

-t, --threads INT         Worker threads (default: 4)
```

---

## Quick Start

1. Prepare `variants.tsv` and `reference.fasta`.
2. Install dependencies.
3. Run the command in **Usage**.
4. Inspect `primers_out.tsv` for primer sequences, Tm, and product size.

---

## Tips & Troubleshooting

* **“TARGET beyond end of sequence”**
  Check that `Position` is within the reference sequence for `RefID` and is 1-based.

* **No primers returned**
  Constraints may be too strict. Widen:

  * `--product-range` (e.g., `360 440`)
  * `--tm-min/--tm-max` (e.g., `58 66`)
  * `--primer-min/--primer-max` (e.g., `18 36`)
  * `--max-tm-diff` (e.g., `3.0`)

* **Performance**
  Increase `-t/--threads` for many variants/large references.

---

## What this tool is (and isn’t)

* ✅ A fast, flexible wrapper around **primer3** to design **HUNDREDS** of **SNP-spanning** amplicons with user-controlled constraints.
* ❌ Not an off-target screening tool (no BLAST/in-silico PCR/primer-dimer checks in this script).

---

## License

**MIT License** (free to use, modify, and distribute).

---
