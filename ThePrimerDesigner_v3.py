#!/usr/bin/env python3
"""
ThePrimerDesigner — SNP-anchored PCR primer designer
Version 3.0 (first public release)

Author: Arif Mohammad Tanmoy

Description
-----------
Given a TSV of variant sites (RefID, Genotype/Lineage/Variant/SNPname, Position, RefBase, AltBase) and a reference FASTA, this script designs one or more PCR primer pairs that amplify a region *containing* the SNP. The tool allows to fix PCR product size range, primer length, primer Tm, and the maximum Tm difference within a pair.

"""

# ---------------------- Args ----------------------
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed

from Bio import SeqIO
from primer3 import bindings as p3

def parse_args():
    p = argparse.ArgumentParser(
        description="Design PCR primer pairs around SNPs (no off-target/dimer screening)."
    )
    p.add_argument("-v", "--variants",  required=True,
                   help="TSV: RefID,Genotype,Position,RefBase,AltBase")
    p.add_argument("-r", "--reference", required=True,
                   help="Reference FASTA (Seq IDs must match RefID)")
    p.add_argument("-O", "--output",    default="primers_out.tsv",
                   help="Output TSV (default: primers_out.tsv)")

    # primer3 design knobs
    p.add_argument("--primer-min", type=int,   default=20, help="Min primer length (default 20)")
    p.add_argument("--primer-opt", type=int,   default=25, help="Opt primer length (default 25)")
    p.add_argument("--primer-max", type=int,   default=35, help="Max primer length (default 35)")

    p.add_argument("--tm-min",     type=float, default=59.0, help="Min primer Tm (°C, default 59)")
    p.add_argument("--tm-opt",     type=float, default=62.0, help="Opt primer Tm (°C, default 62)")
    p.add_argument("--tm-max",     type=float, default=65.0, help="Max primer Tm (°C, default 65)")

    p.add_argument("--product-range", nargs=2, type=int, metavar=("MIN", "MAX"),
                   default=[380, 420],
                   help="PCR product size range (bp), e.g. --product-range 360 440")

    p.add_argument("--num-sets",   type=int,   default=5,
                   help="Max primer pairs per SNP (default 5)")

    p.add_argument("--max-tm-diff", type=float, default=2.0,
                   help="Max allowed ΔTm within a pair (°C, default 2.0)")

    p.add_argument("-t", "--threads", type=int, default=4,
                   help="Worker threads (default 4)")
    return p.parse_args()


# ---------------------- IO ----------------------
def load_ref_fasta(path):
    """Load FASTA to dict of SeqRecord keyed by record.id."""
    return SeqIO.to_dict(SeqIO.parse(path, "fasta"))


def read_variants_tsv(path):
    """Yield rows: (chrom, genotype, pos(int), ref, alt). Skips comments/blank."""
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        for row in rdr:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 5:
                raise ValueError("Each variant line must be: RefID,Genotype,Position,RefBase,AltBase")
            chrom, genotype, pos_str, ref, alt = row[:5]
            yield chrom, genotype, int(pos_str), ref, alt


# ---------------------- Primer design ----------------------
def design_primers_for_site(ref_seqs, chrom, pos, args, global_args):
    """
    Pull a window around the SNP and call primer3 with SEQUENCE_TARGET marking the SNP.
    Window size uses the *max* product size so primer3 has room to place primers.
    """
    if chrom not in ref_seqs:
        raise KeyError(f"RefID '{chrom}' not found in reference FASTA")

    full = str(ref_seqs[chrom].seq)
    L = len(full)

    # Make a symmetric window of size ~ max product around the SNP
    max_prod = args.product_range[1] if hasattr(args, "product_range") else args.product_range[1]
    center = pos - 1  # 0-based index of SNP
    start  = max(0, center - max_prod)
    end    = min(L, center + max_prod)
    template = full[start:end]

    # SNP index within the cropped template
    snp_offset = center - start
    if snp_offset < 0 or snp_offset >= len(template):
        raise ValueError("Computed SNP offset lies outside the cropped template")

    seq_args = {
        "SEQUENCE_ID":       f"{chrom}:{start+1}-{end}",
        "SEQUENCE_TEMPLATE": template,
        "SEQUENCE_TARGET":   [snp_offset, 1],  # force the SNP to be inside the product
    }

    res = p3.design_primers(global_args, seq_args)

    out = []
    for i in range(args.num_sets):
        left  = res.get(f"PRIMER_LEFT_{i}_SEQUENCE")
        right = res.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
        if not left or not right:
            continue
        left_tm  = res[f"PRIMER_LEFT_{i}_TM"]
        right_tm = res[f"PRIMER_RIGHT_{i}_TM"]
        prod     = res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]

        # ΔTm filter
        if abs(left_tm - right_tm) > args.max_tm_diff:
            continue

        out.append({
            "left_seq":     left,
            "right_seq":    right,
            "left_tm":      left_tm,
            "right_tm":     right_tm,
            "product_size": prod
        })
    return out


def process_variant(row, ref_seqs, args, global_args):
    chrom, genotype, pos, _ref, _alt = row
    pairs = design_primers_for_site(ref_seqs, chrom, pos, args, global_args)
    results = []
    for p in pairs:
        results.append({
            "chrom":        chrom,
            "genotype":     genotype,
            "pos":          pos,
            "left_seq":     p["left_seq"],
            "right_seq":    p["right_seq"],
            "left_tm":      round(p["left_tm"], 2),
            "right_tm":     round(p["right_tm"], 2),
            "product_size": int(p["product_size"]),
        })
    return results


# ---------------------- Main ----------------------
def main():
    args = parse_args()

    # primer3 global args
    global_args = {
        "PRIMER_MIN_SIZE": args.primer_min,
        "PRIMER_OPT_SIZE": args.primer_opt,
        "PRIMER_MAX_SIZE": args.primer_max,

        "PRIMER_MIN_TM":   args.tm_min,
        "PRIMER_OPT_TM":   args.tm_opt,
        "PRIMER_MAX_TM":   args.tm_max,

        # product range is a list of [min, max]
        "PRIMER_PRODUCT_SIZE_RANGE": [args.product_range],
        "PRIMER_NUM_RETURN": args.num_sets,
    }

    ref_seqs = load_ref_fasta(args.reference)
    variants = list(read_variants_tsv(args.variants))

    results = []
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = [pool.submit(process_variant, row, ref_seqs, args, global_args)
                   for row in variants]
        for fut in as_completed(futures):
            results.extend(fut.result())

    # write results
    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(
            out,
            fieldnames=["chrom","genotype","pos","left_seq","right_seq","left_tm","right_tm","product_size"],
            delimiter="\t"
        )
        w.writeheader()
        w.writerows(results)


if __name__ == "__main__":
    main()

