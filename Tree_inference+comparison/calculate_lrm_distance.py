#!/usr/bin/env python3
"""
Compute the Lin–Rajan–Moret (LRM) distance between a gene tree and a species tree.

Usage:
  python lrm_distance.py --gene path/to/gene_tree.nwk --species path/to/species_tree.nwk
  # optional flags:
  #   --verbose                 : print details (taxon counts, missing taxa)
  #   --require-exact           : error if the two trees don't have identical tip sets
  #   --min-taxa 4              : require at least N common taxa (default 3)
"""

import argparse
import sys
from cogent3 import load_tree


def load_newick(path: str):
    try:
        return load_tree(filename=path)
    except Exception as e:
        sys.stderr.write(f"[error] failed to load tree '{path}': {e}\n")
        sys.exit(2)


def main():
    ap = argparse.ArgumentParser(
        description="Compute Lin–Rajan–Moret (LRM) tree distance between a gene tree and a species tree."
    )
    ap.add_argument("--gene", required=True, help="Path to the gene tree (Newick).")
    ap.add_argument("--species", required=True, help="Path to the species tree (Newick).")
    ap.add_argument(
        "--require-exact",
        action="store_true",
        help="Fail if the two trees do not have identical tip sets (no pruning).",
    )
    ap.add_argument(
        "--min-taxa",
        type=int,
        default=3,
        help="Minimum number of common taxa required after pruning (default: 3).",
    )
    ap.add_argument(
        "--verbose",
        action="store_true",
        help="Print details about tip sets and pruning.",
    )
    args = ap.parse_args()

    gene_tr = load_newick(args.gene)
    species_tr = load_newick(args.species)

    g_tips = set(gene_tr.get_tip_names())
    s_tips = set(species_tr.get_tip_names())

    if args.verbose:
        print(f"[info] gene tips:    {len(g_tips)}")
        print(f"[info] species tips: {len(s_tips)}")

    if args.require_exact and g_tips != s_tips:
        only_in_gene = sorted(g_tips - s_tips)
        only_in_spec = sorted(s_tips - g_tips)
        sys.stderr.write(
            "[error] tip sets differ and --require-exact was specified.\n"
            f"        only in gene:   {len(only_in_gene)}\n"
            f"        only in species:{len(only_in_spec)}\n"
        )
        if args.verbose:
            if only_in_gene:
                sys.stderr.write("        (examples only-in-gene) " + ", ".join(only_in_gene[:10]) + "\n")
            if only_in_spec:
                sys.stderr.write("        (examples only-in-spec) " + ", ".join(only_in_spec[:10]) + "\n")
        sys.exit(1)

    # Prune to the intersection (standard approach before any topological distance)
    common = sorted(g_tips & s_tips)
    if args.verbose:
        print(f"[info] common taxa:  {len(common)}")

        missing_gene = sorted(s_tips - g_tips)
        missing_spec = sorted(g_tips - s_tips)
        if missing_gene:
            print(f"[info] taxa only in species (dropped): {len(missing_gene)}")
        if missing_spec:
            print(f"[info] taxa only in gene (dropped):    {len(missing_spec)}")

    if len(common) < args.min_taxa:
        sys.stderr.write(
            f"[error] not enough common taxa after pruning: {len(common)} (min required: {args.min_taxa})\n"
        )
        sys.exit(1)

    # Build induced subtrees on the common taxa
    gene_sub    = gene_tr.get_sub_tree(common)
    species_sub = species_tr.get_sub_tree(common)

    # LRM is defined for unrooted trees; convert to unrooted explicitly
    gene_u = gene_sub.unrooted()
    species_u = species_sub.unrooted()

    try:
        lrm = gene_u.tree_distance(species_u, method="lin_rajan_moret")
    except Exception as e:
        sys.stderr.write(f"[error] computing LRM distance failed: {e}\n")
        sys.exit(2)

    # Plain numeric output (easy to parse); add verbose context if requested
    if args.verbose:
        print(f"[result] LRM distance: {lrm}")
    else:
        print(lrm)


if __name__ == "__main__":
    main()
