#!/usr/bin/env python3

import sys
import argparse

from rich.progress import track
from rich.console import Console
from pysam import VariantFile


def main(args):
    ref_paths = {}
    node2var = {}
    variants = {}
    for line in open(args.GFA):
        if line.startswith("P"):
            _, idx, nodes, _ = line.split("\t")
            nodes = [int(x[:-1]) for x in nodes.split(",")]
            if idx.startswith("_alt_"):
                vidx = idx.split("_")[2]
                a = int(idx.split("_")[-1])
                if vidx not in variants:
                    variants[vidx] = {}
                variants[vidx][a] = nodes
            elif idx.startswith(args.tprefix):
                ref_paths[idx] = nodes
            else:
                # reference path
                pass

    for idx, alleles in variants.items():
        assert len(alleles) == 2

    vcf = VariantFile(args.VCF)
    samples = {}
    for rec in track(
        vcf, description="Parsing VCF..", console=Console(file=sys.stderr)
    ):
        if any([a[0] == "<" for a in rec.alts]):
            # Skip symbolic
            continue
        for name, gt in rec.samples.items():
            if name not in samples:
                samples[name] = [{}, {}]
            h1, h2 = gt.allele_indices
            if h1 != 0 or h2 != 0:
                assert rec.id in variants
            # We store for each reference node that is not used by the haplotype, how we have to replace it
            # If the reference allele is made up by multiple nodes, each node will have the same set of alternate nodes
            # but we will fix this later
            if h1 != 0:
                for ref_a in variants[rec.id][0]:
                    samples[name][0][ref_a] = variants[rec.id][h1]
            if h2 != 0:
                for ref_a in variants[rec.id][0]:
                    samples[name][1][ref_a] = variants[rec.id][h2]

    haplotypes = []
    for sample in track(
        samples, description="Parsing samples..", console=Console(file=sys.stderr)
    ):
        hap1, hap2 = samples[sample]
        if hap1 == {}:
            for ref_path_name, ref_path_nodes in ref_paths.items():
                haplotypes.append((f"{sample}_1.{ref_path_name}", ref_path_nodes))
        else:
            for ref_path_name, ref_path_nodes in ref_paths.items():
                new_path = []
                for node in ref_path_nodes:
                    if node in samples[sample][0]:
                        alt_nodes = samples[sample][0][node]
                        if new_path[-1] == alt_nodes[-1]:
                            # reference allele was split in multiple nodes and we already added the alternative nodes
                            continue
                        new_path += alt_nodes
                    else:
                        new_path += [node]
                haplotypes.append((f"{sample}_1.{ref_path_name}", new_path))
        if hap2 == {}:
            for ref_path_name, ref_path_nodes in ref_paths.items():
                haplotypes.append((f"{sample}_2.{ref_path_name}", ref_path_nodes))
        else:
            for ref_path_name, ref_path_nodes in ref_paths.items():
                new_path = []
                for node in ref_path_nodes:
                    if node in samples[sample][1]:
                        alt_nodes = samples[sample][1][node]
                        if new_path[-1] == alt_nodes[-1]:
                            # reference allele was split in multiple nodes and we already added the alternative nodes
                            continue
                        new_path += alt_nodes
                    else:
                        new_path += [node]
                haplotypes.append((f"{sample}_2.{ref_path_name}", new_path))

    for line in open(args.GFA):
        if line.startswith("P"):
            _, idx, nodes, _ = line.split("\t")
            if idx.startswith("_alt_"):
                continue
        print(line, end="")
    for h in haplotypes:
        print("P", h[0], ",".join(f"{x}+" for x in h[1]), "*", sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Haplotype adder",
        description="",
    )
    parser.add_argument("GFA", help="Spliced pangenome in GFA format")
    parser.add_argument("VCF", help="Phased variations in VCF format")
    parser.add_argument(
        "-t",
        dest="tprefix",
        help="Transcript prefix (default: ENST)",
        required=False,
        default="ENST",
    )
    args = parser.parse_args()
    main(args)
