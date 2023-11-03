import sys
import numpy as np
import statistics

ETYPES = ["ES", "IR", "A3", "A5"]
EMAP = {"SE": "ES", "RI": "IR", "A3SS": "A3", "A5SS": "A5"}





def get_interval_s(c, s, e):
    return f"{c}:{s}-{e}"


def parse_rmats_se(fpath, novel, pvalue=0.05):
    EVENTS = list()
    for line in open(fpath):
        if line.startswith("ID"):
            continue
        try:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                ex_s,
                ex_e,
                usex_s,
                usex_e,
                dsex_s,
                dsex_e,
                _idx,
                inc_jc_1,
                sk_jc_1,
                inc_jc_2,
                sk_jc_2,
                inc_len,
                sk_len,
                pv,
                fdr,
                inclvl_1,
                inclvl_2,
                delta_incl,
            ) = line.strip("\n").split("\t")
            pv = float(pv)
            if pv > pvalue:
                continue
        except ValueError:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                ex_s,
                ex_e,
                usex_s,
                usex_e,
                dsex_s,
                dsex_e,
            ) = line.strip("\n").split("\t")

        ex_s, usex_s, dsex_s = int(ex_s), int(usex_s), int(dsex_s)
        ex_e, usex_e, dsex_e = int(ex_e), int(usex_e), int(dsex_e)
        # Starts are 0-based, ends aren't (or they are exclusive)
        # ex_s += 1
        usex_s += 1
        # dsex_s += 1

        chrom = chrom[3:]

        intron1 = (int(usex_e), int(ex_s))
        intron2 = (int(ex_e), int(dsex_s))
        # intron3 = (intron1[0], intron2[1])
        c1 = inclvl_1.split(",")
        c2 = inclvl_2.split(",")
        if "NA" in c1:
            inclvl_1 = "NAN"
        else:
            inclvl_1 = statistics.mean([float(i) for i in c1])
        if "NA" in c2:
            inclvl_2 = "NAN"
        else:
            inclvl_2 = statistics.mean([float(i) for i in c2])
        # k = f"{chrom}:{intron1[0]}-{intron1[1]}-{intron2[0]}-{intron2[1]}"
        k = [
            chrom,
            gene,
            strand,
            get_interval_s(chrom, intron1[0] + 1, intron2[1]),
            get_interval_s(chrom, intron1[0] + 1, intron1[1]),
            get_interval_s(chrom, intron2[0] + 1, intron2[1]),
            "W1",
            "w2",
            inclvl_1,
            inclvl_2,
            delta_incl,
            pv
        ]
     
        EVENTS.append(k)
    return EVENTS


def parse_rmats_a3(fpath, novel, pvalue=0.05):
    EVENTS = list()
    for line in open(fpath):
        if line.startswith("ID"):
            continue
        # Exon position are 0-based
        try:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                lex_s,
                lex_e,
                sex_s,
                sex_e,
                ex_s,
                ex_e,
                _idx,
                inc_jc_1,
                sk_jc_1,
                inc_jc_2,
                sk_jc_2,
                inc_len,
                sk_len,
                pv,
                fdr,
                inclvl_1,
                inclvl_2,
                delta_incl,
            ) = line.strip("\n").split("\t")
            pv = float(pv)
            if pv > pvalue:
                continue
        except ValueError:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                lex_s,
                lex_e,
                sex_s,
                sex_e,
                ex_s,
                ex_e,
            ) = line.strip("\n").split("\t")

        # Starts are 0-based, ends aren't (or they are exclusive)
        ex_s, lex_s, sex_s = int(ex_s), int(lex_s), int(sex_s)
        ex_e, lex_e, sex_e = int(ex_e), int(lex_e), int(sex_e)
        ex_s += 1
        sex_s += 1
        sex_s += 1

        chrom = chrom[3:]

        # strand +
        longer_intron = (int(ex_e), int(sex_s))
        shorter_intron = (int(ex_e), int(lex_s))
        if strand == "-":
            longer_intron = (int(sex_e), int(ex_s))
            shorter_intron = (int(lex_e), int(ex_s))

        c1 = inclvl_1.split(",")
        c2 = inclvl_2.split(",")
        if "NA" in c1:
            inclvl_1 = "NAN"
        else:
            inclvl_1 = statistics.mean([float(i) for i in c1])
        if "NA" in c2:
            inclvl_2 = "NAN"
        else:
            inclvl_2 = statistics.mean([float(i) for i in c2])
        # k = ""
        if strand == "+":

            k = [
                chrom,
                gene,
                strand,
                get_interval_s(chrom, longer_intron[0] + 1, longer_intron[1] - 2),
                get_interval_s(chrom, shorter_intron[0] + 1, shorter_intron[1]),
                ".",
                "W1",
                "w2",
                inclvl_1,
                inclvl_2,
                delta_incl,
            ]
        else:
            k = [
                chrom,
                gene,
                strand,
                get_interval_s(chrom, longer_intron[0] + 1, longer_intron[1] - 1),
                get_interval_s(chrom, shorter_intron[0] + 1, shorter_intron[1] - 1),
                ".",
                "W1",
                "w2",
                inclvl_1,
                inclvl_2,
                delta_incl,
                pv
            ]
      
        EVENTS.append(k)
    return EVENTS


def parse_rmats_a5(fpath, novel, pvalue=0.05):
    EVENTS = list()
    for line in open(fpath):
        if line.startswith("ID"):
            continue
        # Exon position are 0-based
        try:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                lex_s,
                lex_e,
                sex_s,
                sex_e,
                ex_s,
                ex_e,
                _idx,
                inc_jc_1,
                sk_jc_1,
                inc_jc_2,
                sk_jc_2,
                inc_len,
                sk_len,
                pv,
                fdr,
                inclvl_1,
                inclvl_2,
                delta_incl,
            ) = line.strip("\n").split("\t")
            pv = float(pv)
            if pv > pvalue:
                continue
        except ValueError:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                lex_s,
                lex_e,
                sex_s,
                sex_e,
                ex_s,
                ex_e,
            ) = line.strip("\n").split("\t")

        # Starts are 0-based, ends aren't (or they are exclusive)
        ex_s, lex_s, sex_s = int(ex_s), int(lex_s), int(sex_s)
        ex_e, lex_e, sex_e = int(ex_e), int(lex_e), int(sex_e)
        ex_s += 1
        sex_s += 1
        sex_s += 1

        chrom = chrom[3:]

        # strand +
        longer_intron = (int(sex_e), int(ex_s))
        shorter_intron = (int(lex_e), int(ex_s))
        if strand == "-":
            longer_intron = (int(ex_e), int(sex_s))
            shorter_intron = (int(ex_e), int(lex_s))

        c1 = inclvl_1.split(",")
        c2 = inclvl_2.split(",")
        if "NA" in c1:
            inclvl_1 = "NAN"
        else:
            inclvl_1 = statistics.mean([float(i) for i in c1])
        if "NA" in c2:
            inclvl_2 = "NAN"
        else:
            inclvl_2 = statistics.mean([float(i) for i in c2])

        if strand == "+":

            k = [
                chrom,
                gene,
                strand,
                get_interval_s(chrom, longer_intron[0] + 1, longer_intron[1] - 1),
                get_interval_s(chrom, shorter_intron[0] + 1, shorter_intron[1] - 1),
                ".",
                "W1",
                "w2",
                inclvl_1,
                inclvl_2,
                delta_incl,
            ]
        else:

            k = [
                chrom,
                gene,
                strand,
                get_interval_s(chrom, shorter_intron[0] + 1, longer_intron[1] - 2),
                get_interval_s(chrom, longer_intron[0] + 1, shorter_intron[1]),
                ".",
                "W1",
                "w2",
                inclvl_1,
                inclvl_2,
                delta_incl,
                pv
            ]
        EVENTS.append(k)
    return EVENTS


def parse_rmats_ri(fpath, novel, pvalue=0.05):
    EVENTS = list()
    for line in open(fpath):
        if line.startswith("ID"):
            continue
        # Exon position are 0-based
        try:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                ex_s,  # retained exon
                ex_e,
                fex_s,  # first exon
                fex_e,
                sex_s,  # second exon
                sex_e,
                _idx,
                inc_jc_1,
                sk_jc_1,
                inc_jc_2,
                sk_jc_2,
                inc_len,
                sk_len,
                pv,
                fdr,
                inclvl_1,
                inclvl_2,
                delta_incl,
            ) = line.strip("\n").split("\t")
            pv = float(pv)
            if pv > pvalue:
                continue
        except ValueError:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                ex_s,  # retained exon
                ex_e,
                fex_s,  # first exon
                fex_e,
                sex_s,  # second exon
                sex_e,
            ) = line.strip("\n").split("\t")

        # Starts are 0-based, ends aren't (or they are exclusive)
        ex_s, fex_s, sex_s = int(ex_s), int(fex_s), int(sex_s)
        ex_e, fex_e, sex_e = int(ex_e), int(fex_e), int(sex_e)
        ex_s += 1
        # sex_s += 1
        sex_s += 1
        fex_s += 1

        chrom = chrom[3:]

        assert ex_s == fex_s and ex_e == sex_e
        # strand +
        k = f"{chrom}:{ex_s}-{fex_e}-{sex_s-1}-{ex_e}"
        c1 = inclvl_1.split(",")
        c2 = inclvl_2.split(",")
        if "NA" in c1:
            inclvl_1 = "NAN"
        else:
            inclvl_1 = statistics.mean([float(i) for i in c1])
        if "NA" in c2:
            inclvl_2 = "NAN"
        else:
            inclvl_2 = statistics.mean([float(i) for i in c2])
        k = [
            chrom,
            gene,
            strand,
            get_interval_s(chrom, fex_e + 1, sex_s - 1),
            ".",
            ".",
            "W1",
            "w2",
            inclvl_1,
            inclvl_2,
            delta_incl,
            pv
        ]
        EVENTS.append(k)
    return EVENTS


def main(rmats_prefix, min_pvalue, output):


    rmats = {x: set() for x in ETYPES}
    rmats["ES"] = parse_rmats_se(rmats_prefix + "/SE.MATS.JC.txt", False, min_pvalue)
    rmats["A3"] = parse_rmats_a3(rmats_prefix + "/A3SS.MATS.JC.txt", False, min_pvalue)
    rmats["A5"] = parse_rmats_a5(rmats_prefix + "/A5SS.MATS.JC.txt", False, min_pvalue)
    rmats["IR"] = parse_rmats_ri(rmats_prefix + "/RI.MATS.JC.txt", False, min_pvalue)

    with open(output, "w") as f:
        for etype in rmats:
            for e in rmats[etype]:
                print(etype, "annotated", *map(lambda x: str(x).strip('"'), e), sep=",", file=f)



if __name__ == "__main__":
    rmats_path = snakemake.params.rmats
    min_pvalue = snakemake.params.p_value
    output = snakemake.output.csv
    main(rmats_path, min_pvalue, output)

