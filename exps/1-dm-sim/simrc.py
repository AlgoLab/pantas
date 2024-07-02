from collections import defaultdict
from pprint import pprint
import sys
from Bio import SeqIO
from intervaltree import Interval, IntervalTree


def main(args):
    genes = defaultdict(dict)

    exonscount = defaultdict(dict)
    exonsinfo = defaultdict(dict)
    junctionscount = defaultdict(dict)

    header = True
    for line in open(args.JUN):
        if header:
            header = False
            continue
        data = line.strip().split("\t")
        geneid = data[9]
        transcriptid = data[10]
        if not transcriptid in genes[geneid]:
            genes[geneid][transcriptid] = IntervalTree()
            exonscount[geneid][transcriptid] = defaultdict(int)
            exonsinfo[geneid][transcriptid] = defaultdict(dict)
            junctionscount[geneid][transcriptid] = defaultdict(int)
        exon_num = data[12]
        exonsinfo[geneid][transcriptid]["strand"] = data[4]
        if exon_num:
            trstart = int(data[13])
            trend = int(data[14])

            exonsinfo[geneid][transcriptid][exon_num]["genome"] = data[:3]
            exonsinfo[geneid][transcriptid][exon_num]["tr"] = [trstart, trend]
            genes[geneid][transcriptid][trstart : trend + 1] = exon_num
            # print(geneid, transcriptid,exon_num,trstart, trend)

    # print(genes["FBgn0039900"]["FBgn0039900_template"][20:50])

    retainedintrons = defaultdict(dict)
    for line in open(args.ANN):
        if line.startswith("event"):
            continue
        (
            etype,
            transcript,
            template,
            genomic_start,
            genomic_end,
            transcriptomic_start,
            transcriptomic_end,
        ) = line.strip("\n").split("\t")
        if etype != "ir":
            continue
        geneid = template.split("_")[0]
        transcriptid = transcript
        # CHECME: assuming we have only one IR per gene
        # FIXME: I'm not storing all info as per other exons
        retainedintrons[geneid][transcriptid] = [
            int(transcriptomic_start),
            int(transcriptomic_end),
            int(genomic_start),
            int(genomic_end),
            0,
        ]

    for record in SeqIO.parse(args.FQ, "fastq"):
        rname, mate1, mate2 = record.name.split(";")
        transcriptid = rname.split("/")[1]
        geneid = "_".join(transcriptid.split("_")[:-1])
        # try:
        if True:
            for m in [mate1, mate2]:
                # for m in [mate1]:
                # print(rname, mate1, mate2)
                se = m.split(":")[1].split("-")
                if len(se) != 2:
                    print("Read mate with -", file=sys.stderr)
                    continue
                s, e = se
                s = int(s)
                e = int(e)
                qres = sorted(genes[geneid][transcriptid][s:e])
                # print(qres, len(qres), len(qres) == 1)
                if len(qres) == 1:
                    # qres.data["count"] += 1
                    # print('oooo', qres[0].data)
                    exonscount[geneid][transcriptid][qres[0].data] += 1
                else:
                    for i, j in zip(qres, qres[1:]):
                        junctionscount[geneid][transcriptid][(i.data, j.data)] += 1

                # Retained intron
                if transcriptid not in retainedintrons[geneid]:
                    continue
                ts, te, _, _, _ = retainedintrons[geneid][transcriptid]
                if (ts <= s and s <= te) or (ts <= e and e <= te):
                    retainedintrons[geneid][transcriptid][4] += 1
        # except:
        #     # print("skipping", rname, file=sys.stderr)
        #     continue

    # pprint(exonscount)
    # pprint(exonsinfo)
    # pprint(junctionscount)

    print(
        "seqnames,start,end,strand,type,gene_id,transcript_id,gene_exon_number,tr_start,tr_end,read_count"
    )
    for geneid in exonsinfo:
        for transcriptid in exonsinfo[geneid]:
            prevex = None
            for exon in exonsinfo[geneid][transcriptid]:
                if exon == "strand":
                    continue

                if prevex:
                    _seq = exonsinfo[geneid][transcriptid][exon]["genome"][0]
                    if exonsinfo[geneid][transcriptid]["strand"] == "-":
                        _jstart = exonsinfo[geneid][transcriptid][exon]["genome"][2]
                        _jend = exonsinfo[geneid][transcriptid][prevex]["genome"][1]
                    else:
                        # This needs to be checked
                        _jstart = exonsinfo[geneid][transcriptid][prevex]["genome"][2]
                        _jend = exonsinfo[geneid][transcriptid][exon]["genome"][1]
                    print(
                        # *[0,0,0], #*exonsinfo[geneid][transcriptid][exon]["genome"],
                        _seq,
                        _jstart,
                        _jend,
                        exonsinfo[geneid][transcriptid]["strand"],
                        "junction",
                        geneid,
                        transcriptid,
                        f"{prevex}-{exon}",
                        ".",
                        ".",  # *exonsinfo[geneid][transcriptid][exon]["tr"],
                        junctionscount[geneid][transcriptid][(prevex, exon)],
                        sep=",",
                    )
                prevex = exon

                print(
                    *exonsinfo[geneid][transcriptid][exon]["genome"],
                    exonsinfo[geneid][transcriptid]["strand"],
                    "exon",
                    geneid,
                    transcriptid,
                    exon,
                    *exonsinfo[geneid][transcriptid][exon]["tr"],
                    exonscount[geneid][transcriptid][exon],
                    sep=",",
                )

            # retained intron
            if transcriptid not in retainedintrons[geneid]:
                continue
            ts, te, gs, ge, c = retainedintrons[geneid][transcriptid]
            print(
                exonsinfo[geneid][transcriptid]["1"]["genome"][0],
                gs,
                ge,
                exonsinfo[geneid][transcriptid]["strand"],
                "-exon",
                geneid,
                transcriptid,
                0,
                s,
                e,
                c,
                sep=",",
            )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Compute correct read counts on exons and junctions"
    )
    parser.add_argument("FQ", type=str, help="Path to simulated FQ")
    parser.add_argument("JUN", type=str, help="Path to ASimulatoR exon junction tsv")
    parser.add_argument("ANN", type=str, help="Path to ASimulatoR event annotation tsv")
    args = parser.parse_args()

    main(args)
