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
    for line in open(args.TSV):
        if header:
            header = False
            continue
        data = line.strip().split("\t")
        geneid = data[9]
        trascriptid = data[10]
        if not trascriptid in genes[geneid]:
            genes[geneid][trascriptid] = IntervalTree()
            exonscount[geneid][trascriptid] = defaultdict(int)
            exonsinfo[geneid][trascriptid] = defaultdict(dict)
            junctionscount[geneid][trascriptid] = defaultdict(int)
        exon_num = data[12]
        exonsinfo[geneid][trascriptid]["strand"] = data[4]
        if exon_num:
            trstart = int(data[13])
            trend = int(data[14])

            exonsinfo[geneid][trascriptid][exon_num]["genome"] = data[:3]
            exonsinfo[geneid][trascriptid][exon_num]["tr"] = [trstart, trend]

            genes[geneid][trascriptid][trstart:trend] = exon_num
            # print(geneid, trascriptid,exon_num,trstart, trend)

    # print(genes["FBgn0039900"]["FBgn0039900_template"][20:50])

    for record in SeqIO.parse(args.FQ, "fastq"):
        rname, mate1, mate2 = record.name.split(";")
        trascriptid = rname.split("/")[1]
        geneid = "_".join(trascriptid.split("_")[:-1])
        try:
            for m in [mate1, mate2]:
                # for m in [mate1]:
                # print(rname, mate1, mate2)
                s, e = m.split(":")[1].split("-")
                s = int(s)
                e = int(e)
                qres = sorted(genes[geneid][trascriptid][s:e])
                # print(qres, len(qres), len(qres) == 1)
                if len(qres) == 1:
                    # qres.data["count"] += 1
                    # print('oooo', qres[0].data)
                    exonscount[geneid][trascriptid][qres[0].data] += 1
                else:
                    for i, j in zip(qres, qres[1:]):
                        junctionscount[geneid][trascriptid][(i.data, j.data)] += 1

        except:
            # print("skipping", rname, file=sys.stderr)
            continue

    # pprint(exonscount)
    # pprint(exonsinfo)
    # pprint(junctionscount)

    print(
        "seqnames,start,end,strand,type,gene_id,transcript_id,gene_exon_number,tr_start,tr_end,read_count"
    )
    for geneid in exonsinfo:
        for trascriptid in exonsinfo[geneid]:
            prevex = None
            for exon in exonsinfo[geneid][trascriptid]:
                if exon == "strand":
                    continue

                if prevex:
                    _seq = exonsinfo[geneid][trascriptid][exon]["genome"][0]
                    if exonsinfo[geneid][trascriptid]["strand"] == '-':
                        _jstart = exonsinfo[geneid][trascriptid][exon]["genome"][2]
                        _jend = exonsinfo[geneid][trascriptid][prevex]["genome"][1]
                    else:
                        # This needs to be checked
                        _jstart = exonsinfo[geneid][trascriptid][prevex]["genome"][2]
                        _jend = exonsinfo[geneid][trascriptid][exon]["genome"][1]
                    print(
                        # *[0,0,0], #*exonsinfo[geneid][trascriptid][exon]["genome"],
                        _seq, _jstart, _jend,
                        exonsinfo[geneid][trascriptid]["strand"],
                        "junction",
                        geneid,
                        trascriptid,
                        f"{prevex}-{exon}",
                        ".", ".", #*exonsinfo[geneid][trascriptid][exon]["tr"],
                        junctionscount[geneid][trascriptid][(prevex, exon)],
                        sep=","
                    )
                prevex = exon

                print(
                    *exonsinfo[geneid][trascriptid][exon]["genome"],
                    exonsinfo[geneid][trascriptid]["strand"],
                    "exon",
                    geneid,
                    trascriptid,
                    exon,
                    *exonsinfo[geneid][trascriptid][exon]["tr"],
                    exonscount[geneid][trascriptid][exon],
                    sep=","
                )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Check that GFA Path are identical to input FAs"
    )
    parser.add_argument("FQ", type=str, help="Path to simulated FQ")
    parser.add_argument("TSV", type=str, help="Path to ASimulatoR exon junction tsv")
    args = parser.parse_args()

    main(args)
