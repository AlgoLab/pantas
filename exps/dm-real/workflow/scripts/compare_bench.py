#!/usr/bin/env python3

#!/usr/bin/env python3
import sys
import os
import numpy as np
import pandas as pd


def parse_time_verbose(time_file, tool):
    res = {}
    res["tool"] = [tool]
    res["file_name"] = [time_file.split("/")[-1].split(".")[0]]
    res["w"] = ["none"]
    for line in open(time_file):
        line = line[1:-1]
        tokens = line.split(sep=":")
        if tokens[0] == "User time (seconds)":
            res["user_time"] = [float(tokens[1].lstrip())]
        if tokens[0] == "System time (seconds)":
            res["sys_time"] = [float(tokens[1].lstrip())]
        if tokens[0] == "Maximum resident set size (kbytes)":
            res["max_mem"] = [int(tokens[1].lstrip())]
        if tokens[0] == "Elapsed (wall clock) time (h":
            tot = 0.0
            for x in tokens[4:]:
                tot = tot * 60 + float(x.lstrip())
                res["wall_clock"] = [tot]
    return res


def main(argv):
    bench_dir = snakemake.params.bench_dir
    if bench_dir.endswith("/"):
        bench_dir = bench_dir[:-1]
    output_dir = snakemake.params.res_dir
    if output_dir.endswith("/"):
        output_dir = output_dir[:-1]

    # w = argv[2]
    Ws = snakemake.params.Ws
    output = f"{output_dir}"

    df = pd.DataFrame(
        columns=["tool", "file_name", "user_time", "sys_time", "wall_clock", "max_mem"]
    )

    rmats_bench_file = f"{bench_dir}/rmats.time"

    rmats_time = parse_time_verbose(rmats_bench_file, "rmats + STAR")
    df = pd.concat([df, pd.DataFrame(rmats_time)], ignore_index=True, axis=0)

    whippet_bench_dir = f"{bench_dir}/whippet"

    whippet_bench_delta_file = f"{whippet_bench_dir}/delta.time"
    whippet_bench_index_file = f"{whippet_bench_dir}/index.time"

    whippet_bench_quant_files = []
    for file in os.listdir(whippet_bench_dir):
        filename = os.fsdecode(file)
        if filename.endswith(".time") and filename.startswith("quant"):
            whippet_bench_quant_files.append(f"{whippet_bench_dir}/{filename}")
    whippet_bench_quant_files.sort()

    whippet_delta_time = parse_time_verbose(whippet_bench_delta_file, "whippet")
    whippet_index_time = parse_time_verbose(whippet_bench_index_file, "whippet")
    whippet_quant_times = [
        parse_time_verbose(x, "whippet") for x in whippet_bench_quant_files
    ]

    df = pd.concat([df, pd.DataFrame(whippet_index_time)], ignore_index=True, axis=0)
    df = pd.concat([df, pd.DataFrame(whippet_delta_time)], ignore_index=True, axis=0)
    for q in whippet_quant_times:
        df = pd.concat([df, pd.DataFrame(q)], ignore_index=True, axis=0)

    star_bench_dir = f"{bench_dir}/STAR"

    star_bench_index_file = f"{star_bench_dir}/index.time"

    star_bench_map_files = []
    for file in os.listdir(star_bench_dir):
        filename = os.fsdecode(file)
        if filename.endswith(".time") and filename.startswith("map"):
            star_bench_map_files.append(f"{star_bench_dir}/{filename}")
    star_bench_map_files.sort()

    star_index_time = parse_time_verbose(star_bench_index_file, "rmats + STAR")
    star_map_times = [
        parse_time_verbose(x, "rmats + STAR") for x in star_bench_map_files
    ]

    df = pd.concat([df, pd.DataFrame(star_index_time)], ignore_index=True, axis=0)

    for m in star_map_times:
        df = pd.concat([df, pd.DataFrame(m)], ignore_index=True, axis=0)

    pantas_bench_dir = f"{bench_dir}/pantas2"

    pantas_bench_index_file = f"{pantas_bench_dir}/index.time"

    pantas_bench_mpmap_files = []
    for file in os.listdir(pantas_bench_dir):
        filename = os.fsdecode(file)
        if filename.endswith(".time") and filename.startswith("mpmap"):
            pantas_bench_mpmap_files.append(f"{pantas_bench_dir}/{filename}")
    pantas_bench_mpmap_files.sort()

    pantas_bench_weight_files = []
    for file in os.listdir(pantas_bench_dir):
        filename = os.fsdecode(file)
        if filename.endswith(".time") and filename.startswith("weight"):
            pantas_bench_weight_files.append(f"{pantas_bench_dir}/{filename}")
    pantas_bench_weight_files.sort()

    pantas_bench_call_files = {}
    for file in os.listdir(pantas_bench_dir):
        filename = os.fsdecode(file)
        for w in Ws:
            if filename.endswith(f"{w}.time") and filename.startswith("call"):
                # w = filename.split(".")[1][1:]
                if w not in pantas_bench_call_files.keys():
                    pantas_bench_call_files[w] = [f"{pantas_bench_dir}/{filename}"]
                else:
                    pantas_bench_call_files[w].append(f"{pantas_bench_dir}/{filename}")
    for w in pantas_bench_call_files.keys():
        pantas_bench_call_files[w].sort()

    pantas_bench_quant_files = {}
    for file in os.listdir(pantas_bench_dir):
        filename = os.fsdecode(file)
        for w in Ws:
            if filename.endswith(f"{w}.time") and filename.startswith("quant"):
                # w = filename.split(".")[1][1:]
                if w not in pantas_bench_quant_files.keys():
                    pantas_bench_quant_files[w] = [f"{pantas_bench_dir}/{filename}"]
                else:
                    pantas_bench_quant_files[w].append(f"{pantas_bench_dir}/{filename}")
    for w in pantas_bench_quant_files.keys():
        pantas_bench_quant_files[w].sort()

    pantas_index_time = parse_time_verbose(pantas_bench_index_file, "pantas2")
    pantas_mpmap_times = [
        parse_time_verbose(x, "pantas2") for x in pantas_bench_mpmap_files
    ]
    pantas_weight_times = [
        parse_time_verbose(x, "pantas2") for x in pantas_bench_weight_files
    ]
    pantas_call_times = {}
    for w in pantas_bench_call_files.keys():
        pantas_call_times[w] = [
            parse_time_verbose(x, "pantas2") for x in pantas_bench_call_files[w]
        ]
        for x in pantas_call_times[w]:
            x["w"] = [w]
    pantas_quant_times = {}
    for w in pantas_bench_quant_files.keys():
        pantas_quant_times[w] = [
            parse_time_verbose(x, "pantas2") for x in pantas_bench_quant_files[w]
        ]
        for x in pantas_quant_times[w]:
            x["w"] = [w]

    # print(pantas_call_times)
    df = pd.concat([df, pd.DataFrame(pantas_index_time)], ignore_index=True, axis=0)
    for m in pantas_mpmap_times:
        df = pd.concat([df, pd.DataFrame(m)], ignore_index=True, axis=0)
    for w in pantas_weight_times:
        df = pd.concat([df, pd.DataFrame(w)], ignore_index=True, axis=0)
    for w in pantas_call_times.keys():
        for k in pantas_call_times[w]:
            df = pd.concat(
                [df, pd.DataFrame(k)],
                ignore_index=True,
                axis=0,
            )
    for w in pantas_quant_times.keys():
        for k in pantas_quant_times[w]:
            df = pd.concat(
                [df, pd.DataFrame(k)],
                ignore_index=True,
                axis=0,
            )

    df.to_csv(f"{output}/bench.csv", index=False)


if __name__ == "__main__":
    main(sys.argv[1:])
