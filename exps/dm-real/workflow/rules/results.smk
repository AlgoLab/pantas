rule analyze_bench:
    input:
        pjoin(ODIR, "rMATS", "summary.csv"),
        pjoin(ODIR, "whippet", "psi.diff"),
        pjoin(ODIR, "salmon_suppa", "suppa.csv"),
        expand(pjoin(ODIR, "pantas", "quant.w{w}.csv"), w=Ws)
    output:
        pjoin(ODIR, "results", "bench.csv")
    params:
        bench_dir = pjoin(ODIR, "bench"),
        res_dir = pjoin(ODIR, "results"),
        Ws = Ws
    conda: "../envs/plot.yaml"
    script: "../scripts/compare_bench.py"


rule analyze_results:
    input:
        r=pjoin(ODIR, "rMATS", "summary.csv"),
        w=pjoin(ODIR, "whippet", "psi.diff"),
        s=pjoin(ODIR, "salmon_suppa", "suppa.csv"),
        q=expand(pjoin(ODIR, "pantas", "quant.w{w}.csv"), w=Ws)
    output:
        pjoin(ODIR, "results", "res.csv")
    params:
        pantas_dir = pjoin(ODIR, "pantas"),
        res_dir = pjoin(ODIR, "results"),
        Ws = Ws,
        p_value = p_value,
        min_dpsi = min_dpsi,
        min_prob = min_prob,
        min_coverage = min_coverage,
        relax = relax
    conda: "../envs/plot.yaml"
    script: "../scripts/parse_res.py"
