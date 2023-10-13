rule analyze_bench:
    input:
        pjoin(ODIR, "rMATS", "summary.txt"),
        pjoin(ODIR, "whippet", "psi.diff"),
        expand(pjoin(ODIR, "pantas2", "quant.w{w}.csv"), w=Ws)
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
        r=pjoin(ODIR, "rMATS", "summary.txt"),
        w=pjoin(ODIR, "whippet", "psi.diff"),
        q=expand(pjoin(ODIR, "pantas2", "quant.w{w}.csv"), w=Ws)
    output:
        pjoin(ODIR, "results", "res.csv")
    params:
        pantas_dir = pjoin(ODIR, "pantas2"),
        rmats_dir = pjoin(ODIR, "rMATS"),
        res_dir = pjoin(ODIR, "results"),
        Ws = Ws,
        p_value = p_value,
        min_dpsi = min_dpsi,
        min_prob = min_prob,
        relax = relax
    conda: "../envs/plot.yaml"
    script: "../scripts/parse_res.py"
