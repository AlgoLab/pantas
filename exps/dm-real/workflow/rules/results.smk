rule analyze_bench:
    input:
        pjoin(ODIR, "rMATS", "summary.txt"),
        pjoin(ODIR, "whippet", "psi.diff"),
        pjoin(ODIR, "pantas2", "quant.w{w}.csv")
    output:
        pjoin(ODIR, "results", "bench_{w}.csv")
    params:
        bench = pjoin(ODIR, "bench"),
        res_dir = pjoin(ODIR, "results")
    conda: "../envs/plot.yaml"
    shell:
        """
        python3 workflow/scripts/compare_bench.py {params.bench} {params.res_dir} {wildcards.w}
        """


rule analyze_results:
    input:
        r=pjoin(ODIR, "rMATS", "summary.txt"),
        w=pjoin(ODIR, "whippet", "psi.diff"),
        p=pjoin(ODIR, "pantas2", "quant.w{w}.csv")
    output:
        pjoin(ODIR, "results", "res_{w}.csv")
    params:
        rmats = pjoin(ODIR, "rMATS"),
        res_dir = pjoin(ODIR, "results")
    conda: "../envs/plot.yaml"
    shell:
        """
        python3 workflow/scripts/parse_res.py {input.p} {params.rmats} {input.w} {params.res_dir} {wildcards.w}
        touch {output}
        """
