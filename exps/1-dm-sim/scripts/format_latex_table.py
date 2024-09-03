import sys

data = {}
for line in open(sys.argv[1]):
    if line.startswith("p-supp"):
        continue
    _, tool, etype, _, mincov, TP, FN, FP, Prec, Rec, F1, _ = line.strip("\n").split(
        ","
    )
    if mincov not in data:
        data[mincov] = {}
    if etype not in data[mincov]:
        data[mincov][etype] = []
    data[mincov][etype].append([tool, TP, FN, FP, Prec, Rec, F1])

print(
    "True Support (ω)",
    "Event Type",
    "Tool",
    "TP",
    "FN",
    "FP",
    "Precision",
    "Recall",
    "F1",
    sep=" & ",
    end=" \\\\\n",
)

latex = {
    "pantas": "\\pantas",
    "rMATS": "\\rmats",
    "Whippet": "\\whippet",
    "SUPPA2": "\\suppa",
}
for c in data:
    for e in data[c]:
        for x in data[c][e]:
            t = latex[x[0]]
            print(c, e, t, *x[1:], sep=" & ", end=" \\\\\n")
