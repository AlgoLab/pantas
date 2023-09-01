import sys
import gzip

print("Input file: " + sys.argv[1] + "\nOutput file: " + sys.argv[2])
with gzip.open(sys.argv[1], "rt") as fr:
    with gzip.open(sys.argv[2], "wt") as fw:
        for i, line in enumerate(fr):
            if line.startswith("2R"):
                tokens = line.split()
                tokens[2] = str(i)
                fw.write("\t".join(tokens) + "\n")
            else:
                fw.write(line)
