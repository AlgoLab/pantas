import sys
import gzip

i = 1
print("Input file: " + sys.argv[1] + "\nOutput file: stdout", file=sys.stderr)
with gzip.open(sys.argv[1], "rt") as fr:
    for line in fr:
        if line.startswith("#"):
            print(line, end="")
        else:
            tokens = line.split()
            tokens[2] = tokens[0] + "." + str(i)
            print("\t".join(tokens))
            i += 1
