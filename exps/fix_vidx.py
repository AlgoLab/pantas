import sys
import gzip

from pysam import VariantFile


def main():
    in_vcf = sys.argv[1] if len(sys.argv) > 1 else sys.stdin
    vcf = VariantFile(sys.stdin)

    for line in str(vcf.header).split("\n"):
        if line == "":
            continue
        if line.startswith("##contig="):
            line = line.replace("chr", "")
        print(line)

    last_pos = -1
    i = 1
    for record in vcf:
        if record.pos != last_pos:
            i = 1
        record.id = (
            record.contig.replace("chr", "") + "-" + str(record.pos) + "." + str(i)
        )
        if record.contig.startswith("chr"):
            print(str(record)[3:], end="")


if __name__ == "__main__":
    main()
