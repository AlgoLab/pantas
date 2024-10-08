# python alignments_augmentation_from_gaf.py alignment.gaf input.gfa 20 > output.gfa

import sys
import re

# import pickle
# import os


def parse_cigar(cigar):
    patterns = {
        "=": r"[ACGTN]+",
        ":": r"[0-9]+",
        "*": r"[acgtn][acgtn]",
        "+": r"[acgtn]+",
        "-": r"[acgtn]+",
        "~": r"[acgtn]{2}[0-9]+[acgtn]{2}",
    }
    results = []
    regex = re.compile("(=|:|\*|\+|\-|\~)")
    tokens = regex.split(cigar)
    curr_op = None
    for token in tokens:
        if token in patterns:
            curr_op = token
        else:
            if curr_op:
                if curr_op == "*":
                    results.append((curr_op, 1))
                else:
                    if token.isdigit():
                        results.append((curr_op, int(token)))
                    else:
                        results.append((curr_op, len(token)))
                curr_op = None

    return results


def cigar_clipping(cigar_list, start_pos, end_pos):
    cigar_list_new = cigar_list
    new_start = start_pos
    new_end = end_pos
    if cigar_list[0][0] == "+" and cigar_list[1][0] == ":":
        cigar_list_new = [cigar_list[1]]
        new_start = start_pos + cigar_list[0][1]
    elif cigar_list[0][0] == ":" and cigar_list[1][0] == "+":
        cigar_list_new = [cigar_list[0]]

    return cigar_list_new, new_start, new_end


def no_good_cigar(cigar_list, c_thr=3):
    no_good = False
    for c in cigar_list:
        if (c[0] == "+" and c[1] > c_thr) or (c[0] == "-" and c[1] > c_thr):
            no_good = True
            break
    return no_good


## TODO add handle of multiple * at begin
def compact_align(align):
    new_al = []
    for i, a in enumerate(align[1]):
        s_e = False
        if i == 0:
            if a[0] == "*":
                s_e = True
            else:
                new_al.append(a)
        else:
            if len(new_al) == 0:
                # if s_e:
                last_op = a[0]
                new_len = a[1] + 1
                new_al.append((last_op, new_len))
                s_e = False
                # else:
                #     new_al.append(a)
                continue
            if a[0] == new_al[-1][0] or (a[0] == "*"):
                # print(new_al[-1][1], a[1], file=sys.stderr)
                last_op = new_al[-1][0]
                new_len = new_al[-1][1] + a[1]
                if s_e:
                    new_len = new_len + 1
                    s_e = False
                new_al.pop()
                new_al.append((last_op, new_len))
                # new_al[-1][1] = new_al[-1][1] + a[1]
            else:
                new_al.append(a)
    return (align[0], new_al)


def clear_align(align):
    final_align = []
    for i, al in enumerate(align):
        # al = rem_edit_align(al)
        if len(al[1]) == 1 and (al[1][0][0] == "-" or al[1][0][0] == "+"):
            continue
        else:
            ## TODO I can do it when i create align
            cal = compact_align(al)
            final_align.append(cal)
    return final_align


def main(argv):
    gaf_file = argv[0]
    gfa_file = argv[1]
    thr = int(argv[2]) if len(argv) > 2 else 20
    # gfa_file_out = argv[3]
    weights = {}
    nodes_weights = {}
    revs = {}
    nodes_info = {}
    rej = 0
    print("Read GFA", file=sys.stderr)
    with open(gfa_file, "r") as f:
        for line in f:
            if not line.startswith("S"):
                continue
            tokens = line.strip().split()
            nodes_info[tokens[1]] = (len(tokens[2]), [{}, {}])

    # with open(f"test_dump.ser", "wb") as f:
    #     pickle.dump(nodes_info, f)

    # with open(f"test_dump.ser", "rb") as f:
    #     nodes_info = pickle.load(f)

    print("Augmentation by GAF alignments", file=sys.stderr)
    with open(gaf_file) as f:
        last_read = ""
        count = 0
        for line in f:
            # print(line)
            # if count == 100000:
            #    sys.exit()
            tokens = line.strip().split()
            score = int(tokens[11])
            if score < thr:
                rej = rej + 1
                continue
            if tokens[5] == "*":
                continue
            read_name = tokens[0]
            path = tokens[5]
            path_len = int(tokens[6])
            start_pos = int(tokens[7])
            end_pos_rel = path_len - int(tokens[8])
            cigar_re = re.search(
                r"cs:.*?(?=\s|$)", " ".join(item for item in tokens[12:])
            )
            if cigar_re:
                cigar = cigar_re.group(0).replace("cs:Z:", "")
            else:
                cigar = "*"

            cigar_vals = parse_cigar(cigar)
            # handle clipping
            if len(cigar_vals) == 2:
                cigar_vals, start_pos, end_pos_rel = cigar_clipping(
                    cigar_vals, start_pos, end_pos_rel
                )

            # if no_good_cigar(cigar_vals):
            #     rej = rej + 1
            #     continue
            dv_re = re.search(
                r"dv:f:(\d+(\.\d+)?)", " ".join(item for item in tokens[12:])
            )
            if dv_re:
                dv = dv_re.group(0).replace("dv:f:", "")
            else:
                dv = "*"
            if float(dv) > 0.1:
                continue
            # print("score", score, cigar)
            rev = False
            # if score == 0:
            #    continue
            nodes = []
            if path[0] == ">":
                for n in path.split(">")[1:]:
                    if len(nodes) == 0 or nodes[-1] != n:
                        nodes.append(n)
            else:
                for n in path.split("<")[1:]:
                    if len(nodes) == 0 or nodes[-1] != n:
                        nodes.append(n)
                rev = True
                # nodes.reverse()
            # print(nodes)
            assert len(nodes) > 0

            # print(cigar_vals)

            # print(cigar_vals)
            # print(line)
            # print("score", score, cigar, cigar_vals, nodes)
            # print("start", start_pos, "end", end_pos_rel)
            align = []
            # print(nodes)
            # for n in nodes:
            # print(n, nodes_info[n][0])
            n_nodes = len(nodes)

            for i, id_node in enumerate(nodes):
                # print(nodes_info[n])
                # id_node = nodes_info[n][0]
                len_seq_node = nodes_info[id_node][0]
                if i == 0:
                    len_seq_node = len_seq_node - start_pos
                if i == n_nodes - 1:
                    len_seq_node = len_seq_node - end_pos_rel + 1
                tmp_len = len_seq_node
                first = True
                # print("data", id_node, len_seq_node)
                # print("cigar", cigar_vals)
                # print("align", align)
                while tmp_len > 0:
                    # print(tmp_len, cigar_vals, align, file=sys.stderr)

                    curr_cigar_op = cigar_vals[0][0]
                    if curr_cigar_op == ":":
                        curr_cigar_len = cigar_vals[0][1]
                    elif curr_cigar_op == "*":
                        curr_cigar_len = cigar_vals[0][1]
                    elif curr_cigar_op in ["=", "-", "+"]:
                        curr_cigar_len = cigar_vals[0][1]
                    if tmp_len <= curr_cigar_len:
                        if first:
                            align.append((id_node, [(curr_cigar_op, tmp_len)]))
                            first = False
                        else:
                            align[-1][1].append((curr_cigar_op, tmp_len))
                        cigar_vals[0] = (curr_cigar_op, curr_cigar_len - tmp_len)
                        if cigar_vals[0][1] == 0:
                            cigar_vals.pop(0)
                        tmp_len = 0
                    else:
                        if first:
                            align.append((id_node, [(curr_cigar_op, curr_cigar_len)]))
                            first = False
                        else:
                            align[-1][1].append((curr_cigar_op, curr_cigar_len))
                        tmp_len = tmp_len - curr_cigar_len
                        cigar_vals.pop(0)
                    if len(cigar_vals) == 0:
                        if len(align) != len(nodes):
                            pass # print("warning", file=sys.stderr)
                        break

            final_align = clear_align(align)
            # if align != final_align:
            #     print("prer", align, file=sys.stderr)
            #     print("clear", final_align, file=sys.stderr)

            #     print(file=sys.stderr)
            clear_nodes = []
            for n in final_align:
                clear_nodes.append(n[0])
                if n[0] not in nodes_weights.keys():
                    nodes_weights[n[0]] = 1
                else:
                    nodes_weights[n[0]] = nodes_weights[n[0]] + 1
            # if stop:
            #    print(align, final_align, file=sys.stderr)
            #    sys.exit()
            # print("final", final_align, file=sys.stderr)
            ## TODO check semantic of "+" and "-" in cigar and check if all
            ## cigar elements are parsed
            for i, elem in enumerate(final_align):
                node_id = elem[0]
                cigar_values = elem[1]
                for j, c in enumerate(cigar_values):
                    if not rev:
                        if c[0] == "-":
                            if i != 0 and j == 0:
                                seq_len = c[1]
                                if seq_len in nodes_info[node_id][1][0].keys():
                                    nodes_info[node_id][1][0][seq_len] = (
                                        nodes_info[node_id][1][0][seq_len] + 1
                                    )
                                else:
                                    nodes_info[node_id][1][0][seq_len] = 1
                            if i != len(final_align) - 1 and j == len(cigar_values) - 1:
                                seq_len = nodes_info[node_id][0] - c[1] - 1
                                if seq_len in nodes_info[node_id][1][1].keys():
                                    nodes_info[node_id][1][1][seq_len] = (
                                        nodes_info[node_id][1][1][seq_len] + 1
                                    )
                                else:
                                    nodes_info[node_id][1][1][seq_len] = 1
                        elif c[0] != "*":
                            if i != 0:
                                if 0 in nodes_info[node_id][1][0].keys():
                                    nodes_info[node_id][1][0][0] = (
                                        nodes_info[node_id][1][0][0] + 1
                                    )
                                else:
                                    nodes_info[node_id][1][0][0] = 1
                            if i != len(final_align) - 1:
                                seq_len = nodes_info[node_id][0] 
                                if seq_len in nodes_info[node_id][1][1].keys():
                                    nodes_info[node_id][1][1][seq_len] = (
                                        nodes_info[node_id][1][1][seq_len] + 1
                                    )
                                else:
                                    nodes_info[node_id][1][1][seq_len] = 1
                        else:
                            continue
                    else:
                        if c[0] == "-":
                            if i != 0 and j == 0:
                                seq_len = nodes_info[node_id][0] - 1 - c[1]
                                if seq_len in nodes_info[node_id][1][1].keys():
                                    nodes_info[node_id][1][1][seq_len] = (
                                        nodes_info[node_id][1][1][seq_len] + 1
                                    )
                                else:
                                    nodes_info[node_id][1][1][seq_len] = 1
                            if i != len(final_align) - 1 and j == len(cigar_values) - 1:
                                seq_len = c[1]
                                if seq_len in nodes_info[node_id][1][0].keys():
                                    nodes_info[node_id][1][0][seq_len] = (
                                        nodes_info[node_id][1][0][seq_len] + 1
                                    )
                                else:
                                    nodes_info[node_id][1][0][seq_len] = 1

                        elif c[0] != "*":
                            if i != len(final_align) - 1:
                                if 0 in nodes_info[node_id][1][0].keys():
                                    nodes_info[node_id][1][0][0] = (
                                        nodes_info[node_id][1][0][0] + 1
                                    )
                                else:
                                    nodes_info[node_id][1][0][0] = 1

                            if i != 0:
                                seq_len = nodes_info[node_id][0] 
                                if seq_len in nodes_info[node_id][1][1].keys():
                                    nodes_info[node_id][1][1][seq_len] = (
                                        nodes_info[node_id][1][1][seq_len] + 1
                                    )
                                else:
                                    nodes_info[node_id][1][1][seq_len] = 1
                        else:
                            continue
                    # if node_id == "310" and seq_len == 29:
                    #    print(line, align, file=sys.stderr)

            for n1, n2 in zip(clear_nodes, clear_nodes[1:]):
                if rev:
                    (n1, n2) = (n2, n1)
                if (n1, n2) in weights.keys():
                    weights[(n1, n2)] = weights[(n1, n2)] + 1
                else:
                    weights[(n1, n2)] = 1
            # if read_name == last_read:
            #     out.write(f"{path}\t{cigar}\n")
            # else:
            #     out.write(f"R:{read_name}\n{path}\t{cigar}\n")
            last_read = read_name
            # if count == 10:
            # break
            count = count + 1
    # for n, d in nodes_info.items():
    #    if d[1] != [{}, {}]:
    #        print(n, d)
    print(f"Rejected alignments: {rej}", file=sys.stderr)
    print("Annotating GFA", file=sys.stderr)
    with open(gfa_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("S"):
                tokens = line.split()
                node_id = tokens[1]
                node_w = 0
                if node_id in nodes_weights.keys():
                    node_w = nodes_weights[node_id]
                len_seq = nodes_info[node_id][0]
                in_w = nodes_info[node_id][1][0]
                in_l = []
                out_w = nodes_info[node_id][1][1]
                out_l = []
                for s, c in in_w.items():
                    in_l.append(f"{s}.{c}")
                for s, c in out_w.items():
                    out_l.append(f"{s}.{c}")
                if len(tokens) == 3:
                    if len(in_l) != 0 and len(out_l):
                        print(
                            f"{line}\tNC:i:{node_w}\tIL:Z:{','.join(in_l)}\tOL:Z:{','.join(out_l)}"
                        )
                    elif len(in_l) != 0 and len(out_l) == 0:
                        print(f"{line}\tNC:i:{node_w}\tIL:Z:{','.join(in_l)}")
                    elif len(in_l) == 0 and len(out_l) != 0:
                        print(f"{line}\tNC:i:{node_w}\tOL:Z:{','.join(out_l)}")
                    else:
                        print(f"{line}\tNC:i:{node_w}")
                elif len(tokens) > 3:
                    if len(in_l) != 0 and len(out_l):
                        print(
                            f"{line}\tNC:i:{node_w}\tIL:Z:{','.join(in_l)}\tOL:Z:{','.join(out_l)}"
                        )
                    elif len(in_l) != 0 and len(out_l) == 0:
                        print(f"{line}\tNC:i:{node_w}\tIL:Z:{','.join(in_l)}")
                    elif len(in_l) == 0 and len(out_l) != 0:
                        print(f"{line}\tNC:i:{node_w}\tOL:Z:{','.join(out_l)}")
                    else:
                        print(f"{line}\tNC:i:{node_w}")
            elif line.startswith("L"):
                if len(line) == 1:
                    continue
                tokens = line.split()
                w = weights.pop((tokens[1], tokens[3]), 0)
                print(f"{line}\tRC:i:{w}")
            else:
                print(line)

    for k, v in weights.items():
        print(f"L\t{k[0]}\t+\t{k[1]}\t+\t*\tRC:i:{v}\tID:Z:N")


if __name__ == "__main__":
    main(sys.argv[1:])

