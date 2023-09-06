# python alignments_augmentation.py  alignment.json output.path input.gfa > output.gfa


import sys
import json


def get_sources(graph):
    in_nodes = {}
    for n1, adj in enumerate(graph):
        if n1 not in in_nodes:
            in_nodes[n1] = 0
        for n2 in adj:
            if n2 not in in_nodes:
                in_nodes[n2] = 1
            else:
                in_nodes[n2] += 1

    return [n1 for n1, count in in_nodes.items() if count == 0]


def find_paths(graph, curr_node, curr_path, paths):
    curr_path.append(curr_node)
    if len(graph[curr_node]) == 0:
        paths.append(curr_path[:])
    else:
        for child in graph[curr_node]:
            find_paths(graph, child, curr_path, paths)
    curr_path.pop()


def get_full_paths(g):
    paths = []
    curr_path = []
    # print(g)
    sources = get_sources(g)
    # print("sources", sources, file=sys.stderr)
    for s in sources:
        find_paths(g, s, curr_path, paths)
    return paths


def main(argv):
    json_file = argv[0]
    output_file = argv[1]
    gfa_file = argv[2]
    # gfa_file_out = argv[3]
    weights = {}
    revs = {}
    print("Building paths and weights", file=sys.stderr)
    with open(json_file) as f, open(output_file, "w") as out:
        for line in f:
            data = json.loads(line)
            # print("loading...")
            read_name = data["name"]
            # print(read_name, file=sys.stderr)
            sequence = data["sequence"]
            if "subpath" not in data.keys():
                continue
            paths_list = [[]] * len(data["subpath"])
            next_list = [[]] * len(data["subpath"])
            for i, subpaths in enumerate(data["subpath"]):
                if "next" in subpaths.keys():
                    next_list[i] = subpaths["next"]
                nodes_tmp = []
                for elem in subpaths["path"]["mapping"]:
                    id_node = elem["position"]["node_id"]
                    if "is_reverse" in elem["position"].keys():
                        if elem["position"]["is_reverse"]:
                            dir_node = "-"
                        else:
                            dir_node = "+"
                    else:
                        dir_node = "+"
                    node = f"{id_node}{dir_node}"
                    nodes_tmp.append(node)
                    # print(nodes_tmp)
                # print("adding ", i, nodes_tmp)
                paths_list[i] = nodes_tmp

            paths_list_full = get_full_paths(next_list)
            # read226217/FBgn0052000_template;mate1:49-198;mate2:167-316
            # print(next_list, file=sys.stderr)
            # print(paths_list, file=sys.stderr)
            # print(paths_list_full)
            paths = []
            for _path in paths_list_full:
                tmp = []
                for _p in _path:
                    tmp = tmp + paths_list[_p]
                paths.append(tmp)
            # print(read_name, paths, file=sys.stderr)
            paths_final = []
            for _path in paths:
                if _path[0][-1] == "+":
                    tmp = []
                    for p in _path:
                        tmp.append(p[:-1])
                    paths_final.append((tmp, "+"))
                else:
                    tmp = []
                    for p in _path:
                        tmp.append(p[:-1])
                    tmp.reverse()
                    paths_final.append((tmp, "-"))
            # print(paths_final, file=sys.stderr)
            # for p in paths_final:
            #    print(p, file=sys.stderr)
            #    print("------\t", file=sys.stderr)

            for p in paths_final:
                for s, t in zip(p[0], p[0][1:]):
                    # print(s, t)
                    if p[1] == "+":
                        key = (s, t)
                        revs[(s, t)] = False
                    else:
                        key = (t, s)
                        revs[(s, t)] = True
                    if key in weights.keys():
                        weights[key] = weights[key] + 1
                    else:
                        weights[key] = 1
            # for k, v in weights.items():
            #    print(k, v, file=sys.stderr)
            out.write(f">{read_name}\n")
            for p in paths_final:
                if p[1] == "+":
                    d = ">"
                else:
                    d = "<"
                path_str = d.join(p[0])
                out.write(f"{path_str}\n")
    print("Annotating GFA", file=sys.stderr)
    with open(gfa_file, "r") as f:

        for line in f:
            line = line.strip()
            if not line.startswith("L"):
                print(line)
            else:
                if len(line) == 1:
                    # print(line)
                    continue
                tokens = line.split()
                # print("line:", line)
                # print(tokens)
                # print(tokens[1], tokens[3])
                w = weights.pop((tokens[1], tokens[3]), 0)

                # print(w)
                print(f"{line}\tRC:i:{w}")
                # out.write(f"{line}\t{w}\n")

    ## TODO fix the direction
    for k, v in weights.items():
        if revs.get(k[0], k[1]) or not revs.get(k[1], k[0]):
            print(f"L\t{k[1]}\t+\t{k[0]}\t+\t*\tRC:i:{v}")
        else:
            print(f"L\t{k[0]}\t+\t{k[1]}\t+\t*\tRC:i:{v},ID:Z:N")


if __name__ == "__main__":
    main(sys.argv[1:])

