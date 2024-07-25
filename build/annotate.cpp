#include <cinttypes>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <unordered_set>

#include "bdsg/packed_graph.hpp"
#include "gbwt/gbwt.h"

using namespace std;

typedef uint64_t edge_t; // 32 bits for source, 32 for target
typedef unordered_set<edge_t> edges_t;

int main(int argc, char *argv[]) {
  time_t timestamp;

  char *pg_fn = argv[1];
  char *info_fn = argv[2];
  char *hap_gbwt_fn = argv[3];
  char *tr_gbwt_fn = argv[4];

  // Open packed graph
  bdsg::PackedGraph *graph = new bdsg::PackedGraph();
  graph->deserialize(pg_fn);

  // Assuming bidirectional gbwt for reference and haplotypes (seems the case
  // with vg)

  int c = 0; // Assuming single contig (this will be run on single chromosome
             // anyway)

  // maps to link path names to path identifiers in corresponding gBWT
  map<string, int> HAPS; // reference + haplotypes
  map<string, int> HATR; // haplotype-aware transcripts
  // for each haplotype, store all transcripts it "explains"
  map<string, vector<string>> H2HAT;

  time(&timestamp);
  cerr << "\n@ " << ctime(&timestamp) << "Loading " << info_fn << ".." << endl;
  ifstream info_f(info_fn);
  string tpath, len, tidx, haps, hpath;
  while (info_f >> tpath >> len >> tidx >> haps) {
    if (tpath.compare("Name") == 0)
      continue;
    int i = 0;
    while (i < haps.size() && haps.at(i) != ',')
      ++i;
    hpath = haps.substr(0, i);
    // cerr << tpath << " " << hpath << endl;
    HAPS[hpath] = -1;
    HATR[tpath] = -1;
    H2HAT[hpath].push_back(tpath);
  }
  info_f.close();

  // Opening haplotypes gBWT
  time(&timestamp);
  cerr << "\n@ " << ctime(&timestamp) << "Loading " << hap_gbwt_fn << ".."
       << endl;
  gbwt::GBWT hap_gbwt_idx;
  ifstream hap_gbwt_s(hap_gbwt_fn, ifstream::in);
  hap_gbwt_idx.load(hap_gbwt_s);
  hap_gbwt_s.close();
  assert(hap_gbwt_idx.bidirectional());
  string contig = hap_gbwt_idx.metadata.contig(c);
  cerr << contig << endl;
  for (const auto &i : hap_gbwt_idx.metadata.pathsForContig(c)) {
    auto fpn = hap_gbwt_idx.metadata.full_path(i);
    if (fpn.sample_name.compare("_gbwt_ref") == 0)
      hpath = contig;
    else
      hpath = fpn.sample_name + "#" + to_string(fpn.haplotype) + "#" + contig;
    if (HAPS.find(hpath) == HAPS.end())
      continue;
    HAPS[hpath] = gbwt::Path::encode(i, false);
  }
  assert(HAPS.find(contig) != HAPS.end());

  // Opening haplotype-aware transcripts gBWT
  time(&timestamp);
  cerr << "\n@ " << ctime(&timestamp) << "Loading " << tr_gbwt_fn << ".."
       << endl;
  gbwt::GBWT tr_gbwt_idx;
  ifstream tr_gbwt_s(tr_gbwt_fn, ifstream::in);
  tr_gbwt_idx.load(tr_gbwt_s);
  tr_gbwt_s.close();
  assert(!tr_gbwt_idx.bidirectional());
  for (const auto &pi : tr_gbwt_idx.metadata.pathsForContig(c)) {
    tpath = tr_gbwt_idx.metadata.full_path(pi).sample_name;
    // cout << "T: " << tpath << endl;
    HATR[tpath] = pi;
  }

  map<int, string> VTAGS;
  map<edge_t, string> JTAGS;
  int hh = 1, tt = 1, n;
  uint32_t x, y;
  handlegraph::handle_t hx, hy;
  edge_t e;
  vector<string> tpaths;
  gbwt::vector_type path;
  edges_t haplotype;
  edges_t htranscript;

  // TODO: parallelize

  set<int>
      kept_ht; // set of haplotype transcripts that have been kept after pruning
  time(&timestamp);
  cerr << "\n@ " << ctime(&timestamp) << HAPS.size() << " haplotypes and "
       << HATR.size() << " haplotype-aware transcripts" << endl;
  for (const auto &h2hat : H2HAT) {
    hpath = h2hat.first;
    tpaths = h2hat.second;

    time(&timestamp);
    cerr << "\n@ " << ctime(&timestamp) << "Building haplotype " << hpath
         << " (" << hh << "/" << HAPS.size() << ")" << endl;
    ++hh;

    // build set for haplotype edges
    haplotype.clear();
    path = hap_gbwt_idx.extract(HAPS[hpath]);
    // TODO: check path is on + strand
    for (n = 0; n < path.size() - 1; ++n) {
      x = (uint32_t)gbwt::Node::id(path[n]);
      y = (uint32_t)gbwt::Node::id(path[n + 1]);
      haplotype.insert(((uint64_t)x << 32) | y);
    }
    tt = 1;
    for (const auto &t : tpaths) {
      cerr << tt << "/" << tpaths.size() << "\r" << flush;
      ++tt;
      path = tr_gbwt_idx.extract(HATR[t]);

      // CHECKME: this doesn't seem to work. Need time to better read gbwt doc
      // if(gbwt::Path::is_reverse(pi))
      if (path.size() > 1 && gbwt::Node::id(path[0]) > gbwt::Node::id(path[1]))
        gbwt::reversePath(path);

      // Check if haplotype-aware transcript has been pruned
      for (n = 0; n < path.size() - 1; ++n) {
        x = (uint32_t)gbwt::Node::id(path[n]);
        y = (uint32_t)gbwt::Node::id(path[n + 1]);
        if (!graph->has_node(x) or !graph->has_node(y))
          break;
        hx = graph->get_handle(x);
        hy = graph->get_handle(y);
        if (!graph->has_edge(hx, hy))
          break;
      }
      if (n < path.size() - 1)
        break;
      kept_ht.insert(HATR[t]);

      int en = 1; // exon number

      // first exonic vertex is tagged with 1
      x = (uint32_t)gbwt::Node::id(path[0]);
      if (VTAGS.find(x) == VTAGS.end())
        VTAGS[x] = "EX:Z:";
      VTAGS[x] +=
          t + "." + to_string(en) + ","; // TODO: improve how we store this?
      for (int n = 0; n < path.size() - 1; ++n) {
        x = (uint32_t)gbwt::Node::id(path[n]);
        y = (uint32_t)gbwt::Node::id(path[n + 1]);

        // if (y == x+1)
        //     continue;
        e = ((uint64_t)x << 32) | y;
        if (haplotype.find(e) == haplotype.end()) {
          // here we are on a junction since we do not have this edge in the
          // corresponding haplotype path
          if (JTAGS.find(e) == JTAGS.end())
            JTAGS[e] = "JN:Z:";
          ++en;
          JTAGS[e] += t + "." + to_string(en - 1) + "." + to_string(en) +
                      ","; // TODO: improve how we store this?
          // cout << x << "," << y << " is a junction (" << t + "." +
          // to_string(en-1) + "." + to_string(en) << ")" << endl;
        }

        // add tag for "next" exonic vertex
        if (VTAGS.find(y) == VTAGS.end())
          VTAGS[y] = "EX:Z:";
        VTAGS[y] +=
            t + "." + to_string(en) + ","; // TODO: improve how we store this?
      }
    }
    cerr << tpaths.size() << "/" << tpaths.size() << endl;
  }

  // dump annotated GFA
  time(&timestamp);
  cerr << "\n@ " << ctime(&timestamp) << "Dumping gfa to stdout ("
       << VTAGS.size() << " exonic vertices and " << JTAGS.size()
       << " junctions)" << endl;

  // print header
  cout << "H\tVN:Z:1.1" << endl;

  // print segments
  for (x = graph->min_node_id(); x <= graph->max_node_id(); ++x) {
    if (!graph->has_node(x))
      // move to next
      continue;
    hx = graph->get_handle(x);
    if (VTAGS.find(x) != VTAGS.end()) {
      VTAGS[x].pop_back();
      cout << "S"
           << "\t" << x << "\t" << graph->get_sequence(hx) << "\t" << VTAGS[x]
           << endl;
    } else {
      cout << "S"
           << "\t" << x << "\t" << graph->get_sequence(hx) << endl;
    }
  }

  // print links
  for (x = graph->min_node_id(); x <= graph->max_node_id(); ++x) {
    if (!graph->has_node(x))
      // move to next
      continue;
    hx = graph->get_handle(x);
    set<bdsg::handle_t> outs;
    graph->follow_edges(hx, false, [&outs](const bdsg::handle_t &hy) {
      outs.insert(hy);
      return true;
    });
    for (const bdsg::handle_t hy : outs) {
      y = (uint32_t)graph->get_id(hy);
      edge_t e = ((uint64_t)x << 32) | y;
      if (JTAGS.find(e) != JTAGS.end()) {
        JTAGS[e].pop_back();
        cout << "L"
             << "\t" << x << "\t"
             << "+"
             << "\t" << y << "\t"
             << "+"
             << "\t"
             << "*"
             << "\t" << JTAGS[e] << endl;
      } else {
        cout << "L"
             << "\t" << x << "\t"
             << "+"
             << "\t" << y << "\t"
             << "+"
             << "\t"
             << "*" << endl;
      }
    }
  }

  // print reference path
  path = hap_gbwt_idx.extract(HAPS[contig]);
  cout << "P"
         << "\t" << contig << "\t" << gbwt::Node::id(path[0]) << "+";
    for (n = 1; n < path.size(); ++n) {
      cout << "," << gbwt::Node::id(path[n]) << "+";
    }
    cout << "\t"
         << "*" << endl;

  // print transcript paths
  string sep = "+";
  for (const auto &t : kept_ht) {
    tpath = tr_gbwt_idx.metadata.full_path(t).sample_name;
    path = tr_gbwt_idx.extract(t);
    if (path.size() > 1 && gbwt::Node::id(path[0]) > gbwt::Node::id(path[1])) {
      sep = "-";
    }
    cout << "P"
         << "\t" << tpath << "\t" << gbwt::Node::id(path[0]) << sep;
    for (n = 1; n < path.size(); ++n) {
      cout << "," << gbwt::Node::id(path[n]) << sep;
    }
    cout << "\t"
         << "*" << endl;

    // CHECKME: it is not clear to me if mpmap needs both directions or just one
    // gbwt::reversePath(path);
    // sep = sep == "+" ? "-" : "+";
    // cout << "P"
    //      << "\t" << tpath << "'" << "\t" << gbwt::Node::id(path[0]) << sep;
    // for (n = 1; n < path.size(); ++n) {
    //   cout << "," << gbwt::Node::id(path[n]) << sep;
    // }
    // cout << "\t"
    //      << "*" << endl;
  }

  time(&timestamp);
  cerr << "\n@ " << ctime(&timestamp) << "Done." << endl;
  return 0;
}