#include <cinttypes>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <string.h>
#include <string>
#include <unordered_set>

#include "gbwt/gbwt.h"

using namespace std;

typedef uint64_t edge_t; // 32 bits for source, 32 for target
typedef unordered_set<edge_t> edges_t;

int main(int argc, char *argv[]) {
  time_t timestamp;

  char *gfa_fn = argv[1];
  char *info_fn = argv[2];
  char *hap_gbwt_fn = argv[3];
  char *tr_gbwt_fn = argv[4];

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
  int hh = 1, tt = 1;
  uint32_t x, y;
  edge_t e;
  vector<string> tpaths;
  gbwt::vector_type path;
  edges_t haplotype;
  edges_t htranscript;

  // TODO: parallelize

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
    for (int n = 0; n < path.size() - 1; ++n) {
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
  ifstream gfa_f(gfa_fn);
  string line;
  int p, pp;
  while (getline(gfa_f, line)) {
    if (line[0] == 'S') {
      p = 2;
      while (line[p] != '\t')
        ++p;
      x = (uint32_t)stoi(line.substr(2, p - 2));
      if (VTAGS.find(x) != VTAGS.end()) {
        VTAGS[x].pop_back();
        cout << line << "\t" << VTAGS[x] << endl;
      } else {
        cout << line << endl;
      }
    } else if (line[0] == 'L') {
      p = 2;
      while (line[p] != '\t')
        ++p;
      x = (uint32_t)stoi(line.substr(2, p - 2));
      p += 3;
      pp = p;
      while (line[pp] != '\t')
        ++pp;
      y = (uint32_t)stoi(line.substr(p, pp - p));
      edge_t e = ((uint64_t)x << 32) | y;
      if (JTAGS.find(e) != JTAGS.end()) {
        JTAGS[e].pop_back();
        cout << line << "\t" << JTAGS[e] << endl;
      } else {
        cout << line << endl;
      }
    } else {
      cout << line << endl;
    }
  }
  gfa_f.close();

  time(&timestamp);
  cerr << "\n@ " << ctime(&timestamp) << "Done." << endl;
  return 0;
}