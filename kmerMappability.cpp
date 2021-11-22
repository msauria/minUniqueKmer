#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

size_t Nchr=0;
vector<size_t> cindices;
vector<uint8_t> mul, mur;
vector<std::string> cnames;

void print_info(void)
{
  std::cout << "Find mean unique k-mer coverage for a given k-mer size at each genomic position ver. 1.0\n"
        << "\nUsage:\nmeanKmerCoverage <k-mer> <chrom_sz> <mul_wig> <out_wig>\n"
        << "Parameters:\n"
        << "<k-mer>         - Maximum unique k-mer size\n"
        << "<chrom_sz>      - Two-column tab-separated chromosome size file\n"
        << "<mul_wig>       - Left-anchored minimum unique k-mer wiggle file\n"
        << "<out_bed>       - Output k-mer mappability bed file\n";
}

bool help(int argc, char** argv)
{
  const std::string help = "--help";
  for (int i = 1; i < argc; ++i)
  {
    if (argv[i] == help)
      return true;
  }
  return false;
}

int main(int argc, char *argv[])
{
  fprintf(stderr, "\r                                     \r");
  fflush(stderr);
  if (argc <= 4 || help(argc, argv))
  {
    print_info();
    return 1;
  }

  uint8_t kmer = stoi(argv[1]);

  // Load chromosome sizes
  ifstream input(argv[2]);
  cindices.resize(0);
  cnames.resize(0);
  string cur;
  size_t pos=0;
  cindices.push_back(pos);
  fprintf(stderr, "\rReading chromosome sizes ");
  fflush(stderr);
  while (getline(input, cur))
  {
    size_t name_stop = cur.find("\t");
    pos += std::stoi(cur.substr(name_stop + 1));
    Nchr += 1;
    cnames.push_back(cur.substr(0, name_stop));
    cindices.push_back(pos);
  }
  input.close();

  ifstream input2(argv[3]);
  ofstream outfile;
  outfile.open(argv[4]);
  string field, chrom, pchrom="";
  size_t cend, split, val, bstart=0, bend=0, cindex;
  while (getline(input2, cur))
  {
    if (cur.find("f") == 0)
    {
      while (cur.length() > 0)
      {
        split = cur.find(" ");
        if (split == std::string::npos)
          split = cur.length();
        field = cur.substr(0, split);
        if (split == cur.length())
          cur.erase(0, split);
        else
          cur.erase(0, split + 1);
        if (field.find("chrom") == 0)
          chrom = field.substr(field.find("=") + 1, field.length());
          if (chrom != pchrom)
          {
            fprintf(stderr, "\r                         \r");
            std::cerr << "\r                       \r";
            fprintf(stderr, "\rProcessing %s", chrom.c_str());
            fflush(stderr);
            if (bstart != bend)
                outfile << pchrom << "\t" << to_string(bstart) << "\t" << to_string(min(bend, cend)) << "\n";
            bstart = 0;
            bend = 0;
            for (size_t j=0; j<cnames.size(); j++) {
              if (chrom == cnames[j]) {
                cindex = j;
                break;
              }
            }
            cend = cindices[cindex + 1] - cindices[cindex];
            pchrom = chrom;
          }
        if (field.find("start") == 0)
          pos = stoi(field.substr(field.find("=") + 1, field.length())) - 1;
      }
    } else {
      val = stoi(cur);
      if (pos > bend && bstart != bend)
      {
        outfile << chrom << "\t" << to_string(bstart) << "\t" << to_string(bend) << "\n";
        bstart = 0;
        bend = 0;
      }
      if (val <= kmer)
      {
        if (bstart == bend) bstart = pos;
        bend = pos + kmer;
      }
      pos++;
    }
  }
  if (bstart != bend)
    outfile << chrom << "\t" << to_string(bstart) << "\t" << to_string(min(bend, cend)) << "\n";
  input2.close();
  outfile.close();
  fprintf(stderr, "\r                         \r");
  fflush(stderr);
  return EXIT_SUCCESS;
}