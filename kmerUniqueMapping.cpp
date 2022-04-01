#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

void print_info(void)
{
  std::cout << "Find intervals which can be uniquely mapped for a given k-mer ver. 1.0\n"
        << "\nUsage:\nkmerMappability <k-mer> <coverage_wig>\n"
        << "Parameters:\n"
        << "<k-mer>         - k-mer size\n"
        << "<coverage_wig>  - K-mer coverage wiggle file\n";
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
  if (argc <= 2 || help(argc, argv))
  {
    print_info();
    return 1;
  }

  int32_t kmer_len = stoi(argv[1]);
  string cur, field, chrom="", new_chrom="";
  int32_t pos=0, new_pos, split, bufsize, start;
  ifstream wigfile(argv[2]);
  while (getline(wigfile, cur))
  {
    if (cur[0] == 'f')
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
          new_chrom = field.substr(field.find("=") + 1, field.length());
        else if (field.find("start") == 0)
          new_pos = stoi(field.substr(field.find("=") + 1, field.length())) - 1;
        if (chrom != new_chrom)
          fprintf(stderr, "\r                               \r");
      }
      if (chrom == "")
      {
        chrom = new_chrom;
        start = pos;
      }
      else if (chrom != new_chrom || pos != new_pos)
      {
        fprintf(stdout, "%s\t%i\t%i\n", chrom.c_str(), start, pos);
        fflush(stdout);
        chrom = new_chrom;
        pos = new_pos;
        start = new_pos;
      }
    }
    else
    {
      pos++;
      if (pos % 1000000 == 0)
      {
        fprintf(stderr, "\rRead %iMbp of %s", pos/1000000, chrom.c_str());
        fflush(stderr);
      }
    }
  }
  fprintf(stderr, "\r                               \r");
  fflush(stderr);
  wigfile.close();

  return EXIT_SUCCESS;
}