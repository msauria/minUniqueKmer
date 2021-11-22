#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

size_t Nchr=0;
vector<size_t> cindices, mul, mur;
vector<std::string> cnames;

void print_info(void)
{
  std::cout << "Find minimum unique fragment size covering each genomic position ver. 1.0\n"
        << "\nUsage:\nminFragmentCoverage <chrom_sz> <mul_wig> <mur_wig> <out_wig>\n"
        << "Parameters:\n"
        << "<chrom_sz>      - Two-column tab-separated chromosome size file\n"
        << "<mul_wig>       - Left-anchored minimum unique k-mer wiggle file\n"
        << "<mur_wig>       - Right-anchored minimum unique k-mer wiggle file\n"
        << "<out_wig>       - Output minimum fragment size wiggle file\n";
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

void load_wiggle(char* fname, std::vector<size_t>* wig, bool forward)
{
  ifstream input(fname);
  string cur, field, chrom;
  size_t cstart, pos, split;
  while (getline(input, cur))
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
        else if (field.find("start") == 0)
          pos = stoi(field.substr(field.find("=") + 1, field.length())) - 1;
      }
      size_t index;
      for (index=0; index<Nchr; index++)
        if (cnames[index].compare(chrom) == 0 && cnames[index].length() == chrom.length())
          break;
      cstart = cindices[index];
      pos += cstart;
    } else {
      (*wig)[pos] = stoi(cur);
      pos++;
    }
  }
  input.close();
  for (size_t i=0; i<Nchr; i++){
    if (forward){
      for (size_t j=cindices[i+1]-1; j>cindices[i]; j--){
        if ((*wig)[j-1] > 0 && (*wig)[j-1] < (*wig)[j]){
          for (size_t k=j; k<std::min(j + (*wig)[j-1] - 1, cindices[i+1]); k++)
          {
            if ((*wig)[k] > 0) (*wig)[k] = std::min((*wig)[k], (*wig)[j-1]);
            else (*wig)[k] = (*wig)[j-1];
          }
        }
      }
    } else {
      for (size_t j=cindices[i]; j<cindices[i + 1]-1; j++){
        if ((*wig)[j+1] > 0 && (*wig)[j+1] < (*wig)[j]){
          for (size_t k=j; k>=std::max(j + 2 - (*wig)[j+1], cindices[i]); k--)
          {
            if ((*wig)[k] > 0) (*wig)[k] = std::min((*wig)[k], (*wig)[j+1]);
            else (*wig)[k] = (*wig)[j+1];
          }
        }
      }
    }
  }
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

  // Load chromosome sizes
  ifstream input(argv[1]);
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

  // Load wiggle files
  fprintf(stderr, "\r                         \r");
  fprintf(stderr, "\rAllocating MUL array ");
  fflush(stderr);
  mul.resize(cindices[Nchr]);
  fprintf(stderr, "\r                         \r");
  fprintf(stderr, "\rLoading MUL array ");
  fflush(stderr);
  load_wiggle(argv[2], &mul, true);
  fprintf(stderr, "\r                         \r");
  fprintf(stderr, "\rAllocating MUR array ");
  fflush(stderr);
  mur.resize(cindices[Nchr]);
  fprintf(stderr, "\r                         \r");
  fprintf(stderr, "\rLoading MUR array ");
  fflush(stderr);
  load_wiggle(argv[3], &mur, false);

  // Write minimum unique fragment wiggle
  fprintf(stderr, "\r                         \r");
  fprintf(stderr, "\rWriting output ");
  fflush(stderr);
  ofstream outfile;
  outfile.open (argv[4]);
  bool prev;
  size_t minval;
  for (size_t i=0; i<Nchr; i++)
  {
    prev = false;
    for (size_t j=cindices[i]; j<cindices[i+1]; j++)
    {
      if (mul[j] == 0)
      {
        if (mur[j] == 0)
        {
          prev = false;
          continue;
        } else
          minval = mur[j];
      } else if (mur[j] == 0) {
        minval = mul[j];
      } else {
        minval = std::min(mul[j], mur[j]);
      }
      if (!prev) outfile << "fixedStep chrom=" + cnames[i] + " start=" +
                            to_string(j - cindices[i] + 1) + " step=1\n";
      outfile << to_string(minval) + "\n";
      prev = true;
    }
  }
  outfile.close();
  fprintf(stderr, "\r                         \r");
  fflush(stderr);
  return EXIT_SUCCESS;
}