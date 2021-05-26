using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

size_t size, size2, cindex;
vector<size_t> inv, idx, mu, cindices;
vector<char> str;
vector<std::string> cnames;

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}

vector<size_t> getLCP() 
{
  fprintf(stderr, "Inverting suffix array\n");
  fflush(stderr);
  inv.resize(size2);
  for(size_t i=0; i<size2; i++)
  {
    if (i % 1000000 == 0)
    {
      fprintf(stderr, "\rInv SA %luM of %luM", i/1000000, size2/1000000);
      fflush(stderr);
    }
    inv[idx[i]] = i;
    // inv = ref[i] -> sa
    // idx = sa[i] -> ref
  }

  fprintf(stderr, "\r                                     \r");
  fprintf(stderr, "Calculating LCP\n");
  fflush(stderr);
  vector<size_t> lcp(size2, 0);
  size_t curr=0;
  for (size_t i=0; i<size2; i++) 
  {
    size_t k = inv[i];
    if (k > 0)
    {
      size_t j = idx[k - 1];
      if (i % 1000000 == 0)
      {
        fprintf(stderr, "\rLCP %luM of %luM", i/1000000, size2/1000000);
        fflush(stderr);
      }
      //printf("j %d %c\n", j, str[j]);
      while (i + curr < size2 && j + curr < size2 && str[i + curr] == str[j + curr])
      {
        curr++;
      }
      //printf("curr %d\n", curr);
      lcp[k] = curr;
    }
    if (curr > 0)
    {
      curr--;
    }
  }
  fprintf(stderr, "\r                                     \r");
  fflush(stderr);
  return lcp;
}

void print_info(void)
{
  std::cout << "Find minimum unique k-mer bedgraphs ver. 1.0\n"
        << "\nUsage:\nminUniqueKmer <genome_fa> <suffix_array> <mul_bg> <mur_bg>\n"
        << "Parameters:\n"
        << "<genome_fa>    - Genome multi-fasta sequence file\n"
        << "<suffix_array> - Suffix array produced from reverse-complimented genome\n"
        << "<mul_bg>       - Output name for left-anchored minimum unique k-mer bedgraph\n"
        << "<mur_bg>       - Output name for right-anchored minimum unique k-mer bedgraph\n";
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

  // Get total genome size
  FILE* fp = fopen(argv[2], "rb");
  size2 = filesize(argv[2]) / sizeof(size_t);
  size = size2 / 2;
  cout << size << " " << size2 << endl;

  // Load chromosome sequences and sizes
  ifstream input(argv[1]);
  cindices.resize(0);
  cnames.resize(0);
  string cur;
  size_t pos=0;
  str.resize(size2);
  fprintf(stderr, "Reading reference file\n");
  fflush(stderr);
  while (getline(input, cur))
  {
    if(cur[0] == '>')
    {
      cnames.push_back(cur.substr(1));
      cindices.push_back(pos);
    }
    else
    {
      for(size_t i = 0; i<cur.length(); i++)
      {
        if(cur[i] >= 'a' && cur[i] <= 'z') cur[i] += 'A' - 'a';
        if(cur[i] >= 'A' && cur[i] <= 'Z')
        {
          str[pos] = cur[i];
          pos++;
        }
      }
    }
  }
  cindices.push_back(pos);
  if(pos != size)
  {
    fprintf(stderr, "Reference and SA don't match lengths (%lu %lu)\n", pos, size);
    perror(NULL);
    exit(EXIT_FAILURE);
  }
  size_t Nchr = cindices.size() - 1;

  // Add reverse-complimented sequence
  fprintf(stderr, "Reverse complimenting reference\n");
  fflush(stderr);
  for (size_t i=0; i<size; i++)
  {
    if (str[i] == 'A') str[size2 - i - 1] = 'T';
    else if (str[i] == 'T') str[size2 - i - 1] = 'A';
    else if (str[i] == 'C') str[size2 - i - 1] = 'G';
    else if (str[i] == 'G') str[size2 - i - 1] = 'C';
    else str[size2 - i - 1] = 'N';
  }

  // Load suffix array
  fprintf(stderr, "Loading suffix array\n");
  fflush(stderr);
  idx.resize(size2);
  fread(&idx[0], sizeof(size_t), size2, fp);
  fclose(fp);

  // Find LCP
  vector<size_t> lcp = getLCP();

  // Find minimum unique kmer length
  mu.resize(size2);
  fprintf(stderr, "Calculating MUL and MUR\n");
  fflush(stderr);
  for (size_t i=0; i<size2-1; i++)
  {
    if (i % 1000000 == 0)
    {
      fprintf(stderr, "\rMinUnique %luM of %luM", i/1000000, size2/1000000);
      fflush(stderr);
    }
    mu[idx[i]] = 1 + max(lcp[i], lcp[i + 1]);
  }
  mu[idx[size2-1]] = lcp[size2-1] + 1;
  fprintf(stderr, "\r                                     \r");
  fflush(stderr);

  // Trim MUs
  fprintf(stderr, "Trimming MUL and MUR\n");
  fflush(stderr);
  for (size_t i=0; i<Nchr; i++)
  {
    // Remove MUL chromosome end overlaps
    size_t j = cindices[i+1] - 1, k;
    while (true)
    {
      if (j < cindices[i]) break;
      if (mu[j] + j <= cindices[i+1]) break;
      mu[j] = 0;
      j--;
    }

    // Remove MUR chromosome end overlaps
    j = size2 - cindices[i] - 1;
    while (true)
    {
      if (j < size2 - cindices[i + 1]) break;
      if (mu[j] + j <= size2 - cindices[i]) break;
      mu[j] = 0;
      j--;
    }

    // Remove MUL N overlaps
    for (j=cindices[i]; j<cindices[i + 1]; j++)
    {
      if (str[j] == 'N')
      {
        mu[j] = 0;
        k = j;
        while (true)
        {
          if (k <= cindices[i]) break;
          k--;
          if (str[k] == 'N') break;
          if (mu[k] + k <= j) break;
          mu[k] = 0;
        }
      }
    }

    // Remove MUR N overlaps
    for (j=cindices[i]; j<cindices[i + 1]; j++)
    {
      if (str[j] != 'N') continue;
      mu[size2 - j - 1] = 0;
      k = j + 1;
      while (true)
      {
        if (k >= cindices[i + 1]) break;
        if (str[k] == 'N') break;
        if (k - mu[size2 - k - 1] >= j) break;
        mu[size2 - k - 1] = 0;
        k++;
      }
    }
  }

  // Write  minimum unique length, left-anchored wiggle
  fprintf(stderr, "Writing MUL\n");
  fflush(stderr);
  ofstream outfile;
  outfile.open (argv[3]);
  bool prev;
  for (size_t i=0; i<Nchr; i++)
  {
    prev = false;
    for (size_t j=cindices[i]; j<cindices[i+1]; j++)
    {
      if (mu[j] == 0) prev = false;
      else
      {
        if (!prev) outfile << "fixedStep chrom=" + cnames[i] + " start=" +
                              to_string(j - cindices[i] + 1) + " step=1\n";
        outfile << to_string(mu[j]) + "\n";
        prev = true;
      }
    }
  }
  outfile.close();

  // Write  minimum unique length, right-anchored wiggle
  fprintf(stderr, "Writing MUR\n");
  fflush(stderr);
  outfile.open (argv[4]);
  for (size_t i=0; i<Nchr; i++)
  {
    prev = true;
    for (size_t j=cindices[i]; j<cindices[i+1]; j++)
    {
      if (mu[size2 - j - 1] == 0) prev = false;
      else
      {
        if (!prev) outfile << "fixedStep chrom=" + cnames[i] + " start=" +
                              to_string(j - cindices[i] + 1) + " step=1\n";
        outfile << to_string(mu[size2 - j - 1]) + "\n";
        prev = true;
      }
    }
  }
  outfile.close();

  return 0;
}
