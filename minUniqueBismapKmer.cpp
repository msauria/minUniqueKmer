
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

size_t size, size2, cindex, N2chr;
vector<size_t> inv, idx, mu, cindices, lcp;
vector<char> str;
vector<std::string> cnames;
const char Nchar='N';

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}

size_t find_LCP(size_t idx1, size_t idx2, size_t curr)
{
  while (idx1 + curr < size2 && idx2 + curr < size2 \
         && str[idx1 + curr] == str[idx2 + curr] \
         && str[idx1 + curr] != Nchar)
  {
    curr++;
  }
  if (idx1 + curr == size2 || idx2 + curr == size2 || str[idx1 + curr] == Nchar)
    return 0;
  return curr;
}

void fill_LCP() 
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
  lcp.resize(size2);
  fill(lcp.begin(), lcp.end(), 0);
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
      curr = find_LCP(i, j, curr);
      lcp[k] = curr;
      if (curr > 0) curr--;
    }
  }
  fprintf(stderr, "\r                                     \r");
  fflush(stderr);
}

size_t find_chrom(size_t seq_idx, size_t start, size_t end)
{
  if (end - start == 1) return start;
  size_t mid = (end + start) / 2;
  if (seq_idx >= cindices[mid]) return find_chrom(seq_idx, mid, end);
  else return find_chrom(seq_idx, start, mid);
}

size_t resolve_chrom_runoff(size_t seq_idx)
{
  //fprintf(stderr, "%lu %lu\n", chr_idx, seq_idx);
  size_t sa_idx=inv[seq_idx], sa_idx1=sa_idx-2, best_lcp=0;
  size_t seq_idx1, chr, curr;
  while (sa_idx1 > 0)
  {
    seq_idx1 = idx[sa_idx1];
    curr = 0;
    best_lcp = find_LCP(seq_idx, seq_idx1, curr);
    chr = find_chrom(seq_idx1, 0, N2chr);
    if (seq_idx1 + best_lcp < cindices[chr + 1])
        break;
    else
      sa_idx1--;
  }
  if (sa_idx1 == 0) return 0;
  return best_lcp;
}

void print_info(void)
{
  std::cout << "Find minimum unique k-mer bedgraphs ver. 1.0\n"
        << "\nUsage:\nminUniqueKmer <genome_fa> <suffix_array> <mul_wig> <mur_wig>\n"
        << "Parameters:\n"
        << "<genome_fa>    - Genome multi-fasta sequence file\n"
        << "<suffix_array> - Suffix array produced from reverse-complimented genome\n"
        << "<mul_wig>       - Output name for left-anchored minimum unique k-mer wiggle file\n"
        << "<mur_wig>       - Output name for right-anchored minimum unique k-mer wiggle file\n";
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
    fprintf(stderr, "Reference and SA don't match lengths (%lu %lu %lu)\n", pos, size, sizeof(size_t));
    perror(NULL);
    exit(EXIT_FAILURE);
  }
  size_t Nchr = cindices.size() - 1;
  N2chr = Nchr * 2;
  for (size_t i=Nchr; i<N2chr; i++){
    cindices.push_back(cindices[i]
                       + cindices[N2chr - i]
                       - cindices[N2chr - i - 1]);
  }

  // Add reverse-complimented sequence
  fprintf(stderr, "Reverse complimenting reference\n");
  fflush(stderr);
  for (size_t i=0; i<size; i++)
  {
    if (str[i] == 'A') str[size2 - i - 1] = 'T';
    else if (str[i] == 'T') str[size2 - i - 1] = 'A';
    else if (str[i] == 'C')
    {
        str[size2 - i - 1] = 'G';
        str[i] = 'T';
    }
    else if (str[i] == 'G') str[size2 - i - 1] = 'T';
    else str[size2 - i - 1] = 'N';
  }

  // Load suffix array
  fprintf(stderr, "Loading suffix array\n");
  fflush(stderr);
  idx.resize(size2);
  fread(&idx[0], sizeof(size_t), size2, fp);
  fclose(fp);

  // Find LCP
  fill_LCP();

  // Trim LCPs
  fprintf(stderr, "Trimming LCPs\n");
  fflush(stderr);
  vector<size_t> trim, resolve; 
  trim.resize(0);
  resolve.resize(0);
  // Identify chromosome end LCP overlaps
  for (size_t i=0; i<N2chr; i++)
  {
    size_t inv_j;
    size_t cstart=cindices[i], cend=cindices[i+1], j;
    j = cend - 1;
    while (true)
    {
      if (j == cstart) break;
      inv_j = inv[j];
      if (lcp[inv_j] + j <= cend) break;
      trim.push_back(j);
      inv_j++;
      if (inv_j < size2)
        resolve.push_back(idx[inv_j]);
      j--;
    }
  }

  // Resolve chromosome end LCP overlap-adjacent suffices
  for (size_t j=0; j<resolve.size(); j++)
    lcp[inv[resolve[j]]] = resolve_chrom_runoff(resolve[j]);
  // Trim chromosome end LCP overlaps
  for (size_t j=0; j<trim.size(); j++)
    lcp[inv[trim[j]]] = 0;

  ofstream outfile1;
  outfile1.open (argv[5]);
  for (size_t i=0; i<N2chr; i++){
    if (i < Nchr) outfile1 << "chrom " + cnames[i] + "\n";
    else outfile1 << "chrom_RC " + cnames[i-Nchr] + "\n";
    for (size_t j=cindices[i]; j<cindices[i+1]; j++){
      outfile1 << to_string(lcp[inv[j]]) + "\n";
    }
  }
  outfile1.close();

  // Find minimum unique kmer length
  mu.resize(size2);
  size_t seq_idx;
  fprintf(stderr, "Calculating MUL and MUR\n");
  fflush(stderr);
  for (size_t i=0; i<size2-1; i++)
  {
    if (i % 1000000 == 0)
    {
      fprintf(stderr, "\rMinUnique %luM of %luM", i/1000000, size2/1000000);
      fflush(stderr);
    }
    seq_idx = idx[i];
    mu[seq_idx] = 1 + max(lcp[i], lcp[i + 1]);
    if (mu[seq_idx] == 1) mu[seq_idx] = 0;
  }
  mu[idx[size2-1]] = lcp[size2-1] + 1;
  if (mu[idx[size2-1]] == 1) mu[idx[size2-1]] = 0;
  fprintf(stderr, "\r                                     \r");

  // Trim MUs
  fprintf(stderr, "Trimming MUs\n");
  fflush(stderr);
  // Identify chromosome end MU overlaps
  for (size_t i=0; i<N2chr-1; i++)
  {
    size_t cstart=cindices[i], cend=cindices[i+1], j;
    j = cend - 1;
    while (true)
    {
      if (j == cstart) break;
      if (mu[j] + j <= cend && mu[j] > 0) break;
      mu[j] = 0;
      j--;
    }
  }
  size_t cstart=cindices[N2chr-1], cend=cindices[N2chr], j;
  j = cend - 1;
  while (true)
  {
    if (j == cstart) break;
    if (mu[j] + j < cend && mu[j] > 0) break;
    mu[j] = 0;
    j--;
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
    prev = false;
    for (size_t j=cindices[N2chr-i]-1; j>=cindices[N2chr-i-1]; j--)
    {
      if (mu[j] == 0) prev = false;
      else
      {
        if (!prev) outfile << "fixedStep chrom=" + cnames[i] + " start=" +
                              to_string(cindices[N2chr-i] - j) + " step=1\n";
        outfile << to_string(mu[j]) + "\n";
        prev = true;
      }
    }
  }
  outfile.close();

  return 0;
}
