#include <stdlib.h>
#include <stdio.h>
#include <ctpl_stl.h>
#include <algorithm>
#include <thread>
#include <mutex>
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <ctime>

using namespace std;

size_t size, size2, cindex, N2chr, threadN;
vector<size_t> inv, idx, mu, cindices, lcp;
vector<char> str;
vector<std::string> cnames;
size_t finished;
std::mutex finished_mutex;
std::mutex lcp_mutex;
const char Nchar='N';

void print_info(void)
{
  std::cout << "Find minimum unique k-mer bedgraphs ver. 1.0\n"
        << "\nUsage:\nminUniqueKmer [options] <genome_fa> <suffix_array>\n"
        << "Parameters:\n"
        << "<genome_fa>    - Genome multi-fasta sequence file\n"
        << "<suffix_array> - Suffix array produced from reverse-complimented genome\n"
        << "Options:\n"
        << "-b             - Find minimum unique k-mer lengths for a bisulfite-converted genome\n"
        << "-p<value>      - Output prefix for the minimum unique k-mer wiggle files\n"
        << "-t<value>      - Number of threads for concurrent processing\n";
}

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}

void read_fasta(const std::string& fname)
{
  ifstream input(fname);
  cindices.resize(0);
  cnames.resize(0);
  string cur;
  size_t pos=0, cend;
  std::string chrom;
  time_t current_time, new_time;
  current_time -= 1;
  str.resize(size2);
  fprintf(stderr, "Reading reference file\n");
  fflush(stderr);
  while (getline(input, cur))
  {
    if(cur[0] == '>')
    {
      cend = 2;
      while ((cend < cur.length()) & (cur[cend] != ' ')) cend++;
      chrom = cur.substr(1, cend - 1);
      cnames.push_back(chrom);
      cindices.push_back(pos);
      time(&new_time);
      if (current_time != new_time)
      {
        std::cerr << "\r                                         \r";
        std::cerr << "\rLoading " << chrom << " sequence";
        current_time = new_time;
      }
    }
    else
    {
      std::transform(cur.begin(), cur.end(), cur.begin(), ::toupper);
      std:strcpy(&str[pos], cur.c_str());
      pos += cur.length();
    }
  }
  cindices.push_back(pos);
  std::cerr << "\r                                         \r";
}

void reverse_compliment(int id, size_t index)
{
  size_t start=index*1000000, end=start+1000000;
  end = std::min(end, size);
  for (size_t i=start; i<end; i++)
  {
    if (str[i] == 'A') str[size2 - i - 1] = 'T';
    else if (str[i] == 'T') str[size2 - i - 1] = 'A';
    else if (str[i] == 'C') str[size2 - i - 1] = 'G';
    else if (str[i] == 'G') str[size2 - i - 1] = 'C';
    else str[size2 - i - 1] = 'N';
  }
  const std::lock_guard<std::mutex> lock(finished_mutex);
  finished++;
}

void reverse_compliment_bs(int id, size_t index)
{
  size_t start=index*1000000, end=start+1000000;
  end = std::min(end, size);
  for (size_t i=start; i<end; i++)
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
  const std::lock_guard<std::mutex> lock(finished_mutex);
  finished++;
}

size_t find_LCP(size_t idx1, size_t idx2, size_t curr)
{
  while (idx1 + curr < size2 && idx2 + curr < size2 \
         && str[idx1 + curr] == str[idx2 + curr] \
         && str[idx1 + curr] != Nchar)
      {
        curr++;
      }
  if (str[idx1 + curr] == Nchar)
    return 0;
  return curr;
}

void fill_LCP(int id, size_t index) 
{
  size_t start=index*1000000, end=start+1000000;
  end = std::min(end, size2);
  size_t curr=0;
  for (size_t i=start; i<end; i++) 
  {
    size_t k = inv[i];
    if (k > 0)
    {
      size_t j = idx[k - 1];
      curr = find_LCP(i, j, curr);
      lcp[k] = curr;
      if (curr > 0) curr--;
    }
  }
  const std::lock_guard<std::mutex> lock(finished_mutex);
  finished++;
}

void invert_sa(int id, size_t index)
{
  size_t start=index*1000000, end=start+1000000;
  end = std::min(end, size2);
  size_t length=end-start;
  for(size_t i=start; i<end; i++)
  {
    inv[idx[i]] = i;
    // inv = ref[i] -> sa
    // idx = sa[i] -> ref
  }
  const std::lock_guard<std::mutex> lock(finished_mutex);
  finished++;
}

void trim_LCP(int id, size_t sa_idx1)
{
  size_t sa_idx0, seq_idx0, seq_idx1, new_lcp;
  sa_idx0 = sa_idx1 - 2;
  seq_idx1 = idx[sa_idx1];
  while ((sa_idx0 > 0) & (lcp[sa_idx0] == 0)) sa_idx0--;
  if (lcp[sa_idx0] == 0) return;
  seq_idx0 = idx[sa_idx0];
  new_lcp = find_LCP(seq_idx0, seq_idx1, 0);
  lcp[sa_idx1] = new_lcp;
  const std::lock_guard<std::mutex> lock1(finished_mutex);
  finished++;
}

void find_mu(int id, size_t index)
{
  size_t start=index*1000000, end=start+1000000;
  end = std::min(end, size2 - 1);
  size_t seq_idx;
  for (size_t i=start; i<end; i++)
  {
    seq_idx = idx[i];
    inv[seq_idx] = 1 + max(lcp[i], lcp[i + 1]);
    if (inv[seq_idx] == 1) inv[seq_idx] = 0;
  }
  if (index == 0)
  {
    inv[idx[size2-1]] = lcp[size2-1] + 1;
    if (inv[idx[size2-1]] == 1) inv[idx[size2-1]] = 0;
  }
  const std::lock_guard<std::mutex> lock(finished_mutex);
  finished++;
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

    //------------------------------------------------------------
  // Parse input parameters
  //------------------------------------------------------------
  bool bisulfite=false;
  threadN = 1;
  std::string prefix="";
  if(argc < 2)
  {
    print_info();
    return 0;
  }

  int32_t apos;
  for(apos = 1; apos < argc; ++apos)
    if(argv[apos][0] == '-')
    {
      if(std::strncmp(argv[apos], "-p", 2) == 0)
        prefix = &argv[apos][2];
      else if (std::strncmp(argv[apos], "-t", 2) == 0)
        threadN = atoi(&argv[apos][2]);
      else if (std::strncmp(argv[apos], "-b", 2) == 0)
        bisulfite = true;
    }
    else
      break;

  if(argc - apos < 2)
  { 
    print_info();
    return 0;
  }

  std::string fasta_fname="", sa_fname="", mul_fname="", mur_fname="";
  fasta_fname.assign(argv[apos++]);
  sa_fname.assign(argv[apos]);
  if (prefix == "")
    prefix.assign(fasta_fname);
  mul_fname.assign(prefix + ".mul.wig");
  mur_fname.assign(prefix + ".mur.wig");

  threadN = std::min(std::thread::hardware_concurrency(), \
                       static_cast<unsigned>(threadN));
  ctpl::thread_pool pool(threadN);

  // Get total genome size
  FILE* fp = std::fopen(sa_fname.c_str(), "rb");
  size2 = filesize(sa_fname.c_str()) / sizeof(size_t);
  size = size2 / 2;

  // Load chromosome sequences and sizes
  read_fasta(fasta_fname);
  size_t Nchr = cindices.size() - 1;
  if(cindices[Nchr] != size)
  {
    fprintf(stderr, "Reference and SA don't match lengths\n");
    perror(NULL);
    exit(EXIT_FAILURE);
  }

  // Add reverse-compliments to chromosome indices
  N2chr = Nchr * 2;
  for (size_t i=Nchr; i<N2chr; i++){
    cindices.push_back(cindices[i]
                       + cindices[N2chr - i]
                       - cindices[N2chr - i - 1]);
  }

  // Add reverse-complimented sequence
  fprintf(stderr, "Reverse complimenting reference\n");
  fflush(stderr);
  finished = 0;
  size_t jobN = (size - 1) / 1000000 + 1;
  fprintf(stderr, "\rRevcomp 0M of %luM", jobN);
  fflush(stderr);
  for (size_t i=0; i<jobN; i++) {
    if (bisulfite)
      pool.push(reverse_compliment_bs, i);
    else
      pool.push(reverse_compliment, i);
  }
  while (finished < jobN)
  {
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    fprintf(stderr, "\rRevcomp %luM of %luM", finished, jobN);
    fflush(stderr);
  }
  fprintf(stderr, "\r                                  \r");
  fflush(stderr);

  // Load suffix array
  fprintf(stderr, "Loading suffix array\n");
  fflush(stderr);
  idx.resize(size2);
  fread(&idx[0], sizeof(size_t), size2, fp);
  fclose(fp);

  // Inverting suffix array
  fprintf(stderr, "Inverting suffix array\n");
  fflush(stderr);
  inv.resize(size2);
  jobN = (size2 - 1) / 1000000 + 1;
  finished = 0;
  for (size_t i=0; i<jobN; i++)
  {
    pool.push(invert_sa, i);
  }
  fprintf(stderr, "\rInv SA 0M of %luM", jobN);
  fflush(stderr);
  while (finished < jobN)
  {
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    fprintf(stderr, "\rInv SA %luM of %luM", finished, jobN);
    fflush(stderr);
  }
  fprintf(stderr, "\r                                  \r");

  // Find LCP
  fprintf(stderr, "\r                                     \r");
  fprintf(stderr, "Calculating LCP\n");
  fflush(stderr);
  lcp.resize(size2);
  fill(lcp.begin(), lcp.end(), 0);
  jobN = (size2 - 1) / 1000000 + 1;
  for (size_t i=0; i<jobN; i++)
  {
    pool.push(fill_LCP, i);
  }
  finished=0;
  fprintf(stderr, "\rLCP 0M of %luM", jobN);
  fflush(stderr);
  while (finished < jobN)
  {
    fprintf(stderr, "\rLCP %luM of %luM", finished, jobN);
    fflush(stderr);
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
  }
  fprintf(stderr, "\r                                     \r");
  fflush(stderr);

  // Trim LCPs
  fprintf(stderr, "Trimming LCPs\n");
  fflush(stderr);
  // Identify chromosome end LCP overlaps
  for (size_t i=0; i<N2chr; i++)
  {
    fprintf(stderr, "\r                                     \r");
    if (i < Nchr)
      fprintf(stderr, "\rTrimming %s", cnames[i].c_str());
    else
      fprintf(stderr, "\rTrimming rev-%s", cnames[N2chr - 1 - i].c_str());
    fflush(stderr);
    size_t inv_j;
    size_t cstart=cindices[i], cend=cindices[i+1], j=cend-1;
    while (true)
    {
      if (j == cstart) break; // Reached start of chromosome
      inv_j = inv[j];
      // If lcp fits on chromosome, then finished trimming
      if (lcp[inv_j] + j + 1 < cend) break;
      lcp[inv_j] = 0;
      j--;
    }
  }
  fprintf(stderr, "\r                                     \r");
  fprintf(stderr, "Fixing adjacent LCPs\n");
  fflush(stderr);
  finished = 0;
  size_t changed=0;
  time_t current_time, new_time;
  current_time -= 1;
  for (size_t i=0; i<N2chr; i++)
  {
    size_t inv_j;
    size_t cstart=cindices[i], cend=cindices[i+1], j=cend-1;
    while (true)
    {
      if (j == cstart) break; // Reached start of chromosome
      inv_j = inv[j];
      // If lcp fits on chromosome, then finished trimming
      if (lcp[inv_j] + j + 1 < cend) break;
      if ((lcp[inv_j] == 0) & (inv_j + 1 < size2) & (lcp[inv_j + 1] > 0))
      {
        // Trim chromosome end LCP overlaps
        pool.push(trim_LCP, inv_j + 1);
        time(&new_time);
        if (current_time != new_time)
        {
          fprintf(stderr, "\rFixing adjacent LCP %lu of %lu", \
                  finished, changed);
          fflush(stderr);
          current_time = new_time;
        }
        changed++;
      }
      j--;
    }
  }
  while (finished < changed)
  {
    fprintf(stderr, "\rFixing adjacent LCP %lu of %lu", \
            finished, changed);
    fflush(stderr);
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
  }
  fprintf(stderr, "\r                                                 \r");
  fflush(stderr);

  // Find minimum unique kmer length
  fill(inv.begin(), inv.end(), 0); // Reusing inv for memory savings
  fprintf(stderr, "Calculating MUL and MUR\n");
  fflush(stderr);
  finished=0;
  jobN = (size2 - 1) / 1000000 + 1;
  for (size_t i=0; i<jobN; i++)
  {
    pool.push(find_mu, i);
  }
  fprintf(stderr, "\rMU 0M of %luM", jobN);
  fflush(stderr);
  while (finished < jobN)
  {
    fprintf(stderr, "\rMU %luM of %luM", finished, jobN);
    fflush(stderr);
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
  }
  fprintf(stderr, "\r                                                 \r");
  fflush(stderr);

  // Trim MUs
  fprintf(stderr, "Trimming MUs\n");
  fflush(stderr);
  // Identify chromosome end MU overlaps
  for (size_t i=0; i<N2chr-1; i++)
  {
    fprintf(stderr, "\r                                     \r");
    if (i < Nchr)
      fprintf(stderr, "\rTrimming %s", cnames[i].c_str());
    else
      fprintf(stderr, "\rTrimming rev-%s", cnames[N2chr - 1 - i].c_str());
    fflush(stderr);
    size_t cstart=cindices[i], cend=cindices[i+1], j;
    j = cend - 1;
    while (true)
    {
      if (j == cstart) break; // Reached start of chromosome
      if (inv[j] + j <= cend && inv[j] > 0) break; // Found valid mu
      inv[j] = 0;
      j--;
    }
  }
  fprintf(stderr, "\rTrimming rev-%s", cnames[N2chr - 1].c_str());
  fflush(stderr);
  size_t cstart=cindices[N2chr-1], cend=cindices[N2chr], j;
  j = cend - 1;
  while (true)
  {
    if (j == cstart) break;
    if (inv[j] + j < cend && inv[j] > 0) break;
    inv[j] = 0;
    j--;
  }
  fprintf(stderr, "\r                                                 \r");
  fflush(stderr);

  // Write  minimum unique length, left-anchored wiggle
  fprintf(stderr, "Writing MUL\n");
  fflush(stderr);
  ofstream outfile;
  outfile.open (mul_fname);
  bool prev;
  size_t counter=1000000, sizeM=size/1000000;
  for (size_t i=0; i<Nchr; i++)
  {
    prev = false;
    for (size_t j=cindices[i]; j<cindices[i+1]; j++)
    {
      if (j >= counter)
      {
        fprintf(stderr, "\rWriting MUL %lu of %luM", counter/1000000, sizeM);
        fflush(stderr);
        counter += 1000000;
      }
      if (inv[j] == 0) prev = false;
      else
      {
        if (!prev) outfile << "fixedStep chrom=" + cnames[i] + " start=" +
                              to_string(j - cindices[i] + 1) + " step=1\n";
        outfile << to_string(inv[j]) + "\n";
        prev = true;
      }
    }
  }
  outfile.close();

  // Write  minimum unique length, right-anchored wiggle
  fprintf(stderr, "\r                                                 \r");
  fprintf(stderr, "Writing MUR\n");
  fflush(stderr);
  outfile.open (mur_fname);
  counter = 1000000;
  for (size_t i=0; i<Nchr; i++)
  {
    prev = false;
    for (size_t j=cindices[i]; j<cindices[i+1]; j++)
    {
      if (j >= counter)
      {
        fprintf(stderr, "\rWriting MUR %lu of %luM", counter/1000000, sizeM);
        fflush(stderr);
        counter += 1000000;
      }
      if (inv[size2 - j - 1] == 0) prev = false;
      else
      {
        if (!prev) outfile << "fixedStep chrom=" + cnames[i] + " start=" +
                              to_string(j - cindices[i] + 1) + " step=1\n";
        outfile << to_string(inv[size2 - j - 1]) + "\n";
        prev = true;
      }
    }
  }
  outfile.close();
  fprintf(stderr, "\r                                                 \r");
  fflush(stderr);
  return 0;
}
