#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <cstring>
#include <string.h>
#include <ctpl_stl.h>
#include <thread>

using namespace std;

int32_t kmer_len, threadN;
vector<int32_t> mul;
std::vector<std::string> wigout;
ifstream wigfile;

void print_info(void)
{
  std::cout << "Find mean unique k-mer coverage for a given k-mer size at each genomic position ver. 1.0\n"
        << "\nUsage:\nmeanKmerCoverage [options] <k-mer> <mapping_wig>\n"
        << "Parameters:\n"
        << "<k-mer>     - k-mer size\n"
        << "<mul_wig>   - Minimum unique k-mer length (left-anchored) wiggle file\n"
        << "Options:\n"
        << "-t<value>   - Number of additional threads to use\n";
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

bool load_wiggle_chrom(string &chrom, int32_t &end)
{
  string cur, cur0, field, new_chrom="";
  int32_t pos=0, split, bufsize, start, counter=1000000;
  mul.clear();
  while (getline(wigfile, cur))
  {
    if (cur[0] == 'f')
    {
      end = pos;
      cur0 = cur;
      bufsize = cur0.size() + 1;
      while (cur0.length() > 0)
      {
        split = cur0.find(" ");
        if (split == std::string::npos)
          split = cur0.length();
        field = cur0.substr(0, split);
        if (split == cur.length())
          cur0.erase(0, split);
        else
          cur0.erase(0, split + 1);
        if (field.find("chrom") == 0)
        {
          new_chrom = field.substr(field.find("=") + 1, field.length());
          if (chrom == "") chrom = new_chrom;
        }
        else if (field.find("start") == 0)
        {
          start = stoi(field.substr(field.find("=") + 1, field.length())) - 1;
          if (start > pos)
            for (int32_t i=pos; i<start;i++) mul.push_back(0);
          pos = start;
        }
      }
      if (new_chrom != chrom)
      {
        wigfile.seekg(-bufsize, ios_base::cur);
        break;
      }
    } else {
      mul.push_back(stoi(cur));
      pos++;
      if (pos >= counter)
      {
        fprintf(stderr, "\rRead %iMbp of %s", counter / 1000000, chrom.c_str());
        fflush(stderr);
        counter += 1000000;
      }
    }
  }
  fprintf(stderr, "\r                               \r");
  fflush(stderr);
  if (new_chrom == "")
    return false;
  else
    return true;
}

std::string rtrim(const std::string &s)
{
    int32_t end = s.find_last_not_of("0");
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

void find_coverage(const uint id, const string chrom, const int32_t index,
  const int32_t maxlen)
{
  int32_t start=index*1000000, end=(index+1)*1000000;
  end = std::min(end, maxlen);
  int32_t count=0, upos=start-kmer_len;
  bool header=true;
  float fkmer_len = kmer_len;
  std::string outline="", val; 
  for (int32_t i=std::max(0, upos); i<start; i++)
    if ((mul[i] > 0) & (mul[i] <= kmer_len)) count++;
  for (int32_t i=start; i<end; i++)
  {
    if (upos >= 0)
      if ((mul[upos] > 0) & (mul[upos] <= kmer_len)) count -= 1;
    if ((mul[i] > 0) & (mul[i] <= kmer_len)) count += 1;
    if (count == 0) header = true;
    else
    {
      if (header)
      {
        outline += "fixedStep chrom=";
        outline += chrom;
        outline += " start=";
        outline += std::to_string(i + 1);
        outline += " step=1\n";
        header = false;
      }
      val = std::to_string(count / fkmer_len);
      val = rtrim(val);
      outline += val;
      outline += "\n";
    }
    upos++;
  }
  if (outline.length() == 0) outline += "x";
  wigout[index] = outline;
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

  int32_t apos;
  threadN = 1;
  for(apos = 1; apos < argc; ++apos)
    if(argv[apos][0] == '-')
    {
      if (strncmp(argv[apos], "-t", 2) == 0)
        threadN = atoi(&argv[apos][2]);
    } else
      break;

  if(argc - apos < 2)
  { 
    print_info();
    return 0;
  }

  kmer_len = stoi(argv[apos++]);
  // Create thread pool
  threadN = std::min(std::thread::hardware_concurrency(), \
                     static_cast<unsigned>(threadN));
  ctpl::thread_pool pool(threadN);

  // Begin loading wiggle and submitting jobs
  wigfile.open(argv[apos]);
  string chrom="";
  int32_t chrom_end;
  while (load_wiggle_chrom(chrom, chrom_end))
  {
    // Submit jobs
    int32_t num_jobs = (chrom_end - 1) / 1000000 + 1;
    wigout.clear();
    for (int32_t j=0; j<num_jobs; j++)
    {
      wigout.push_back("");
      pool.push(find_coverage, chrom, j, chrom_end);
    }
    // Wait for jobs to finish
    time_t current_time, new_time;
    time(&current_time);
    int32_t outpos=0;
    current_time -= 1;
    while (outpos < num_jobs)
    {
      if (wigout[outpos].length() > 0)
      {
        if (wigout[outpos] != "x")
        {
          fputs(wigout[outpos].c_str(), stdout);
          fflush(stdout);
        }
        wigout[outpos] = "";
        outpos ++;
      }
      time(&new_time);
      if (new_time != current_time)
      {
        std::cerr << "\rScoring " << chrom << " job " << outpos
                  << " of " << num_jobs;
        current_time = new_time;
      }
    }
    chrom = "";
    std::cerr << "\r                                     \r";
  }
  wigfile.close();

  return EXIT_SUCCESS;
}