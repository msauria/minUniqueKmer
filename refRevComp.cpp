using namespace std;

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>


void print_info(void)
{
  std::cout << "Create concatenated reverse-complimented genome character file ver. 1.0\n"
        << "\nUsage:\nrefRevComp <genome_fa> <revcomp_ref>\n"
        << "Parameters:\n"
        << "<genome_fa>   - Genome multi-fasta sequence file\n"
        << "<revcomp_ref> - Output file for reverse-complimented genome\n";
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

int main(int argc, char **argv)
{
  printf("\r                                     \r");
  fflush(stdout);
  if (argc <= 2 || help(argc, argv))
  {
    print_info();
    return 1;
  }

  ifstream input(argv[1]);
  string cur;
  ostringstream out("");
  cout << "Reading reference genome" << endl;
  while (getline(input, cur))
  {
    if(cur[0] != '>')
    {
      for(int i = 0; i<cur.length(); i++)
      {
        if(cur[i] >= 'a' && cur[i] <= 'z')
          cur[i] += 'A' - 'a';
        if(cur[i] >= 'A' && cur[i] <= 'Z')
          out << cur[i];
      }
    }
  }
  string reference = out.str();
  size_t n = reference.length();
  ostringstream out2("");
  for(size_t i=0; i < n; i++)
  {
    if(reference[n-i-1] == 'A')
      out2 << 'T';
    else
    {
      if(reference[n-i-1] == 'T')
        out2 << 'A';
      else 
      {
        if(reference[n-i-1] == 'C')
          out2 << 'G';
        else
        {
          if(reference[n-i-1] == 'G')
            out2 << 'C';
          else
            out2 << 'N';
        }
      }
    }
  }
  string reference2 = out2.str();

  cout << n << endl;
  cout << reference.length() << " " << reference2.length() << "\n";

  FILE *refoutfile = fopen (argv[2], "wb");
  fwrite(&reference[0], sizeof(char), n, refoutfile);
  fwrite(&reference2[0], sizeof(char), n, refoutfile);
}
