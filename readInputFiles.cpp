#include <math.h>
#include <cstdlib>
#include <random>
#include <iostream>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <boost/filesystem/operations.hpp>
//#include "boost/filesystem.hpp"
#include <sys/types.h>
#include <dirent.h>
#include "utilities.hpp"
using namespace std;
using namespace boost;
using namespace boost::filesystem;

void read_SMUA_file(string filename, vector<vector<string> > &UA_alignments_and_markup)
{
  ifstream file(filename.c_str(), std::ios::in );
  if (!file.is_open()) {cerr << "could not open " << filename << " ...exiting...\n"; exit(1);}
  
  vector<string> UA_markup_file_contents;
  
  string file_str;
  while (!getline(file, file_str).eof())
    {
      chomp(file_str);
      UA_markup_file_contents.push_back(file_str);
    }

  //cout << UA_markup_file_contents.size()<<"\n";
  for(int i=0; i<UA_markup_file_contents.size(); i+=6)
    {
      string sequence_name, trimmed_sequence, uca_name, uca_sequence, markup_name, markup_sequence;
      chomp(UA_markup_file_contents[i]);
      chomp(UA_markup_file_contents[i+1]);
      chomp(UA_markup_file_contents[i+2]);
      chomp(UA_markup_file_contents[i+3]);
      chomp(UA_markup_file_contents[i+4]);
      chomp(UA_markup_file_contents[i+5]);

      sequence_name=UA_markup_file_contents[i].substr(1);

      trimmed_sequence=UA_markup_file_contents[i+1];
      uca_name=UA_markup_file_contents[i+2].substr(1);
      uca_sequence=UA_markup_file_contents[i+3];
      markup_name=UA_markup_file_contents[i+4].substr(1);
      markup_sequence=UA_markup_file_contents[i+5];

      //store all in UA_alignments vector
      vector<string> temp;
      temp.push_back(sequence_name);
		      
      temp.push_back(trimmed_sequence);
      temp.push_back(uca_name);
      temp.push_back(uca_sequence);
      temp.push_back(markup_name);
      temp.push_back(markup_sequence);
      UA_alignments_and_markup.push_back(temp);
    }
  return;
}

void read_PARTIS_file(string filename, vector<vector<string> > &UA_alignments_and_markup)
{
  ifstream file(filename.c_str(), std::ios::in );
  if (!file.is_open()) {cerr << "could not open " << filename << " ...exiting...\n"; exit(1);}
  
  vector<string> UA_markup_file_contents;
  
  string file_str;
  while (!getline(file, file_str).eof())
    {
      chomp(file_str);
      UA_markup_file_contents.push_back(file_str);
    }
  for(int i=0;i<UA_markup_file_contents.size();i++)
    {
      string sequence_name, trimmed_sequence, uca_name, uca_sequence, markup_name;
      //cout << UA_markup_file_contents[i]<<"\n";
      size_t found = UA_markup_file_contents[i].find("unique_ids");
      if(found!=std::string::npos)
	{
	  continue;
	}

      vector<string> parts;
      boost::split(parts,UA_markup_file_contents[i],is_any_of(","));

      sequence_name=parts[0];
      trimmed_sequence=parts[8];
      uca_name=parts[0]+"|UCA";
      uca_sequence=parts[13];
      markup_name=parts[0]+"|"+parts[2]+"|"+parts[3]+"|"+parts[4];
      
      if(trimmed_sequence.back()=='N')
	{
	  trimmed_sequence.pop_back();
	  uca_sequence.pop_back();
	}
      string markup_sequence(trimmed_sequence.length(), 'U');
      //cout << uca_name<<endl;
      //getchar();
      
      vector<string> temp;
      temp.push_back(sequence_name);	      
      temp.push_back(trimmed_sequence);
      temp.push_back(uca_name);
      temp.push_back(uca_sequence);
      temp.push_back(markup_name);
      temp.push_back(markup_sequence);
      UA_alignments_and_markup.push_back(temp);
    }
  //getchar();
  return;
}
