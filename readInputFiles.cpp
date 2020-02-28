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
#include <boost/algorithm/string/replace.hpp>
#include <algorithm>
#include <boost/filesystem/operations.hpp>
//#include "boost/filesystem.hpp"
#include <sys/types.h>
#include <dirent.h>
#include "utilities.hpp"
#include "readInputFiles.hpp"
using namespace std;
using namespace boost;
using namespace boost::filesystem;

void read_SMUA_file(string filename, vector<vector<string> > &UA_alignments_and_markup)//function to read cloanalyst files
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

void read_PARTIScsv_file(string filename, vector<vector<string> > &UA_alignments_and_markup)//function to read partis csv files
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
      size_t found = UA_markup_file_contents[i].find("unique_ids");
      if(found!=std::string::npos)
	{
	  continue;
	}

      vector<string> parts;
      boost::split(parts,UA_markup_file_contents[i],is_any_of(","));
      
      sequence_name=parts[0];
      if(parts[8].length()>parts[11].length())
	trimmed_sequence=parts[8];
      else
	trimmed_sequence=parts[11];
      replace_all(trimmed_sequence,".","-");
      
      uca_name=parts[0]+"|UCA";
      if(parts[13].length()>parts[12].length())
	uca_sequence=parts[13];
      else
	uca_sequence=parts[12];
      replace_all(uca_sequence,".","-");
      
      markup_name=parts[0]+"|"+parts[2]+"|"+parts[3]+"|"+parts[4];

      cleanSeqs(trimmed_sequence,uca_sequence);
      string markup_sequence(trimmed_sequence.length(), 'U');
      
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

void read_PARTISyaml_file(string filename, vector<vector<string> > &UA_alignments_and_markup)//function to read partis yaml files
{
  ifstream file(filename.c_str(), std::ios::in );
  if (!file.is_open()) {cerr << "could not open " << filename << " ...exiting...\n"; exit(1);}
  
  vector<string> UA_markup_file_contents;
  vector<string> parts;
  int n=0;
  string file_str;
  getline(file, file_str);
  boost::split(parts,file_str,is_any_of(","));

  string inputSeq=" ",naiveSeq=" ",glSeq=" ",qrSeq=" ";
  string sequence_name=" ",trimmed_sequence=" ",uca_sequence=" ";
  string vGene,dGene,jGene;

  for(int i=0;i<parts.size();i++)
    {
      if(parts[i].find("input_seqs")!= std::string::npos)
	{
	  inputSeq=cleanYAMLline(parts[i]);
	}
      else if(parts[i].find("gl_gap_seqs")!=std::string::npos)
	{
	  glSeq=cleanYAMLline(parts[i]);
	}
      else if(parts[i].find("qr_gap_seqs")!=std::string::npos)
	{
	  qrSeq=cleanYAMLline(parts[i]);
	}
      else if(parts[i].find("naive_seq")!=std::string::npos)
	{
	  naiveSeq=cleanYAMLline(parts[i]);
	}
      else if(parts[i].find("unique_ids")!=std::string::npos)
	{
	  sequence_name=cleanYAMLline(parts[i]);
	}
      else if(parts[i].find("v_gene")!=std::string::npos)
	{
	  vGene=cleanYAMLline(parts[i]);
	}
      else if(parts[i].find("j_gene")!=std::string::npos)
	{
	  jGene=cleanYAMLline(parts[i]);
	}
      else if(parts[i].find("d_gene")!=std::string::npos)
	{
	  dGene=cleanYAMLline(parts[i]);
	  if(inputSeq.length()>qrSeq.length())
	    trimmed_sequence=inputSeq;
	  else
	    trimmed_sequence=qrSeq;

	  if(naiveSeq.length()>glSeq.length())
	    uca_sequence=naiveSeq;
	  else
	    uca_sequence=glSeq;
	  vector<string> temp;

	  cleanSeqs(trimmed_sequence,uca_sequence);
	  
	  string markup_sequence(trimmed_sequence.length(),'U');
	  
	  temp.push_back(sequence_name);	      
	  temp.push_back(trimmed_sequence);
	  temp.push_back(sequence_name+"|UCA");
	  temp.push_back(uca_sequence);
	  temp.push_back(sequence_name+"|"+vGene+"|"+dGene+"|"+jGene);
	  temp.push_back(markup_sequence);
	  UA_alignments_and_markup.push_back(temp);
	  
	}
    }
  return;
}

string cleanYAMLline(string line)
{
  string output="";
  int pos=line.find_last_of(":");
  output=line.substr(pos+1);
  replace_all(output,"\"","");
  replace_all(output,"]","");
  replace_all(output,"[","");
  replace_all(output," ","");
  replace_all(output,"}","");
  return output;
}

void cleanSeqs(string &seq1, string &seq2)
{
 
  replace_all(seq1,".","-");
  replace_all(seq2,".","-");
  if(seq1.back()=='N' || seq2.back()=='N')
    {
      seq1.pop_back();
      seq2.pop_back();
    }
  string seq1tmp="",seq2tmp="";
  for (int i = 0; i < seq1.size(); i++)
    {
      if(seq1[i]!='N' && seq2[i]!='N')
	{
	  seq1tmp+=seq1[i];
	  seq2tmp+=seq2[i];
	}
    }

  seq1=seq1tmp;
  seq2=seq2tmp;
}
