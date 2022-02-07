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
#include "HTML.hpp"
#include "utilities.hpp"

using namespace std;

//classes
class S5F_mut
{
public:
  string fivemer, group, subst_group;
  double score,score25,score75;
  map<char, double> substitutions;

  ~S5F_mut(){};
  S5F_mut(){};
  S5F_mut(string _fivemer,double _score, string _group, double _score25, double _score75)
  {
    fivemer=_fivemer;
    score=_score;
    group=_group;
    score25=_score25;
    score75=_score75;
  }
  
};

class Seq //nucleotide level sequence object (later could have each base point to a aa object)
{
public:
  int aa_num;
  string aa;
  string base;
  double S5F_mut_score;
  double simulated_aa_positional_frequency;
  map<char,double> all_simulated_aa_positional_frequencies_map;
  string SMUA_code;
  string CDR_markup;
  bool isMut;
  ~Seq(){};
  Seq(){};
  Seq(string _base, int _aa_num,string _aa,double _S5F_mut_score)
  {
    base=_base;
    aa_num=_aa_num;
    aa=_aa;
    S5F_mut_score=_S5F_mut_score;
  }
};

///functions
void read_SMUA_file(string, vector<vector<string> > &,bool);
void load_S5F_files(string,string, map<string,S5F_mut> &);
void process_fasta_sequence_to_seq_vector(string &,vector<Seq> &);
void process_SMUA_sequence_to_seq_vector(string &, string &, vector<Seq> &, map<string,string> &, map<string, S5F_mut> &);
void convert_seq_vector_to_HTML_table(vector<Seq> &, string , HTML::Table &, double &, vector<double>  &, vector<double> &);
void number_of_mutations_two_seqs(string &, string &, int &);
void simulate_S5F_mutation(string , int &, map<string,S5F_mut> &, mt19937 &, uniform_real_distribution<double> &, bool, vector<string> &, bool,  vector<bool> &);
vector<pair<char,double> > sort_map_into_pair_vctr(map<char,double> &);
bool mycompare(pair<char,double> A, pair<char,double> B){return A.second > B.second;}
void correct_for_fivemer_with_gap(int, string, string &);
void print_output(string, vector<vector<Seq> > &, vector<string>, int, double, vector<double> &);
void print_pct_progress(int, int, int);
void get_mutability_scores(map<string,S5F_mut> &, string, int, bool, vector<bool> &, vector<double> &, vector<double> &, double &, double &);
void print_freq_table_to_file(string,  map<int, map<char,double> > &);
void cleanup_SMUA_sequences(string, string, string , string , string , string &, string &, string &, string , string , int &, bool &);
bool sequence_has_ambiguities(string);
void convert_2D_seq_vector_to_HTML_table(vector<vector<Seq> >&, vector<string> &, HTML::Table &, double &, vector<double> &, int &);

///templated functions
template <typename Type>
void vector2D_to_3D(vector<vector<Type> > &, int, vector<vector<vector<Type> > > &);
template <typename Type>
void vector1D_to_2D(vector<Type> &, int, vector<vector<Type> > &);
template <typename Type>
string convert_to_string(Type);
  
///TODO:
//0. The gap bases should have n/a for mutability score in HTML 
//1. Use the markup string to highlight CDRs in the HTML
//2. Make sure last cell in ladder is 1 always

int main(int argc, char *argv[])
{  
  if (argc <2){cout << "USAGE: analyze_mutations -SMUA [SMUA file] -w [line wrap length (60)] -m [S5F mutability file] -s [S5F substitution file] -max_iter [cycles of B cell maturation(100)] -c [cutoff for highlighting low prob (1=1%)] -replace_J_upto [number of replacements in J allowed] -chain [chain type (heavy=default|kappa|lambda)] -species [(human=default|rhesus)] -clean_first [clean the SMUA prior to running] -random_seed [provide a random seed] -color_probables\n"; exit(1);}
 
  ///get cmdline args
  int i=0, line_wrap_length=60, max_iter=100, mutation_count_from_cmdline=-1, replace_J_upto=0, random_seed=0;
  string fasta_filename="", mutability_filename="", substitution_filename="", SMUA_filename="", species="human", chain_type="heavy", freq_table_dir="", output_filename="";
  int SMUA_start=0, SMUA_end=-1;
  double low_prob_cutoff=.02;
  bool ignore_CDR3=false, clean_SMUA_first=false, user_provided_random_seed=false;
   while(i<argc)
     {
       string arg=argv[i];
       string next_arg;
       if (i<argc-1){next_arg=argv[i+1];}else{next_arg="";}
       if (arg == "-SMUA")
	 {
	   SMUA_filename=next_arg;
	 }
       if (arg == "-freq_table_dir")
	{
	  freq_table_dir=next_arg;
	}
       if (arg == "-output_file")
	 {
	   output_filename=next_arg;
	 }
       if (arg == "-w")
	 {
	   line_wrap_length=atoi(next_arg.c_str());
	 }
       if (arg == "-ignore_cdr3")
	 {
	   ignore_CDR3=true;
	 }
      
       i++;
     }

   cerr << "highlighting residues with less than " << low_prob_cutoff << " probability for mutation\n"; 

   ///amino acids vector
   vector<char> amino_acids={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
  
   ///setup random num generation
   std::random_device rd;
   int seed;
   if (user_provided_random_seed)
     {seed=random_seed;}
   else
     {seed=rd();}
   
   std::mt19937 gen(seed);
   std::uniform_real_distribution<double> dis(0, 1);

   ///load dna_to_aa map
   map<string,string> dna_to_aa_map;
   get_aa_tranx_map(dna_to_aa_map);

   ///read input sequence alignment
   map <string, string> sequences;
   vector <string> sequence_names;   

   //define color ladder
   vector<double> color_ladder{0.0001, 0.001, 0.01, 0.02, 0.10, 0.20, 0.5, 1};

   //load fasta
   //read_fasta_file(fasta_filename, sequences, sequence_names);

   //load SMUA file
   vector<vector<string> > SMUA_alignments_and_markup;
   read_SMUA_file(SMUA_filename, SMUA_alignments_and_markup);

   vector<Seq> seq_vector;
   
   //go through each sequence and print HTML
   vector<vector<Seq> > all_sequences;
   vector<string> all_sequences_names;
   for(int i=0; i<SMUA_alignments_and_markup.size(); i++)
     {
       string sequence_name=SMUA_alignments_and_markup[i][0];
       string dna_sequence=SMUA_alignments_and_markup[i][1];
       string UCA_sequence_name=SMUA_alignments_and_markup[i][2];
       string dna_UCA_sequence=SMUA_alignments_and_markup[i][3];
       string markup_header=SMUA_alignments_and_markup[i][4];
       string dna_markup_string=SMUA_alignments_and_markup[i][5];
       string sequence, UCA_sequence;
       translate_dna_to_aa(dna_UCA_sequence, UCA_sequence, 1, dna_to_aa_map);
       translate_dna_to_aa(dna_sequence, sequence, 1, dna_to_aa_map);

       string cdr_dna_markup_string, cdr_markup_string;
       for(int j=0; j<dna_markup_string.length(); j++)
	 {
	   if ((dna_markup_string[j]=='V')||(dna_markup_string[j]=='n')||(dna_markup_string[j]=='D')||(dna_markup_string[j]=='J')){cdr_dna_markup_string+="3";}
	   else if (dna_markup_string[j]=='B'){cdr_dna_markup_string+="2";}
	   else if (dna_markup_string[j]=='A'){cdr_dna_markup_string+="1";}
	   else {cdr_dna_markup_string+="0";}
	 }
       for(int j=0; j<cdr_dna_markup_string.length(); j+=3)
	 {
	   cdr_markup_string+=cdr_dna_markup_string[j];
	 }
       cout << sequence << "\n" << cdr_markup_string << "\n"; 
       process_fasta_sequence_to_seq_vector(sequence, seq_vector);
      
       //open corresponding freq table file for sequence
       bool error_status=false;
       vector<vector<string> > freq_table_str_vector=read_delimited_file(freq_table_dir+sequence_name+".freq_table.txt", ",", error_status);
       if (error_status){continue;}
   
       if (freq_table_str_vector.size()==0){cerr << "ERROR: could not read file " << freq_table_dir+sequence_name+".freq_table.txt\n"; exit(1);}
       vector<string> AAs{"A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"};
       map<string, map<string, double> > freq_table_map;
       
       if (sequence.length() != freq_table_str_vector.size()-1){cerr << "lengths differ\n"; 
	 cerr << sequence_name << "\t" << sequence.length() << "\t" << freq_table_str_vector.size() << "\n"; }
       
       for(int j=1; j<freq_table_str_vector.size(); j++)
	 {
	   for(int k=1; k<freq_table_str_vector[j].size(); k++)
	     {
	       freq_table_map[freq_table_str_vector[j][0]][AAs[k-1]]=atof(freq_table_str_vector[j][k].c_str());
	     }
	 }

       //assign aa probabilities to seq vector
       for(int j=0; j<seq_vector.size(); j++)
	 {
	   string pos=freq_table_str_vector[j+1][0];
	   string aa(1,sequence[j]);
	   double aa_freq=freq_table_map[pos][aa];
	   if ((sequence[j] == '-')|| (UCA_sequence[j]=='-')){aa_freq=1;}
	   if ((ignore_CDR3) && (cdr_markup_string[j]=='3')){aa_freq=1;}
	   seq_vector[j].simulated_aa_positional_frequency=aa_freq;
	   seq_vector[j].CDR_markup=cdr_markup_string[j];
	   if (sequence[j] != UCA_sequence[j]){seq_vector[j].isMut=true;}else{seq_vector[j].isMut=false;}
	   //  cout << pos << "\t" << aa << "\t" << aa_freq << "\n"; 
	 }
       
       all_sequences.push_back(seq_vector);
       all_sequences_names.push_back(sequence_name);
      
     }
   
   //equalize lengths
   int max_seq_length=0;
   for(int i=0; i<all_sequences.size(); i++)
     {
       if (all_sequences[i].size()>max_seq_length){max_seq_length=all_sequences[i].size();}
     }
   for(int i=0; i<all_sequences.size(); i++)
     {
       for(int j=all_sequences[i].size(); j<max_seq_length; j++)
	 {
	   Seq temp;
	   temp.aa="Z";
	   temp.simulated_aa_positional_frequency=1;
	   all_sequences[i].push_back(temp);
	 }
     }
   

   HTML::Table html_table;
   html_table.hclass="results";
   //convert seq vector to html table
   print_output(output_filename, all_sequences, all_sequences_names, line_wrap_length, low_prob_cutoff, color_ladder);
   
   //convert_seq_vector_to_HTML_table(seq_vector,sequence_name,html_table, low_prob_cutoff, color_ladder);
   
   //print HTML tables
   /*
   html_table.print(file_string);
   file_string+="<p></p>\n"; 
   
   file_string+="<body>\n</html>\n"; 
   //print HTML string to file
   ofstream file_out;
   file_out.open(output_filename.c_str());
   file_out << file_string;
   file_out.close();
   */
   
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///                FUNCTION DEFINITIONS
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_output(string filename, vector<vector<Seq> > &all_sequences, vector<string> sequence_names, int line_wrap_length, double low_prob_cutoff, vector<double> &color_ladder)
{
  vector<vector<vector<Seq> > > split_all_sequences;
  vector2D_to_3D(all_sequences,line_wrap_length, split_all_sequences);
  ///output html header
  string file_string="";
  file_string+="<html xmlns='http://www.w3.org/1999/xhtml' xml:lang='en'>\n"; 
  file_string+="<head>\n"; 
  file_string+=" <meta http-equiv='Content-Type' content='text/html; charset=utf-8' />\n";
  file_string+=" <title>Antibody Mutation Analysis</title>\n"; 
  file_string+=" <link rel='stylesheet' href='sequence_color.css' />\n"; 
  file_string+="</head>\n"; 
  file_string+="<body>\n";
  file_string+="<p></p><br>\n"; 
  //cerr << "number of splits: " << split_all_sequences.size() << "\n"; 
  //for each split, make a HTML table and print it
  int counter=1;
  for(int i=0; i<split_all_sequences.size(); i++)
    {
      HTML::Table html_table;
      html_table.hclass="results";
      //convert seq vector to html table
      convert_2D_seq_vector_to_HTML_table(split_all_sequences[i],sequence_names,html_table, low_prob_cutoff, color_ladder, counter);
      //missing step: stylize the table
      
      //print HTML tables
      html_table.print(file_string);
      file_string+="<p></p>\n"; 
    }

  file_string+="<body>\n</html>\n"; 

  ofstream file_out;
  file_out.open(filename.c_str());
  file_out << file_string;
  file_out.close();
}

void convert_2D_seq_vector_to_HTML_table(vector<vector<Seq> >&v2, vector<string> &names, HTML::Table &html_table, double &low_prob_cutoff, vector<double> &color_ladder, int &counter)
{

  //RULER
  HTML::Tr ruler1_row, ruler2_row;
  HTML::Td ruler1("ruler","","","","");
  HTML::Td ruler2("ruler","","","","");
  ruler1_row.cols.push_back(ruler1);
  ruler2_row.cols.push_back(ruler2);

  if (v2.size()==0){cerr << "No data for HTML table. Exiting\n"; exit(1);}
  for(int j=0; j<v2[0].size(); j++) ///iterate over cols
    {
      ruler1.value=""; 
      ruler2.value="";
      if (((j+1)%10==0)||(j==0))
	{
	  ostringstream s;
	  s<<counter+j;
	  ruler1.value=s.str();
	  ruler2.value="|";
	}
      else if ((j+1)%5==0){ruler2.value="|";}
      ruler1_row.cols.push_back(ruler1);
      ruler2_row.cols.push_back(ruler2);
    }      
  counter+=v2[0].size();

  html_table.rows.push_back(ruler1_row);
  html_table.rows.push_back(ruler2_row);

  //REST OF ROWS
  for(int i=0; i<v2.size(); i++) ///iterate over rows
    {
      HTML::Tr row1, cdr_row, spacer_row; 
    
      //string str0="this";
      HTML::Td td0("seq_name","","","",names[i]);
      HTML::Td blank("noborder","","","","");
   
      row1.cols.push_back(td0);
      cdr_row.cols.push_back(blank);
      spacer_row.cols.push_back(blank);

      for(int j=0; j<v2[i].size(); j++) ///iterate over cols
	{
	  string str1="<div class=\"sm\">"+v2[i][j].aa+"</div>"; //
	  HTML::Td td1("","","","",str1);
	  HTML::Td td_cdr("noborder","","","","");
	  HTML::Td spacer("spacer","","","","");
	  //assign category for coloring
	  for(int k=0; k<color_ladder.size(); k++)
	    {
	      if (v2[i][j].simulated_aa_positional_frequency <= color_ladder[k])
		{
		  ostringstream ss;
		  ss << "color_cat" << k+1;
		  td1.hclass=ss.str();
		  break;
		}
	    }
	  if (v2[i][j].aa=="Z"){td1.hclass="noborder"; td1.value="";}
	  if ((v2[i][j].CDR_markup=="1")||(v2[i][j].CDR_markup=="2")||(v2[i][j].CDR_markup=="3")){td_cdr.hclass="CDR";}
	  if (v2[i][j].isMut){td1.hclass+=" mut";}
	
	  row1.cols.push_back(td1); //row 1: amino acid
	  cdr_row.cols.push_back(td_cdr); //row 0 CDR line
	  spacer_row.cols.push_back(spacer);//row 2 spacer line
	}

      //collapse cdr row elements
      for(int j=0; j<cdr_row.cols.size(); j++)
	{
	  if ( (((j==0)&&cdr_row.cols[j].hclass=="CDR")) || ((j>0)&&(cdr_row.cols[j-1].hclass!="CDR")&&(cdr_row.cols[j].hclass=="CDR"))    ) //start of CDR
	    {
	      int cdr_start=j;
	      int k=j;
	      while ((cdr_row.cols[k].hclass=="CDR")&&(k<cdr_row.cols.size()))
		{
		  k++;
		}
	      int cdr_end=k-1;
	      //cerr << cdr_start << "\t" << cdr_end << "\n"; 
	      // int d; cin >> d;
	      ostringstream s;
	      s<<cdr_end-cdr_start+1;
	      cdr_row.cols[cdr_start].colspan=s.str();
	      if (cdr_end-cdr_start>=1){
		cdr_row.cols.erase(cdr_row.cols.begin()+cdr_start+1, cdr_row.cols.begin()+cdr_end+1);}
	     
	    }
	}

      html_table.rows.push_back(cdr_row);
      html_table.rows.push_back(row1);
      html_table.rows.push_back(spacer_row);
    }

  return;
}



void print_freq_table_to_file(string filename,  map<int, map<char,double> > &positional_aa_freqs)
{
  ofstream file_out;
  file_out.open(filename.c_str());

  ///amino acids vector
  vector<char> amino_acids={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};

  file_out << "pos";
  for(int i=0; i<amino_acids.size(); i++){file_out << "," << amino_acids[i];}
  file_out << "\n"; 
  for(int j=0; j<positional_aa_freqs.size(); j++)
    {
      file_out << j+1;
      for(int i=0; i<amino_acids.size(); i++)
	{
	  file_out << "," <<positional_aa_freqs[j][amino_acids[i]];
	}
      file_out << "\n";
    }
  file_out.close();
 
}





void number_of_mutations_two_seqs(string &s1, string &s2, int &mutation_count)
{
  ///Assumes sequences are already properly aligned
  assert(s1.length()==s2.length());

  mutation_count=0;
  for(int i=0; i<s1.size(); i++)
    {
      if ((s1[i] == '-') || (s2[i] == '-')){continue;} ///Not counting gaps as mutations currently
      if (s1[i] != s2[i]) {mutation_count++;}
    }
  return;
}

void convert_seq_vector_to_HTML_table(vector<Seq> &v, string name, HTML::Table &html_table, double &low_prob_cutoff, vector<double> &color_ladder)
{
  vector<char> amino_acids={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
   
 
  HTML::Tr row1; ///equal to the number of attributes in seq obj we want to print
  string str0="<div class=\"seq_name\">"+name+"</div>";
  HTML::Td td0("seq_name","","","",str0);
  row1.cols.push_back(td0);
  for(int j=0; j<v.size(); j++) ///iterate over cols
    {
      string str1="<div class=\"sm\">"+v[j].aa+"</div>"; //
      HTML::Td td1("","","","",str1);

      //assign category for coloring
      for(int k=0; k<color_ladder.size(); k++)
	{
	  if (v[j].simulated_aa_positional_frequency <= color_ladder[k])
	    {
	      ostringstream ss;
	      ss << "color_cat" << k+1;
	      td1.hclass=ss.str();
	      break;
	    }
	}
      
      row1.cols.push_back(td1); //row 1: amino acid
    }
  html_table.rows.push_back(row1);
 
  return;
}

void process_SMUA_sequence_to_seq_vector(string &sequence, string &markup_string, vector<Seq> &seq_vector, map<string,string> &dna_to_aa_map, map<string, S5F_mut> &S5F_5mers)
{
  seq_vector.clear();
  string aa="X";
  int aa_counter=0;
  if (sequence.length() != markup_string.length()){cerr << "WARNING: sequence and markup string lengths are not identical\n"; }

  for(int i=0; i<sequence.length(); i++)
    {
      Seq temp;
      temp.base=sequence.substr(i,1);
      //get amino acid
      if (i%3==0)
	{
	  //how to deal with gaps -> gaps should have to be 3mers in order for alignment to be valid between two fxnl sequences, check for this
	  string codon=sequence.substr(i,3);
	  if (dna_to_aa_map.find(codon) != dna_to_aa_map.end())
	    aa=dna_to_aa_map[codon];
	  else
	    aa="X";
	  aa_counter++;
	}
      //get S5F fivemer mutability score
      double mut_score;
      if ((i<2)|| (i>=sequence.length()-2))
	{
	  mut_score=-1;
	}
      else
	{
	   string fivemer=sequence.substr(i-2,5);
	   string fivemer_copy=fivemer;
	     if (fivemer.find('-') != std::string::npos)
	     {
	       string new_fivemer="";
	       correct_for_fivemer_with_gap(i,sequence,new_fivemer);
	       fivemer=new_fivemer;
	     }
	     if (fivemer == "NO_SCORE"){ mut_score=-1; cerr << "WARNING: gap detected in central position of 5mer: " << fivemer_copy << ". Mutability zeroed out for this position. Make sure insertion relative to UCA is known.\n";}
	     else if (S5F_5mers.find(fivemer) == S5F_5mers.end()){ mut_score=-1; cerr << "ERROR 47: no 5mer for " << fivemer_copy << " found in S5F. Trying to resolve to " << fivemer << " but no lstill no bueno\n"; }
	     else{mut_score=S5F_5mers[fivemer].score;}
	}
      temp.aa=aa;
      temp.aa_num=aa_counter;
      temp.S5F_mut_score=mut_score;
      temp.SMUA_code=markup_string.substr(i,1);
      seq_vector.push_back(temp);
    }
  return;

}

void process_fasta_sequence_to_seq_vector(string &sequence,vector<Seq> &seq_vector)
{
  seq_vector.clear();
  int aa_counter=0;
  for(int i=0; i<sequence.length(); i++)
    {
      Seq temp;
      temp.base="";
      temp.aa=sequence[i];
      temp.aa_num=++aa_counter;
      temp.S5F_mut_score=0;
      seq_vector.push_back(temp);
    }
  return;
}


template <typename Type>
void vector1D_to_2D(vector<Type> &vector1D, int interval, vector<vector<Type> > &vector2D)
{
  vector2D.clear();///start fresh
  assert(interval<=vector1D.size());
  vector<Type> row;
  for(int i=0; i<vector1D.size(); i++)
    {
      row.push_back(vector1D[i]);
      if (((i+1)%interval==0)||(i==vector1D.size()-1))
	{
	  vector2D.push_back(row);
	  row.clear();
	}
    }
  return;
}

template <typename Type>
void vector2D_to_3D(vector<vector<Type> > &vector2D, int interval, vector<vector<vector<Type> > > &vector3D)
{
  vector3D.clear();
  //   cerr << vector2D.size() << " by " << vector2D[0].size() << "\n"; 
  //slice each row by interval, and stack into one big 2D vector
  vector<vector<Type> > all_sliced_rows;
  int max_num_splits=0;
  for(int i=0; i<vector2D.size(); i++)
    {
      vector<vector<Type> > dim2;
      vector1D_to_2D(vector2D[i],interval,dim2);
      for(int j=0; j<dim2.size(); j++)
	{
	  all_sliced_rows.push_back(dim2[j]);
	}
      if (dim2.size() > max_num_splits){max_num_splits=dim2.size();}
    }

  //cerr << "asr: " <<  all_sliced_rows.size() << "\n"; 
  //re-apportion the slices in the right order to a 3D vector
  for(int i=0; i<max_num_splits; i++)
    {

      vector<vector<Type> > temp;
      for(int j=0; j<vector2D.size(); j++)
	{
	  temp.push_back(all_sliced_rows[i+(j*max_num_splits)]);

	}
      vector3D.push_back(temp);
    }

  //cerr << vector3D.size() << " by " << vector3D[0].size() << " by " << vector3D[0][0].size() << "\n"; 
  return;
}



void load_S5F_files(string mutability_filename, string substitution_filename, map<string, S5F_mut> &S5F_5mers)
{
  //Open mutability CSV and read in
  ifstream file1(mutability_filename.c_str(), std::ios::in );
  if (!file1.is_open()) {cerr << "could not open " << mutability_filename << " ...exiting...\n"; exit(1);}

  string file1_str;
  int counter=0;
  while (!getline(file1, file1_str).eof())
    {
      if (counter==0){counter++; continue;}//skip header line

      chomp(file1_str);
      vector<string> tokens; 
      tokenize(file1_str, tokens," ");
      for(int i=0; i<tokens.size(); i++){boost::replace_all(tokens[i],"\"","");}
      S5F_mut temp(tokens[0],atof(tokens[1].c_str()),tokens[2],atof(tokens[3].c_str()),atof(tokens[4].c_str()));
      S5F_5mers[tokens[0]]=temp;
      counter++;
    }
  
  //Open substitution CSV and read in
  ifstream file2(substitution_filename.c_str(), std::ios::in );
  if (!file2.is_open()) {cerr << "could not open " << mutability_filename << " ...exiting...\n"; exit(1);}

  string file2_str;
  counter=0;
  while (!getline(file2, file2_str).eof())
    {
      if (counter==0){counter++; continue;}//skip header line

      chomp(file2_str);
      vector<string> tokens; 
      tokenize(file2_str, tokens," ");
      for(int i=0; i<tokens.size(); i++){boost::replace_all(tokens[i],"\"","");}
      if (S5F_5mers.find(tokens[0]) != S5F_5mers.end())
	{
	  S5F_5mers[tokens[0]].substitutions['A']=atof(tokens[1].c_str());
	  S5F_5mers[tokens[0]].substitutions['C']=atof(tokens[2].c_str());
	  S5F_5mers[tokens[0]].substitutions['G']=atof(tokens[3].c_str());
	  S5F_5mers[tokens[0]].substitutions['T']=atof(tokens[4].c_str());
	  S5F_5mers[tokens[0]].subst_group=tokens[5];
	}
      else{cerr << "ERROR parsing substitution.csv.  No fivemer " << tokens[0] << " found in mutability scores\n"; exit(1);}
      counter++;
    }
  
  return;
}


void read_SMUA_file(string filename, vector<vector<string> > &UA_alignments_and_markup,bool fred)
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


template <typename Type>
string convert_to_string(Type t)
{
  ostringstream convert;
  convert << t;
  return convert.str();
}



vector<pair<char, double > > sort_map_into_pair_vctr(map<char,double> &M)
{
  
  vector<pair<char,double> > V;
  for(map<char,double>::iterator it=M.begin(); it !=M.end(); ++it)
    {
      pair<char,double> p(it->first, it->second);
      V.push_back(p);
    }
  sort(V.begin(), V.end(), mycompare);
  return V;
}

void correct_for_fivemer_with_gap(int i, string sequence, string &new_fivemer)
{
  //case 0: gap in the middle position
  if (sequence[i] =='-')
    {
      new_fivemer="NO_SCORE";
    }
  else //traverse left and right through gaps until two bases are found on each side, or end of string is reached
    {
      //5' (go left)
      int five_prime_count=0;
      int j=i-1;
      string five_prime_bases="";
      while(five_prime_count<2 && j>=0)
	{
	  if (sequence[j] !='-')
	    {
	      five_prime_bases=sequence[j]+five_prime_bases;
	      five_prime_count++;
	    }
	  j--;
	}
      //3 ' (go right)
      int three_prime_count=0;
      j=i+1;
      string three_prime_bases="";
      while(three_prime_count<2 && j<sequence.length())
	{
	  if (sequence[j] !='-')
	    {
	      three_prime_bases+=sequence[j];
	      three_prime_count++;
	    }
	  j++;
	}
      
      new_fivemer=five_prime_bases+sequence[i]+three_prime_bases;
      if (new_fivemer.length() !=5){new_fivemer="NO_SCORE";}
    }
  
  return;
}

void print_pct_progress(int i, int size, int level)
{
  if (size<100){return;}
  double a=size/(100*pow(10,level));
  int b=1;
  if (a>1){b=(int) a;}
  if ((i%b)==0){cerr << setw(3) << fixed << setprecision(level) << (i/(double)size)*100.0 << "%\r" << flush;}
  
}


void cleanup_SMUA_sequences(string sequence_name, string markup_header, string UCA_sequence, string sequence, string markup, string &new_UCA_sequence, string &new_sequence, string &new_markup, string species, string chain_type, int &number_of_replacements, bool &error_status)
{
  if ( ((species =="human") && (chain_type != "heavy")) || ((species == "rhesus")&&(chain_type == "lambda")) ){ cerr << "CLEANUP ONLY WORKS FOR: Human/Rhesus Heavy and Rhesus Kappa for now. Sorry. Exiting.\n"; exit(1);}

  new_UCA_sequence=UCA_sequence;
  new_sequence=sequence;
  new_markup=markup;
  number_of_replacements=0;
  if ((UCA_sequence.length() != sequence.length())||(sequence.length()!=markup.length())||(UCA_sequence.length()!=markup.length())){cerr << "FATAL ERROR: UCA sequence, sequence and markup strings are not same length\n"; exit(1);}

  string H_J4_1="ACTACTTTGACTACTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCAG";
  string H_J4_2="ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG";
  string H_J4_3="GCTACTTTGACTACTGGGGCCAAGGGACCCTGGTCACCGTCTCCTCAG";
  string H_J3_1="TGATGCTTTTGATGTCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG";
  string H_J3_2="TGATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG";
  string H_J5_1="ACAACTGGTTCGACTCCTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCAG";
  string H_J5_2="ACAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG";
  string H_J1_1="GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG";
  string H_J2_1="CTACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG";
  string H_J6_2="ATTACTACTACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCA";
  string H_J6_3="ATTACTACTACTACTACTACATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCCTCA";
  string H_J6_1="ATTACTACTACTACTACGGTATGGACGTCTGGGGGCAAGGGACCACGGTCACCGTCTCCTCAG";
  string H_J6_4="ATTACTACTACTACTACGGTATGGACGTCTGGGGCAAAGGGACCACGGTCACCGTCTCCTCAG";

  string K_J1_1="GTGGACGTTCGGCCAAGGGACCAAGGTGGAAATCAAAC";
  string K_J3_1="ATTCACTTTCGGCCCTGGGACCAAAGTGGATATCAAAC";
  string K_J4_1="GCTCACTTTCGGCGGAGGGACCAAGGTGGAGATCAAAC";
  string K_J5_1="GATCACCTTCGGCCAAGGGACACGACTGGAGATTAAAC";
  string K_J2_1="TGTACACTTTTGGCCAGGGGACCAAGCTGGAGATCAAAC";

  string L_J1_1="TTATGTCTTCGGAACTGGGACCAAGGTCACCGTCCTAG";
  string L_J2_1="TGTGGTATTCGGCGGAGGGACCAAGCTGACCGTCCTAG";
  string L_J3_2="TTGGGTGTTCGGCGGAGGGACCAAGCTGACCGTCCTAG";
  string L_J6_1="TAATGTGTTCGGCAGTGGCACCAAGGTGACCGTCCTCG";
  string L_J7_1="TGCTGTGTTCGGAGGAGGCACCCAGCTGACCGTCCTCG";
  string L_J7_2="TGCTGTGTTCGGAGGAGGCACCCAGCTGACCGCCCTCG";
  
  string H_J4_1_rm="ACTACTTTGACTACTGGGGCCAGGGAGTCCTGGTCACCGTCTCCTCAG";
  string H_J3_1_rm="TGATGCTTTTGATTTCTGGGGCCAAGGGCTCAGGGTCACCGTCTCTTCAG";
  string H_J5_1_1_rm="ACAACCGGTTCGATGTCTGGGGCCCGGGAGTCCTGGTCACCGTCTCCTCAG";
  string H_J5_1_2_rm="ACAACTCATTGGATGTCTGGGGCCAGGGAGTTCTGGTCACCGTCTCCTCAG";
  string H_J5_2_1_rm="ACAACTCATTGGATGTCTGGGGCCGGGGAGTTCTGGTCACCGTCTCCTCAG";
  string H_J1_1_rm="GCTGAATACTTCGAGTTCTGGGGCCAGGGCGCCCTGGTCACCGTCTCCTCAG";
  string H_J2_1_rm="CTACTGGTACTTCGATCTCTGGGGCCCTGGCACCCCAATCACCATCTCCTCAG";
  string H_J6_1_rm="ATTACTACGGTTTGGATTCCTGGGGCCAAGGGGTCGTCGTCACCGTCTCCTCAG";
  string K_J1_1_1_rm="GTGGACGTTCGGCCAAGGGACCAAGGTGGAAATCAAAC";
  string K_J3_1_1_rm="ATTCACTTTCGGCCCCGGGACCAAACTGGATATCAAAC";
  string K_J4_1_1_rm="GCTCACTTTCGGCGGAGGGACCAAGGTGGAGATCAAAC";
  string K_J5_1_1_rm="GATCACCTTCGGCCAAGGGACACGACTGGAGATTAAAC";
  string K_J2_1_1_rm="TGTACAGTTTTGGCCAGGGGACCAAAGTGGAGATCAAAC";

  map<string, map<string, string> > J_genes;
  J_genes["human"]["IGHJ1*01"]=H_J1_1;
  J_genes["human"]["IGHJ2*01"]=H_J2_1;
  J_genes["human"]["IGHJ3*01"]=H_J3_1;
  J_genes["human"]["IGHJ3*02"]=H_J3_2;
  J_genes["human"]["IGHJ4*01"]=H_J4_1;
  J_genes["human"]["IGHJ4*02"]=H_J4_2;
  J_genes["human"]["IGHJ4*03"]=H_J4_3;
  J_genes["human"]["IGHJ5*01"]=H_J5_1;
  J_genes["human"]["IGHJ5*02"]=H_J5_2;
  J_genes["human"]["IGHJ6*01"]=H_J6_1;
  J_genes["human"]["IGHJ6*02"]=H_J6_2;
  J_genes["human"]["IGHJ6*03"]=H_J6_3;
  J_genes["human"]["IGHJ6*04"]=H_J6_4;

  J_genes["human"]["IGKJ1*01"]=K_J1_1;
  J_genes["human"]["IGKJ3*01"]=K_J3_1;
  J_genes["human"]["IGKJ4*01"]=K_J4_1;
  J_genes["human"]["IGKJ5*01"]=K_J5_1;
  J_genes["human"]["IGKJ2*01"]=K_J2_1;
  
  J_genes["human"]["IGLJ1*01"]=L_J1_1;
  J_genes["human"]["IGLJ2*01"]=L_J2_1;
  J_genes["human"]["IGLJ3*02"]=L_J3_2;
  J_genes["human"]["IGLJ6*01"]=L_J6_1;
  J_genes["human"]["IGLJ7*01"]=L_J7_1;
  J_genes["human"]["IGLJ7*02"]=L_J7_2;

  
  //rhesus
  J_genes["rhesus"]["IGHJ4*01"]=H_J4_1_rm;
  J_genes["rhesus"]["IGHJ3*01"]=H_J3_1_rm;
  J_genes["rhesus"]["IGHJ5-1*01"]=H_J5_1_1_rm;
  J_genes["rhesus"]["IGHJ5-1*02"]=H_J5_1_2_rm;
  J_genes["rhesus"]["IGHJ5-2*01"]=H_J5_2_1_rm;
  J_genes["rhesus"]["IGHJ1*01"]=H_J1_1_rm;
  J_genes["rhesus"]["IGHJ2*01"]=H_J2_1_rm;
  J_genes["rhesus"]["IGHJ6*01"]=H_J6_1_rm;
  
  J_genes["rhesus"]["IGKJ1-1*01"]=K_J1_1_1_rm;
  J_genes["rhesus"]["IGKJ2-1*01"]=K_J2_1_1_rm;
  J_genes["rhesus"]["IGKJ3-1*01"]=K_J3_1_1_rm;
  J_genes["rhesus"]["IGKJ4-1*01"]=K_J4_1_1_rm;
  J_genes["rhesus"]["IGKJ5-1*01"]=K_J5_1_1_rm;
  
  map<string, map<string, bool> > J_genes_too_long;
  J_genes_too_long["human"]["IGHJ1*01"]=true;
  J_genes_too_long["human"]["IGHJ2*01"]=true;
  J_genes_too_long["human"]["IGHJ3*01"]=true;
  J_genes_too_long["human"]["IGHJ3*02"]=true;
  J_genes_too_long["human"]["IGHJ4*01"]=true;
  J_genes_too_long["human"]["IGHJ4*02"]=true;
  J_genes_too_long["human"]["IGHJ4*03"]=true;
  J_genes_too_long["human"]["IGHJ5*01"]=true;
  J_genes_too_long["human"]["IGHJ5*02"]=true;
  J_genes_too_long["human"]["IGHJ6*01"]=true;
  J_genes_too_long["human"]["IGHJ6*02"]=false;
  J_genes_too_long["human"]["IGHJ6*03"]=false;
  J_genes_too_long["human"]["IGHJ6*04"]=true;
  
  J_genes_too_long["rhesus"]["IGHJ4*01"]=true;
  J_genes_too_long["rhesus"]["IGHJ3*01"]=true;
  J_genes_too_long["rhesus"]["IGHJ5-1*01"]=true;
  J_genes_too_long["rhesus"]["IGHJ5-1*02"]=true;
  J_genes_too_long["rhesus"]["IGHJ5-2*01"]=true;
  J_genes_too_long["rhesus"]["IGHJ1*01"]=true;
  J_genes_too_long["rhesus"]["IGHJ2*01"]=true;
  J_genes_too_long["rhesus"]["IGHJ6*01"]=true;
  
  J_genes_too_long["rhesus"]["IGKJ1-1*01"]=true;
  J_genes_too_long["rhesus"]["IGKJ2-1*01"]=true;
  J_genes_too_long["rhesus"]["IGKJ3-1*01"]=true;
  J_genes_too_long["rhesus"]["IGKJ4-1*01"]=true; 
  J_genes_too_long["rhesus"]["IGKJ5-1*01"]=true;

  J_genes_too_long["human"]["IGKJ1*01"]=true;
  J_genes_too_long["human"]["IGKJ3*01"]=true;
  J_genes_too_long["human"]["IGKJ4*01"]=true;
  J_genes_too_long["human"]["IGKJ5*01"]=true;
  J_genes_too_long["human"]["IGKJ2*01"]=true;
  
  J_genes_too_long["human"]["IGLJ1*01"]=true;
  J_genes_too_long["human"]["IGLJ2*01"]=true;
  J_genes_too_long["human"]["IGLJ3*02"]=true;
  J_genes_too_long["human"]["IGLJ6*01"]=true;
  J_genes_too_long["human"]["IGLJ7*01"]=true;
  J_genes_too_long["human"]["IGLJ7*02"]=true;

  vector<string> tokens;
  tokenize(markup_header,tokens,"|");

  if ((chain_type=="heavy")&&(tokens[tokens.size()-1].substr(0,4) != "IGHJ")){cerr << "SMUA file is incorrectly formatted, missing J gene name in markup header\n"; error_status=true; return;}
  if ((chain_type=="kappa")&&(tokens[tokens.size()-1].substr(0,4) != "IGKJ")){cerr << "SMUA file is incorrectly formatted, missing J gene name in markup header\n"; error_status=true; return;}
  if ((chain_type=="lambda")&&(tokens[tokens.size()-1].substr(0,4) != "IGLJ")){cerr << "SMUA file is incorrectly formatted, missing J gene name in markup header\n"; error_status=true; return;}
 
  string J_gene_name=tokens[tokens.size()-1];

  string J_gene_sequence=J_genes[species][J_gene_name];
  //make sure the last 10 chars in UCA sequence match the last 10 chars in J gene 
  bool UCA_matches_J10=true;
  int UCA_len=UCA_sequence.length();
  int j=J_gene_sequence.length()-1;
  for(int i=UCA_len-1; i>UCA_len-1-10; i--)
    {
      if (UCA_sequence[i] != J_gene_sequence[j]){UCA_matches_J10=false;}
      j--;
    }
  if (!UCA_matches_J10){cerr << "UCA does not match J gene sequence for " <<  J_gene_name << " in last 10 bases!\n";}
  
  // trim off the one base overhang when J overextends past last codon
  int trim_end=UCA_len-1;
  int seq_end=sequence.length()-1;
  int markup_end=markup.length()-1;

  if (J_genes_too_long[species][J_gene_name]){trim_end--; } 
  new_UCA_sequence=UCA_sequence.substr(0,trim_end+1);
  new_sequence=sequence.substr(0,trim_end+1);
  new_markup=markup.substr(0,trim_end+1);

  //fill in end gaps with UCA sequence
  for(int i=new_sequence.length()-1; i>=0; i--)
    {
      if (new_sequence[i]!='-'){break;}
      new_sequence[i]=new_UCA_sequence[i];
      if (new_markup[i]=='U'){new_markup[i]='4';}
      number_of_replacements++;
    }

  //deal with indels if there are any
  bool indels_present=false;
  for(int i=0; i<new_sequence.length(); i++)
    {
      if ((new_sequence[i] == '-') || (new_UCA_sequence[i]=='-')){indels_present=true;}
    }

  if (sequence_has_ambiguities(new_sequence)){cerr << "ERROR: " << sequence_name << " has ambiguities, cannot proceed\n"; error_status=true; return;}
  if (indels_present == false){return;} //no indels, we're done
  
  //for now trust that Cloanalyst makes accurate alignments
  //the rest assumes that the UCA sequence starts in frame 1
  //INSERTION IN OBS
  for(int i=1; i<new_UCA_sequence.length()-1; i++)
    {
      //Fix such that gaps cover codons without overhang
      if ((new_UCA_sequence[i]=='-')&&(new_UCA_sequence[i-1]!='-'))
	{
	  int gap_start=i;
	  int j=i;
	  while(j<new_UCA_sequence.length())
	    {
	      if (new_UCA_sequence[j]=='-'){j++;}
	      else{break;}
	    }
	  int gap_end=j-1;
	  int gap_length=gap_end-gap_start+1;
	  if (gap_length%3!=0){cerr << "ERROR1: " << sequence_name << " has gap at " << gap_start << " which is not multiple of 3\n"; error_status=true; return;}

	  int frame_of_gap=gap_start%3;
	  if (frame_of_gap==0){continue;}
	  else if (frame_of_gap==1) //shift 3' base to right
	    {
	      char NU_hold=new_UCA_sequence[gap_start-1];
	      char NM_hold=new_markup[gap_start-1];
	      new_UCA_sequence[gap_start-1]='-';
	      new_markup[gap_start-1]='-';
	      new_UCA_sequence[gap_end]=NU_hold;
	      new_markup[gap_end]=NM_hold;
	    }
	  else if (frame_of_gap==2) //shift 5' base to left
	    {
	      char NU_hold=new_UCA_sequence[gap_end+1];
	      char NM_hold=new_markup[gap_end+1];
	      new_UCA_sequence[gap_end+1]='-';
	      new_markup[gap_end+1]='-';
	      new_UCA_sequence[gap_start]=NU_hold;
	      new_markup[gap_start]=NM_hold;
	    }
	}
    }

   //DELETION IN OBS
  for(int i=1; i<new_sequence.length()-1; i++)
    {
      //Fix such that gaps cover codons without overhang
      if ((new_sequence[i]=='-')&&(new_sequence[i-1]!='-'))
	{
	  int gap_start=i;
	  int j=i;
	  while(j<new_sequence.length())
	    {
	      if (new_sequence[j]=='-'){j++;}
	      else{break;}
	    }
	  int gap_end=j-1;
	  int gap_length=gap_end-gap_start+1;
	  if (gap_length%3!=0){cerr << "ERROR2: " << sequence_name << " has gap at " << gap_start << " which is not multiple of 3\n"; error_status=true; return;}

	  int frame_of_gap=gap_start%3;
	  if (frame_of_gap==0){continue;}
	  else if (frame_of_gap==1) //shift 3' base to right
	    {
	      char NU_hold=new_sequence[gap_start-1];
	      new_sequence[gap_start-1]='-';
	      new_sequence[gap_end]=NU_hold;
	    }
	  else if (frame_of_gap==2) //shift 5' base to left
	    {
	      char NU_hold=new_sequence[gap_end+1];
	      char NM_hold=new_markup[gap_end+1];
	      new_sequence[gap_end+1]='-';
	      new_sequence[gap_start]=NU_hold;
	    }
	}
    }
  
  return;
  

}

bool sequence_has_ambiguities(string sequence)
{
  for(int i=0; i<sequence.length(); i++)
    {
      if ((sequence[i]!='A')&&(sequence[i]!='C')&&(sequence[i]!='T')&&(sequence[i]!='G')&&(sequence[i]!='-')){ return true;}
    }
  return false;
}
