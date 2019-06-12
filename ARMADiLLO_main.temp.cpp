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
  string CDR_markup;
  bool isMut;
  map<char,double> all_simulated_aa_positional_frequencies_map;
  string SMUA_code;

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

///GLOBALS
int stop_codon_count;

///functions
void read_SMUA_file(string, vector<vector<string> > &);
void load_S5F_files(string,string, map<string,S5F_mut> &);
void process_fasta_sequence_to_seq_vector(string &,vector<Seq> &, map<string,string> &, map<string,S5F_mut> &);
void process_SMUA_sequence_to_seq_vector(string &, string &, vector<Seq> &, map<string,string> &, map<string, S5F_mut> &);
void convert_2D_seq_vector_to_HTML_table(vector<vector<Seq> >&, vector<string> &, HTML::Table &, double &);
void number_of_mutations_two_seqs(string &, string &, int &);
void simulate_S5F_mutation(string , int &, map<string,S5F_mut> &, mt19937 &, uniform_real_distribution<double> &, bool, vector<string> &, bool,  vector<bool> &);
vector<pair<char,double> > sort_map_into_pair_vctr(map<char,double> &);
bool mycompare(pair<char,double> A, pair<char,double> B){return A.second > B.second;}
void correct_for_fivemer_with_gap(int, string, string &);
void print_output(string, vector<vector<Seq> > &, vector<string>, int, double);
void print_tile_view(string,vector<vector<Seq> > &, vector<string>, int, double, vector<double> &);
void print_pct_progress(int, int, int);
void get_mutability_scores(map<string,S5F_mut> &, string, int, bool, vector<bool> &, vector<double> &, vector<double> &, double &, double &);
void print_freq_table_to_file(string,  map<int, map<char,double> > &);
void cleanup_SMUA_sequences(string, string, string , string , string , string &, string &, string &, string , string , int &, bool &);
bool sequence_has_ambiguities(string);
void convert_2D_seq_vector_to_HTML_table_for_tiles_view(vector<vector<Seq> >&, vector<string> &, HTML::Table &, double &, vector<double> &, int &);
void print_output_for_tiles_view(string, vector<vector<Seq> > &, vector<string>, int, double, vector<double> &);
void replace_UCA_sequence_in_SMUA(string, string, string, string, string &, string &, string &, bool);

///templated functions
template <typename Type>
void vector2D_to_3D(vector<vector<Type> > &, int, vector<vector<vector<Type> > > &);
template <typename Type>
void vector1D_to_2D(vector<Type> &, int, vector<vector<Type> > &);
template <typename Type>
string convert_to_string(Type);

///OLD TODO:
//0. The gap bases should have n/a for mutability score in HTML 
//1. Use the markup string to highlight CDRs in the HTML
//2. Make sure last cell in ladder is 1 always


int main(int argc, char *argv[])
{  
  if (argc <2){cout << "USAGE: analyze_mutations -SMUA [SMUA file] -w [line wrap length (60)] -m [S5F mutability file] -s [S5F substitution file] -max_iter [cycles of B cell maturation(100)] -c [cutoff for highlighting low prob (1=1%)] -replace_J_upto [number of replacements in J allowed] -chain [chain type (heavy=default|kappa|lambda)] -species [(human=default|rhesus)] -clean_first [clean the SMUA prior to running] -output_seqs [output sim seqs] -random_seed [provide a random seed]\n"; exit(1);}
 
  ///get cmdline args
  int i=0, line_wrap_length=60, max_iter=100, mutation_count_from_cmdline=-1, replace_J_upto=0, random_seed=0;
  string fasta_filename="", mutability_filename="", substitution_filename="", SMUA_filename="", species="human", chain_type="heavy", input_UCA_sequence="";
  int SMUA_start=0, SMUA_end=-1;
  double low_prob_cutoff=.02;
  bool ignore_CDR3=false, clean_SMUA_first=false, user_provided_random_seed=false, remutate=false, output_seqs=false, ignore_warnings=false;
   while(i<argc)
     {
       string arg=argv[i];
       string next_arg;
       if (i<argc-1){next_arg=argv[i+1];}else{next_arg="";}
       
        //  if ((arg.substr(0,1)=="-")&&(next_arg.substr(0,1)=="-")){cerr << "incorrectly formatted cmdline\n"; exit(1);}
       if (arg == "-SMUA")
	 {
	   SMUA_filename=next_arg;
	 }
       if (arg == "-m")
	 {
	   mutability_filename=next_arg;
	 }
       if (arg == "-s") 
	 {
	   substitution_filename=next_arg;
	 }
       if (arg == "-w")
	 {
	   line_wrap_length=atoi(next_arg.c_str());
	 }
       if (arg == "-max_iter")
	 {
	   max_iter=atoi(next_arg.c_str());
	 }
       if (arg == "-mut_count")
	 {
	   mutation_count_from_cmdline=atoi(next_arg.c_str());
	 }
       if (arg == "-c")
	 {
	   low_prob_cutoff=atof(next_arg.c_str())/100.0;
	 }
       if (arg == "-ignore_CDR3")
	 {
	   ignore_CDR3=true;
	 }
       if (arg == "-start")
	 {
	   SMUA_start=atoi(next_arg.c_str());
	 }
       if (arg == "-end")
	 {
	   SMUA_end=atoi(next_arg.c_str());
	 }
       if (arg == "-replace_J_upto")
	{
	  replace_J_upto=atoi(next_arg.c_str());
	}
      if (arg == "-species")
	{
	  species=next_arg;
	}
      if (arg == "-chain")
	{
	  chain_type=next_arg; 
	}
      if (arg == "-clean_first")
	{
	  clean_SMUA_first=true;
	}
      if (arg == "-random_seed")
	{
	  random_seed=atoi(next_arg.c_str());
	  user_provided_random_seed=true;
	}
      if (arg == "-remutate")
	{
	  bool remutate=true;
	  //cerr << "REMUTATE HAS NOT YET BEEN IMPLEMENTED. FEATURE COMING SOON!\n";
	}
      if (arg == "-output_seqs")
	{
	  output_seqs=true;
	}
      if (arg == "-UCA_sequence")
	{
	  input_UCA_sequence=next_arg;
	}
      if (arg == "-ignore_warnings")
	{
	  ignore_warnings=true;
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

   //define color ladder
   vector<double> color_ladder{0.0001, 0.001, 0.01, 0.02, 0.10, 0.20, 0.5, 1};

   ///load dna_to_aa map
   map<string,string> dna_to_aa_map;
   get_aa_tranx_map(dna_to_aa_map);

   ///read input sequence alignment
   map <string, string> sequences;
   vector <string> sequence_names;   
   vector<vector<string> > SMUA_alignments_and_markup;

   read_SMUA_file(SMUA_filename, SMUA_alignments_and_markup);
   
   ///load S5F files
   map <string, S5F_mut> S5F_5mers;
   load_S5F_files(mutability_filename,substitution_filename, S5F_5mers);

   cout << "NAME\t#AA_MUTS\t#MUTS\t<.02\t<.01\t<.001\t<.0001\t#INS\t#DEL\t#INDELS/3\tCDR3_LEN\n";
   //iterate through the SMUA file and perform mutation analysis for each sequence
   double total_elapsed_time=0;
   if (SMUA_end==-1){SMUA_end=SMUA_alignments_and_markup.size();}
   if (SMUA_end>SMUA_alignments_and_markup.size()){SMUA_end=SMUA_alignments_and_markup.size();}
   //for(int i=0; i<SMUA_alignments_and_markup.size(); i++)
   for(int i=SMUA_start; i<SMUA_end; i++)
     {
       clock_t begin=clock();
       string sequence_name=SMUA_alignments_and_markup[i][0];
       string sequence=SMUA_alignments_and_markup[i][1];
       string UCA_sequence_name=SMUA_alignments_and_markup[i][2];
       string UCA_sequence=SMUA_alignments_and_markup[i][3];
       string markup_header=SMUA_alignments_and_markup[i][4];
       string markup_string=SMUA_alignments_and_markup[i][5];

       string new_sequence="", new_UCA_sequence="", new_markup_string="";
       int number_of_replacements=0;
       bool error_status=false;

      
       if (clean_SMUA_first)
	 {
	   cerr << "Cleaning SMUA step\n"; 
	   cleanup_SMUA_sequences(sequence_name, markup_header, UCA_sequence, sequence, markup_string, new_UCA_sequence, new_sequence, new_markup_string, species, chain_type, number_of_replacements, error_status);
	   if (number_of_replacements>replace_J_upto)
	     {
	       cerr << "Sequence: " << sequence_name << " required " << number_of_replacements << " replacements at the end of the J where " << replace_J_upto << " is allowed (user defined). Skipping that sequence\n"; 
	       cout << sequence_name << "\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\n"; 
	       continue;
	     }
	   if (error_status){
	     cerr << "SMUA incorrectly formatted for sequence " << sequence_name << ". Skipping that sequence\n"; 
	     cout << sequence_name << "\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\n"; 
	     continue;
	   }
	   //cerr << new_UCA_sequence << "\n" << new_sequence << "\n"; 
	   UCA_sequence=new_UCA_sequence;
	   sequence=new_sequence;
	   markup_string=new_markup_string;
	 }

        if (input_UCA_sequence!="")//replace all UCAs with input UCA
	 {
	   new_UCA_sequence="", new_sequence="", new_markup_string="";
	   replace_UCA_sequence_in_SMUA(sequence, UCA_sequence, markup_string, input_UCA_sequence, new_sequence, new_UCA_sequence, new_markup_string, ignore_warnings);
	   string new_UCA_aa_sequence="",input_UCA_aa_sequence="",UCA_aa_sequence="";
	   translate_dna_to_aa(new_UCA_sequence, new_UCA_aa_sequence, 1, dna_to_aa_map);
	   translate_dna_to_aa(input_UCA_sequence, input_UCA_aa_sequence, 1, dna_to_aa_map);
	   translate_dna_to_aa(UCA_sequence, UCA_aa_sequence, 1, dna_to_aa_map);
	   cout << "inp_UCA: " << input_UCA_aa_sequence << "\n    UCA: " << UCA_aa_sequence << "\nnew_UCA: " << new_UCA_aa_sequence << "\n"; 
	   cout << ">inp_UCA\n" << input_UCA_sequence << "\n>UCA\n" << UCA_sequence << "\n>new_UCA\n" << new_UCA_sequence << "\n";
	   UCA_sequence=new_UCA_sequence;
	   sequence=new_sequence;
	   markup_string=new_markup_string;
	   //	   cout << UCA_sequence << "\n" << sequence << "\n" << markup_string << "\n"; 
	 } 
       

       string output_filename=sequence_name+".ARMADiLLO.html";
       string tiles_output_filename=sequence_name+".tiles.html";
 
       //if CDR3 is to be ignored, get CDR3 bases from markup string and store in bool vector
       vector<bool> shield_mutations(markup_string.length(), false);
       int shield_counter=0;
       string UCA_CDR3="", seq_CDR3="";
       for(int j=0; j<markup_string.length(); j++)
	 {
	   if ((markup_string[j]=='V')||(markup_string[j]=='n')||(markup_string[j]=='D')||(markup_string[j]=='J'))
	     {
	       shield_mutations[j]=true; 
	       shield_counter++;
	       UCA_CDR3+=UCA_sequence[j];
	       seq_CDR3+=sequence[j];
	     }
	   else
	     { 
	       shield_mutations[j]=false; 
	     }
	 }
       //convert dna markup string to aa frame
       //for(int j=0; j<cdr_dna_markup_string.length(); j+=3)
       // {
       //   cdr_markup_string+=cdr_dna_markup_string[j];
       //	 }
      
       if (ignore_CDR3){cerr << "Shielding CDR3 from mutation. CDR3 has " << shield_counter << " bases, mutability zeroed out for these\n"; }
       int CDR3_length=shield_counter;

       vector<Seq> seq_vector, UCA_seq_vector, aa_seq_vector, aa_UCA_seq_vector;
       process_SMUA_sequence_to_seq_vector(sequence, markup_string, seq_vector, dna_to_aa_map, S5F_5mers);
       
       process_SMUA_sequence_to_seq_vector(UCA_sequence, markup_string, UCA_seq_vector, dna_to_aa_map, S5F_5mers);
       
       //calc number of dna mutations
       int mut_count=0, CDR3_mut_count=0;
       number_of_mutations_two_seqs(UCA_sequence, sequence, mut_count);
       number_of_mutations_two_seqs(UCA_CDR3, seq_CDR3, CDR3_mut_count); //get mut count for CDR3

       //calc number of amino acid mutations
       int aa_mut_count=0, CDR3_aa_mut_count=0;
       string UCA_aa_sequence="", aa_sequence="", UCA_aa_CDR3="", seq_aa_CDR3="";
       
       translate_dna_to_aa(UCA_sequence, UCA_aa_sequence, 1, dna_to_aa_map);
       cout << "AA: " << UCA_aa_sequence << "\n"; 
       translate_dna_to_aa(sequence, aa_sequence, 1, dna_to_aa_map);
       translate_dna_to_aa(seq_CDR3, seq_aa_CDR3, 1, dna_to_aa_map);
       translate_dna_to_aa(UCA_CDR3, UCA_aa_CDR3, 1, dna_to_aa_map);
       number_of_mutations_two_seqs(UCA_aa_sequence, aa_sequence, aa_mut_count);
       number_of_mutations_two_seqs(UCA_aa_CDR3, seq_aa_CDR3, CDR3_aa_mut_count);
       
       //       cerr << "CDRs:\n" << aa_sequence << "\n" << cdr_markup_string << "\n"; 
       

       if (ignore_CDR3)
	 {
	   mut_count-=CDR3_mut_count;
	   aa_mut_count-=CDR3_aa_mut_count;
	   cerr << "ignoring " <<  CDR3_mut_count << " mutations\n" << UCA_CDR3 << "\n" << seq_CDR3 << "\n";
	 }
       cerr << "processing " << sequence_name << " which has  " << mut_count << " mutations\n"; 

       //check for any UCA weirdness
       //cout << sequence_name << "\n" << sequence << "\n" << UCA_sequence << "\n" << markup_string << "\n"; 
       if (dna_sequence_has_stop_codon_in_reading_frame(UCA_sequence))
	 {
	   cerr << "germline has stop codon...skipping this sequence\n"; 
	   cout << sequence_name << "\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\n"; 
	   continue;
	 }

       //number of indels 
       int insertion_count=0, deletion_count=0; 
       for(int j=0; j<UCA_sequence.length(); j++)
	 {
	   if (UCA_sequence[j] == '-'){insertion_count++;}
	   if (sequence[j] == '-'){deletion_count++;}
	 }

     
       
       //simulate maturation at mutation frequency = to observed
       cerr << "Simulating maturation...\n"; 
       vector<string> mature_mutant_sequences(max_iter);
       stop_codon_count=0;

       for(int j=1; j<=max_iter; j++)
	 {
	   print_pct_progress(j, max_iter, 1);
	   vector<string> mutant_sequences;
	   simulate_S5F_mutation(UCA_sequence, mut_count, S5F_5mers, gen, dis,true, mutant_sequences, ignore_CDR3, shield_mutations);
	   string aa_sequence;
	   translate_dna_to_aa(mutant_sequences[mutant_sequences.size()-1], aa_sequence, 1, dna_to_aa_map);
	   mature_mutant_sequences[j-1]=aa_sequence;
	 }
       cerr << "STOP CODON #: " << stop_codon_count << "\n"; 
       cerr << "done\n"; 

       //output simulated seqs if necessary
       if (output_seqs)
	 {
	   string seqs_fasta_file=sequence_name+".ARMADiLLO.simulated_seqs.fasta";
	    ofstream file_out;
	    file_out.open(seqs_fasta_file.c_str());
	    for(int j=0; j<mature_mutant_sequences.size(); j++)
	      {
		file_out << ">seq_" << j+1 << "\n"; 
		file_out << mature_mutant_sequences[j] << "\n"; 
	      }
	    file_out.close();
	 }

       ///get positional frequency of aa from simulated sequences (with same num maturation mutations)
       map<int, map<char,double> > mature_mutant_positional_aa_freqs;
       ///init map to 0
       for(int j=0; j<mature_mutant_sequences[0].length(); j++)
	 {
	   for(int k=0; k<amino_acids.size(); k++)
	     {
	       mature_mutant_positional_aa_freqs[j][amino_acids[k]]=0;
	     }
	 }

       for(int j=0; j<mature_mutant_sequences.size(); j++)
	 {
	   //cout << ">mutant" << i+1 << "\n" << mature_mutant_sequences[i] << "\n"; 
	   for(int k=0; k<mature_mutant_sequences[j].length(); k++)
	     {
	       mature_mutant_positional_aa_freqs[k][mature_mutant_sequences[j][k]]+=(1/(double)mature_mutant_sequences.size());
	     }
	 }

       ///annotate seq vector with positional aa freq
       for(int j=0; j<seq_vector.size(); j++) //per position
	 {
	   if ((UCA_seq_vector[j].base == "-") || (seq_vector[j].base== "-")) //indel special case
	     {
	       seq_vector[j].simulated_aa_positional_frequency=-99.99;
	     } 
	   else
	     {
	       seq_vector[j].simulated_aa_positional_frequency=mature_mutant_positional_aa_freqs[seq_vector[j].aa_num-1][seq_vector[j].aa[0]];
	     }
	   seq_vector[j].all_simulated_aa_positional_frequencies_map=mature_mutant_positional_aa_freqs[seq_vector[j].aa_num-1];
	 }
       
       ///tabulate 
       int p02_count=0, p01_count=0, p001_count=0, p0001_count=0;
       map<string, int> region_counts;
       cerr << seq_vector.size() << "\t" << UCA_seq_vector.size() << "\n"; 
       for(int j=0; j<seq_vector.size(); j+=3)
	 {
	   // cerr << "j: " << j << "\n"; 
	   
	   if ((UCA_seq_vector[j].aa=="-")||(UCA_seq_vector[j].base=="-")){continue;} //skip if insertion i.e. gap in UCA_seq sequence
	   if (j+1<seq_vector.size()){if(UCA_seq_vector[j+1].base=="-"){continue;}}
	   if (j+2<seq_vector.size()){if(UCA_seq_vector[j+2].base=="-"){continue;}}

	   if ((seq_vector[j].aa=="-")||(seq_vector[j].base=="-")){continue;} //skip if insertion i.e. gap in obs sequence  
	   if (j+1<seq_vector.size()){if(seq_vector[j+1].base=="-"){continue;}}
	   if (j+2<seq_vector.size()){if(seq_vector[j+2].base=="-"){continue;}}

	   if (shield_mutations[j]){continue;}
	   if (seq_vector[j].simulated_aa_positional_frequency<.02){p02_count++;}
	   if (seq_vector[j].simulated_aa_positional_frequency<.01){p01_count++;}
	   if (seq_vector[j].simulated_aa_positional_frequency<.001){p001_count++;}
	   if (seq_vector[j].simulated_aa_positional_frequency<.0001){p0001_count++;}
	 }

       if (mut_count==0) 
	 { 
	   cerr << "0 mutations found\n";   
	   cout << sequence_name << "\t" << aa_mut_count << "\t" << mut_count << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << insertion_count << "\t" << deletion_count << "\t" << (insertion_count+deletion_count)/3 << "\t" << CDR3_length << "\n"; 
	  
	 }
       else
	 {
	   cout << sequence_name << "\t" << aa_mut_count << "\t" << mut_count << "\t" << p02_count << "\t" << p01_count << "\t" << p001_count << "\t" << p0001_count << "\t" << insertion_count << "\t" << deletion_count << "\t" << (insertion_count+deletion_count)/3 << "\t" << CDR3_length <<  "\n"; 
	   ///print detailed ARMADiLLO output as HTML
	   vector<vector<Seq> > all_sequences;
	   all_sequences.push_back(UCA_seq_vector);
	   all_sequences.push_back(seq_vector);
	   vector<string> sequence_names;
	   sequence_names.push_back(UCA_sequence_name);
	   sequence_names.push_back(sequence_name);
	   print_output(output_filename, all_sequences, sequence_names, line_wrap_length, low_prob_cutoff);

	   ///print tiles as HTML
	   for(int j=0; j<seq_vector.size(); j+=3)
	     {
	       if (seq_vector[j].aa != UCA_seq_vector[j].aa){seq_vector[j].isMut=true;}else{seq_vector[j].isMut=false;}
	       UCA_seq_vector[j].isMut=false;
	       UCA_seq_vector[j].simulated_aa_positional_frequency=1;
	       aa_UCA_seq_vector.push_back(UCA_seq_vector[j]);
	       aa_seq_vector.push_back(seq_vector[j]);
	     }
	   vector<vector<Seq> > all_aa_sequences;
	   vector<string> aa_sequence_names;
	   all_aa_sequences.push_back(aa_UCA_seq_vector);
	   all_aa_sequences.push_back(aa_seq_vector);
	   aa_sequence_names.push_back("UCA");
	   aa_sequence_names.push_back(sequence_name);
	   print_output_for_tiles_view(tiles_output_filename, all_aa_sequences, aa_sequence_names, line_wrap_length, low_prob_cutoff, color_ladder);
	 }
       string output_freq_table=sequence_name+".freq_table.txt";
       print_freq_table_to_file(output_freq_table,mature_mutant_positional_aa_freqs);

       clock_t end=clock();
       double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
       cerr << "TIME: " << sequence_name << " took " << elapsed_secs << " to process\n"; 
       total_elapsed_time+=elapsed_secs;
     }
   cerr << "TOTAL ELAPSED TIME: " << total_elapsed_time << "\n"; 
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///                FUNCTION DEFINITIONS
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_output_for_tiles_view(string filename, vector<vector<Seq> > &all_sequences, vector<string> sequence_names, int line_wrap_length, double low_prob_cutoff, vector<double> &color_ladder)
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
      convert_2D_seq_vector_to_HTML_table_for_tiles_view(split_all_sequences[i],sequence_names,html_table, low_prob_cutoff, color_ladder, counter);
      //missing step: stylize the table
      
      //print HTML tables
      html_table.print(file_string);
      file_string+="<p></p>\n"; 
    }

  file_string+="<p><br></p><p align=\"center\"><img src=\"Mutation_Probability_legend.png\" alt=\"Mutation Probability Legend\" height=\"25\"></p>\n";
  file_string+="</body>\n</html>\n"; 

  ofstream file_out;
  file_out.open(filename.c_str());
  file_out << file_string;
  file_out.close();
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

void print_output(string filename, vector<vector<Seq> > &all_sequences, vector<string> sequence_names, int line_wrap_length, double low_prob_cutoff)
{
  vector<vector<vector<Seq> > > split_all_sequences;
  vector2D_to_3D(all_sequences,line_wrap_length, split_all_sequences);
  ///output html header
  string file_string="";
  file_string+="<html xmlns='http://www.w3.org/1999/xhtml' xml:lang='en'>\n"; 
  file_string+="<head>\n"; 
  file_string+=" <meta http-equiv='Content-Type' content='text/html; charset=utf-8' />\n";
  file_string+=" <title>Antibody Mutation Analysis</title>\n"; 
  file_string+=" <link rel='stylesheet' href='AMA.css' />\n"; 
  file_string+="</head>\n"; 
  file_string+="<body>\n";
  file_string+="<p></p><br>\n"; 
  //cerr << "number of splits: " << split_all_sequences.size() << "\n"; 
  //for each split, make a HTML table and print it
  for(int i=0; i<split_all_sequences.size(); i++)
    {
      HTML::Table html_table;
      html_table.hclass="results";
      //convert seq vector to html table
      convert_2D_seq_vector_to_HTML_table(split_all_sequences[i],sequence_names,html_table, low_prob_cutoff);
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

void get_mutability_scores(map<string,S5F_mut> &S5F_model, string sequence, int last_mutate_position, bool is_shielded, vector<bool> &shield_mutations, vector<double> &last_mut_scores, vector<double> &mut_scores, double &last_sum_mut_scores, double &sum_mut_scores)
{
  mut_scores.clear();
  
  if (last_mutate_position==-2) //start from scratch
    {
      mut_scores.push_back(1);///first two positions set to neutral
      mut_scores.push_back(1);///
      sum_mut_scores=2.0;
      for(int i=2; i<sequence.length()-2; i++)
	{
	  string fivemer=sequence.substr(i-2,5);
	  if (fivemer.find('-') != std::string::npos) //if there is a gap in fivemer, correct for it if not in middle pos
	    {
	      string new_fivemer="";
	      correct_for_fivemer_with_gap(i,sequence,new_fivemer);
	      fivemer=new_fivemer;
	    }
	  double mut_score;
	  
	  if (fivemer == "NO_SCORE"){mut_score=0;}
	  else if (S5F_model.find(fivemer)==S5F_model.end()){cerr << "ERROR: inside simulation, can't find fivemer " << fivemer << " in S5F model\n"; mut_score=0;}
	  else if ((is_shielded) && (shield_mutations[i]))//shield from mutation
	    {
	      mut_score=0;
	    }
	  else
	    {
	      mut_score=S5F_model[fivemer].score;
	    }
	  mut_scores.push_back(mut_score);
	  sum_mut_scores+=mut_score;
	}
      
      mut_scores.push_back(1);///last two positions set to neutral
      mut_scores.push_back(1);///
      sum_mut_scores+=2.0;
    }
  else //just recalculate in 5mer region around last mutate position
    {
      if (last_mutate_position<2 || last_mutate_position>=sequence.length()-2)
	{
	  
	}
      mut_scores=last_mut_scores;
      sum_mut_scores=last_sum_mut_scores;
      //only recompute the -2 to +2 windows
      //get indices of non-gap -2 +2, in case we're in the middle of a gapped region
     
      int five_prime_count=0, three_prime_count=0;
      string fivemer="";

      // cerr << "last mutate pos: " << last_mutate_position << "\n";
      vector<int> window_indices;
      window_indices.push_back(last_mutate_position);
      int j=last_mutate_position-1;
      while(five_prime_count<2 && j>=0) //5 prime side
	{
	   if (sequence[j] !='-')
	    {
	      //cerr << "5 prime: " << j << "\n"; 
	      window_indices.push_back(j);
	      five_prime_count++;
	    }
	   j--;
	}
      j=last_mutate_position+1;
      while(three_prime_count<2 && j<sequence.length()) //3 prime side
	{
	  if (sequence[j] !='-')
	    {
	      // cerr << "3 prime: " << j << "\n"; 
	      window_indices.push_back(j);
	      three_prime_count++;
	    }
	  j++;
	}
      sort(window_indices.begin(), window_indices.end()); //OPTIMIZE: does this really need to be sorted?

      //iterate through indices stored from previous step and recompute the fivemer score for these
      for(int i=0; i<window_indices.size(); i++)
	{
	  // cerr << "window indices " << i << "\t" << window_indices[i] << "\n"; 
	  int index=window_indices[i];
	  //if (index == -1){cerr << "FATAL ERROR: Fivemer fail in gappy UCA region\n"; cerr << "last mutate pos: " << last_mutate_position << "\n"; exit(1);}

	  if (index<2 || index >= sequence.length()-2){sum_mut_scores-=mut_scores[index]; mut_scores[index]=1; sum_mut_scores++; continue;}//if at edges because of following gaps

	  string fivemer=sequence.substr(index-2,5);
	  if (fivemer.find('-') != std::string::npos) //if there is a gap in fivemer, correct for it if not in middle pos
	    {
	      string new_fivemer="";
	      correct_for_fivemer_with_gap(index,sequence,new_fivemer);
	      fivemer=new_fivemer;
	    }
	  double mut_score;
	  
	  if (fivemer == "NO_SCORE"){mut_score=0;}
	  else if (S5F_model.find(fivemer)==S5F_model.end()){cerr << "ERROR: inside simulation, can't find fivemer " << fivemer << " in S5F model\n"; mut_score=0;}
	  else if ((is_shielded) && (shield_mutations[index]))//shield from mutation
	    {
	      mut_score=0;
	    }
	  else
	    {
	      mut_score=S5F_model[fivemer].score;
	    }
	  sum_mut_scores-=mut_scores[index]; //subtract out old score from +/- 2 window from sum
	  mut_scores[index]=mut_score; //new score
	  sum_mut_scores+=mut_score; //add back new score to sum
	}
    }
}

void simulate_S5F_mutation(string sequence, int &num_mutations, map<string,S5F_mut> &S5F_model, mt19937 &gen, uniform_real_distribution<double> &dis, bool kill_stop_seqs, vector<string> &mutant_sequences, bool is_shielded, vector<bool> &shield_mutations)
{
  if (num_mutations==0){mutant_sequences.push_back(sequence); return;}
  ///iterate num_mutations times
  int last_mutate_position=-2;
  vector<double> last_mut_scores(sequence.length(), 0.0);
  double last_sum_mut_scores=0;

  for(int j=1; j<=num_mutations; j++)
    {
      ///get mutability scores  
      //clock_t begin, end, total_start; double elapsed;
      //begin=clock();
      //total_start=begin;

      vector<double> mut_scores(sequence.length(), 0.0);
      double sum_mut_scores=0;
      get_mutability_scores(S5F_model, sequence, last_mutate_position, is_shielded, shield_mutations, last_mut_scores, mut_scores, last_sum_mut_scores, sum_mut_scores);
      last_mut_scores=mut_scores;
      last_sum_mut_scores=sum_mut_scores;

      //end=clock();
      //elapsed = double(end - begin);// / CLOCKS_PER_SEC;
      //cerr << "S5F SCORING: " << elapsed << "\n";
      //begin=clock();

       //check
      if(mut_scores.size() != shield_mutations.size()){cerr << "FATAL ERROR: shield mutations vector is not same size as mut_scores. Exiting\n"; exit(1);} 

      // cerr << j << "\tsum of mut scores: " << sum_mut_scores << "\n"; 
      
      ///convert mutability scores to probability of mutating position  
      vector<double> mut_probability_ladder(sequence.length(),0); ///cumulative distribution of probabilities
      mut_probability_ladder[0]=mut_scores[0]/(double) sum_mut_scores;
      // cerr << "0" << "\t" << sequence[0] << "\t" << (mut_scores[0]/(double) sum_mut_scores) << "\t" << mut_probability_ladder[0] << "\n"; 
      for(int i=1; i<sequence.length(); i++)
	{
	  mut_probability_ladder[i]=mut_probability_ladder[i-1]+(mut_scores[i]/(double) sum_mut_scores); 
	  //	    cerr << i << "\t" << sequence[i] << "\t" << (mut_scores[i]/(double) sum_mut_scores) << "\t" << mut_probability_ladder[i] << "\n"; 
	}
      //clock_t mut_ladder=clock();
      ///draw position randomly according to probability ladder
      double R=dis(gen);
      //   cerr << "R: " << R << "\n"; 
      int mutate_position_i=-1;
      for(int i=0; i<mut_probability_ladder.size(); i++)  ///OPTIMIZE: combine the two loops into one
	{
	  if (R<mut_probability_ladder[i])
	    {
	      mutate_position_i=i;
	      break;
	    }
	}
      //    int d; cin >> d;
      //end=clock();
      //elapsed = double(end - begin);// / CLOCKS_PER_SEC;
      //cerr << "LADDER BUILDING: " << elapsed << "\n";
      //catch when mutate_position_i is not set 
      if (mutate_position_i == -1)
	{
	  cerr << "FATAL ERROR: Setting mutation_position_i failed\n"; exit(1);
	}

      last_mutate_position=mutate_position_i;

      if (sequence.substr(mutate_position_i,1) == "-"){cerr << "mutate position should never be a gap\n"; int d; cin >> d; }

      //cerr << j << "\t" << R << "\t" << mutate_position_i << "\t" << sequence[mutate_position_i] << "\t" << mut_scores[mutate_position_i] << "\t"; 

      ///mutate position according to substitution model
      map<char, double> substitution_probs;
      map<char, double> substitution_probs_uniform;

      //make uniform sub model when we need it
      substitution_probs_uniform['A']=1.0/3.0; substitution_probs_uniform['C']=1.0/3.0; substitution_probs_uniform['G']=1.0/3.0; substitution_probs_uniform['T']=1.0/3.0; 
      substitution_probs_uniform[sequence[mutate_position_i]]=0;

      if ((mutate_position_i<=1) || (mutate_position_i>=sequence.length()-2))///edge cases
	{
	  substitution_probs=substitution_probs_uniform;
	}
      else
	{
	  string fivemer_to_mutate=sequence.substr(mutate_position_i-2,5);
	  if (fivemer_to_mutate.find('-') != std::string::npos)
	    {
	      string new_fivemer="";
	      correct_for_fivemer_with_gap(mutate_position_i,sequence,new_fivemer);
	      fivemer_to_mutate=new_fivemer;
	    }
	  if ((S5F_model.find(fivemer_to_mutate)==S5F_model.end()) || (fivemer_to_mutate=="NO_SCORE"))
	    {
	      cerr << "ERROR: could not find " << fivemer_to_mutate << " in S5F substitution model\n";
	      substitution_probs=substitution_probs_uniform;
	    }
	  else
	    {
	      substitution_probs=S5F_model[fivemer_to_mutate].substitutions;
	    }
	}
      
      // cerr << "A: " << substitution_probs['A'] << " C: " << substitution_probs['C'] << " G: " << substitution_probs['G'] << " T: " <<  substitution_probs['T'] << "\t"; 
      double R2=dis(gen);
      double cuml=0;
      char base_to_mutate_to='X';
      for(map<char,double>::iterator it = substitution_probs.begin(); it != substitution_probs.end(); ++it)
	{
	  if (R2<(cuml+it->second))
	    {
	      base_to_mutate_to=it->first;
	      break;
	    }
	  cuml+=it->second;
	}
      if (base_to_mutate_to == 'X'){
	cerr << "should not get an X ever, paused\n"; 
	cerr << "A: " << substitution_probs['A'] << " C: " << substitution_probs['C'] << " G: " << substitution_probs['G'] << " T: " <<  substitution_probs['T'] << "\n"; 
	cerr << j << "\t" << mutate_position_i << "\t" << sequence[mutate_position_i] << "\n"; 
	int d; cin >> d; 
      }
      string sequence_copy=sequence;
      sequence[mutate_position_i]=base_to_mutate_to;

      ///store 
      if (dna_sequence_has_stop_codon_in_reading_frame(sequence))
	{sequence=sequence_copy; j--; stop_codon_count++;}///discard if hits a stop codon and start from prev sequence
      else
	{
	  mutant_sequences.push_back(sequence);
	}
      //end=clock();
      //elapsed = double(end - begin);// / CLOCKS_PER_SEC;
      //cerr << "THE REST: " << elapsed << "\n";
      //cerr << "TOTAL: " << double(end - total_start); // / CLOCKS_PER_SEC << "\n"; 
      //int d; cin >> d; 
    }
  return;
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

void convert_2D_seq_vector_to_HTML_table(vector<vector<Seq> >&v2, vector<string> &names, HTML::Table &html_table, double &low_prob_cutoff)
{
  vector<char> amino_acids={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};

  HTML::Tr final_row, penultimate_row;
  final_row.hclass="arrow_row";
 
  final_row.cols.push_back(Td("absent","","","",""));
   
  for(int i=0; i<v2.size(); i++) ///iterate over rows
    {
      HTML::Tr row1,row2,row3; ///equal to the number of attributes in seq obj we want to print
      HTML::Td td("seq_name","","","",names[i]);
      td.rowspan="3";
      if (i==v2.size()-1){td.rowspan="4";}
      row1.cols.push_back(td);
      
      for(int j=0; j<v2[i].size(); j++) ///iterate over cols
	{
	  HTML::Td td1("","","","3",v2[i][j].aa);
	  HTML::Td td2("","","","3",convert_to_string(v2[i][j].aa_num));

	  char c_str[10];
	  sprintf(c_str,"%.2f",v2[i][j].S5F_mut_score);
	  string mut_score_str(c_str);
	  if (v2[i][j].S5F_mut_score == -1){mut_score_str="N/A";}

	  string str3="<div class=\"sm\">"+v2[i][j].base+"<br/>"+mut_score_str+"</div>"; //
	  HTML::Td td3("","","","",str3);
	  if (v2[i][j].S5F_mut_score>2){td3.hclass="highlighthotspot";}
	  if ((v2[i][j].S5F_mut_score<.3)&&(v2[i][j].S5F_mut_score!=-1)){td3.hclass="highlightcoldspot";}

	  char d_str[10];
	  sprintf(d_str,"%.5f",v2[i][j].simulated_aa_positional_frequency);
	  string pos_freq_num_str(d_str);
	  if (v2[i][j].simulated_aa_positional_frequency == -99.99){pos_freq_num_str=" N/A ";}
	  //get frequencies of all aa at this position and print for tooltip hover table
	  
	  string tooltip_table_str="";
	  vector<pair<char,double> > sorted_map_as_pair_vctr=sort_map_into_pair_vctr(v2[i][j].all_simulated_aa_positional_frequencies_map);
	  for(int k=0; k<sorted_map_as_pair_vctr.size(); k++)
	    { 
	      char aa_freq_str[20];
	      sprintf(aa_freq_str,"%c: %.3f<br>",sorted_map_as_pair_vctr[k].first, sorted_map_as_pair_vctr[k].second); 
	      tooltip_table_str+=aa_freq_str;
	    }
	  string pos_freq_str="<div class=\"tooltip\">"+pos_freq_num_str+"<span class=\"tooltiptext\">"+tooltip_table_str+"</span></div>";
	  HTML::Td tdP1("","","","3",pos_freq_str);
	  if ((v2[i][j].simulated_aa_positional_frequency<low_prob_cutoff)&&(v2[i][j].simulated_aa_positional_frequency!=-99.99)){tdP1.hclass="highlightaalowprob";}
	  //if (v2[i][j].simulated_aa_positional_frequency<low_prob_cutoff){tdP1.hclass="highlightaalowprob";}
	  HTML::Td tdF("absent","","","","");
	  if (i>0) ///highlight aa/dna mutation
	    {
	       if (v2[i][j].aa!=v2[0][j].aa){td1.hclass="highlightaamutation";}
	       if (v2[i][j].base!=v2[0][j].base)//dna mutation
		 {
		   td3.style="color:#4d0000";
		   if (i==v2.size()-1)
		     {
		       if ((v2[0][j].S5F_mut_score<.3)&&(v2[0][j].S5F_mut_score!=-1)){tdF.value="<img style='vertical-align:bottom' src='unusual.png' alt='Unsusual'/>";}
		       
		     }
		 }
	    }
	  if (j%3==0){row1.cols.push_back(td1);} //row 1: aa 
	  if (j%3==0){row2.cols.push_back(td2);} //row 2: aa number
	  row3.cols.push_back(td3); //row 3: base and mut score
	  if (i==v2.size()-1)
	    {
	      if (j%3==0){penultimate_row.cols.push_back(tdP1);}//pen row: pos aa freq
	      final_row.cols.push_back(tdF);
	    } //final rows
	}
      html_table.rows.push_back(row1);
      html_table.rows.push_back(row2);
      html_table.rows.push_back(row3);
    }
  html_table.rows.push_back(penultimate_row);
  html_table.rows.push_back(final_row);
  return;
}

void process_SMUA_sequence_to_seq_vector(string &sequence, string &markup_string, vector<Seq> &seq_vector, map<string,string> &dna_to_aa_map, map<string, S5F_mut> &S5F_5mers)
{
  seq_vector.clear();
  string aa="X";
  int aa_counter=0;
  if (sequence.length() != markup_string.length()){cerr << "WARNING: sequence and markup string lengths are not identical\n"; cerr << sequence << "\n" << markup_string << "\n"; }

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
      if ((markup_string[i]=='V')||(markup_string[i]=='n')||(markup_string[i]=='D')||(markup_string[i]=='J')){temp.CDR_markup="3";}
      else if (markup_string[i]=='B'){temp.CDR_markup="2";}
      else if (markup_string[i]=='A'){temp.CDR_markup="1";}
      else {temp.CDR_markup="0";}
      seq_vector.push_back(temp);
    }
  return;

}

void process_fasta_sequence_to_seq_vector(string &sequence,vector<Seq> &seq_vector, map<string,string> &dna_to_aa_map, map<string,S5F_mut> &S5F_5mers)
{
  seq_vector.clear();
  string aa="X";
  int aa_counter=0;
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
	   if (fivemer.find('-') != std::string::npos)
	     {
	       string new_fivemer="";
	       correct_for_fivemer_with_gap(i,sequence,new_fivemer);
	       fivemer=new_fivemer;
	     }
	   if (S5F_5mers.find(fivemer) == S5F_5mers.end()){ mut_score=-1; cerr << "ERROR 48: no 5mer for " << fivemer << " found in S5F\n"; }
	   else{mut_score=S5F_5mers[fivemer].score;}
	}
      temp.aa=aa;
      temp.aa_num=aa_counter;
      temp.S5F_mut_score=mut_score;
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
  // cerr << vector2D.size() << " by " << vector2D[0].size() << "\n"; 
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
      //cerr << i << "\n";
      vector<vector<Type> > temp;
      for(int j=0; j<vector2D.size(); j++)
	{
	  temp.push_back(all_sliced_rows[i+(j*max_num_splits)]);
	  //cerr << " " << j << "\n"; 
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
  if ((UCA_sequence.length() != sequence.length())||(sequence.length()!=markup.length())||(UCA_sequence.length()!=markup.length()))
    {
      cerr << "FATAL ERROR: UCA sequence, sequence and markup strings are not same length\n";
      cerr << "UCA: " << UCA_sequence << "\nSEQ: " << sequence << "\nMRK: " << markup << "\n"; 
      exit(1);
    }

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

void convert_2D_seq_vector_to_HTML_table_for_tiles_view(vector<vector<Seq> >&v2, vector<string> &names, HTML::Table &html_table, double &low_prob_cutoff, vector<double> &color_ladder, int &counter)
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
      if ((counter+j)%10==0)
	{
	  ostringstream s;
	  s<<counter+j;
	  ruler1.value=s.str();
	  ruler2.value="|";
	}
      else if ((counter+j)%5==0){ruler2.value="|";}
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

void replace_UCA_sequence_in_SMUA(string old_obs_sequence, string old_uca_sequence, string old_markup_string, string input_uca_sequence, string &new_obs_sequence, string &new_uca_sequence, string &new_markup_string, bool ignore_warnings)
{
  new_obs_sequence="", new_uca_sequence="", new_markup_string="";
  map<char,map<char,int> > scoring_matrix;
  load_EDNAFULL_matrix(scoring_matrix);
  double gap_open=-20, gap_extend=-2, score=0;
  string alignment_input_uca_sequence, alignment_old_uca_sequence;
  pairwise_align_sequences_semiglobal_w_affine_gap(scoring_matrix, input_uca_sequence, old_uca_sequence, gap_open, gap_extend, alignment_input_uca_sequence, alignment_old_uca_sequence, score);

  string status="no_indel";
  for(int i=0; i<alignment_input_uca_sequence.length(); i++)
    {
      if (alignment_input_uca_sequence[i]=='-'){status="new_del";}
      if (alignment_old_uca_sequence[i]=='-'){status="new_ins";}
    }

  //vector of the three strings to keep track of indels
  vector<vector<string> > SMUA_alignment, new_SMUA_alignment;
  for(int i=0; i<old_obs_sequence.length(); i++)
    {
      vector<string> temp;
      temp.push_back(string(1,old_obs_sequence[i]));
      temp.push_back(string(1,old_uca_sequence[i]));
      temp.push_back(string(1,old_markup_string[i]));
      SMUA_alignment.push_back(temp);
    }
  cout << "ALIGN\nNEW_UCA: " << alignment_input_uca_sequence << "\nOLD_UCA: " << alignment_old_uca_sequence << "\nEND\n";
  
  if (status=="no_indel")
    {
      new_uca_sequence = input_uca_sequence;
      new_obs_sequence = old_obs_sequence;
      new_markup_string = old_markup_string;
    }
  else //indels are present
    {
      if (ignore_warnings==false){
      cout << "NEW UCA ALERT: INDEL DETECTED RELATIVE TO OLD UCA. MANUAL CHECK IS ADVISED!\n"; 
      cout << "PRESS ANY KEY TO CONTINUE\n";
      int d; cin >>d;
      }
      int SMUA_counter=0;
      for(int i=0; i<alignment_old_uca_sequence.length(); i++)
	{
	  if (alignment_old_uca_sequence[i]=='-')//ins in new uca rel to old uca
	    {
	      int five_prime_non_gap=i-1;
	      if (five_prime_non_gap<0){cerr << "ERROR GAP TO START, CAN'T REPLACE UCA\n";}
	      int j=i;
	      while(alignment_old_uca_sequence[j]=='-') //within ins
		{
		  vector<string> temp;
		  temp.push_back("-");
		  temp.push_back(string(1,alignment_input_uca_sequence[j])); 
		  temp.push_back(SMUA_alignment[five_prime_non_gap][2]); //just continue region from 5p call prior to insertion
		  new_SMUA_alignment.push_back(temp);
		  j++;
		}
	      i=j;
	    }
	  else if (alignment_input_uca_sequence[i]=='-')//del in new uca rel to old uca
	    {
	      SMUA_counter++;
	      continue;
	    }
	  else
	    {
	      new_SMUA_alignment.push_back(SMUA_alignment[SMUA_counter]); //pushing old SMUA alignment which is wrong
	      //THIS IS WRONG!!!! Need to do a better job traversing alignment
	
	      SMUA_counter++;
	    }
	}
      for(int i=0; i<new_SMUA_alignment.size(); i++)
	{
	  new_obs_sequence+=new_SMUA_alignment[i][0];
	  new_uca_sequence+=new_SMUA_alignment[i][1];
	  new_markup_string+=new_SMUA_alignment[i][2];
	}
      cout << ">new_obs\n" << new_obs_sequence << "\n>new_uca\n" << new_uca_sequence << "\n>new_markup\n" << new_markup_string << "\n"; 
    }
  return;
}
