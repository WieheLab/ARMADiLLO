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
#include <algorithm>
#include <boost/filesystem/operations.hpp>
//#include "boost/filesystem.hpp"
#include <sys/types.h>
#include <dirent.h>
#include <thread>
#include <mutex>
//#include <boost/filesystem/fstream.hpp>
//#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/optional.hpp>
// include headers that implement a archive in simple text format
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
using namespace std;
using namespace boost;
using namespace boost::filesystem;

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
  double rank;
  double percentile;
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

struct Arguments
{
  bool ignore_CDR3=false;
  bool ignoreV=false;
  bool ignoreJ=false;
  bool quick=false;
  string species="human";
  string chain_type="heavy";
  string input_UCA_sequence="";
  string aaMuts="";
  string outputMode="HTML";
  string outDirectory="";
  vector<double> color_ladder{0.0001, 0.001, 0.01, 0.02, 0.10, 0.20, 0.5, 1};
  vector<double> color_rank_ladder{0.001, 0.01, 0.05, 0.1, 0.25, 0.50, 0.75, 1};
  int numbMutations=-1;
  double low_prob_cutoff=.02;
  bool clean_SMUA_first=false, remutate=false, output_seqs=false, ignore_warnings=false;
  bool lineage=false;
  int branches=1;
  int line_wrap_length=60, max_iter=1000,replace_J_upto=0;
  bool annotateFlag=false;
  bool rank=false;
  std::mt19937 gen;
  std::uniform_real_distribution<double> dis;
};

void process_fasta_sequence_to_seq_vector(string &,vector<Seq> &, map<string,string> &, map<string,S5F_mut> &);
void process_SMUA_sequence_to_seq_vector(string &, string &, vector<Seq> &, map<string,string> &, map<string, S5F_mut> &);
int simulate_S5F_mutation(string , int &, map<string,S5F_mut> &, mt19937 &, uniform_real_distribution<double> &, bool, vector<string> &, bool,  vector<bool> &);
void simulate_S5F_lineage(string , int, int &, map<string,S5F_mut> &, mt19937 &, uniform_real_distribution<double> &, bool, vector<string> &, bool,  vector<bool> &);
void print_output(string, vector<vector<Seq> > &, vector<string>, int, double);
void print_output_for_tiles_view(string, vector<vector<Seq> > &, vector<string>, int, double, vector<double> &,bool);
void print_freq_table_to_file(string,  map<int, map<char,double> > &);
void print_HTML_freq_table_to_file(string,  map<int, map<char,double> > &, string, vector<double> &);
///functions
//void read_SMUA_file(string, vector<vector<string> > &);
void load_S5F_files(string,string, map<string,S5F_mut> &);
void convert_2D_seq_vector_to_HTML_table(vector<vector<Seq> >&, vector<string> &, HTML::Table &, double &);
vector<double> estimate_S5F_mutation(string , int &, map<string,S5F_mut> &, mt19937 &, uniform_real_distribution<double> &, bool, vector<string> &, bool,  vector<bool> &);
vector<pair<char,double> > sort_map_into_pair_vctr(map<char,double> &);
bool mycompare(pair<char,double> A, pair<char,double> B){return A.second > B.second;}
void correct_for_fivemer_with_gap(int, string, string &);
void print_tile_view(string,vector<vector<Seq> > &, vector<string>, int, double, vector<double> &);
void get_mutability_scores(map<string,S5F_mut> &, string, int, bool, vector<bool> &, vector<double> &, vector<double> &, double &, double &);
void cleanup_SMUA_sequences(string, string, string , string , string , string &, string &, string &, string , string , int &, bool &);
bool sequence_has_ambiguities(string);
void convert_2D_seq_vector_to_HTML_table_for_tiles_view(vector<vector<Seq> >&, vector<string> &, HTML::Table &, double &, vector<double> &, int &,bool);
void replace_UCA_sequence_in_SMUA(string, string, string, string, string &, string &, string &, bool);
void helpMenu();
void read_V(const std::string &, map<string,map<int, map<char,double> >> &);
///templated functions
template <typename Type>
void vector2D_to_3D(vector<vector<Type> > &, int, vector<vector<vector<Type> > > &);
template <typename Type>
void vector1D_to_2D(vector<Type> &, int, vector<vector<Type> > &);
template <typename Type>
string convert_to_string(Type);
//map<string, map<string, string> >
string J_genes_list(string species,string IGname);

void run_entry(map<string,S5F_mut> &, map<string,string> &, vector<string> , map<string,map<int,map<char,double>>> &, Arguments &, map<string,vector<string>> &,map<string,vector<Seq>> &);
//function to run the entry object to generate simulations
void printTileStack(string filename, map<string,vector<Seq>> &seq_map, vector<string> sequence_names, int line_wrap_length, double low_prob_cutoff, vector<double> &color_ladder);
void convert_2D_seq_vector_to_HTML_table(vector<vector<Seq> >&v2, vector<string> &names, HTML::Table &html_table, double &low_prob_cutoff, vector<double> &color_ladder, int &counter);
