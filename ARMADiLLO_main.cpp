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

//#include <algorithm>
#include <boost/filesystem/operations.hpp>
//#include "boost/filesystem.hpp"
#include <sys/types.h>
//#include <dirent.h>
#include <thread>
#include <mutex>
//#include <boost/filesystem/fstream.hpp>
//#define BOOST_FILESYSTEM_NO_DEPRECATED
//#include <boost/optional.hpp>
// include headers that implement a archive in simple text format
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include "ARMADiLLO_main.hpp"
#include "HTML.hpp"
#include "utilities.hpp"
#include "nab.hpp"
#include "readInputFiles.hpp"

using namespace std;
using namespace boost;
using namespace boost::filesystem;

///GLOBALS
int stop_codon_count;

///OLD TODO:
//0. The gap bases should have n/a for mutability score in HTML 
//1. Use the markup string to highlight CDRs in the HTML
//2. Make sure last cell in ladder is 1 always

void helpMenu()
{
  cout << "ARMADiLLO\n\n";
  cout << "USAGE:\n\t ARMADiLLO [seq file options] -m [S5F mutability file] -s [S5F substitution file] <opt arguments>\n\n";
  cout << "required arguments:\n";
  cout << "\t -m [S5F mutability file]\n";
  cout << "\t -s [S5F substitution file]\n";
  cout << "Sequence Files options - either SMUA, partis or seq argument required:\n";
  cout << "\t -SMUA [SMUA file] : argument for SMUA file from Cloanalyst\n";
  cout << "\t -partis [partis file] : file from partis either yaml or cvs\n";
  cout << "\t -seq [seq fasta file] : fasta file containing sequences to process requires a uca file\n";
  cout << "\t -uca [uca fasta file] : UCA fasta can contain either 1 seq or matching sequences to the seq file\n";
  cout << "\t -markup [markup fasta file] : optional fasta for seq and uca sequence files\n";
  cout << "output arguments\n";
  //cout << "\t -no_text\n";
  cout << "\t -simple_text : flag to print out simple text files\n";
  cout << "\t -text : flag to print out all text files\n";
  cout << "\t -HTML : (default) flag to print out HTML files\n";
  cout << "\t -fulloutput : flag to print out all text and HTML files\n";
  cout << "\t -annotate : flag to print out annotation of the sequences\n";
  cout << "optional arguments:\n";
  //cout << "\t -(d)ir : sets the output directory of the files\n";
  cout << "\t -freq_dir [V, J Frequency file directory] : directory to pull the frequency tables for quick analysis\n";
  cout << "\t -amofile [amo file] : sets the amo file to use for the quick analysis\n";
  cout << "\t -resetamo   : flag to reset the amo file associated\n";
  cout <<"\t -w [line wrap length (60)]\n\t -max_iter [cycles of B cell maturation(1000)]\n\t -c [cutoff for highlighting low prob (1=1%)]\n\t -replace_J_upto [number of replacements in J allowed]\n\t -chain [chain type (heavy=default|kappa|lambda)]\n\t -species [(human=default|rhesus)]\n\t -(l)ineage [number of trees] : argument to generate the mutations through a lineage generation instead of linear generation\n\t -(n)umber [number of mutations] : arguemnt to set number of mutations to generate instead of taking from mutant sequence\n\t -clean_first : flag to turn on cleaning the SMUA prior to running\n\t -output_seqs : flag to turn on printing out simulated seq]\n";
  cout <<"\t -ignore_CDR3 : flag to ignore CDR3, default is false\n";
  cout <<"\t -ignore_V    : flag to ignore V, default is false\n";
  cout <<"\t -ignore_J    : lag to ignore J, default is false\n";
  cout << "\t -threads [number] : sets the number of threads to use during processing - default is number of processors\n";
  cout <<"\t -random_seed [provide a random seed]\n";
  exit(1);

  return;
}

int main(int argc, char *argv[])
{  
  if (argc <2){
    helpMenu();
  }
  ///get cmdline args
  int i=0, mutation_count_from_cmdline=-1, random_seed=0;
  string fasta_filename="", mutability_filename="", substitution_filename="", SMUA_filename="";
  string treefile="";
  int SMUA_start=0, SMUA_end=-1;
  string freq_dir="";
  //  bool ignore_CDR3=false,ignoreV=false,ignoreJ=false;
  bool resetARMOfile=false;
  bool user_provided_random_seed=false;
  string amoFile="";
  string UCAtype="cloanalyst";
  string ucaname="";
  string markupFile="";
  string  annotationFile;
  int num_threads=0, max_num_threads=thread::hardware_concurrency();
  bool estimate=false;
  bool reverse=false;
  string htmlStackName="";
 
  map<string,vector<string>> annotation;
  Arguments arguments;
  while(i<argc)
    {
      string arg=argv[i];
      string next_arg;
      if (arg=="-h" or arg=="-help")
	{
	  helpMenu();
	}
      if (i<argc-1){next_arg=argv[i+1];}else{next_arg="";}
      
      //  if ((arg.substr(0,1)=="-")&&(next_arg.substr(0,1)=="-")){cerr << "incorrectly formatted cmdline\n"; exit(1);}
      if (arg == "-SMUA")
	{
	 UCAtype="cloanalyst"; 
	  SMUA_filename=next_arg;
	 }
      if (arg=="-partis")
	{
	  UCAtype="partis";
	  SMUA_filename=next_arg;
	}
      if(arg=="-seq" or arg=="-seqs")
	{
	  UCAtype="individual";
	  SMUA_filename=next_arg;

	}
      if(arg=="-uca")
	{
	  ucaname=next_arg;	  
	}
      if(arg=="-markup")
	{
	  markupFile=next_arg;
	}
      if (arg == "-m")
	{
	  mutability_filename=next_arg;
	}
      if (arg == "-s") 
	 {
	   substitution_filename=next_arg;
	 }
      if (arg == "-freq_dir")
	{
	  freq_dir=next_arg;
	  arguments.quick=true;

	  char sep = '/';
          #ifdef _WIN32
	  sep = '\\';
          #endif
          amoFile=freq_dir+sep+"Freq_Table.amo";
	}
      if(arg=="-amofile" or arg=="-amo_file")
	{
	  amoFile=next_arg;
	  arguments.quick=true;
	  if(!fexists(amoFile))
	    {
	      map<string,map<int, map<char,double> >>  v_input;
	      ofstream outfs(amoFile);
	      boost::archive::binary_oarchive outa(outfs);
	      outa << v_input;
	      //cout << amoFile<< " can not be found\n";
	      //exit(1);
	    }
	}
      if (arg == "-w")
	{
	  arguments.line_wrap_length=atoi(next_arg.c_str());
	}
      if (arg == "-max_iter")
	{
	  arguments.max_iter=atoi(next_arg.c_str());
	}
      if (arg == "-number" or arg=="-n")
	{
	  arguments.numbMutations=atoi(next_arg.c_str());
	}
      if (arg == "-dir" or arg=="-d")
	{
	  arguments.outDirectory=next_arg.c_str();
	  
	}
      if (arg == "-lineage" or arg == "-l")
	{
	  arguments.lineage=true;
	  arguments.branches=atoi(next_arg.c_str());
	}
      if (arg == "-mut_count")
	{
	  mutation_count_from_cmdline=atoi(next_arg.c_str());
	  //where did this argument come from?
	}
      if(arg == "-no_text" || arg=="-notext")
	{
	  cout << "no output!"<<endl;
	  getchar();
	  arguments.outputMode="none";
	}
      if(arg == "-simple_text" || arg=="-simpletext")
	{
	  arguments.outputMode="simple";
	}
      if(arg == "-text" || arg=="-fulltext")
	{
	  arguments.outputMode="fulltext";
	}
      if(arg== "-HTML" || arg=="-html" || arg=="-web")
	{
	  arguments.outputMode="HTML";
	}
      if(arg== "-fulloutput" || arg=="-alloutput")
	{
	  arguments.outputMode="all";
	}
      if(arg=="-annotate")
	{
	  cout << "turning annotation on"<<endl;
	  arguments.annotateFlag=true;
	}
      if (arg == "-c")
	{
	  arguments.low_prob_cutoff=atof(next_arg.c_str())/100.0;
	}
      if (arg == "-ignore_CDR3" || arg =="-ignoreCDR3" || arg =="-ignorecdr3")
        {
          arguments.ignore_CDR3=true;
        }
      if(arg == "-ignoreJ" || arg=="-ignore_J" || arg=="-ignore_j"|| arg=="-ignorej")
        {
          arguments.ignoreJ=true;
        }
      if(arg == "-ignoreV" || arg =="-ignore_V"|| arg=="-ignore_v"|| arg=="-ignorev")
        {
          arguments.ignoreV=true;
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
	  arguments.replace_J_upto=atoi(next_arg.c_str());
	}
      if (arg == "-species")
	{
	  arguments.species=next_arg;
	}
      if (arg == "-chain")
	{
	  arguments.chain_type=next_arg; 
	}
      if (arg == "-clean_first")
	{
	  arguments.clean_SMUA_first=true;
	}
      if (arg == "-random_seed" && arg == "-seed")
	{
	  random_seed=atoi(next_arg.c_str());
	  user_provided_random_seed=true;
	}
      if(arg == "-threads" || arg == "-num_threads")
	{
	  num_threads=atoi(next_arg.c_str());
	}
      if (arg == "-remutate")
	{
	  bool remutate=true;
	  //cerr << "REMUTATE HAS NOT YET BEEN IMPLEMENTED. FEATURE COMING SOON!\n";
	}
      if (arg == "-output_seqs")
	{
	  arguments.output_seqs=true;
	}
      if (arg == "-UCA_sequence")
	{
	  arguments.input_UCA_sequence=next_arg;
	}
      if (arg == "-ignore_warnings")
	{
	  arguments.ignore_warnings=true;
	}
      if(arg == "-resetamo" or arg == "-reset_amo" or arg == "-reload_amo")
	{
	  cout << "reseting AMO file\n";
	  resetARMOfile=true;
	}
      if(arg=="-findpairs")
	{
	  arguments.aaMuts=next_arg;
	}
      if(arg=="-rank")
	{
	  arguments.rank=true;
	}
      if(arg=="-reverse")
	{
	  reverse=true;
	}
      if(arg=="-stack")
	{
	  htmlStackName=next_arg;
	}
      i++;
    }
  
  if( SMUA_filename.size()<1 || !fexists(SMUA_filename))
     {
       cout << "Error in sequence file\n\n";
       helpMenu();
     }
   if(mutability_filename.size()<1 || !fexists(mutability_filename))
     {
       cout << "Error in mutability file\n\n";
       helpMenu();
     }
   if(substitution_filename.size()<1 || !fexists(substitution_filename))
     {
       cout << "Error in substitution file\n\n";
       helpMenu();
     }
   if ((num_threads==0)||(num_threads>max_num_threads)){num_threads=max_num_threads;}
  
   ///setup random num generation
   std::random_device rd;
   int seed;
   if (user_provided_random_seed)
     {seed=random_seed;}
   else
     {seed=rd();}
   
   std::mt19937 gen(seed);
   arguments.gen=gen;
   std::uniform_real_distribution<double> dis(0, 1);
   arguments.dis=dis;

   //define color ladder
   //vector<double> color_ladder{0.0001, 0.001, 0.01, 0.02, 0.10, 0.20, 0.5, 1};
   ///load dna_to_aa map
   map<string,string> dna_to_aa_map;
   get_aa_tranx_map(dna_to_aa_map);

   if(!fexists("AMA.css") && (arguments.outputMode=="HTML" || arguments.outputMode=="all"))
     writeAMA();
   if(!fexists("sequence_color.css") && (arguments.outputMode=="HTML" || arguments.outputMode=="all"))
     writeColor();
   
   ///read input sequence alignment
   map <string, string> sequences;
   vector <string> sequence_names;   
   vector<vector<string> > SMUA_alignments_and_markup;

   if (UCAtype.compare("cloanalyst")==0)
     {
       read_SMUA_file(SMUA_filename, SMUA_alignments_and_markup);
     }
   else if (UCAtype.compare("partis")==0)
     {
       int pos =SMUA_filename.find_last_of(".");
       string ext =SMUA_filename.substr(pos+1);
       if(string("csv").compare(ext)==0)
	 {
	   read_PARTIScsv_file(SMUA_filename,SMUA_alignments_and_markup);
	 }
       else if(string("yaml").compare(ext)==0)
	 {
	   read_PARTISyaml_file(SMUA_filename,SMUA_alignments_and_markup);
	 }
       else
	 {
	   cerr << "Unidentified partis file type.\n ARMADiLLO supports yaml or csv file types from partis"<<endl;
	   exit(1);
	 }
     }
   else if(UCAtype.compare("individual")==0)
     {
       if( ucaname.size()<1 || !fexists(ucaname))
	 {
	   cout << "Error in uca file\n\n";
	   helpMenu();
	 }
       cout << "mark up: "<<markupFile<<endl;
       if(markupFile.size()<1 || !fexists(markupFile))
	 readIndividualFiles(SMUA_filename,ucaname,SMUA_alignments_and_markup);
       else
	 readIndividualFiles(SMUA_filename,ucaname,markupFile,SMUA_alignments_and_markup);
     }
   else
     {
       cout << "unknown UCA collection sequence"<<endl;
       cout << "exiting"<<endl;
       exit(1);
       return 0;
     }

   if(reverse)//records name either reverse or forward
     {
       for(int j=SMUA_alignments_and_markup.size()-1;  j>-1; j--)
	 {
	   sequence_names.push_back(SMUA_alignments_and_markup[j][0]);
	 }
     }
   else
     {
       for(int j=0;  j<SMUA_alignments_and_markup.size(); j++)
	 {
	   sequence_names.push_back(SMUA_alignments_and_markup[j][0]);
	 }
     }
   
   cerr << "highlighting residues with less than " << arguments.low_prob_cutoff << " probability for mutation\n"; 
   ///amino acids vector
   vector<char> amino_acids={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
   
   ///load S5F files
   map <string, S5F_mut> S5F_5mers;
   load_S5F_files(mutability_filename,substitution_filename, S5F_5mers);
   const std::string freq_directory=freq_dir;
   map<string,map<int, map<char,double> >>  v_input;
   if (arguments.quick==true)
     {
       arguments.ignore_CDR3=true;

       if(fexists(amoFile) && !resetARMOfile)
	 {
	   clock_t begin=clock();
	   std::ifstream infs(amoFile,std::ios::binary);
	   boost::archive::binary_iarchive ina(infs);
	   ina >> v_input;
	   clock_t end=clock();
	   double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	 }
       else
	 {     
	   std::mutex mtx;
	   mtx.lock();
	   read_V(freq_directory, v_input);
	   mtx.unlock();
	 }
     }
   cout << "NAME\t#AA_MUTS\t#MUTS\t<.02\t<.01\t<.001\t<.0001\t#INS\t#DEL\t#INDELS/3\tCDR3_LEN\tsum(log(P))\n";
   //map<string, map<string, string> > J_genes=J_genes_list();
   //iterate through the SMUA file and perform mutation analysis for each sequence
   //double total_elapsed_time=0;
   if (SMUA_end==-1)
     {
       SMUA_end=SMUA_alignments_and_markup.size();
     }
   if (SMUA_end>SMUA_alignments_and_markup.size())
     {
       SMUA_end=SMUA_alignments_and_markup.size();
     }
   //for(int i=0; i<SMUA_alignments_and_markup.size(); i++)
   int MAX_THREADS=num_threads;
   int size=SMUA_end;
   if(arguments.annotateFlag)
     {
       for(int i=SMUA_start;i<SMUA_end;i++)
	 {
	   vector<string> tmp;
	   annotation[SMUA_alignments_and_markup[i][0]]=tmp;
	 }
     }

   map<string, vector<Seq> > seq_map;
   for(int i=SMUA_start;i<SMUA_end;i+=MAX_THREADS)
     {
       int number_of_threads=min(MAX_THREADS,size-i);
       vector<thread> thread_list;
       for(int j=0;j<number_of_threads;j++)
	 {
	   thread_list.push_back(thread(run_entry,std::ref(S5F_5mers),std::ref(dna_to_aa_map),SMUA_alignments_and_markup[i+j],std::ref(v_input),std::ref(arguments),std::ref(annotation),std::ref(seq_map)));
	 }
       for(int j=0;j<number_of_threads;j++)
	 {
	   thread_list[j].join();
	 }
       //getchar();
     }

   if (arguments.quick==true)
     {
       ofstream outfs(amoFile);
       boost::archive::binary_oarchive outa(outfs);
       outa << v_input;
     }

   if(arguments.annotateFlag)
     {
       int pos =SMUA_filename.find_last_of(".");
       annotationFile=SMUA_filename.substr(0,pos)+".annotation.txt";

       ofstream file_out;
       file_out.open(annotationFile.c_str());
       
       for(int i=SMUA_start;i<SMUA_end;i++)
	 {
	   file_out << SMUA_alignments_and_markup[i][0]<<";";
	   for(int j=0;j<annotation[SMUA_alignments_and_markup[i][0]].size();j++)
	     {
	       if(j==0)
		 file_out << annotation[SMUA_alignments_and_markup[i][0]][j];
	       else
		 file_out << ","<<annotation[SMUA_alignments_and_markup[i][0]][j];
	     }
	   file_out << endl;
	 }
     }

   //sequence color code - need to build the inputs
   //HTML::Table html_table;
   //html_table.hclass="results";
   if(htmlStackName.size()>1)//test if html stack name is given
     {
       int len=htmlStackName.size();//if name doesn't end in html, add it
       if(len > 3 && htmlStackName.substr(len - 4)!="html")
	{
	  htmlStackName=htmlStackName+".html";
	}
       //pass to function to take the seq_map and create the html
      printTileStack(htmlStackName, seq_map, sequence_names, arguments.line_wrap_length, arguments.low_prob_cutoff, arguments.color_ladder);
     }

   //cerr << "TOTAL ELAPSED TIME: " << total_elapsed_time << "\n"; 
   return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///                FUNCTION DEFINITIONS
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void run_entry(map<string,S5F_mut> &S5F_5mers,map<string,string> &dna_to_aa_map, vector<string>  SMUA_entries, map<string,map<int, map<char,double> >> &v_input, Arguments &arg, map<string,vector<string>> &annotation,map<string,vector<Seq>> &seq_map)
{
  //double total_elapsed_time=0;
  //clock_t begin=clock();
  
  string seq=SMUA_entries[1];
  for (auto & c: seq) c = toupper(c);
    
  int Acount = countChar(seq,'A');
  int Gcount = countChar(seq,'G');
  int Ccount = countChar(seq,'C');
  int Tcount = countChar(seq,'T');
  int dashcount = countChar(seq,'-');

  if(seq.size()>Acount+Gcount+Ccount+Tcount+dashcount)
    {
      cerr << "Unidentified nucleotide in sequence: "<<SMUA_entries[0]<<endl;
      cerr <<"ARMADiLLO cannot accept amino acid sequences or nucleotide sequences with ambuguity codes"<<endl;
      cerr << "Skipping sequence: "<<SMUA_entries[0]<<endl;
      writeError(SMUA_entries[0]+".tiles.html",SMUA_entries[0]);
      writeError(SMUA_entries[0]+".ARMADiLLO.html",SMUA_entries[0]);
      return;
    }
    
  NabEntry nab(SMUA_entries,arg);
  nab.rank=arg.rank;
  
  if(arg.clean_SMUA_first)
    {
      if (nab.cleanSMUA(arg.species,arg.chain_type,arg.replace_J_upto)<0)
	return;
      vector<bool> _shield_mutations(nab.markup_string.length(),false);
      nab.shield_mutations=_shield_mutations;
    }
  if(arg.input_UCA_sequence!="")
    {
      nab.replaceUCA(dna_to_aa_map,arg.input_UCA_sequence,arg.ignore_warnings);
      nab.markup_mask=nab.parseMarkup();
    }
  //cout << arg.ignoreV<<"\t"<<arg.ignore_CDR3<<"\t"<<arg.ignoreJ<<endl;
  nab.createShield(arg.ignore_CDR3,arg.ignoreV,arg.ignoreJ);
  if(arg.quick)
    {
      //nab.Jgene_sequence=J_genes_list(arg.species,nab.Jgene);
      nab.replaceTable(S5F_5mers,dna_to_aa_map, v_input,arg);
      if(arg.output_seqs)//write out simulated sequences if argument is true
	cerr << "No sequences generated during quick generations\n";
    }
  else
    {
      bool goodSim=nab.SimulateSequences(S5F_5mers,dna_to_aa_map,arg.gen,arg.dis,arg.max_iter,arg.branches,arg.lineage);
      if(!goodSim)
	{
	  nab.printlog();
	  return;
	}
      if(arg.output_seqs)//write out simulated sequences if argument is true
	{
	  nab.outputSimSeqs(arg.max_iter,arg.branches);
	}
    }
  if(arg.aaMuts!="")
    {
      nab.countAAPairs(arg.aaMuts);
    }
  vector<string> SNPs;
  if(arg.rank)
    SNPs=nab.printResults(S5F_5mers,dna_to_aa_map,arg.line_wrap_length,arg.low_prob_cutoff,arg.color_rank_ladder);
  else
    SNPs=nab.printResults(S5F_5mers,dna_to_aa_map,arg.line_wrap_length,arg.low_prob_cutoff,arg.color_ladder);

  if(arg.annotateFlag)
    {
      annotation[nab.sequence_name]=SNPs;
    }
  nab.printlog();
  
  seq_map[nab.sequence_name]=nab.aa_out;//creates the aa_out map to pass to the tile out function


  //clock_t end=clock();
  //double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  //cerr << "TIME: " << nab.sequence_name << " took " << elapsed_secs << " to process\n"; 
  //total_elapsed_time+=elapsed_secs;
  return;
}

void print_output_for_tiles_view(string filename, vector<vector<Seq> > &all_sequences, vector<string> sequence_names, int line_wrap_length, double low_prob_cutoff, vector<double> &color_ladder,bool rank)
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
      convert_2D_seq_vector_to_HTML_table_for_tiles_view(split_all_sequences[i],sequence_names,html_table, low_prob_cutoff, color_ladder, counter,rank);
      //missing step: stylize the table
      
      //print HTML tables
      html_table.print(file_string);
      file_string+="<p></p>\n"; 
    }

  file_string+="<br><br><br>\n";
  /*file_string+="<p><table class=\"results\" align=center><tr><td class=\"noborder\"><font align=center size=\"4\"><b>Mutation Probability:&nbsp;</b></font></td>";
  file_string+="<td class=\"color_cat7 mut\"><div class=\"mm\">&nbsp;&nbsp;&nbsp;&nbsp;</div></td><td class=\"noborder\">&ge;20\%&nbsp;&nbsp;</td>";
  file_string+="<td class=\"color_cat6 mut\"><div class=\"mm\">&nbsp;&nbsp;&nbsp;&nbsp;</div></td><td class=\"noborder\">20-10\%&nbsp;&nbsp;</td>";
  file_string+="<td class=\"color_cat5 mut\"><div class=\"mm\">&nbsp;&nbsp;&nbsp;&nbsp;</div></td><td class=\"noborder\">10-2\%&nbsp;&nbsp;</td>";
  file_string+="<td class=\"color_cat4 mut\"><div class=\"mm\">&nbsp;&nbsp;&nbsp;&nbsp;</div></td><td class=\"noborder\">2-1\%&nbsp;&nbsp;</td>";
  file_string+="<td class=\"color_cat3 mut\"><div class=\"mm\">&nbsp;&nbsp;&nbsp;&nbsp;</div></td><td class=\"noborder\">1.0-0.1\%&nbsp;&nbsp;</td>";
  file_string+="<td class=\"color_cat2 mut\"><div class=\"mm\">&nbsp;&nbsp;&nbsp;&nbsp;</div></td><td class=\"noborder\">0.10-0.01\%&nbsp;&nbsp;</td> ";
  file_string+="<td class=\"color_cat1 mut\"><div class=\"mm\">&nbsp;&nbsp;&nbsp;&nbsp;</div></td><td class=\"noborder\">&lt;0.01\%</td></tr></table></p>\n";*/

  
  //file_string+="<p><br></p><p align=\"center\"><img src=\"Mutation_Probability_legend.png\" alt=\"Mutation Probability Legend\" height=\"25\"></p>\n";
  if(rank)
    file_string+="<p><table class=\"results\" align=center><tr><td class=\"noborder\"><font align=center size=\"4\"><b>Probability Rank:&nbsp;</b></font></td>";
  else
    file_string+="<p><table class=\"results\" align=center><tr><td class=\"noborder\"><font align=center size=\"4\"><b>Mutation Probability:&nbsp;</b></font></td>";
  for(int i=6;i>-1;i--)
    {
      if(i==6)
	file_string+="<td class=\"color_cat7 mut\"><div class=\"mm\">&nbsp;&nbsp;&nbsp;&nbsp;</div></td><td class=\"noborder\">&ge;"+formatDouble(color_ladder[i-1]*100)+"\%&nbsp;&nbsp;</td>";
      else if(i==0)
	file_string+="<td class=\"color_cat"+to_string(i+1)+ " mut\"><div class=\"mm\">&nbsp;&nbsp;&nbsp;&nbsp;</div></td><td class=\"noborder\">&lt;"+formatDouble(color_ladder[i]*100) +"\%</td>";
      else
	file_string+="<td class=\"color_cat"+to_string(i+1)+" mut\"><div class=\"mm\">&nbsp;&nbsp;&nbsp;&nbsp;</div></td><td class=\"noborder\">"+formatDouble(color_ladder[i]*100)+"-"+formatDouble(color_ladder[i-1]*100)+"\%&nbsp;&nbsp;</td> ";
    }
  file_string+="</tr></table></p>\n";
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

  /*  vector<string> mutation_list;
  vector<float> values_list;
  
  for(int j=0; j<positional_aa_freqs.size(); j++)
    {
      string mut;
      for(int i=0; i<amino_acids.size(); i++)
	{
	  mut=to_string(positional_aa_freqs[j][amino_acids[i]])+"_"+amino_acids[i]+to_string(j);
	  mutation_list.push_back(mut);
	  values_list.push_back(positional_aa_freqs[j][amino_acids[i]]);
	}
    }

  sort(mutation_list.begin(),mutation_list.end());
  for(int i=0; i<mutation_list.size();i++)
    {
      //cout << mutation_list[i]<<endl;
    }

  sort(values_list.begin(),values_list.end());
  vector<float>::iterator ip;
  ip=std::unique(values_list.begin(),values_list.begin()+values_list.size());
  values_list.resize(std::distance(values_list.begin(),ip));
  
  for(int i=0;i<values_list.size();i++)
    cout << values_list[i]<<endl;
  */
  
}


void print_HTML_freq_table_to_file(string filename,  map<int, map<char,double> > &positional_aa_freqs,string aa_sequence,vector<double> &color_ladder)
{
  vector<char> amino_acids={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
  string file_string="";
  file_string+="<html xmlns='http://www.w3.org/1999/xhtml' xml:lang='en'>\n"; 
  file_string+="<head>\n"; 
  file_string+=" <meta http-equiv='Content-Type' content='text/html; charset=utf-8' />\n";
  file_string+=" <title>Antibody Mutation Analysis</title>\n"; 
  file_string+=" <link rel='stylesheet' href='sequence_color.css' />\n";
  //file_string+=" <link rel='stylesheet' href='AMA.css' />\n";
  file_string+="<style>\n";
  file_string+=".tooltip{\n position: relative;\ndisplay: inline-block;\n}";
  file_string+=".tooltiptext {\n visibility: hidden;\n  width: 100px;\n  background-color: #fffff5;\n color: #000000;\n  font-size: 12px;\n font-family:'Monospace';\n  text-align: center;\n  padding: 5px 0;\n position: absolute;\n  top: 100%;\n  left: 50%;\n   margin-left: -30px;\n   z-index: 1;\n}";
  file_string+=".tooltiptext::after {\n   content: \" \";\n   position: absolute;\n    bottom: 100%;\n    left: 50%;\n    margin-left: -5px;\n    border-width: 5px;\n    border-style: solid;\n   border-color: transparent transparent black transparent;\n}";
  file_string+=".tooltip:hover .tooltiptext {\n    visibility: visible;\n}\n";
  file_string+=" </style>\n";
  file_string+="</head>\n"; 
  file_string+="<body>\n";
  

  file_string+="<table cellspacing=\"0\" border=\"1\">\n<colgroup span=\"1\" width=\"100\"></colgroup><colgroup span=\"20\" width=\"60\"></colgroup>\n<tr>\n";
  file_string+="<td height=\"20\" align=\"left\"></td>\n";
  for(int i=0; i<amino_acids.size(); i++)
    {
      string s;
      file_string+= "<td align=\"center\"><b>" +s.replace(0,1,1,amino_acids[i])+"</b></td>\n";
      }
  file_string+= "</tr>\n";
  
  for(int j=0; j<positional_aa_freqs.size(); j++)
    {
      //html_file_out << j+1;
     
      file_string+= "<tr><td height=\"20\" align=\"center\" sdval=\"3\"><b>"+to_string(j+1)+":"+aa_sequence[j] +"</b></td>\n";
      for(int i=0; i<amino_acids.size(); i++)
	{
	  char aa_char[10];
	  sprintf(aa_char,"%.2f%%",100*positional_aa_freqs[j][amino_acids[i]]);
	  string aa_str(aa_char);
	  
	  for(int k=0; k<color_ladder.size(); k++)
	    {
	      if (positional_aa_freqs[j][amino_acids[i]] <= color_ladder[k])
		{
		  ostringstream ss;
		  if(aa_sequence[j]=='X')
		    file_string+= "<td height=\"20\" align=\"center\" sdval=\"3\" class=\"color_cat8\"><div class=\"tooltip\">"+aa_str;
		  else
		    file_string+= "<td height=\"20\" align=\"center\" sdval=\"3\" class=\"color_cat"+to_string(k+1)+"\"><div class=\"tooltip\">"+aa_str;
		  break;
		}
	    }

	  char tooltip_char[75];
	  sprintf(tooltip_char,"pos:%d<br>%c->%c : %.2f%%",j+1,aa_sequence[j],amino_acids[i],100*positional_aa_freqs[j][amino_acids[i]]);
	  string tooltip(tooltip_char);
	  
	  file_string+="<span class=\"tooltiptext\"><b>"+tooltip+"</b><br></span></div></td>\n";

	}
      file_string+= "</tr>\n";
    }
  file_string+="<body>\n</html>\n"; 
  ofstream html_file_out;
  replace_all(filename,"txt","html");
  html_file_out.open(filename);
  html_file_out<<file_string;
  html_file_out.close();
  
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

int simulate_S5F_mutation(string sequence, int &num_mutations, map<string,S5F_mut> &S5F_model, mt19937 &gen, uniform_real_distribution<double> &dis, bool kill_stop_seqs, vector<string> &mutant_sequences, bool is_shielded, vector<bool> &shield_mutations)
{
  if (num_mutations==0){mutant_sequences.push_back(sequence); return 0;}
  ///iterate num_mutations times
  int last_mutate_position=-2;
  vector<double> last_mut_scores(sequence.length(), 0.0);
  double last_sum_mut_scores=0;
  int stop_codon_count=0;
  for(int j=1; j<=num_mutations; j++)
    {
      vector<double> mut_scores(sequence.length(), 0.0);
      double sum_mut_scores=0;
      get_mutability_scores(S5F_model, sequence, last_mutate_position, is_shielded, shield_mutations, last_mut_scores, mut_scores, last_sum_mut_scores, sum_mut_scores);     
      last_mut_scores=mut_scores;
      last_sum_mut_scores=sum_mut_scores;
      if(mut_scores.size() != shield_mutations.size())
	{
	  cerr << "FATAL ERROR: shield mutations vector is not same size as mut_scores. Exiting\n";
	  exit(1);
	} 
      // cerr << j << "\tsum of mut scores: " << sum_mut_scores << "\n"; 
      ///convert mutability scores to probability of mutating position  
      vector<double> mut_probability_ladder(sequence.length(),0); ///cumulative distribution of probabilities
      mut_probability_ladder[0]=mut_scores[0]/(double) sum_mut_scores;
      // cerr << "0" << "\t" << sequence[0] << "\t" << (mut_scores[0]/(double) sum_mut_scores) << "\t" << mut_probability_ladder[0] << "\n"; 
      for(int i=1; i<sequence.length(); i++)
	{
	  mut_probability_ladder[i]=mut_probability_ladder[i-1]+(mut_scores[i]/(double) sum_mut_scores); 
	}

      //cout << "ladder"<<endl;
      //for(int x=0;x<mut_probability_ladder.size();x++)
      //cout << mut_probability_ladder[x]<<" ";
      //cout << endl;

      //clock_t mut_ladder=clock();
      ///draw position randomly according to probability ladder
      double R=dis(gen);
      //cerr << "R: " << R << "\t"; 
      int mutate_position_i=-1;

      for(int i=0; i<mut_probability_ladder.size(); i++)  ///OPTIMIZE: combine the two loops into one
	{
	  //cout << i<<" ";
	  
	  if (R<mut_probability_ladder[i] && !shield_mutations[i])
	    {
	      mutate_position_i=i;
	      break;
	    }
	  if(mut_probability_ladder.size()-1==i)
	    {
	      while(shield_mutations[i])
		{
		  i--;
		}
	      mutate_position_i=i;
	      break;
	    }
	}
      //end=clock();
      //elapsed = double(end - begin);// / CLOCKS_PER_SEC;
      //cerr << "LADDER BUILDING: " << elapsed << "\n";
      //catch when mutate_position_i is not set 
      if (mutate_position_i == -1)
	{
	  cerr << "FATAL ERROR: Setting mutation_position_i failed-684\n";
	  exit(1);
	}

      last_mutate_position=mutate_position_i;
      
      if (sequence.substr(mutate_position_i,1) == "-")
	{
	  cerr << "mutate position should never be a gap\n";
	  //int d; cin >> d;
	}
      
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
	cerr << "should not be an X ever, paused\n"; 
	cerr << "A: " << substitution_probs['A'] << " C: " << substitution_probs['C'] << " G: " << substitution_probs['G'] << " T: " <<  substitution_probs['T'] << "\n"; 
	cerr << j << "\t" << mutate_position_i << "\t" << sequence[mutate_position_i] << "\n"; 
	//int d; cin >> d; 
      }
      string sequence_copy=sequence;
      sequence[mutate_position_i]=base_to_mutate_to;

      ///store
      if (dna_sequence_has_stop_codon_in_reading_frame(sequence,mutate_position_i,base_to_mutate_to))
	//if (dna_sequence_has_stop_codon_in_reading_frame(sequence))
	{ sequence=sequence_copy; j--; stop_codon_count++;}///discard if hits a stop codon and start from prev sequence
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

  return stop_codon_count;
}


void simulate_S5F_lineage(string sequence, int number_branches, int &num_mutations, map<string,S5F_mut> &S5F_model, mt19937 &gen, uniform_real_distribution<double> &dis, bool kill_stop_seqs, vector<string> &clonal_sequences, bool is_shielded, vector<bool> &shield_mutations)
{
  if (num_mutations==0){clonal_sequences.push_back(sequence); return;}
  vector<string> clonal_lineage;
  int mutation_count;
  int new_num_mutations;
  simulate_S5F_mutation(sequence, num_mutations, S5F_model, gen, dis,kill_stop_seqs, clonal_lineage, is_shielded, shield_mutations);
  clonal_sequences.push_back(clonal_lineage.back());
  for(int k=1;k<=number_branches-1;k++)
    {
      vector<int> mutation_count_vect;
      vector<int> rand_index;

      for(int i=0;i<clonal_lineage.size();i++)
	{
	  number_of_mutations_two_seqs(sequence,clonal_lineage[i],mutation_count);
	  mutation_count_vect.push_back(mutation_count);
	}
      uniform_int_distribution<int> distribution(2,mutation_count-1);
      int vect=distribution(gen);
      for(int i=0;i<mutation_count_vect.size();i++)
	{
	  if(mutation_count_vect[i]==vect)
	    {
	      rand_index.push_back(i);
	    }
	}
      random_shuffle(rand_index.begin(),rand_index.end());
      new_num_mutations=num_mutations-vect;
      string new_seq=clonal_lineage[rand_index[0]];
      number_of_mutations_two_seqs(sequence,new_seq,mutation_count);
      simulate_S5F_mutation(new_seq, new_num_mutations, S5F_model, gen, dis,kill_stop_seqs, clonal_lineage, is_shielded, shield_mutations);
      clonal_sequences.push_back(clonal_lineage.back());
    }
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
	  if (v2[i][j].aa=="X")
	    {
	      v2[i][j].aa="-";
	    }

	  HTML::Td td1("","","","3",v2[i][j].aa);

	  HTML::Td td2("","","","3",convert_to_string(v2[i][j].aa_num));

	  char c_str[10];
	  sprintf(c_str,"%.2f",v2[i][j].S5F_mut_score);
	  string mut_score_str(c_str);
	  if (v2[i][j].S5F_mut_score == -1)
	    {
	      mut_score_str="N/A";
	    }
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
	  if ((v2[i][j].simulated_aa_positional_frequency<low_prob_cutoff)&&(v2[i][j].simulated_aa_positional_frequency!=-99.99))
	    {
	      tdP1.hclass="highlightaalowprob";
	    }
	  //if (v2[i][j].simulated_aa_positional_frequency<low_prob_cutoff){tdP1.hclass="highlightaalowprob";}
	  HTML::Td tdF("absent","","","","");
	  if (i>0) ///highlight aa/dna mutation
	    {
	       if (v2[i][j].aa!=v2[0][j].aa)
		 {
		   td1.hclass="highlightaamutation";
		 }
	       if (v2[i][j].base!=v2[0][j].base)//dna mutation
		 {
		   td3.style="color:#4d0000";
		   if (i==v2.size()-1)
		     {
		       if ((v2[0][j].S5F_mut_score<.3)&&(v2[0][j].S5F_mut_score!=-1))
			 {
			   tdF.value="<img style='vertical-align:bottom' src='unusual.png' alt='^'/>";
			 }
		     }
		 }
	    }
	  if (j%3==0)
	    {
	      row1.cols.push_back(td1);
	    } //row 1: aa 
	  if (j%3==0)
	    {
	      row2.cols.push_back(td2);
	    } //row 2: aa number
	  row3.cols.push_back(td3); //row 3: base and mut score
	  if (i==v2.size()-1)
	    {
	      if (j%3==0)
		{
		  penultimate_row.cols.push_back(tdP1);
		}//pen row: pos aa freq
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
  if (sequence.length() != markup_string.length()){
    cerr << "WARNING: sequence and markup string lengths are not identical\n";
    cerr << sequence << "\n" << markup_string << "\n";
  }

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

void cleanup_SMUA_sequences(string sequence_name, string markup_header, string UCA_sequence, string sequence, string markup, string &new_UCA_sequence, string &new_sequence, string &new_markup, string species, string chain_type, int &number_of_replacements, bool &error_status)
{

  if ( string("human").compare(species)){ cerr << "CLEANUP ONLY WORKS FOR: Human sepeices. Sorry. Exiting.\n"; exit(1);}

  //if ( species =="human" ){ cerr << "CLEANUP ONLY WORKS FOR: Human sepeices. Sorry. Exiting.\n"; exit(1);}
  //if ( ((species =="human") && (chain_type != "heavy")) || ((species == "rhesus")&&(chain_type == "lambda")) ){ cerr << "CLEANUP ONLY WORKS FOR: Human/Rhesus Heavy and Rhesus Kappa for now. Sorry. Exiting.\n"; exit(1);}

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
      if ((new_sequence[i] == '-') || (new_UCA_sequence[i]=='-'))
	{
	  indels_present=true;
	}
    }

  if (sequence_has_ambiguities(new_sequence))
    {
      cerr << "ERROR: " << sequence_name << " has ambiguities, cannot proceed\n";
      error_status=true;
      return;
    }
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

void printTileStack(string filename, map<string,vector<Seq>> &seq_map, vector<string> sequence_names, int line_wrap_length, double low_prob_cutoff, vector<double> &color_ladder)
{
  vector<vector<Seq> > all_sequences;

  for(int i=0; i<sequence_names.size(); i++)
    {
      all_sequences.push_back(seq_map[sequence_names[i]]);
    }

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


file_string+="<p><table class=\"results\" align=center><tr><td class=\"noborder\"><font align=center size=\"4\"><b>Mutation Probability:&nbsp;</b></font></td>";
  for(int i=6;i>-1;i--)
    {
      if(i==6)
	file_string+="<td class=\"color_cat7 mut\"><div class=\"mm\">&nbsp;&nbsp;&nbsp;&nbsp;</div></td><td class=\"noborder\">&ge;"+formatDouble(color_ladder[i-1]*100)+"\%&nbsp;&nbsp;</td>";
      else if(i==0)
	file_string+="<td class=\"color_cat"+to_string(i+1)+ " mut\"><div class=\"mm\">&nbsp;&nbsp;&nbsp;&nbsp;</div></td><td class=\"noborder\">&lt;"+formatDouble(color_ladder[i]*100) +"\%</td>";
      else
	file_string+="<td class=\"color_cat"+to_string(i+1)+" mut\"><div class=\"mm\">&nbsp;&nbsp;&nbsp;&nbsp;</div></td><td class=\"noborder\">"+formatDouble(color_ladder[i]*100)+"-"+formatDouble(color_ladder[i-1]*100)+"\%&nbsp;&nbsp;</td> ";
    }
  file_string+="</tr></table></p>\n";

  
  file_string+="<body>\n</html>\n"; 

  ofstream file_out;
  file_out.open(filename.c_str());
  file_out << file_string;
  file_out.close();
}

void convert_2D_seq_vector_to_HTML_table_for_tiles_view(vector<vector<Seq> >&v2, vector<string> &names, HTML::Table &html_table, double &low_prob_cutoff, vector<double> &color_ladder, int &counter, bool RANK)
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
	  string str1;
	  if(v2[i][j].aa=="X")
	    str1="<div class=\"sm\">&mdash;</div>"; //
	  else
	    str1="<div class=\"sm\">"+v2[i][j].aa+"</div>"; //
	  HTML::Td td1("","","","",str1);
	  HTML::Td td_cdr("noborder","","","","");
	  HTML::Td spacer("spacer","","","","");
	  //assign category for coloring
	  for(int k=0; k<color_ladder.size(); k++)
	    {
	      if(RANK)
		{
		  if (v2[i][j].rank <= color_ladder[k])
		    {
		      ostringstream ss;
		      if(v2[i][j].aa=="X")
			ss << "color_cat" << 8;
		      else
			ss << "color_cat" << k+1;
		      td1.hclass=ss.str();
		      break;
		    }
		}
	      else
		{
		  if (v2[i][j].simulated_aa_positional_frequency <= color_ladder[k])
		    {
		      ostringstream ss;
		  if(v2[i][j].aa=="X")
		    ss << "color_cat" << 8;
		  else
		    ss << "color_cat" << k+1;
		  td1.hclass=ss.str();
		  break;
		    }
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
	  if ( (((j==0)&&cdr_row.cols[j].hclass=="CDR")) || ((j>0)&&(cdr_row.cols[j-1].hclass!="CDR")&&(cdr_row.cols[j].hclass=="CDR")) ) //start of CDR
	    {
	      int cdr_start=j;
	      int k=j;
	      while ((k<cdr_row.cols.size())&&(cdr_row.cols[k].hclass=="CDR"))
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
  string alignment_input_uca_sequence, alignment_old_obs_sequence;
  pairwise_align_sequences_semiglobal_w_affine_gap(scoring_matrix, input_uca_sequence, old_obs_sequence, gap_open, gap_extend, alignment_input_uca_sequence, alignment_old_obs_sequence, score);

  int oms_counter=0;
  for(int i=0; i<alignment_input_uca_sequence.length(); i++)
    {
      if (alignment_input_uca_sequence[i]=='-') //del
	{
	  new_markup_string+=old_markup_string[oms_counter];
	  oms_counter++;
	}
      else if (alignment_old_obs_sequence[i]=='-') //ins 
	{
	  new_markup_string+="-";
	}
      else //no indel
	{
	  new_markup_string+=old_markup_string[oms_counter];
	  oms_counter++;
	}
    }

  new_obs_sequence=alignment_old_obs_sequence;
  new_uca_sequence=alignment_input_uca_sequence;

  cout << ">new_obs\n" << new_obs_sequence << "\n>new_uca\n" << new_uca_sequence << "\n>new_markup\n" << new_markup_string << "\n";
  return;
}

void read_V(const std::string & dirname, map<string,map<int, map<char,double> >>  & v_input)
{
    vector<char> amino_acids={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
  DIR* dpdf= opendir(dirname.c_str());
  struct dirent * epdf;

  //map<int, map<char,double> > mature_mutant_positional_aa_freqs1;
  
  clock_t begin=clock();
  cerr <<"Loading, Please wait ......";
  vector<thread> thread_list;
  while ((epdf=readdir(dpdf))!=NULL){
    string filename= epdf->d_name;
    if(filename.length()<3)
      continue;
        map<int, map<char,double> > mature_mutant_positional_aa_freqs;
    //       cout << name << "\t" << filename << "\n";
    ifstream file;
    string Vgen = dirname+'/'+ filename;
    file.open(Vgen, std::ios::in );
    if(!file)
      {
	return;
      }
    string dummyline;
    int pos1; int pos=-1; double A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y;
    while((!getline(file,dummyline).eof())){
      if(pos==-1){pos++;continue;}
      istringstream ss(dummyline);
      getline(ss,dummyline,',')>>A;
      getline(ss,dummyline,',')>>C;
      getline(ss,dummyline,',')>>D;
      getline(ss,dummyline,',')>>E;
      getline(ss,dummyline,',')>>F;
      getline(ss,dummyline,',')>>G;
      getline(ss,dummyline,',')>>H;
      getline(ss,dummyline,',')>>I;
      getline(ss,dummyline,',')>>K;
      getline(ss,dummyline,',')>>L;
      getline(ss,dummyline,',')>>M;
      getline(ss,dummyline,',')>>N;
      getline(ss,dummyline,',')>>P;
      getline(ss,dummyline,',')>>Q;
      getline(ss,dummyline,',')>>R;
      getline(ss,dummyline,',')>>S;
      getline(ss,dummyline,',')>>T;
      getline(ss,dummyline,',')>>V;
      getline(ss,dummyline,',')>>W;
      getline(ss,dummyline,',')>>Y;
      
      mature_mutant_positional_aa_freqs[pos][amino_acids[0]]=A;
      mature_mutant_positional_aa_freqs[pos][amino_acids[1]]=C;
      mature_mutant_positional_aa_freqs[pos][amino_acids[2]]=D;
      mature_mutant_positional_aa_freqs[pos][amino_acids[3]]=E;
      mature_mutant_positional_aa_freqs[pos][amino_acids[4]]=F;
      mature_mutant_positional_aa_freqs[pos][amino_acids[5]]=G;
      mature_mutant_positional_aa_freqs[pos][amino_acids[6]]=H;
      mature_mutant_positional_aa_freqs[pos][amino_acids[7]]=I;
      mature_mutant_positional_aa_freqs[pos][amino_acids[8]]=K;
      mature_mutant_positional_aa_freqs[pos][amino_acids[9]]=L;
      mature_mutant_positional_aa_freqs[pos][amino_acids[10]]=M;
      mature_mutant_positional_aa_freqs[pos][amino_acids[11]]=N;
      mature_mutant_positional_aa_freqs[pos][amino_acids[12]]=P;
      mature_mutant_positional_aa_freqs[pos][amino_acids[13]]=Q;
      mature_mutant_positional_aa_freqs[pos][amino_acids[14]]=R;
      mature_mutant_positional_aa_freqs[pos][amino_acids[15]]=S;
      mature_mutant_positional_aa_freqs[pos][amino_acids[16]]=T;
      mature_mutant_positional_aa_freqs[pos][amino_acids[17]]=V;
      mature_mutant_positional_aa_freqs[pos][amino_acids[18]]=W;
      mature_mutant_positional_aa_freqs[pos][amino_acids[19]]=Y;
      
      pos++;
    }
    file.close();
    v_input[filename]=mature_mutant_positional_aa_freqs;  
  }
  
  clock_t end=clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  cerr << "Reading in files took " << elapsed_secs << " to process\n";
  closedir(dpdf);
  return ;
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
	  // cout << names[i] << "\t" << v2[i][j].aa_num << "\t" << v2[i][j].aa <<  "\t" <<  v2[i][j].base << "\t" << v2[i][j].simulated_aa_positional_frequency << "\t" << td1.hclass << "\n";
	  
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


//map<string, map<string, string> >
string J_genes_list(string species,string IGname)
{
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
  return J_genes[species][IGname];
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

template <typename Type>
string convert_to_string(Type t)
{
  ostringstream convert;
  convert << t;
  return convert.str();
}
