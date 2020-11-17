#include <math.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <boost/algorithm/string.hpp>
#include "utilities.hpp"

using namespace std;

std::vector < std::vector <string> > read_delimited_file(string filename, string delimiter)
{
  std::vector < std::vector <string> > file_contents;

    ifstream file(filename.c_str(), ios::in );

    string line_str;
    while (!getline(file, line_str).eof())
    {
      chomp(line_str);
      std::vector <string> strs;
        tokenize(line_str, strs, delimiter);
        file_contents.push_back(strs);
    }
    file.close();
    if (file_contents.size()==0){cerr << "Problem reading file: " << filename << " ...exiting\n"; exit(1);}
    return file_contents;
}

void translate_dna_to_aa(string &dna, string &aa, int reading_frame, map<string,string> &dna_to_aa_tranx_map)
{
  for(int i=reading_frame-1; i<dna.size(); i+=3)
    {
      string codon=dna.substr(i,3);
      boost::to_upper(codon);
      if (dna_to_aa_tranx_map.find(codon)!=dna_to_aa_tranx_map.end())
	{aa+=dna_to_aa_tranx_map[codon];}
      else
	{aa+="X";}

    }
  return;
}

void tokenize(const string& str, std::vector<string>& tokens, const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}


void chomp(string &str)
{
  std::string whitespaces ("\n\r");
  
  std::size_t found = str.find_last_not_of(whitespaces);
  if (found!=std::string::npos)
    str.erase(found+1);
  else
    str.clear();            // str is all whitespace
}


void get_aa_tranx_map(map<string,string> &dna_to_aa_tranx_map)
{
  ///ALA
  dna_to_aa_tranx_map["GCT"]="A";  
  dna_to_aa_tranx_map["GCC"]="A";  
  dna_to_aa_tranx_map["GCA"]="A";  
  dna_to_aa_tranx_map["GCG"]="A";  
  ///ARG
  dna_to_aa_tranx_map["CGT"]="R";  
  dna_to_aa_tranx_map["CGC"]="R";  
  dna_to_aa_tranx_map["CGA"]="R";  
  dna_to_aa_tranx_map["CGG"]="R"; 
  dna_to_aa_tranx_map["AGA"]="R"; 
  dna_to_aa_tranx_map["AGG"]="R"; 
  ///ASN
  dna_to_aa_tranx_map["AAT"]="N"; 
  dna_to_aa_tranx_map["AAC"]="N"; 
  ///ASP
  dna_to_aa_tranx_map["GAT"]="D"; 
  dna_to_aa_tranx_map["GAC"]="D"; 
  ///CYS
  dna_to_aa_tranx_map["TGT"]="C"; 
  dna_to_aa_tranx_map["TGC"]="C"; 
  ///GLN
  dna_to_aa_tranx_map["CAA"]="Q"; 
  dna_to_aa_tranx_map["CAG"]="Q"; 
  ///GLU
  dna_to_aa_tranx_map["GAA"]="E"; 
  dna_to_aa_tranx_map["GAG"]="E"; 
  ///GLY
  dna_to_aa_tranx_map["GGT"]="G"; 
  dna_to_aa_tranx_map["GGC"]="G"; 
  dna_to_aa_tranx_map["GGA"]="G"; 
  dna_to_aa_tranx_map["GGG"]="G"; 
  ///HIS
  dna_to_aa_tranx_map["CAT"]="H"; 
  dna_to_aa_tranx_map["CAC"]="H"; 
  ///ILE
  dna_to_aa_tranx_map["ATT"]="I"; 
  dna_to_aa_tranx_map["ATC"]="I"; 
  dna_to_aa_tranx_map["ATA"]="I"; 
  ///LEU
  dna_to_aa_tranx_map["TTA"]="L"; 
  dna_to_aa_tranx_map["TTG"]="L"; 
  dna_to_aa_tranx_map["CTT"]="L"; 
  dna_to_aa_tranx_map["CTC"]="L"; 
  dna_to_aa_tranx_map["CTA"]="L"; 
  dna_to_aa_tranx_map["CTG"]="L"; 
  ///LYS
  dna_to_aa_tranx_map["AAA"]="K"; 
  dna_to_aa_tranx_map["AAG"]="K"; 
  ///MET
  dna_to_aa_tranx_map["ATG"]="M"; 
  ///PHE
  dna_to_aa_tranx_map["TTT"]="F"; 
  dna_to_aa_tranx_map["TTC"]="F"; 
  ///PRO
  dna_to_aa_tranx_map["CCT"]="P"; 
  dna_to_aa_tranx_map["CCC"]="P"; 
  dna_to_aa_tranx_map["CCA"]="P"; 
  dna_to_aa_tranx_map["CCG"]="P"; 
  //SER
  dna_to_aa_tranx_map["TCT"]="S"; 
  dna_to_aa_tranx_map["TCC"]="S"; 
  dna_to_aa_tranx_map["TCA"]="S"; 
  dna_to_aa_tranx_map["TCG"]="S"; 
  dna_to_aa_tranx_map["AGT"]="S"; 
  dna_to_aa_tranx_map["AGC"]="S"; 
  //THR
  dna_to_aa_tranx_map["ACT"]="T"; 
  dna_to_aa_tranx_map["ACC"]="T"; 
  dna_to_aa_tranx_map["ACA"]="T"; 
  dna_to_aa_tranx_map["ACG"]="T"; 
  ///TRP
  dna_to_aa_tranx_map["TGG"]="W"; 
  //TYR
  dna_to_aa_tranx_map["TAT"]="Y"; 
  dna_to_aa_tranx_map["TAC"]="Y"; 
  ///VAL
  dna_to_aa_tranx_map["GTT"]="V"; 
  dna_to_aa_tranx_map["GTC"]="V"; 
  dna_to_aa_tranx_map["GTA"]="V"; 
  dna_to_aa_tranx_map["GTG"]="V"; 
  ///STOP
  dna_to_aa_tranx_map["TAA"]="*"; 
  dna_to_aa_tranx_map["TGA"]="*"; 
  dna_to_aa_tranx_map["TAG"]="*"; 
 
}

bool dna_sequence_has_stop_codon_in_reading_frame(string inSequence)
{
  string sequence=findAndReplaceAll(inSequence,"-","");
  for(int i=0; i<sequence.length(); i+=3)
    {
      if (sequence[i]=='T')
	{
	  if (sequence[i+1] == 'G')
	    {
	      if (sequence[i+2] == 'A'){ return true;}
	    }
	  if (sequence[i+1] == 'A')
	    {
	      if ((sequence[i+2] == 'A') || (sequence[i+2] == 'G')){return true;}
	    }
	}
    }
  return false;
}

bool dna_sequence_has_stop_codon_in_reading_frame(string sequence,int position,char newBase)
{

  if (newBase == 'C')
    {
      return false;
    }
  
  int codonStart=position/3;
  string codon=sequence.substr(3*codonStart,3);
  if (codon=="TAG" || codon=="TAA" || codon=="TGA")
    {
      return true;
    }
  return false;
  
}



void read_fasta_file(string filename, map<string, string> &sequence_hash, vector<string> &sequence_names)
{
   ///OPEN SEQUENCE FILE AND PUT INTO MAP
  ifstream file(filename.c_str(), std::ios::in );
  if (!file.is_open()) {cerr << "ERROR reading fasta file: could not open \"" << filename << "\"...exiting...\n"; exit(1);}

  string file_str;
  string sequence_string, name;
  int name_count=0, seq_count=0;

  while (!getline(file, file_str).eof())
    {
      chomp(file_str);
      if (file_str.substr(0,1) == ">")
	{
	  name=file_str.substr(1,file_str.size()-1);
	  sequence_names.push_back(name);
	  name_count++;
	}
      else
	{
	  sequence_string=file_str.substr(0,file_str.size());
	  sequence_hash[name]+=sequence_string;
	  seq_count++;
	}
    }
}

void read_fastq_file(string filename, vector<string> &names, map<string, vector<int> > &q_score_map, map<string,string> &fastq_sequences)
{
   ///OPEN SEQUENCE FILE AND PUT INTO MAP
  ifstream file(filename.c_str(), std::ios::in );
  if (!file.is_open()) {cerr << "ERROR2: could not open '" << filename << "' ...exiting...\n"; exit(1);}

  string file_str;
  vector<string> file_contents, duplicate_names;
  string sequence_string, name;
  int name_count=0, seq_count=0;
  map <string,int> observed_names;
  while (!getline(file, file_str).eof())
    {
      chomp(file_str);
      file_contents.push_back(file_str);
    }
  for(int i=0; i<file_contents.size(); i+=4)
    {
      //check for duplicate names
      vector <string> tokens;
      tokenize(file_contents[i],tokens," ");
      name=tokens[0].substr(1,file_contents[i].size());
      //      boost::replace_all(name,":",""); ///remove colons
      names.push_back(name);
      //SKIP DUPLICATES 
      if (observed_names[name]==1){duplicate_names.push_back(name);continue;}else{observed_names[name]=1;}//cerr << "Duplicate names in sequence file: " << filename << ". Name " << name << " at line " << i+1 << " appears twice...exiting\n"; exit(1);}else{observed_names[name]=1;}
      
      //THROW ERROR ON DUPLICATES
      //if (observed_names[name]==1){duplicate_names.push_back(name); cerr << "Duplicate names in sequence file: " << filename << ". Name " << name << " at line " << i+1 << " appears twice...exiting\n"; exit(1);}else{observed_names[name]=1;}

      vector<char> fastq_vctr(file_contents[i+3].begin(),file_contents[i+3].end());
      vector<int> q_score_vctr(file_contents[i+3].length(),0);
      for(int j=0; j<fastq_vctr.size(); j++)
	{
	  int q_score=((int) fastq_vctr[j])-33;
	  q_score_vctr[j]=q_score;
	}
      q_score_map[name]=q_score_vctr;
      fastq_sequences[name]=file_contents[i+1];
    }
  
  //remove all instances of duplicate names; usually occurs due to colon removal, the punctuation not the organ
  
  for(int i=0; i<duplicate_names.size(); i++)
    {
      std::map<string,vector<int> >::iterator it;
      it=q_score_map.find(duplicate_names[i]);
      q_score_map.erase(it);         
      std::map<string,string >::iterator it2;
      it2=fastq_sequences.find(duplicate_names[i]);
      fastq_sequences.erase(it2); 
      //!!!!need to remove from names vector too...unfinished here!!!!!
    }
  
}

void load_EDNAFULL_matrix(map<char, map<char,int> > &scoring_matrix)
{
  map<char, int> A = {{'A',5},{'T',-4},{'G',-4},{'C',-4},{'S',-4},{'W',1},{'R',1},{'Y',-4},{'K',-4},{'M',1},{'B',-4},{'V',-1},{'H',-1},{'D',-1},{'N',-2}};
  map<char, int> T = {{'A',-4},{'T',5},{'G',-4},{'C',-4},{'S',-4},{'W',1},{'R',-4},{'Y',1},{'K',1},{'M',-4},{'B',-1},{'V',-4},{'H',-1},{'D',-1},{'N',-2}};
  map<char, int> G = {{'A',-4},{'T',-4},{'G',5},{'C',-4},{'S',1},{'W',-4},{'R',1},{'Y',-4},{'K',1},{'M',-4},{'B',-1},{'V',-1},{'H',-4},{'D',-1},{'N',-2}};
  map<char, int> C = {{'A',-4},{'T',-4},{'G',-4},{'C',5},{'S',1},{'W',-4},{'R',-4},{'Y',1},{'K',-4},{'M',1},{'B',-1},{'V',-1},{'H',-1},{'D',-4},{'N',-2}};
  map<char, int> S = {{'A',-4},{'T',-4},{'G',1},{'C',1},{'S',-1},{'W',-4},{'R',-2},{'Y',-2},{'K',-2},{'M',-2},{'B',-1},{'V',-1},{'H',-3},{'D',-3},{'N',-1}};
  map<char, int> W = {{'A',1},{'T',1},{'G',-4},{'C',-4},{'S',-4},{'W',-1},{'R',-2},{'Y',-2},{'K',-2},{'M',-2},{'B',-3},{'V',-3},{'H',-1},{'D',-1},{'N',-1}};
  map<char, int> R = {{'A',1},{'T',-4},{'G',1},{'C',-4},{'S',-2},{'W',-2},{'R',-1},{'Y',-4},{'K',-2},{'M',-2},{'B',-3},{'V',-1},{'H',-3},{'D',-1},{'N',-1}};
  map<char, int> Y = {{'A',-4},{'T',1},{'G',-4},{'C',1},{'S',-2},{'W',-2},{'R',-4},{'Y',-1},{'K',-2},{'M',-2},{'B',-1},{'V',-3},{'H',-1},{'D',-3},{'N',-1}};
  map<char, int> K = {{'A',-4},{'T',1},{'G',1},{'C',-4},{'S',-2},{'W',-2},{'R',-2},{'Y',-2},{'K',-1},{'M',-4},{'B',-1},{'V',-3},{'H',-3},{'D',-1},{'N',-1}};
  map<char, int> M = {{'A',1},{'T',-4},{'G',-4},{'C',1},{'S',-2},{'W',-2},{'R',-2},{'Y',-2},{'K',-4},{'M',-1},{'B',-3},{'V',-1},{'H',-1},{'D',-3},{'N',-1}};
  map<char, int> B = {{'A',-4},{'T',-1},{'G',-1},{'C',-1},{'S',-1},{'W',-3},{'R',-3},{'Y',-1},{'K',-1},{'M',-3},{'B',-1},{'V',-2},{'H',-2},{'D',-2},{'N',-1}};
  map<char, int> V = {{'A',-1},{'T',-4},{'G',-1},{'C',-1},{'S',-1},{'W',-3},{'R',-1},{'Y',-3},{'K',-3},{'M',-1},{'B',-2},{'V',-1},{'H',-2},{'D',-2},{'N',-1}};
  map<char, int> H = {{'A',-1},{'T',-1},{'G',-4},{'C',-1},{'S',-3},{'W',-1},{'R',-3},{'Y',-1},{'K',-3},{'M',-1},{'B',-2},{'V',-2},{'H',-1},{'D',-2},{'N',-1}};
  map<char, int> D = {{'A',-1},{'T',-1},{'G',-1},{'C',-4},{'S',-3},{'W',-1},{'R',-1},{'Y',-3},{'K',-1},{'M',-3},{'B',-2},{'V',-2},{'H',-2},{'D',-1},{'N',-1}};
  map<char, int> N = {{'A',-2},{'T',-2},{'G',-2},{'C',-2},{'S',-1},{'W',-1},{'R',-1},{'Y',-1},{'K',-1},{'M',-1},{'B',-1},{'V',-1},{'H',-1},{'D',-1},{'N',-1}};
 
  scoring_matrix['A']=A;
  scoring_matrix['T']=T;
  scoring_matrix['G']=G;
  scoring_matrix['C']=C;
  scoring_matrix['S']=S;
  scoring_matrix['W']=W;
  scoring_matrix['R']=R;
  scoring_matrix['Y']=Y;
  scoring_matrix['K']=K;
  scoring_matrix['M']=M;
  scoring_matrix['B']=B;
  scoring_matrix['V']=V;
  scoring_matrix['H']=H;
  scoring_matrix['D']=D;
  scoring_matrix['N']=N;
}

void pairwise_align_sequences_semiglobal_w_affine_gap(map<char, map<char, int> > &scoring_matrix, string sequence1,string sequence2,double gap_open,double gap_extend,string &alignment_sequence1, string &alignment_sequence2,double &score)
{
  //convert all to uppercase
  transform(sequence1.begin(), sequence1.end(), sequence1.begin(), ::toupper);
  transform(sequence2.begin(), sequence2.end(), sequence2.begin(), ::toupper);

  ///initialize
  alignment_sequence1="";
  alignment_sequence2="";

  ///convert sequences to vectors of chars
  vector<char> sequence1_vctr(sequence1.begin(),sequence1.end());
  vector<char> sequence2_vctr(sequence2.begin(),sequence2.end());
  sequence1_vctr.insert(sequence1_vctr.begin(),'-'); ///start with gap
  sequence2_vctr.insert(sequence2_vctr.begin(),'-'); ///start with gap
  
  ///initialize matrices
  vector<vector<double> > M(sequence1_vctr.size(), vector<double>(sequence2_vctr.size(),0));
  vector<vector<double> > A(sequence1_vctr.size(), vector<double>(sequence2_vctr.size(),0));
  vector<vector<double> > B(sequence1_vctr.size(), vector<double>(sequence2_vctr.size(),0));
  vector<vector<char> > M_ptr(sequence1_vctr.size(), vector<char>(sequence2_vctr.size(),'x'));
  vector<vector<char> > A_ptr(sequence1_vctr.size(), vector<char>(sequence2_vctr.size(),'x'));
  vector<vector<char> > B_ptr(sequence1_vctr.size(), vector<char>(sequence2_vctr.size(),'x'));

  ///set boundaries
  double boundary=-INT_MAX;
  for(int i=1; i<sequence1_vctr.size(); i++)
    {
      M[i][0]=0;
      B[i][0]=boundary;
    }
  for(int j=1; j<sequence2_vctr.size(); j++)
    {
      M[0][j]=0;
      A[0][j]=boundary; 
    }
  
  ///fill in pointer matrices
  M[0][0]=0; 
  for(int i=1; i<sequence1_vctr.size(); i++)
    {
      for(int j=1; j<sequence2_vctr.size(); j++)
	{
	  
	   int num;
	   A[i][j]=Max(M[i-1][j]+gap_open+gap_extend, A[i-1][j]+gap_extend, B[i-1][j]+gap_open+gap_extend,num);
	   if (num==1){A_ptr[i][j]='M';}else if (num==2){A_ptr[i][j]='A';}else if (num==3){A_ptr[i][j]='B';}

	   B[i][j]=Max(M[i][j-1]+gap_open+gap_extend,B[i][j-1]+gap_extend,A[i][j-1]+gap_open+gap_extend,num);
	   if (num==1){B_ptr[i][j]='M';}else if (num==2){B_ptr[i][j]='B';}else if (num==3){B_ptr[i][j]='A';}

	   M[i][j]=scoring_matrix[sequence1_vctr[i]][sequence2_vctr[j]]+Max(M[i-1][j-1],A[i-1][j-1],B[i-1][j-1], num);
	   if (num==1){M_ptr[i][j]='M';}else if (num==2){M_ptr[i][j]='A';}else if (num==3){M_ptr[i][j]='B';}

	 }
     }
  // cout << "M:\n"; 
  //  print_vector_of_vector(M);
  
  for(int i=1; i<sequence1_vctr.size(); i++)
    {
      M_ptr[i][0]='A';
      A_ptr[i][0]='A';//new
      B_ptr[i][0]='A';//new
      
    }
  for(int i=1; i<sequence2_vctr.size(); i++)
    {
      M_ptr[0][i]='B';
      B_ptr[0][i]='B';//new
      A_ptr[0][i]='B';//new
    }

  // cout << "A:\n";
  //print_vector_of_vector(A);

  // cout << "B:\n";
  // print_vector_of_vector(B);
  
   string status="continue";
   int s1=sequence1_vctr.size()-1, s2=sequence2_vctr.size()-1;
   vector<vector< char> > align(sequence1_vctr.size()+sequence2_vctr.size(),vector<char>(2));
   int counter=0;
   int num;
  
   ///semi global alignment starts with max of last row and last column, not corner cell
   double max_last_val=0;
   int maxi,maxj;
   for(int i=0; i<sequence1_vctr.size(); i++)
     {
       if (M[i][sequence2_vctr.size()-1]>max_last_val){max_last_val=M[i][sequence2_vctr.size()-1]; maxi=i;maxj=sequence2_vctr.size()-1;}
     }
   for(int j=0; j<sequence2_vctr.size(); j++)
     {
       if (M[sequence1_vctr.size()-1][j]>max_last_val){max_last_val=M[sequence1_vctr.size()-1][j]; maxi=sequence1_vctr.size()-1;maxj=j;}
     }
   ///record score as maximum of last row and last col
   score=max_last_val;
    
   /*
   cout << "M_ptr:\n"; 
   print_vector_of_vector(M_ptr);
   cout << "A_ptr:\n"; 
   print_vector_of_vector(A_ptr);
   cout << "B_ptr:\n"; 
   print_vector_of_vector(B_ptr);
   */
   // exit(1);
   ///go from last cell to maxi/maxj
   if (s1==maxi){for(int j=s2; j>maxj; j--){align[counter][0]='-';align[counter][1]=sequence2_vctr[j];counter++;}}
   if (s2==maxj){for(int i=s1; i>maxi; i--){align[counter][0]=sequence1_vctr[i];align[counter][1]='-';counter++;}}
   s1=maxi;
   s2=maxj;
   //cout << s1 << "\t" << s2 << "\n"; 
   char current_matrix='M';
   int iter=0;

   while(status=="continue")
     {
       if (++iter>2000){cerr << "infinite loop\n"; exit(1);}
       if (s1==0)
	 {
	   for(int i=s2; i>0; i--){align[counter][0]='-'; align[counter][1]=sequence2_vctr[i];counter++;}
	   status="end";
	   break;
	 }
       if (s2==0)
	 {
	   for(int i=s1; i>0; i--){align[counter][1]='-'; align[counter][0]=sequence1_vctr[i];counter++;}
	   status="end";
	   break;
	 }

       /*
       if ((s1==0) && (s2!=0))//base condition, top row, go from here to 0,0 filling with gaps
	 {
	   for(int j=s2; j>=0; j--){cout << s1 << "," << j << "\n"; align[counter][0]='-'; counter++;}
	   break;
	 }
     
       if ((s1!=0) && (s2==0))//base condition, first col, go from here to 0,0 filling with gaps
	 {
	   for(int j=s1; j>=0; j--){align[counter][1]='-'; counter++;}
	   break;
	   }
       */
       // cout << s1 << "," << s2 << " in " << current_matrix;
       //    cout << s1 << "\t" << s2 << "\t" << current_matrix; 
       if (current_matrix=='M')
	 {
	   current_matrix=M_ptr[s1][s2]; //go to this matrix next
	   ///if in M, go diagonal, if going diagonal it's a match
	   align[counter][0]=sequence1_vctr[s1]; 
	   align[counter][1]=sequence2_vctr[s2];
	   s1--; 
	   s2--; 
	 }
       else if (current_matrix=='A')
	 {
	   current_matrix=A_ptr[s1][s2];
	   ///if in A, go up, if going up, it's a del relative to first seq
	   align[counter][0]=sequence1_vctr[s1]; 
	   align[counter][1]='-';
	   s1--; 
	 }
       else if (current_matrix=='B')
	 {
	   current_matrix=B_ptr[s1][s2];
	   //if in B, go left, if going left, it's an ins relative to first seq 
	   align[counter][0]='-'; 
	   align[counter][1]=sequence2_vctr[s2];
	   s2--; 
	 }
       else {cerr << "should not reach here\n"; exit(1);}
       // cout << "--> " << s1 << "," << s2 << " in " << current_matrix << "\n"; 
       ///<< align[counter][0] << "|" << align[counter][1] << "\n"; 
       if (current_matrix=='x'){cerr << "should not reach x\n"; exit(1);}
       counter++;
       
     }

   for(int i=counter-1; i>=0; i--)
     {
       alignment_sequence1+=align[i][0];
     }
   for(int i=counter-1; i>=0; i--)
     {
       alignment_sequence2+=align[i][1];
     }

   return;
}

double Max(double a, double b, int &dir)
{
  if (a>=b){dir=1;return a;}else{dir=2;return b;}
}
double Max(double a, double b, double c, int &dir)
{
  //cout << "MAX: " << a  << " " << b << " " << c << "\n"; 
  if ((a>=b) && (a>=c)){dir=1;return a;}
  if ((b>=a) && (b>=c)){dir=2;return b;}
  else{dir=3; return c;}
}

void expected_number_of_errors_from_fastq(vector< vector<int> > &vctr_of_q_scores, vector<double> &expected_errors)
{
  for(int i=0; i<vctr_of_q_scores.size(); i++)
    {
      double error_sum=0;
      for(int j=0; j<vctr_of_q_scores[i].size(); j++)
	{
	  error_sum+=pow(10,(-1*vctr_of_q_scores[i][j]/(double) 10));
	}
      expected_errors.push_back(error_sum);
    }
  return;
}

void parse_recombination_summaries_file(string RS_file, map<string, map<string, string> > &recombination_summaries)
{
  vector<vector<string> > file_contents=read_delimited_file(RS_file,"\t");
  if (file_contents[0][0] != "UID"){cerr << "RecombinationSummaries.txt is in incorrect format. Exiting...\n"; exit(1);}

  vector<string> headers=file_contents[0];
  for(int i=1; i<file_contents.size(); i++)
    {
      string read_name=file_contents[i][0];
      if (read_name==""){continue;} //skip a blank line (does this happen?)
      for(int j=1; j<file_contents[i].size(); j++)
	{
	  recombination_summaries[read_name][headers[j]]=file_contents[i][j];
	}
    }
  
  return;
}

void sequence_identity(string alignment_seq1, string alignment_seq2, double &identity, double &gapless_identity)
{
  int match_count=0, total_count=0, gapless_match_count=0, gapless_total_count=0;
  if (alignment_seq1.length() != alignment_seq2.length()){cerr << "ERROR: identity % calculation on two sequences that are unequal. Something went horribly wrong. All hope is lost. Exiting...\n"; exit(1);}

  
  for(int j=0; j<alignment_seq1.size(); j++)
    {
      if ((alignment_seq1[j]!='-')&&(alignment_seq2[j]!='-'))
	{
	  if(alignment_seq1[j]==alignment_seq2[j]){gapless_match_count++;}//skip positions with gaps
	  gapless_total_count++;
	}
    
      if (alignment_seq1[j]==alignment_seq2[j]){match_count++;}
      total_count++;
    }
  gapless_identity=gapless_match_count/(double) gapless_total_count;
  identity=match_count/(double) total_count;
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


void print_pct_progress(int i, int size, int level)
{
  if (size<100){return;}
  double a=size/(100*pow(10,level));
  int b=1;
  if (a>1){b=(int) a;}
  if ((i%b)==0){cerr << setw(3) << fixed << setprecision(level) << (i/(double)size)*100.0 << "%\r" << flush;}
  
}

bool sequence_has_ambiguities(string sequence)
{
  for(int i=0; i<sequence.length(); i++)
    {
      if ((sequence[i]!='A')&&(sequence[i]!='C')&&(sequence[i]!='T')&&(sequence[i]!='G')&&(sequence[i]!='-')){ return true;}
    }
  return false;
}

bool fexists(const std::string& filename)
{
  std::ifstream ifile(filename.c_str());
  return (bool)ifile;
}



int countChar(string sample, char findIt)
{
    vector<int> characterLocations;
    for(int i =0; i < sample.size(); i++)
        if(sample[i] == findIt)
            characterLocations.push_back(i);

    return characterLocations.size();
}


string findAndReplaceAll(std::string & data, std::string toSearch, std::string replaceStr)
{
	// Get the first occurrence
	size_t pos = data.find(toSearch);
	string newString = data;
 
	// Repeat till end is reached
	while( pos != std::string::npos)
	{
		// Replace this occurrence of Sub String
		newString.replace(pos, toSearch.size(), replaceStr);
		// Get the next occurrence from the current position
		pos =newString.find(toSearch, pos + replaceStr.size());
	}
	return newString;
}

string formatDouble(double value)
{
  char buffer[50];
  
  sprintf(buffer,"%2.2f",value);

  string str(buffer);
  return str;
}
