#ifndef NabEntry_H
#define NabEntry_H

#include <iostream>
//#include <boost/filesystem.hpp>

using namespace std;
//namespace fs=std::filesystem;
using namespace boost;

//using namespace boost::filesystem;

class NabEntry//class for collecting and doing the sequences
{
public:
  string sequence_name;
  string sequence;
  string UCA_sequence_name;
  string UCA_sequence;
  string markup_header;
  string markup_string;

  ///amino acids vector
  vector<char> amino_acids={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
  
  string output_filename,tiles_output_filename,seqsDNA_fasta_file,seqs_fasta_file,output_freq_table;
  string UCA_seq_name,Vgene,Dgene,Jgene;
  string Jgene_sequence;
  int cdr_length=0,Vgene_length=0,Jgene_length=0;
  map<int, map<char,double> > mature_mutant_positional_aa_freqs;
  int mut_count=0,Vgene_mut_count=0,Jgene_mut_count=0,CDR3_mut_count=0;//mutations counts for all and for the pieces
  int insertion_count=0, deletion_count=0;//insertion and deletion counts
  vector<string> mature_mutant_sequences;
  vector<string> DNA_mutant_sequences;
  vector<string> markup_mask;
  bool ignore_CDR3=false,ignoreJ=false,ignoreV=false,setMutcount=false;
  string log_cout="";
  string log_cerr="";
  string outputMode="HTML";//options:HTML,none,simple,fulltext,all
  vector<bool> shield_mutations;
  vector<Seq> aa_out;

  bool rank=false;
  
  ~NabEntry() {};
  NabEntry(vector<string>  Entry, Arguments arg)//constructor with number of mutations
  {
    sequence_name=Entry[0];
    sequence=Entry[1];
    UCA_sequence_name=Entry[2];
    UCA_sequence=Entry[3];
    markup_header=Entry[4];
    markup_string=Entry[5];
    to_upper(sequence);
    to_upper(UCA_sequence);
    mut_count=arg.numbMutations;
    ignore_CDR3=arg.ignore_CDR3;
    ignoreJ=arg.ignoreJ;
    ignoreV=arg.ignoreV;
    outputMode=arg.outputMode;

	
    if(UCA_sequence.length()!=sequence.length())
      {
	if(UCA_sequence.length()==sequence.length()+1)//if UCA is just one longer
	  {
	    UCA_sequence=UCA_sequence.substr(0,sequence.length());
	    markup_string=markup_string.substr(0,sequence.length());
	  }
	else//alignment if sequences are different sizes
	  {
	  string newUCA,newSeq;
	  map<char,map<char,int> > scoring_matrix;
	  load_EDNAFULL_matrix(scoring_matrix);
	  double gap_open=-20, gap_extend=-2, score=0;
    
	  pairwise_align_sequences_semiglobal_w_affine_gap(scoring_matrix, UCA_sequence,sequence ,gap_open,gap_extend,newUCA, newSeq,score);
	  UCA_sequence=newUCA;
	  sequence=newSeq;
	  }
      }
    
    for(int j=max(UCA_sequence.length(),sequence.length())-1;j>0;j--)
      {
	if(sequence[j]=='-' && sequence.length()%3>0)
	  {
	    UCA_sequence=UCA_sequence.substr(0,j);
	    sequence=sequence.substr(0,j);
	    markup_string=markup_string.substr(0,j);
	  }
	else
	  {
	    break;
	  }
      }
    
    for(int j=0; j<UCA_sequence.length(); j++)
      {
	if (UCA_sequence[j] == '-'){insertion_count++;}
	if (sequence[j] == '-'){deletion_count++;}
      }

    markup_mask=parseMarkup();
    vector<bool> _shield_mutations(markup_string.length(), false);
    shield_mutations=_shield_mutations;
    //cout << "outdir\t"<<arg.outDirectory<<"\n";
    //if(!arg.outDirectory.empty())
    //  {
    //	boost::filesystem::create_directory(arg.outDirectory);
    //	cout << "outdir not empty";
    //  }

    if(arg.numbMutations<0)
      {
	output_filename=sequence_name+".ARMADiLLO.html";
	tiles_output_filename=sequence_name+".tiles.html";
	seqs_fasta_file=sequence_name+".ARMADiLLO.simulated_seqs.fasta";
	seqsDNA_fasta_file=sequence_name+".DNA.ARMADiLLO.simulated_seqs.fasta";
	output_freq_table=sequence_name+".freq_table.txt";
	//number_of_mutations_two_seqs(UCA_sequence, sequence, mut_count);
	getNumberMutations();
      }
    else
      {
	setMutcount=true;
	output_filename=sequence_name+".N"+to_string(arg.numbMutations)+".ARMADiLLO.html";
	tiles_output_filename=sequence_name+".N"+to_string(arg.numbMutations)+".tiles.html";
	seqs_fasta_file=sequence_name+".N"+to_string(arg.numbMutations)+".ARMADiLLO.simulated_seqs.fasta";
	seqsDNA_fasta_file=sequence_name+".N"+to_string(arg.numbMutations)+".DNA.ARMADiLLO.simulated_seqs.fasta";
	output_freq_table=sequence_name+".N"+to_string(arg.numbMutations)+".freq_table.txt";
      }
  }
  
  void printlog()
  {
    cerr << log_cerr;
    cout << log_cout;
  }

  map<int, map<char,double>> generateVInputTable(map<string,S5F_mut> &S5F_5mers, map<string,string> &dna_to_aa_map,mt19937 &gen, uniform_real_distribution<double> &dis,int mutation_count, string gene)
  {
    int iter=200000;
    //map<string,map<string, string>> J_genes_list=J_genes_list();
    map<int, map<char,double> > sim_aa_freqs;
    vector<string> sub_mature_mutant_sequences;
    vector<bool> sub_shield_mutations;
    string sub_sequence;
    int start=0;

    if(gene.compare("J")==0)
      {
	//sub_sequence=Jgene_sequence;
	sub_sequence=UCA_sequence.substr(cdr_length+Vgene_length);
      }
    else
      {
	sub_sequence=UCA_sequence.substr(0,Vgene_length+2);
      }
    //boost::replace_all(sub_sequence,"-","");

    for(int i=0;i<sub_sequence.length();i++)
      sub_shield_mutations.push_back(0);
    
    for(int j=1; j<=iter; j++)
      {
	vector<string> mutant_sequences;
	string _aa_sequence;
	simulate_S5F_mutation(sub_sequence, mutation_count, S5F_5mers, gen, dis,true, mutant_sequences, true, sub_shield_mutations);
	DNA_mutant_sequences.push_back(mutant_sequences[mutant_sequences.size()-1]);
	translate_dna_to_aa(mutant_sequences[mutant_sequences.size()-1], _aa_sequence, 1, dna_to_aa_map);
	//	   mature_mutant_sequences[j-1]=aa_sequence2;
	sub_mature_mutant_sequences.push_back(_aa_sequence);
      }

    for(int j=start; j<sub_mature_mutant_sequences[0].length(); j++)
      {
	for(int k=0; k<amino_acids.size(); k++)
	  {
	    sim_aa_freqs[j][amino_acids[k]]=0;
	  }
      }
    for(int j=0; j<sub_mature_mutant_sequences.size(); j++)
      {
	//cout << ">mutant" << i+1 << "\n" << mature_mutant_sequences[i] << "\n";
	for(int k=start; k<sub_mature_mutant_sequences[j].size(); k++)
	  {
	    sim_aa_freqs[k][sub_mature_mutant_sequences[j][k]]+=(1/((double)sub_mature_mutant_sequences.size()));
	  }
      }
    
    return sim_aa_freqs;
  }
  
  void getNumberMutations()
  {
    assert(UCA_sequence.length()==sequence.length());

    mut_count=0;
    Vgene_mut_count=0;
    Jgene_mut_count=0;
    CDR3_mut_count=0;

    for(int i=0; i<UCA_sequence.size(); i++)
    {
      if ((UCA_sequence[i] == '-') || (sequence[i] == '-')){continue;} ///Not counting gaps as mutations currently
      if (UCA_sequence[i] != sequence[i]) {mut_count++;}
      if (UCA_sequence[i] != sequence[i] && markup_mask[i]=="V") {Vgene_mut_count++;}
      if (UCA_sequence[i] != sequence[i] && markup_mask[i]=="J") {Jgene_mut_count++;}
      if (UCA_sequence[i] != sequence[i] && markup_mask[i]=="C") {CDR3_mut_count++;}
    }
    return;
  }
  
  vector<string> parseMarkup()
  {
    //parses the markup Header
    stringstream ss(markup_header);
    string tmp;
    vector<string> markup_pieces;
    while (getline(ss,tmp,'|'))
      {
	if(tmp.substr(0,2)=="IG" && tmp[3]=='V')
	  {
	    Vgene=tmp;
	  }
	else if(tmp.substr(0,2)=="IG" && tmp[3]=='D')
	  {
	    Dgene=tmp;
	  }
	else if(tmp.substr(0,2)=="IG" && tmp[3]=='J')
	  {
	  Jgene=tmp;
	  }
	else
	  UCA_seq_name=tmp;
      }
    
    //parses the markup string
    vector<string> _markup_mask(markup_string.length(),"-");
    string lastChar="U";
    for(int j=0;j<markup_string.length();j++)
      {
	if ((markup_string[j]=='V')||(markup_string[j]=='n')||(markup_string[j]=='D')||(markup_string[j]=='J'))
	  {
	    _markup_mask[j]="C";
	    cdr_length++;
	    lastChar="C";
	  }
	else if((markup_string[j]=='1')||(markup_string[j]=='A')||(markup_string[j]=='2')||(markup_string[j]=='B')||(markup_string[j]=='3'))
	  {
	    _markup_mask[j]="V";
	    Vgene_length++;
	    lastChar="V";
	  }
	else if(markup_string[j]=='4')
	  {
	    _markup_mask[j]="J";
	    Jgene_length++;
	    lastChar="J";
	  }
	else if(markup_string[j]=='U' and (lastChar=="U" || lastChar=="V"))
	  {
	    _markup_mask[j]="V";
	    Vgene_length++;
	    lastChar="V";
	  }
	else if(markup_string[j]=='U' and lastChar=="J")
	  {
	    _markup_mask[j]="J";
	    Jgene_length++;
	    lastChar="J";
	  }
	else if(markup_string[j]=='-')
	  {
	    _markup_mask[j]="G";
	  }
      }
    return _markup_mask;
  }

  void createShield(bool ignore_CDR3,bool ignoreV, bool ignoreJ)
  {
    if(ignore_CDR3 && ignoreV && ignoreJ)
      {
	log_cerr+="Whole sequence ignored, nothing to process\n";
	exit(-1);
      }
    for(int j=0;j<markup_string.length();j++)
      {
	if(UCA_sequence.substr(j,1)=="-")
	  shield_mutations[j]=true;
	else if(ignore_CDR3 && markup_mask[j]=="C")
	   shield_mutations[j]=true;
	else if(ignoreV && markup_mask[j]=="V")
	  shield_mutations[j]=true;
	else if(ignoreJ && markup_mask[j]=="J")
	  shield_mutations[j]=true;
      }
    
  }

  int cleanSMUA(string species, string chain_type,int replace_J_upto)
  {
    string new_sequence="", new_UCA_sequence="", new_markup_string="";
    int number_of_replacements=0;
    bool error_status=false;

    log_cerr += "Cleaning SMUA step\n"; 
    cleanup_SMUA_sequences(sequence_name, markup_header, UCA_sequence, sequence, markup_string, new_UCA_sequence, new_sequence, new_markup_string, species, chain_type, number_of_replacements, error_status);
    if (number_of_replacements>replace_J_upto)
      {
	log_cerr+= "Sequence: "+ sequence_name + " required " + to_string(number_of_replacements) + " replacements at the end of the J where " + to_string(replace_J_upto) + " is allowed (user defined). Skipping that sequence\n"; 
	log_cout+= sequence_name + "\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\n"; 
	return -1;
      }
    if (error_status)
      {
	log_cerr+= "SMUA incorrectly formatted for sequence " + sequence_name + ". Skipping that sequence\n"; 
	log_cout+= sequence_name + "\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\n"; 
	return -1;
      }
    //cerr << new_UCA_sequence << "\n" << new_sequence << "\n"; 
    UCA_sequence=new_UCA_sequence;
    sequence=new_sequence;
    markup_string=new_markup_string;

    return 0;
  }

  void replaceUCA(map<string,string> &dna_to_aa_map,string input_UCA_sequence, bool ignore_warnings)
  {
    string new_sequence="", new_UCA_sequence="", new_markup_string="";
    boost::erase_all(input_UCA_sequence,"-");
    boost::erase_all(UCA_sequence,"-");
    
    replace_UCA_sequence_in_SMUA(sequence, UCA_sequence, markup_string, input_UCA_sequence, new_sequence, new_UCA_sequence, new_markup_string, ignore_warnings);
    
    string new_UCA_aa_sequence="",input_UCA_aa_sequence="",UCA_aa_sequence="";
    translate_dna_to_aa(new_UCA_sequence, new_UCA_aa_sequence, 1, dna_to_aa_map);
    translate_dna_to_aa(input_UCA_sequence, input_UCA_aa_sequence, 1, dna_to_aa_map);
    translate_dna_to_aa(UCA_sequence, UCA_aa_sequence, 1, dna_to_aa_map);
    //cout << "inp_UCA: " << input_UCA_aa_sequence << "\n    UCA: " << UCA_aa_sequence << "\nnew_UCA: " << new_UCA_aa_sequence << "\n"; 
    //cout << ">inp_UCA\n" << input_UCA_sequence << "\n>UCA\n" << UCA_sequence << "\n>new_UCA\n" << new_UCA_sequence << "\n";
    log_cout += "inp_UCA: " + input_UCA_aa_sequence + "\n    UCA: " + UCA_aa_sequence + "\nnew_UCA: " + new_UCA_aa_sequence + "\n"; 
    log_cout += ">inp_UCA\n" + input_UCA_sequence + "\n>UCA\n" + UCA_sequence + "\n>new_UCA\n" + new_UCA_sequence + "\n";
    log_cout +=">old markup\n"+markup_string+"\n";
    log_cout +=">new markup\n"+new_markup_string+"\n";
    UCA_sequence=new_UCA_sequence;
    sequence=new_sequence;
    markup_string=new_markup_string;
    
    ///output new SMUA per record
    ofstream file_out;
    string filename=sequence_name+".new_SMUA.fasta";
    file_out.open(filename.c_str());
    file_out << ">" << sequence_name << "\n" << sequence << "\n>" << UCA_sequence_name << "\n" << UCA_sequence << "\n>" << markup_header <<  "\n" << markup_string << "\n";
    file_out.close();

    vector<bool> _shield_mutations(markup_string.length(), false);
    shield_mutations=_shield_mutations;
  }

  vector<string> printResults(map<string,S5F_mut> &S5F_5mers, map<string,string> &dna_to_aa_map, int line_wrap_length,double low_prob_cutoff,vector<double> &color_ladder)
  {
    int aa_mut_count=0;
    vector<string> SNPs;
    string UCA_aa_sequence="", aa_sequence="", UCA_aa_CDR3="", seq_aa_CDR3="";
    vector<Seq> seq_vector, UCA_seq_vector, aa_seq_vector, aa_UCA_seq_vector;
    
    int CDR3_length=0;
    
    process_SMUA_sequence_to_seq_vector(sequence, markup_string, seq_vector, dna_to_aa_map, S5F_5mers);
    process_SMUA_sequence_to_seq_vector(UCA_sequence, markup_string, UCA_seq_vector, dna_to_aa_map, S5F_5mers);
    translate_dna_to_aa(sequence, aa_sequence, 1, dna_to_aa_map);
    translate_dna_to_aa(UCA_sequence, UCA_aa_sequence, 1, dna_to_aa_map);
    number_of_mutations_two_seqs(UCA_aa_sequence, aa_sequence, aa_mut_count);
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
    
    int p02_count=0, p01_count=0, p001_count=0, p0001_count=0;
    float PPvalue=0;
    map<string, int> region_counts;
    for(int j=0; j<seq_vector.size(); j+=3)
      {
	UCA_seq_vector[j].rank=1;
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
	if (seq_vector[j].simulated_aa_positional_frequency>0)
	  {
	    PPvalue+=log(seq_vector[j].simulated_aa_positional_frequency);
	  }
      }
    generateRanking(seq_vector,mature_mutant_positional_aa_freqs);
    
    int AAseqLen= UCA_seq_vector.size()/3;
    if (mut_count==0) 
      { 
	log_cerr += "0 mutations found\n";   
	log_cout += sequence_name + "\t" + to_string(aa_mut_count) + "\t" + to_string(mut_count) + "\t0\t 0\t0\t0\t"+ to_string(insertion_count)+ "\t"+to_string(deletion_count)+"\t"+to_string((insertion_count+deletion_count)/3) + "\t" +to_string(CDR3_length)+"\t"+to_string(PPvalue/AAseqLen)+ "\n";

	//for no mutation sequences
	vector<vector<Seq> > all_sequences;
	all_sequences.push_back(UCA_seq_vector);
	all_sequences.push_back(seq_vector);
	for(int j=0; j<seq_vector.size(); j+=3)
	  {
	    if (seq_vector[j].aa != UCA_seq_vector[j].aa)
	      {
		char snp[50];
		//
		if(rank)
		  sprintf(snp,"%s%d%s:%0.6f",UCA_seq_vector[j].aa.c_str(),(j)/3+1,seq_vector[j].aa.c_str(),seq_vector[j].rank);
		else
		  sprintf(snp,"%s%d%s:%0.6f",UCA_seq_vector[j].aa.c_str(),(j)/3+1,seq_vector[j].aa.c_str(),seq_vector[j].simulated_aa_positional_frequency);

		//sprintf(snp,"%s%d%s:%0.6f:%0.5f",UCA_seq_vector[j].aa.c_str(),j,seq_vector[j].aa.c_str(),seq_vector[j].rank,seq_vector[j].simulated_aa_positional_frequency);
		string snpStr(snp);
		seq_vector[j].isMut=true;
		SNPs.push_back(snpStr.c_str());
	      }
	    else
	      {
		seq_vector[j].isMut=false;
		seq_vector[j].simulated_aa_positional_frequency=1;
	      }
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
	aa_sequence_names.push_back(sequence_name.substr(0,min(20,int(sequence_name.length()))));

	aa_out=all_aa_sequences[1];
	
      }
    else
      {
	log_cout += sequence_name + "\t" + to_string(aa_mut_count) + "\t" + to_string(mut_count) + "\t" + to_string(p02_count)+"\t"+to_string(p01_count)+ "\t" +to_string(p001_count)+ "\t"+to_string(p0001_count)+"\t"+to_string(insertion_count)+"\t"+to_string(deletion_count)+"\t"+to_string((insertion_count+deletion_count)/3) + "\t" + to_string(CDR3_length)+"\t"+to_string( PPvalue/AAseqLen)+"\n"; 
	///print detailed ARMADiLLO output as HTML
	vector<vector<Seq> > all_sequences;
	all_sequences.push_back(UCA_seq_vector);
	all_sequences.push_back(seq_vector);
	vector<string> sequence_names;
	sequence_names.push_back(sequence_name.substr(0,min(20,int(sequence_name.length())))+"|UCA");//UCA_sequence_name
	sequence_names.push_back(sequence_name.substr(0,min(20,int(sequence_name.length()))));
	if(outputMode=="HTML"||outputMode=="all")
	  print_output(output_filename, all_sequences, sequence_names, line_wrap_length, low_prob_cutoff);
	
	///print tiles as HTML
	for(int j=0; j<seq_vector.size(); j+=3)
	  {
	    if (seq_vector[j].aa != UCA_seq_vector[j].aa)
	      {
		char snp[20];
		//
		if(rank)
		  sprintf(snp,"%s%d%s:%0.6f",UCA_seq_vector[j].aa.c_str(),(j)/3+1,seq_vector[j].aa.c_str(),seq_vector[j].rank);
		else
		  sprintf(snp,"%s%d%s:%0.6f",UCA_seq_vector[j].aa.c_str(),(j)/3+1,seq_vector[j].aa.c_str(),seq_vector[j].simulated_aa_positional_frequency);

		//sprintf(snp,"%s%d%s:%0.6f:%0.5f",UCA_seq_vector[j].aa.c_str(),j,seq_vector[j].aa.c_str(),seq_vector[j].rank,seq_vector[j].simulated_aa_positional_frequency);
		string snpStr(snp);
		seq_vector[j].isMut=true;
		SNPs.push_back(snpStr.c_str());
	      }
	    else
	      {
		seq_vector[j].isMut=false;
	      }
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
	aa_sequence_names.push_back(sequence_name.substr(0,min(20,int(sequence_name.length()))));

	aa_out=all_aa_sequences[1];
	
	if(outputMode=="HTML"||outputMode=="all")
	  print_output_for_tiles_view(tiles_output_filename, all_aa_sequences, aa_sequence_names, line_wrap_length, low_prob_cutoff, color_ladder,rank);
	if (outputMode=="fulltext" || outputMode=="all")
	  {
	    detailedTextPrintOut(sequence_name+".ARMADiLLO.Detailed.text",aa_sequence,UCA_aa_sequence,seq_vector,all_sequences);
	  }
	
      }
    if(outputMode=="all" || outputMode=="fulltext")
      print_freq_table_to_file(output_freq_table,mature_mutant_positional_aa_freqs);
    if(outputMode=="HTML" || outputMode=="all")
      print_HTML_freq_table_to_file(output_freq_table,mature_mutant_positional_aa_freqs,UCA_aa_sequence,color_ladder);
    if(outputMode=="simple" || outputMode=="fulltext" || outputMode=="all")
      {
	simpleTextPrintOut(sequence_name+".ARMADiLLO.fasta",aa_sequence,UCA_aa_sequence,seq_vector);
      }
    return SNPs;
  }
  
  bool SimulateSequences(map<string,S5F_mut> &S5F_5mers, map<string,string> &dna_to_aa_map,mt19937 &gen, uniform_real_distribution<double> &dis,int max_iter, int branches, bool lineage)
  {
    if (dna_sequence_has_stop_codon_in_reading_frame(UCA_sequence))
      {
    	log_cerr+= "germline has stop codon...skipping this sequence\n"; 
    	log_cout+= sequence_name + "\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A - UCA has stop codon\n"; 
    	return false;
      }
    int stop_codon_count=0;
    int _mut_count=mut_count;
    vector<int> countStopCodons;

    if(ignore_CDR3 && !setMutcount)
      {
      _mut_count-=CDR3_mut_count;
      }
    if(ignoreJ && !setMutcount)
      _mut_count-=Jgene_mut_count;
    if(ignoreV && !setMutcount)
      _mut_count-=Vgene_mut_count;
    //log_cerr += "Simulating maturation of "+sequence_name+"\n";
    for(int j=1; j<=max_iter; j++)
      {
	//print_pct_progress(j, max_iter, 1);
	if(lineage)
	  {
	    vector<string> mutant_sequences;
	    simulate_S5F_lineage(UCA_sequence, branches,_mut_count, S5F_5mers, gen, dis,true, mutant_sequences, true, shield_mutations);
	    for(int i=0;i<mutant_sequences.size();i++)
	      {
		string _aa_sequence;
		DNA_mutant_sequences.push_back(mutant_sequences[i]);
		translate_dna_to_aa(mutant_sequences[i], _aa_sequence, 1, dna_to_aa_map);
		//	   mature_mutant_sequences[j-1]=aa_sequence2;
		mature_mutant_sequences.push_back(_aa_sequence);
	      }
	  }
	else
	  {
	    vector<string> mutant_sequences;
	    string _aa_sequence;
	    int stopCounts=0;

	    stopCounts=simulate_S5F_mutation(UCA_sequence,_mut_count, S5F_5mers, gen, dis,true, mutant_sequences, true, shield_mutations);
	    countStopCodons.push_back(stopCounts);
	    DNA_mutant_sequences.push_back(mutant_sequences[mutant_sequences.size()-1]);
	    translate_dna_to_aa(mutant_sequences[mutant_sequences.size()-1], _aa_sequence, 1, dna_to_aa_map);
	    //	   mature_mutant_sequences[j-1]=aa_sequence2;
	    mature_mutant_sequences.push_back(_aa_sequence);
	  }
      }
    //cerr << "STOP CODON #: " << stop_codon_count << "\n"; 
    //cerr << "done\n";
    ///get positional frequency of aa from simulated sequences (with same num maturation mutations)
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
	for(int k=0; k<mature_mutant_sequences[j].size(); k++)
	  {
	    mature_mutant_positional_aa_freqs[k][mature_mutant_sequences[j][k]]+=(1/((double)mature_mutant_sequences.size()));
	  }
      }
    return true;
  }

  void generateRanking(vector<Seq> &seq_vector,  map<int, map<char,double> > &positional_aa_freqs)
  {
    vector<double> values_list;
    for(int j=0; j<positional_aa_freqs.size(); j++)
      {
	for(int i=0; i<amino_acids.size(); i++)
	  {
	    values_list.push_back(positional_aa_freqs[j][amino_acids[i]]);
	  }
      }
    
    sort(values_list.begin(),values_list.end());//sorting the values
    vector<double>::iterator ip;
    ip=std::unique(values_list.begin(),values_list.begin()+values_list.size());//takes only unique values
    values_list.resize(std::distance(values_list.begin(),ip));

    for(int i=0;i<seq_vector.size();i++)
      {
	int p=0;
	double v=seq_vector[i].simulated_aa_positional_frequency;
	for(int j=0;j<values_list.size();j++)
	  {
	    if(v==values_list[j])
	      {
		p=values_list.size()-j;//counts back
		continue;
	      }
	  }
	
	seq_vector[i].rank=p;
	seq_vector[i].percentile=(double)p/(double)values_list.size();
      }
  }
  
  bool replaceTable(map<string,S5F_mut> &S5F_5mers,map<string,string> &dna_to_aa_map, map<string,map<int, map<char,double> >>  &v_input, Arguments &arg)
  {
    string UCA_aa_sequence="";
    int V_mut_count=round(((mut_count-CDR3_mut_count)*Vgene_length/(sequence.length()-cdr_length)+Vgene_mut_count*3)/4);
    int J_mut_count=round(((mut_count-CDR3_mut_count)*(Jgene_length)/(sequence.length()-cdr_length)+Jgene_mut_count*3)/4);
    translate_dna_to_aa(UCA_sequence, UCA_aa_sequence, 1, dna_to_aa_map);
    string freq_table=Vgene+"_"+std::to_string(V_mut_count)+".freq_table.txt";
    log_cerr+="Sequence name "+markup_header+'\n';
    log_cerr+= "Chosen V frequency table "+freq_table+"\n";
    if (v_input.find(freq_table)==v_input.end())
      {
	log_cout += freq_table+ " : not found in database, reverting to simulating sequence\n";
	map<int, map<char,double>> freqTable_aa=generateVInputTable(S5F_5mers,dna_to_aa_map,arg.gen,arg.dis,V_mut_count,"V");
	v_input[freq_table]=freqTable_aa;
      }
    map<int, map<char,double>> positional_aa_freqs=v_input.find(freq_table)->second;
    string shift=UCA_aa_sequence;
    for (int i=0; i < UCA_aa_sequence.length(); i++)
      {
      if (UCA_sequence[3*i]=='-' || UCA_sequence[3*i+1]=='-' || UCA_sequence[3*i+2]=='-')
	{
	  shift[i]='1';
	}
      else
	{
	  shift[i]='0';
	}
      
      if (sequence[3*i]=='-' || sequence[3*i+1]=='-' || sequence[3*i+2]=='-')
	{
	  shift[i]='0';
	}
      }
    ///get positional frequency of aa from simulated sequences (with same num maturation mutations)
    ///dealing with insertion and deletion
    int cnt=0 ;
    for(int j=0; j<Vgene_length/3; j++)
      {
	if (j>0 && shift[j-1]=='1' )
	  {
	    cnt+=1;
	  }
	for(int k=0; k<amino_acids.size(); k++)
	  {
	    mature_mutant_positional_aa_freqs[j][amino_acids[k]]=positional_aa_freqs[j-cnt][amino_acids[k]];
	  }
      }
    for(int j=0; j<cdr_length/3; j++)
      {
	for(int k=0; k<amino_acids.size(); k++)
	  {
	    if (amino_acids[k]==UCA_aa_sequence[j+Vgene_length/3])
	      {
		mature_mutant_positional_aa_freqs[j+Vgene_length/3][amino_acids[k]]=1.0;
	      }
	    else
	      {
		mature_mutant_positional_aa_freqs[j+Vgene_length/3][amino_acids[k]]=0.0;
	      }
	  }
      }

    freq_table=Jgene+"_"+std::to_string(J_mut_count)+".freq_table.txt";
    log_cerr+= "Chosen J frequency table:" + freq_table+"\n";
    if (v_input.find(freq_table)==v_input.end())
      {
	log_cout += freq_table+ " : not found in database, reverting to simulating sequence\n";
	map<int, map<char,double>> freqTable_aa=generateVInputTable(S5F_5mers,dna_to_aa_map,arg.gen,arg.dis,J_mut_count,"J");
	v_input[freq_table]=freqTable_aa;
      }
    if (Jgene_mut_count==0)
      {
	for(int j=0; j<Jgene_length/3; j++)
	  {
	    for(int k=0; k<amino_acids.size(); k++)
	      {
		if (amino_acids[k]==UCA_aa_sequence[j+Vgene_length/3+cdr_length/3])
		  {
		    mature_mutant_positional_aa_freqs[j+Vgene_length/3+cdr_length/3][amino_acids[k]]=1.0;
		  }
		else
		  {
		    mature_mutant_positional_aa_freqs[j+Vgene_length/3+cdr_length/3][amino_acids[k]]=0.0;
		  }
	      }
	  }
      }
    else
      {
	positional_aa_freqs=v_input.find(freq_table)->second;
	int cnt=0 ;
	for(int j=0; j<Jgene_length/3; j++)
	  {
	    if (j>0 && shift[j+Vgene_length/3+cdr_length/3]=='1' )
	      {
		cnt+=1;
	      }
	    for(int k=0; k<amino_acids.size(); k++)
	      {
		mature_mutant_positional_aa_freqs[j+Vgene_length/3+cdr_length/3][amino_acids[k]]=positional_aa_freqs[j-cnt][amino_acids[k]];
	      }
	  }
      }
    return true;
  }

  void countAAPairs(string aaMuts)
  {
    vector<string> aaMutsVector;
    split(aaMutsVector,aaMuts,boost::is_any_of(","));
    int pairCount=0;
    for(int j=0;j<mature_mutant_sequences.size();j++)
    {
      int subcount=0;
	for(int i=0;i<aaMutsVector.size();i++)
	  {
	    string startNT=aaMutsVector[i].substr(0,1);
	    int NTpos=stoi(aaMutsVector[i].substr(1,aaMutsVector[i].size()-2));
	    string endNT=aaMutsVector[i].substr(aaMutsVector[i].size()-1,1);
	    if(mature_mutant_sequences[0].substr(NTpos-1,1)==endNT)
	      {
		cout << endNT<<" found ";
		subcount++;
	      }
	    else
	      cout << endNT<<" not found ";
	  }
	cout << subcount <<"\t"<<pairCount<<"\t"<<aaMutsVector.size()<<"\n";
	if (subcount==aaMutsVector.size())
	  pairCount++;
	}
    log_cout+=sequence_name+" Percent of sequences with mutants:  "+ to_string(pairCount/mature_mutant_sequences.size()) +"\n";
    return;
  }
  
  void outputSimSeqs(int max_iter, int branches)//function to write out the simulated sequences
  {
    ofstream file_out;
    file_out.open(seqs_fasta_file.c_str());
    ofstream fileDNA_out;
    fileDNA_out.open(seqsDNA_fasta_file.c_str());
    //for loop to build branches
    for(int j=0; j<max_iter*branches; j++)
      {
	file_out << ">seq_" << j+1 << "\n";
	file_out << mature_mutant_sequences[j] << "\n";
	fileDNA_out << ">seq_" << j+1 << "\n";
	fileDNA_out << DNA_mutant_sequences[j] << "\n"; 
      }
    file_out.close();
    fileDNA_out.close();
   return;
  }


  void simpleTextPrintOut(string filename,string &aa_sequence,string &UCA_aa_sequence,vector<Seq> &seq_vector)
  {
    //cout << "out ln654"<<endl;
    ofstream file_out;
    file_out.open(filename.c_str());
    file_out<<">"<<sequence_name<<endl;
    file_out<<UCA_aa_sequence<<endl;
    file_out<<aa_sequence<<endl;
//cout << "lin663"<<endl;
    for(int j=0; j<aa_sequence.size(); j+=1)
      {
	//cout << "line666"<<endl;
	if (aa_sequence[j]=='X'||UCA_aa_sequence[j]=='X')
	  {
	    file_out<<"X";
	  }
	else if (aa_sequence[j]==UCA_aa_sequence[j])
	  {
	    //cout << "same"<<endl;
	    file_out<<"-";
	  }
	else if (seq_vector[3*j].simulated_aa_positional_frequency<.0001){file_out<<"6";}
	else if (seq_vector[3*j].simulated_aa_positional_frequency<.001){file_out<<"5";}
	else if (seq_vector[3*j].simulated_aa_positional_frequency<.01){file_out<<"4";}
	else if (seq_vector[3*j].simulated_aa_positional_frequency<.02){file_out<<"3";}
	else if (seq_vector[3*j].simulated_aa_positional_frequency<.10){file_out<<"2";}
	else
	  {
	    file_out<<"1";
	  }
      }
    file_out<<endl;
    //getchar();
    file_out.close();
//x	cout << "ln 688"<<endl;
return;
  }

  void detailedTextPrintOut(string filename,string &aa_sequence,string &UCA_aa_sequence,vector<Seq> &seq_vector,vector<vector<Seq>> &all_sequences)
  {
    vector<vector<vector<Seq> > > split_all_sequences;
    vector2D_to_3D(all_sequences,sequence.length(), split_all_sequences);
    
    ofstream file_out;
    file_out.open(filename.c_str());
    //file_out <<"\tUCA\t\t"<<sequence_name<<endl;
    file_out <<"sequence name:\t"<<sequence_name<<endl;
    file_out <<"markup header:\t"<<markup_header<<endl;
    file_out<<endl;

    file_out <<"pos\tUCA AA\tUCA NT\tseq AA\tseq NT\t\tP(NT)\tP(AA)"<<endl;
    for(int j=0;j<sequence.length();j++)
      {
	char c_str[10];
	sprintf(c_str,"%.3f",split_all_sequences[0][1][j].S5F_mut_score);
	string mut_score_str(c_str);
	if (split_all_sequences[0][1][j].S5F_mut_score == -1)
	  {
	    mut_score_str="N/A";
	  }
	
	string symbol="";
	if(UCA_sequence[j]!=sequence[j])
	  symbol="*";
	if(j%3==0)
	  {
	    sprintf(c_str,"%.3f",split_all_sequences[0][1][j].simulated_aa_positional_frequency);
	    file_out <<j<<"\t"<<UCA_aa_sequence[j/3]<<"\t"<<UCA_sequence[j]<<"\t"<<split_all_sequences[0][1][j].aa<<"\t"<<sequence[j]<<"\t"<<symbol<<"\t"<<mut_score_str<<"\t"<<c_str;
	    if(split_all_sequences[0][1][j].simulated_aa_positional_frequency<0.02)
	      file_out<<"<<";
	    
	    file_out<<endl;
	  }
	else
	  file_out <<j<<"\t\t"<<UCA_sequence[j]<<"\t\t"<<sequence[j]<<"\t"<<symbol<<"\t"<<mut_score_str<<"\t "<<endl;

      }
    file_out<<endl;
    file_out.close();
    //cout << seq_vector[1]<<endl;
  }

  
};

#endif
