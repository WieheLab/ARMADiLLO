![Alt text](logos/logoD.png?raw=true "ARMADiLLO log")

ARMADiLLO code base

To compile update Makefile for correct pathing and run make. This will compile ARMADiLLO and should results in an executable called ARMADiLLO. Included in the package is a tar file containing the pregenerated frequency tables of mutations.

A containerized version using Apptainer is available upon request from Wiehelab@duke.edu.

[ARMADiLLO website](https://armadillo.dhvi.duke.edu/)

[ARMADiLLO read-the-docs](https://armadillo-docs.readthedocs.io/en/latest/index.html)

[![DOI](https://zenodo.org/badge/455586018.svg)](https://zenodo.org/badge/latestdoi/455586018)

# Compiling
ARMADiLLo requires the Boost library (https://www.boost.org/).
A makefile has been included for compiling ARMADiLLO. If there are problems with compiling, try fixing the links and libs variables on line 13 & 14 to properly point to the required libraries.

# Running

To run ARMADillo, us the syntax:

```
ARMADiLLO -SMUA <SMUA_fastq seq> -m Mutability.csv -s Substitution.csv -max_iter <NumbIterations> -chain <heavy/light> -species <species>
```
where the SMUA_fastaq seq file in the format output by Cloanalyst. (Cloanalyst can be downloaded from http://www.bu.edu/computationalimmunology/research/software/)

For example:
```
ARMADiLLO -SMUA DH270_natural_pair_clone.VH.SimpleMarkedUAs.fasta -m Mutability.csv -s Substitution.csv -max_iter 100000 -chain heavy -species human -random _seed 12345
```
ARMADiLLO generates 2 HTML files \<name\>.ARMADiLLO.html, \<name\>.tiles.html, and \<name\>.freq_table.txt. If the -output_seqs flag is used, two additional fasta files containing the generated nucleotide and amino acid sequences. The generated HTML files can be viewed in any browser; however, it will require the css files sequence_color.css and AMA.css for proper formatting.

# Arguments
usage: ARMADiLLO -SMUA \<SMUA file\> -m Mutability.csv -s Substitution.csv <optional arguments>

**required arguments:**

-m \[S5F mutability file\]

-s \[S5F substitution file\]
	 
**Sequence Files options - either SMUA, partis or seq argument required:**

-SMUA \[SMUA file\] : argument for SMUA file from Cloanalyst
	 
-partis \[partis file\] : file from partis either yaml or cvs
	 
-seq \[seq fasta file\] : fasta file containing sequences to process requires a uca file
	 
-uca \[uca fasta file\] : UCA fasta can contain either 1 seq or matching sequences to the seq file
	 
-markup \[markup fasta file\] : optional fasta for seq and uca sequence files
	 
**output arguments\:**

-simple_text : flag to print out simple text files
	 
-text : flag to print out all text files
	 
-HTML : \(default\) flag to print out HTML files
	 
-fulloutput : flag to print out all text and HTML files
	 
-annotate : flag to print out annotation of the sequences
	 
*Frequency Table Lookup options\:*

-freq_dir \[V, J Frequency file directory\] : directory to pull the frequency tables for quick analysis
	 
-amofile \[amo file\] : sets the amo file to use for the quick analysis
	 
-resetamo   : flag to reset the amo file associated
	 
-generateTable : option to generate full frequency table using given UCA
	 
**optional arguments\:**

-\(h\)elp prints help menu
	 
-w \[line wrap length \(60)\]
	 
-max_iter \[cycles of B cell maturation\(1000)\]
	 
-c \[cutoff for highlighting low prob \(1=1%)\]
	 
-replace_J_upto \[number of replacements in J allowed\]
	 
-chain \[chain type \(heavy=default|kappa|lambda)\]
	 
-species \[\(human=default|rhesus)\]
	 
-\(l)ineage \[number of trees\] : argument to generate the mutations through a lineage generation instead of linear generation
	 
-\(p)ercent \[percent mutation\] : argument to set percent of mutations to generate instead of taking from mutant sequence
	 
-\(n)umber \[number of mutations\] : argument to set number of mutations to generate instead of taking from mutant sequence
	 
-clean_first : flag to turn on cleaning the SMUA prior to running
	 
-output_seqs : flag to turn on printing out simulated seq\]
	 
-only_V      : flag to only do V gene, default is false
	 
-ignore_CDR3 : flag to ignore CDR3, default is false
	 
-ignore_V    : flag to ignore V, default is false
	 
-ignore_J    : flag to ignore J, default is false
	 
-igoreStopCodon : set flag to ignore stop codons during generation
	 
-threads \[number\] : sets the number of threads to use during processing - default is number of processors
	 
-random_seed \[provide a random seed\]

## Copyright and Licensing

The copyrights of this software are owned by Duke University. As such, two licenses for this software are offered:
1. An open-source license under the CC BY-NC-ND 4.0 license for non-commercial academic use.
2. A custom license with Duke University, for commercial use or uses without the CC BY-NC-ND 4.0 license restrictions.

As a recipient of this software, you may choose which license to receive the code under. Outside contributions to the Duke-owned code base cannot be accepted unless the contributor transfers the copyright to those changes over to Duke University. To enter a custom license agreement without the the CC BY-NC-ND 4.0 license restrictions, please contact the Digital Innovations department at the Duke Office for Translation & Commercialization (OTC) at otcquestions@duke.edu with reference to “OTC File No. 5124” in your email.

Please note that this software is distributed AS IS, WITHOUT ANY WARRANTY; and without the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

**To request a license for non-commercial academic use, please write Wiehelab@duke.edu.**

## Authors

* Kevin Wiehe
* Joshua Martin Beem
* Yunfei Wang

## Citation Papers

*Wiehe K., Bradley T., Meyerhoff R.R., Hart C. Williams W.B., Easterhoff D., Faison W.J., Kepler T.B., Saunders K.O., Alam S.M., Bonsignori M. and Haynes B.F.* (2018) Functional Relevance of Improbable Antibody Mutations for HIV Broadly Neutralizing Antibody Development. *Cell Host & Microbe*. 23(6):759-765.
[https://doi.org/10.1016/j.chom.2018.04.018]

*Martin Beem JS, Venkatayogi S, Haynes BF, Wiehe K.* (2023) ARMADiLLO: a web server for analyzing antibody mutation probabilities. *Nucleic Acids Res*. 2023;gkad398. doi:10.1093/nar/gkad398. [https://academic.oup.com/nar/article/51/W1/W51/7187705]

*Wiehe, Kevin et al.* (2024) Mutation-guided vaccine design: A process for developing boosting immunogens for HIV broadly neutralizing antibody induction. *Cell Host & Microbe.* 32(5):693-709. [https://www.cell.com/cell-host-microbe/fulltext/S1931-3128(24)00126-4]