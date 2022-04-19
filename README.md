![Alt text](logo.png?raw=true "ARMADiLLO log")

ARMADiLLO code base

To compile update Makefile for correct pathing and run make. This will compile ARMADiLLO and should results in an executable called ARMADiLLO. Included in the package is a tar file containing the pregenerated frequency tables of mutations.


[ARMADiLLo website](https://armadillo.dhvi.duke.edu/)

For documentation about armadillo please see:

[ARMADiLLO docs](https://armadillo-docs.readthedocs.io/en/latest/index.html)

# Compiling
ARMADiLLO requires the Boost library (https://www.boost.org/).
A makefile has been included for compiling ARMADiLLO.  Common problems during compiling is missing the proper links and libraries to the Boost libraries. This can be fix in the make file on lines 13 and 14.

# Running

To run ARMADillo, us the syntax:

```
ARMADiLLO -SMUA <SMUA_fastq seq> -m Mutability.csv -s Substitution.csv -max_iter <NumbIterations> -chain <heavy/light> -species <species>
```
where the SMUA_fastaq seq file in the format output by Cloanalyst. (Cloanalyst can be downloaded from [here](http://www.bu.edu/computationalimmunology/research/software/))

For example:
```
ARMADiLLO -SMUA DH270.SimpleMarkedUAs.fasta -m Mutability.csv -s Substitution.csv -max_iter 100000 -chain heavy -species human -random _seed 12345
```
ARMADiLLO can generate a variety of outputs as controlled by the optional flags. The default is three HTML files for each sequence being analyzed \<name\>.ARMADiLLO.html, \<name\>.tiles.html, and \<name\>.freq_table.html. The generated HTML files can be viewed in any browser; however, it will require the css files sequence_color.css and AMA.css for proper formatting which will are included and will be generated with the HTML files. The "-fulloutput" will print out the HTML files and additional text and fasta files containing the results.

ARMADiLLO can also accept files generated by the [Partis] (https://github.com/psathyrella/Partis) program to generate UCA and individual sequences and UCA sequences. Fasta files of the sequences and UCA's can also be used as input to ARMADiLLO using the "-seq" and "-uca" flags.


# Arguments

USAGE:

	 ARMADiLLO [seq file options] -m [S5F mutability file] -s [S5F substitution file] <opt arguments>

required arguments:

	 -m [S5F mutability file]
	 
	 -s [S5F substitution file]
	 
Sequence Files options - either SMUA, partis or seq argument required:

	 -SMUA [SMUA file] : argument for SMUA file from Cloanalyst
	 
	 -partis [partis file] : file from partis either yaml or cvs
	 
	 -seq [seq fasta file] : fasta file containing sequences to process requires a uca file
	 
	 -uca [uca fasta file] : UCA fasta can contain either 1 seq or matching sequences to the seq file
	 
	 -markup [markup fasta file] : optional fasta for seq and uca sequence files
	 
output arguments

	 -simple_text : flag to print out simple text files
	 
	 -text : flag to print out all text files
	 
	 -HTML : (default) flag to print out HTML files
	 
	 -fulloutput : flag to print out all text and HTML files
	 
	 -annotate : flag to print out annotation of the sequences
	 
optional arguments:

	 -freq_dir [V, J Frequency file directory] : directory to pull the frequency tables for quick analysis
	 
	 -amofile [amo file] : sets the amo file to use for the quick analysis
	 
	 -resetamo   : flag to reset the amo file associated
	 
	 -w [line wrap length (60)]
	 
	 -max_iter [cycles of B cell maturation(1000)]
	 
	 -c [cutoff for highlighting low prob (1=1%)]
	 
	 -replace_J_upto [number of replacements in J allowed]
	 
	 -chain [chain type (heavy=default|kappa|lambda)]
	 
	 -species [(human=default|rhesus)]
	 
	 -(l)ineage [number of trees] : argument to generate the mutations through a lineage generation instead of linear generation
	 
	 -(n)umber [number of mutations] : arguemnt to set number of mutations to generate instead of taking from mutant sequence
	 
	 -clean_first : flag to turn on cleaning the SMUA prior to running
	 
	 -output_seqs : flag to turn on printing out simulated seq]
	 
	 -ignore_CDR3 : flag to ignore CDR3, default is false
	 
	 -ignore_V    : flag to ignore V, default is false
	 
	 -ignore_J    : flag to ignore J, default is false
	 
	 -threads [number] : sets the number of threads to use during processing - default is number of processors
	 
	 -random_seed [provide a random seed]
	 


## Authors

* Kevin Wiehe
* Joshua Martin Beem
* Yunfei Wang

## Citation Paper

*Wiehe K., Bradley T., Meyerhoff R.R., Hart C. Williams W.B., Easterhoff D., Faison W.J., Kepler T.B., Saunders K.O., Alam S.M., Bonsignori M. and Haynes B.F.* (2018) Functional Relevance of Improbable Antibody Mutations for HIV Broadly Neutralizing Antibody Development. *Cell Host & Microbe*. 23(6):759-765.
[https://doi.org/10.1016/j.chom.2018.04.018]
