This is changing the file

![Alt text](logos/logoD.png?raw=true "ARMADiLLO log")

ARMADiLLO code base

To compile update Makefile for correct pathing and run make. This will compile ARMADiLLO and should results in an executable called ARMADiLLO. Included in the package is a tar file containing the pregenerated frequency tables of mutations.


[ARMADiLLo website](https://armadillo.dhvi.duke.edu/)

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

-help/-h

-SMUA \[SMUA file\]

-freq_dir \[dir containing freq tables\]

-w \[line wrap length (60)\]

-m \[S5F mutability file - provided with ARMADiLLO\]

-s \[S5F substitution file - provided with ARMADiLLO\]

-max_iter \[cycles of B cell maturation(1000)\]

-c \[cutoff for highlighting low prob (1=1%)\]

-ignore_CDR3 \[flag to turn on ignoring the CDR3\]

-ignore_J \[flag to turn on ignoring the J region\]

-ignore_V \[flag to turn on ignoring the V region\]

-freq_dir \[directory containing frequency tables\]

-amo_file \[AMO file that contains the binary data ofthe frequency tables\]

-replace_J_upto [number of replacements in J allowed\]

-chain \[chain type (heavy=default|kappa|lambda)\]

-species \[(human=default|rhesus)\]

-\(n\)umber \[number of mutations\]

-lineage/-l \[integer number of end branches for lineage generation\]

-clean_first \[clean the SMUA prior to running\]

-output_seqs \[flag to output sim seqs\]

-random_seed \[provide a random seed\]

-resetamo \[flag to reset AMO file in freq_dir\]

## Authors

* Kevin Wiehe
* Joshua Martin Beem
* Yunfei Wang

## Citation Paper

*Wiehe K., Bradley T., Meyerhoff R.R., Hart C. Williams W.B., Easterhoff D., Faison W.J., Kepler T.B., Saunders K.O., Alam S.M., Bonsignori M. and Haynes B.F.* (2018) Functional Relevance of Improbable Antibody Mutations for HIV Broadly Neutralizing Antibody Development. *Cell Host & Microbe*. 23(6):759-765.
[https://doi.org/10.1016/j.chom.2018.04.018]
