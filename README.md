# ARMADiLLO

ARMADiLLO code base

To compile update Makefile for correct pathing and run make. This will compile ARMADiLLO and should results in an executable called ARMADiLLO. Included in the package is a tar file containing the pregenerated frequency tables of mutations.





## Running
To run ARMADillo us the syntax:

```
ARMADiLLO -SMUA <SMUA_fastqseq> -m Mutability.csv -s Substitution.csv -max_iter <NumbIterations> -chain <heavy/light> -species <species> -random _seed <randSeed:optional>
```

For example:
```
ARMADiLLO -SMUA DH270_natural_pair_clone.VH.SimpleMarkedUAs.fasta -m Mutability.csv -s Substitution.csv -max_iter 100000 -chain heavy -species human -random _seed 12345
```

The generated HTML files can be viewed in any browser; however, it will require the css files sequence_color.css and AMA.css for proper formatting. It is also advisable to move unusual.png as well for labeling the mutations.


## Arguments
*Arguments*

-help/-h

-SMUA \[SMUA file\]

-freq_dir \[dir containing freq tables\]

-w \[line wrap length (60)\]

-m \[S5F mutability file\]

-s \[S5F substitution file\]

-max_iter \[cycles of B cell maturation(100)\]

-c \[cutoff for highlighting low prob (1=1%)\]

-replace_J_upto [number of replacements in J allowed\]

-chain \[chain type (heavy=default|kappa|lambda)\]

-species \[(human=default|rhesus)\]

-\(n\)umber \[number of mutations\]

-lineage/-l \[integer number of end branches for lineage generation\]

-clean_first \[clean the SMUA prior to running\]

-output_seqs \[output sim seqs\]

-random_seed \[provide a random seed\]

-ignore_CDR3 \[flag to turn on ignoring the CDR3\]

## Authors

* Kevin Wiehe
* Joshua Martin Beem
* Yunfei ----


## Citation Paper

*Wiehe K., Bradley T., Meyerhoff R.R., Hart C. Williams W.B., Easterhoff D., Faison W.J., Kepler T.B., Saunders K.O., Alam S.M., Bonsignori M. and Haynes B.F.* (2018) Functional Relevance of Improbable Antibody Mutations for HIV Broadly Neutralizing Antibody Development. *Cell Host & Microbe*. 23(6):759-765.
[https://doi.org/10.1016/j.chom.2018.04.018]
