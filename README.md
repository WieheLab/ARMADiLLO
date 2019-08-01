# ARMADiLLO

ARMADiLLO code base

To compile update Makefile and run make




## Running
To run ARMADillo us the syntax:

```
ARMADiLLO -SMUA <SMUA_fastqseq> -m Mutability.csv -s Substitution.csv -max_iter <NumbIterations> -chain <heavy/light> -species <species> -random _seed <randSeed:optional>
```

For example:
```
ARMADiLLO -SMUA DH270_natural_pair_clone.VH.SimpleMarkedUAs.fasta -m Mutability.csv -s Substitution.csv -max_iter 100000 -chain heavy -species human -random _seed 12345
```

lineage

## Arguments



## Authors

* Kevin Wiehe
* Joshua Martin Beem


## Citation Paper

Wiehe K., Bradley T., Meyerhoff R.R., Hart C. Williams W.B., Easterhoff D., Faison W.J., Kepler T.B., Saunders K.O., Alam S.M., Bonsignori M. and Haynes B.F. (2018) Functional Relevance of Improbable Antibody Mutations for HIV Broadly Neutralizing Antibody Development. Cell Host & Microbe. 23(6):759-765.
https://doi.org/10.1016/j.chom.2018.04.018