# ARMADiLLO

ARMADiLLO code base

To compile update Makefile and run make

To run ARMADillo us the syntax:

ARMADiLLO -SMUA <fastqseq> -m Mutability.csv -s Substitution.csv -max_iter <NumbIterations> -chain <heavy/light> -species <species> -random _seed <randSeed:optional>

ARMADiLLO -SMUA DH270_natural_pair_clone.VH.SimpleMarkedUAs.fasta -m Mutability.csv -s Substitution.csv -max_iter 100000 -chain heavy -species human -random _seed 12345