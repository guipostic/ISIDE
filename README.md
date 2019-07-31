## Ab initio protein folding

The ISIDE algorithm can predict the threedimensional structure of proteins through a metaheuristic exploration (tabu search<sup>1,2</sup>) of a conformational space discretized by the use of a structural alphabet (Protein Blocks<sup>3</sup>).  
&nbsp;  
The method is entirely implemented in C++ and can be run in parallel on UNIX-like operating systems.

## Install

Compile ISIDE:
```markdown
$ make
```
Then compile TMscore:
```markdown
$ cd TMscore
$ gfortran -static -O3 -ffast-math -lm -o TMscore TMscore.f
```
or
```markdown
$ g77 -static -O3 -lm -o TMscore TMscore.f
```
To get help:
```markdown
$ ./ISIDE -h
```


### Examples

```markdown
$ ./ISIDE -t 38 -c examples/folds/a_9.pdb -j job01 -n 4 -f examples/pbmtx/a_9.pbmtx -s 1000
```
```markdown
$ ./ISIDE -t 38 -j job01 -n 4 -f examples/pbmtx/a_9.pbmtx -s 1000 -i RVIAMPSVRKYAREKGVDIRLVQGTGKNGRVLKEDIDAFLAG
```
```markdown
$ ./ISIDE -t 38 -c examples/folds/a_9.pdb -j job01 -n 4 -f examples/pbmtx/a_9.pbmtx -s 1000 -y 0.20
```
```markdown
$ ./ISIDE -t 38 -c examples/folds/a_9.pdb -j job01 -n 15 -f examples/pbmtx/a_9.pbmtx -s 1000 -r 0.08
```


Note: The tabu list size should be selected based on the protein sequence length, so that:  
*(seq length) × 0.5 < [-t arg] < (seq length) × 2*

### References
1. [Glover, F. (1989). Tabu search—part I. ORSA Journal on computing, 1(3), 190-206.](https://pubsonline.informs.org/doi/abs/10.1287/ijoc.1.3.190)
2. [Glover, F. (1990). Tabu search—part II. ORSA Journal on computing, 2(1), 4-32.](https://pubsonline.informs.org/doi/abs/10.1287/ijoc.2.1.4)
3. [De Brevern, A. G., Etchebest, C., & Hazout, S. (2000). Bayesian probabilistic approach for predicting backbone structures in terms of protein blocks. Proteins: Structure, Function, and Bioinformatics, 41(3), 271-287.](https://onlinelibrary.wiley.com/doi/full/10.1002/1097-0134%2820001115%2941%3A3%3C271%3A%3AAID-PROT10%3E3.0.CO%3B2-Z)


##### © 2019 Postic G, Santuz H, Deniau R, Gelly JC.
###### Contact: [guillaume.postic@univ-paris-diderot.fr](mailto:guillaume.postic@univ-paris-diderot.fr); [jean-christophe.gelly@univ-paris-diderot.fr](mailto:jean-christophe.gelly@univ-paris-diderot.fr)
