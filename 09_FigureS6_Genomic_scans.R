## Figure S6
## Genomic scans
################

## Figure S6 is composed of a panel obtained with BayeScan and another by Lositan :
## - BayeScan (Foll & Gaggiotti 2008) : http://cmpg.unibe.ch/software/BayeScan/index.html
## - Lositan  (Beaumont & Nichols, 1996; Antao et al., 2008) : https://github.com/tiagoantao/lositan
## (cf. article for complete references)

## The Lositan software has a GUI that output graphs directly. Data must be in the Genepop format, use the freeware PGDSpider to handle them.
## Hereafter are listed the log files of the 3 successive runs made with Lositan :

RUN#1 (ATL only)
on 1012 loci / 827 ATL within 20 pops
100 000 sims
"Neutral" mean Fst (x)
Confidence interval 0.95
False Disc. Rate 0.1
Attempted Fst 0.001
Expected total pops 20
Mutation model: Infinite Alleles
Subsample size: 50


RUN#2 (ATL only)
on 1012 loci / 827 ATL within 20 pops
1 000 000 sims
"Neutral" mean Fst (x)
Confidence interval 0.95
False Disc. Rate 0.1
Attempted Fst 0.001
Expected total pops 20
Mutation model: Infinite Alleles
Subsample size: 50

RUN#3 (ATL and MED)
on 1012 loci / 827 ATL + 10 MED within 22 pops
1 000 000 sims
"Neutral" mean Fst (x)
Force mean Fst (x)
Confidence interval 0.995
False Disc. Rate 0.1
Attempted Fst 0.072
Expected total pops 22
Mutation model: Infinite Alleles
Subsample size: 50


## BayeScan graph was obtained with the function "plot_bayescan", located in the plot_R.r file downloaded with the BayeScan software.
## in a bash terminal
$ ./BayeScan2.1_linux64bits ./Bayescan_827DLAB_1012loci.txt  -o bayescanLongRunOD10_827ATL_1012loci -n 5000 -thin 10 -pilot 5000 -burn 50000 -pr_odds 10
bayescan_longrun <- read.table("./bayescanLongRunOD10_827ATL_1012loci")
plot_bayescan(bayescan_longrun)