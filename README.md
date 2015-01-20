# PCHCluster
This software is used to automatically inspect the hierarchical population structure from genotype data

Usage: PCHCluster accepts 2 or 3 arguments:
a.  the first arguments is the path of the folder,including 3 files:
1.".geno"(genotype file) 2.".ind" (individual file) 3.".snp" (SNP file). Any other file is not permit to be in this folder

b.  the second argument is the path to executable binary file, which usually is "bin/" including in this software. It contains
several program of "EIGENSTRAT" (Price et al. 2006), and the R script PCHCluster will call smartpca and twstats in this folder

c.  the thrid argument is an option sepecifying the output file name. If you don't give the thrid argument, then the output
file will in the folder sepcified with "_res.txt" in the folder named with the first argument

To run the software: the command line should like the following:

Rscript PCHCluster arg1 arg2 arg3 >R_pchcluster.log 2>&1 

In addition, we provide a example testing set named "example_set", then if you download the software in your home directory,
e.g. /home/fuopen/pchcluster/
then in Linux platform, the command line is:
>cd /home/fuopen/pchcluster/

>Rscript PCHCluster example_set bin example_test >R_pchcluster.log 2>&1

Further, if you have any questions or suggestions, or if you find some bugs, please contact to husile@picb.ac.cn.
