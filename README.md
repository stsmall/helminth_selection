# helminth selection
* Scripts for exploring selection on helminths in multiple populations. Specifically in regards to Redman 2015 in PLoS NTDs an H. contortus.

## fasta2pylibset.py  
* Script accepts a fasta file and parses it for stats using pylibseq modules  

## helminth_figs1AB.py  
* Main script for creating figures 1A and 1B in the manuscript. You will have to change the path to msms in the code.  

## fixplotHelminth1B.py  
* A mistake in how I was coding resistant haplotypes. This script should be run on the output csv of helminth_figs1AB.py  

## ms2stats.py  
* Script accepts a ms-formatted file and parses it for use with pylibseq modules  

## permuteFst.py  
* Simple script that takes as input a ms-formatted file and outputs FST calculations as well as p-values based on permutation test. I implemented this after some issues with the current version of msstats and FST calculation from stdin. Probably become deprecated when pylibset and pymsstats is updated.  

## helminth_models.txt  
* Text version of some of the models.  

## helminthSelection.Rmd  
* Markdown for plotting output data 
