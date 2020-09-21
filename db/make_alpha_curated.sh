#! /usr/bin/env bash

echo -n > alpha_curated.fa
while read line 
do 
    sleep 10
    F=($line)
    echo ${F[1]}
    curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${F[1]}&rettype=fasta&retmode=text" >> alpha_curated.fa
    echo ${lines[1]}
done < HOR.list.txt 

