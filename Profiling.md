# Taxonomic and Functional Profiling


## Taxonomic profiling

Login to server:

```
ssh -X ubuntu@xxx.xxx.xxx.xxx
```

We will start by creating a new sub-directory in Projects:

Move into Projects with the `cd` change directory command:
```
cd Projects
```

Create a new directory for our analysis with `mkdir`:
```
mkdir InfantGut
```

Now move into the directory InfantGut (command not supplied) ...

We are going to process a subset of the classic Sharon et al. (2011) Infant Gut data set:
```
ls ~/Data/InfantGut/ReadsSub
```

How many reads are in each sample and how many samples are there?

Lets link `ln` in the Data directory:
```
ln -s ~/Data/InfantGut/ReadsSub .
```

Now we will profile the Infant Gut reads with Kraken. We will use forward reads only:

```
mkdir Kraken
for file in ReadsSub/*R1*fastq
do
    base=${file##*/}
    stub=${base%_R1.fastq}
    echo $stub
    kraken --db ~/Databases/minikraken_20141208/ --threads 8 --preload --output Kraken/${stub}_R1.kraken $file
done
```

Try to understand the anatomy of the above commands.

We match against the 'minikraken' database which corresponds to RefSeq 2014.

Would we expect the profile to differ between R1 and R2?

Can you edit the above to run the R2 reads?

Look at percentage of reads classified. Infant guts are well studied communities.


The output is just a text file:

```
head Kraken/sample1_R1.kraken
```

And we can generate a report:

```
kraken-report --db ~/Databases/minikraken_20141208/  Kraken/sample1_R1.kraken >  Kraken/sample1_R1.kraken.report
```

We can get a report of the predicted genera:
```
cat  Kraken/sample1_R1.kraken.report | awk '$4=="G"'
```

Some people prefer a different format:
```
kraken-mpa-report --db ~/Databases/minikraken_20141208/ Kraken/sample1_R1.kraken  > Kraken/sample1_R1.kraken.mpa.report
```

Compare with reverse reads what is the major difference?

Will combine the two for analysis below.

```
for file in Kraken/*R1.kraken
do
    stub=${file%_R1.kraken}
    rfile=${stub}_R2.kraken
    cat $file $rfile > ${stub}_R12.kraken
done
```


What is awk?

Now lets get reports on all combined samples:
```
for file in Kraken/*_R12.kraken
do
    stub=${file%.kraken}
    echo $stub
    kraken-report --db ~/Databases/minikraken_20141208/ $file >  ${stub}.kraken.report
done
```

Having done this we want to get one table of annotations at the genera level for community comparisons:

```
for file in Kraken/*_R12.kraken.report
do
    stub=${file%_R12.kraken.report}
    cat  $file | awk '$4=="G"' > $stub.genera
done
```

And then run associated script:
```
CollateK.pl Kraken > GeneraKraken.csv
```
There is a clear shift in genera level structure over time.

![Kraken Genera NMDS](./Figures/GeneraKNMDS.png)


We can generate this plot either locally or on the server by:

```
cp ~/Data/InfantGut/Meta.csv .
Rscript ~/bin/GeneraKNMDS.R 
```

Discussion points:
1. Non-metric multidimensional scaling
2. Multivariate permutational ANOVA

<a name="functionalprofiling"/>


## Running Metaphlan2 on the Infant Gut


```
python ~/Installation/metaphlan2/metaphlan2.py ReadsSub/sample1_R1.fastq,ReadsSub/sample1_R2.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 8 --input_type fastq
```



Note these files are not true interleaved fasta.

```
mkdir MetaphlanResults
for file in ReadsSub/*_R1.fastq
do
    base=${file##*/}
    stub=${base%_R1.fastq}
    rfile=ReadsSub/${stub}_R2.fastq
    echo $stub
    
    python ~/Installation/metaphlan2/metaphlan2.py $file,$rfile --bowtie2out MetaphlanResults/${stub}.bowtie2.bz2 --input_type fastq --nproc 8 > MetaphlanResults/${stub}_pm.txt
done
```



Then when we are done we merge these tables:
```
python ~/Installation/metaphlan2/utils/merge_metaphlan_tables.py MetaphlanResults/*_pm.txt > MetaphlanResults/merged_abundance_table.txt
python ~/Installation/metaphlan2/merge_metaphlan_tables.py MetaphlanResults/*_pm.txt > MetaphlanMerged/merged_abundance_table.txt
```

and select species:
```
SelectSpecies.pl < MetaphlanResults/merged_abundance_table.txt > MetaphlanResults/Species.tsv
```

#and generate a heatmap:
#```
#python ~/Installation/metaphlan2/utils/metaphlan_hclust_heatmap.py -c bbcry --top 25 --minv 0.1 -s log --in MetaphlanResults/merged_abundance_table.txt --out MetaphlanResults/abundance_heatmap.png
#```


## Functional gene profiling of Infant Gut reads

To perform functional gene profiling we will use Diamond to map against the KEGG database. 
First we will set an environmental variable to point to our copy of the Kegg:
```
export KEGG_DB=~/Databases/keggs_database/KeggUpdate/
```
```
cd ~/Projects/InfantGut
mkdir KeggD
for file in ReadsSub/*R1.fastq
do 
   
   stub=${file%_R1.fastq}
   stub=${stub#ReadsSub\/}
   echo $stub
   if [ ! -f KeggD/${stub}.m8 ]; then
    echo "KeggD/${stub}.m8"
    diamond blastx -d $KEGG_DB/genes/fasta/genes.dmnd -q $file -p 8 -o KeggD/${stub}.m8
   fi
done
```

Having mapped reads to the KEGG genes we can collate these into ortholog coverages:
```
for file in KeggD/*.m8
do
    stub=${file%.m8}

    echo $stub
    
    python ~/bin/CalcKOCov.py $file $KEGG_DB/ko_genes_length.csv $KEGG_DB/genes/ko/ko_genes.list > ${stub}_ko_cov.csv

done
```

Note this script uses a hard coded read length of 150 nt or 50 aa.

Discussion points:

1. What is coverage?

2. What pipelines exist for doing this, HumanN? Could we use kmers for functional profiling?

3. What is the [KEGG](http://www.genome.jp/kegg/pathway.html)

We collate these into a sample table:
```
mkdir FuncResults
CollateKO.pl KeggD > FuncResults/ko_cov.csv
```

and also KEGG modules:
```
for file in KeggD/*ko_cov.csv
do
    stub=${file%_ko_cov.csv}

    echo $stub
    python ~/bin/MapKO.py $KEGG_DB/genes/ko/ko_module.list $file > ${stub}_mod_cov.csv 
done
```

Collate those across samples:
```
CollateMod.pl KeggD > CollateMod.csv
mv CollateMod.csv FuncResults
```

What about module names? My former PDRA (Umer Ijaz) has a nice one liner for this:

```
cd FuncResults
awk -F"," 'NR>1{print $1}' CollateMod.csv | xargs -I {} curl -s http://rest.kegg.jp/find/module/{} > ModNames.txt
cd ..
```

We can view modules as multivariate data just like the genera relative frequencies. Is there a stronger or weaker relationship between time and module abundance than there was 
for the genera abundances?

![Modules](../Figures/Modules.png)