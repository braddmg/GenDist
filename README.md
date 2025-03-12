# GenDist
GenDist is a pipeline design to measure the co-occurrence and genetic distance of gene pairs in a set of múltiple plasmids, genomes or contigs.

# Installing conda environment 
Copy and activate the conda environment 

```
git clone https://github.com/braddmg/GenDist
cd GenDist
conda env create -f GenDist.yml
conda activate GenDist
```
# Generteting a co-occurrence matrix 
To create a co-occurrence analysis we need a presence absence matrix with contigs as columns and genes as rows. See the file matrix.csv as example.
Then we can run the R script with the matrix. First run the -h function to see the necesarry options.

```
Rscript Co-occurrence.R -h 
```
Options:
        -f FILE, --file=FILE
                Input CSV file (gene presence/absence matrix)

        -t VALUE, --threshold=VALUE
                Jaccard distance threshold (default: 0.5)

        -p, --prevalence
                Scale node size by gene prevalence

        -h, --help
                Show this help message and exit
Lets try to run the matrix.csv with a Jaccard distance of 0.5. This file incorporates the presence of different Antimicrobial genes across a set of IncH plasmids.
```
# -p option allow to scale node size by gene prevalence
Rscript Co-occurrence.R -f matrix.csv -t 0.5 -p -o network.pdf
```
We can evaluate co-occurence genes with the next plot. 
![Co-occurrence plot](https://github.com/braddmg/images/blob/main/network-1.png)
In your folder you also will get a file with the jaccard distances if you need it. 

## Mearusing Genetic Distance 
We can measure the genetic distance between two genes using the next python script. 
To see the option use the next command:
```
python GenDist.py -h
```
```
usage: GenDist.py [-h] -f FILE -g1 GENE1 -g2 GENE2 [--output OUTPUT]

Calculate genomic distances between two genes.

options:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  The input file with genomic data (CSV format).
  -g1 GENE1, --gene1 GENE1
                        The first gene to compare.
  -g2 GENE2, --gene2 GENE2
                        The second gene to compare.
  --output OUTPUT       Base name for output files (CSV and PDF).

- [ ] [Create](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#create-a-file) or [upload](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#upload-a-file) files
- [ ] [Add files using the command line](https://docs.gitlab.com/topics/git/add_files/#add-a-file-using-the-command-line) or push an existing Git repository with the following command:
```
Lets run the script with the resistance.csv file, which was obtained with ABRicate using the CARD database. 

cd existing_repo
git remote add origin https://gitlab.com/INISA/gendist.git
git branch -M main
git push -uf origin main
```
