# GenDist Pipeline: Co-occurrence & Genetic Distance Analysis

GenDist is a pipeline designed to measure the co-occurrence and genetic distance of gene pairs in plasmids, genomes, or contigs.

---

## 1. Install Conda Environment

To get started, clone the repository and create the conda environment:

```bash
git clone https://github.com/braddmg/GenDist
cd GenDist
conda env create -f GenDist.yml
conda activate GenDist
```
## Co-occurrence Matrix Generation
To perform co-occurrence analysis, you'll need a presence-absence matrix where:
- Contigs are columns.
- Genes are rows.
Check the example file matrix.csv.

## Running the Co-occurrence Analysis
Start by reviewing the available options with:
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
## Example Run
To run the co-occurrence analysis using the matrix.csv file and a Jaccard distance threshold of 0.5:
``` bash
# -p option allow to scale node size by gene prevalence
Rscript Co-occurrence.R -f matrix.csv -t 0.5 -p -o network.pdf
```
The plot will evaluate gene co-occurrence. You’ll get a PDF file (network.pdf) with the network and a CSV file with Jaccard distances for further analysis. 
![Co-occurrence plot](https://github.com/braddmg/images/blob/main/network-1.png)

## Measuring Genetic Distance Between Genes 
To measure the genetic distance between two genes, use the following Python script.
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
                        The second gene to compare.Give me
  --output OUTPUT       Base name for output files (CSV and PDF).
```
## Example Run
Run the genetic distance analysis on resistance.csv, which contains resistance gene data. You’ll need the columns: SEQUENCE, GENE, START, and END. 
```
python GenDist.py -f resistance.csv -g1 sul3 -g2 qacL --output sul3-qacL
```
![Co-occurrence plot](https://github.com/braddmg/images/blob/main/sul3-qacL_histogram.pdf)
This command will output a histogram of genetic distances between the sul3 and qacL genes.
The histogram suggests a genetic distance of approximately 1.2 kb, explaining the co-occurrence of these genes across plasmids.
