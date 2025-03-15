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
The pipeline was designed to work with the ouutput of [ABRicate](https://phagescope.deepomics.org) containing the annotation of Antibiotic Resistance Genes. 
Howeer, you can use it with any file containing the columns SEQUENCE, GENE, START and END.
Check the example file resistance.txt.

## Running the Co-occurrence Analysis
Start by reviewing the available options with:
```
Rscript Network.R -h 
```
```
Options:
        -f FILE, --file=FILE
                Input TSV file containing gene presence data

        -t THRESHOLD, --threshold=THRESHOLD
                Jaccard distance threshold for network edges (default: 0.5)

        -s, --scale_node
                Scale node size by gene prevalence

        -o OUTPUT, --output=OUTPUT
                Output filename for the network plot (default: network_plot.pdf)

        -p PREVALENCE, --prevalence=PREVALENCE
                Prevalence threshold (default: 0, no filtering). Genes must be present in at least X fraction of sequences

        -h, --help
                Show this help message and exit
```
## Example Run
To run the co-occurrence analysis using theresistance.txt file and a Jaccard distance threshold of 0.5:
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
usage: GenDist.py [-h] -f FILE -g1 GENE1 -g2 GENE2 [--output OUTPUT] [--circular CIRCULAR] [--all-circular]

Calculate genomic distances between two genes.

options:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  The input file with genomic data (TSV format).
  -g1 GENE1, --gene1 GENE1
                        The first gene to compare.
  -g2 GENE2, --gene2 GENE2
                        The second gene to compare.
  --output OUTPUT       Base name for output files (CSV and PDF).
  --circular CIRCULAR   File containing sequence names, topology (0 or 1), and length for circular sequences.
```
## Example Run
Run the genetic distance analysis on resistance.csv, which contains resistance gene data. The required columns are SEQUENCE, GENE, START, and END.

If your dataset includes circular contigs, such as plasmids, you can provide a .txt file specifying the sequence name, circularity status, and total length. This information will be used to accurately calculate the genetic distance between genes while accounting for circularity. If no such file is provided, all sequences will be treated as linear by default.  
```
python GenDist.py -f resistance.csv -g1 sul3 -g2 qacL --output sul3-qacL
```
```
--output sul3-qacL
Results saved to 'sul3-qacL_results.csv'
Min: 1147.00 bp, Median: 1180.00 bp, Max: 63309.00 bp, Mode: 1180
```
![Co-occurrence plot](https://github.com/braddmg/images/blob/main/sul3-qacL_histogram-1.png)
This command will output a histogram of genetic distances between the sul3 and qacL genes.
The histogram suggests a genetic distance of approximately 1.2 kb, explaining the co-occurrence of these genes across plasmids and suggesting desequilibirum linkage.
