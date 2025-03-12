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
Then we can run the R script with the matrix.

```
git clone https://github.com/braddmg/GenDist
cd GenDist
conda env create -f GenDist.yml
conda activate GenDist
```
- [ ] [Create](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#create-a-file) or [upload](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#upload-a-file) files
- [ ] [Add files using the command line](https://docs.gitlab.com/topics/git/add_files/#add-a-file-using-the-command-line) or push an existing Git repository with the following command:

```
cd existing_repo
git remote add origin https://gitlab.com/INISA/gendist.git
git branch -M main
git push -uf origin main
```
