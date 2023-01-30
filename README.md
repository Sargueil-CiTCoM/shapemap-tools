Shapemap Tools
==============
A set of tools to manage shapemap on RNA families.

## Installation 

Using git repository :

```bash
git clone https://cristal-204.pharmacie.univ-paris5.fr:8090/citcom-trad/shapemap-tools.git
# or
git clone [via koda]

cd shapemap-tools
pip install -e shapemap-tools
```

Via conda (not yet available)

```bash
conda env create shapemap
conda activate shapemap
conda install shapemap-tools -c conda-forge -c bioconda
```

## Usage 

All script are prefixed by `shpm_*`

If you installed shapemap-tools through conda, you might have to launch a virtual environment in order to use scripts .

```bash
conda activate shapemap
```

## Library preparation
## GenTag 

separate repository - Generate Tag/Barcode RNA family performing a basic RNAfold folding to minimize risk of interaction between the Barcode and the RNA.


## Treatment
### Demultiplexing - Fastq-multx

```

```

### ShapeMapper Launcher




## Normalize By conditions

Wrapper around Shapemapper normalization script in order to easy normalize a RNA for all condition

## Aggregate replicates

Aggregate multiple replicate of a given condition for all RNA of a family, with divers statistics

## Generate structure

Wrapper around IPANEMAP in order to generate structure for each RNA of a family


