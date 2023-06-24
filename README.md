# Genomic Annotations
Genomic Annotations is a project aimed at annotating genomic data using data from various databases. The primary objective is to ensure fast annotation to keep up with AI model training speed.

## Setup and Installation
To set up and run the project, follow these steps:

1. Clone the repository:

```bash
git clone <repo_url>
```

2. Install the required packages. Make sure you have Python 3.6 installed. Use the following command to install the dependencies:
```bash
pip install -r requirements.txt
```
3. Download the necessary annotation files. The following files are required for each type of annotation:
   1. 164 cell type regulation annotations:
   - [hg37](https://noble.gs.washington.edu/proj/encyclopedia/segway_encyclopedia.bed.gz)
   - need to add for hg38
   2. Classifications:
   - [hg37](https://ftp.ensembl.org/pub/grch37/current/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz)
   - [hg38](https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/Homo_sapiens.GRCh38.109.gff3.gz)
   3. Regulation regions:
   - [hg37](https://ftp.ensembl.org/pub/grch37/current/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz)
   - [hg38](https://ftp.ensembl.org/pub/current_regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz)

## Building the annotations databases
Run the following command to initiate the database building process:
#### To build cell type regulation DB:
   ```bash
   python3 build_cell_type_regulation_db.py <inputpath> <outputpath> <hg>
   ```
Parameters

- `inputpath`: The local path to cell_type_regulation.bed.gz.

- `outputpath`: The local path to save the DB.

- `hg`: {37, 38}. The desired reference genome. Use `37` for hg37 and `38` for hg38.

#### To build classifications DB:
   ```bash
   to be added
   ```
#### To build regulation regions DB:
   ```bash
   to be added
   ```

This commands will process the files and create the necessary databases for annotations in the path you provided.

## Usage
#### To test the 164 cell type regulation annotation speed, run:
   ```bash
   python3 runtime_test_cell_type_annotation.py <path> <hg> <numberofsamples=1> <outputforamt=flat> <sample>
   ``` 
Parameters

- `path`: The local path to the cell type regulation DB.

- `hg`: {37, 38}. The desired reference genome. Use `37` for hg37 and `38` for hg38.

- `numberofsamples`: The desired number of randomly generated samples on which to test the runtime. The default value is 1.

- `outputforamt`: {'flat', 'matrix'}. The desired output format of the annotation. Use `'flat'` for a one-dimensional feature vector and `'matrix'` for a matrix with features for each nucleotide. The default value is `'flat'`.

- `sample`: `-s chromosome start_pos end_pos flag`. (Optional) Specify a specific sample to annotate.




