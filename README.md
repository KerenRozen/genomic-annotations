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
4. Build the annotation database. Run the following command to initiate the database building process:
   1. For 164 cell type regulation:
   ```bash
   python3 build_cell_type_regulation_db.py <path to the 164_cell_type_regulation.bed.gz from 3i> <path to save DB file> <37/38 for hg37/hg38 respectively> 
   ```
   
   2. For classifications:
   ```bash
   to be added
   ```
   3. For regulation regions:
   ```bash
   to be added
   ```

This command will process the annotation files and create the necessary database for annotations in the path you provided.

## Usage
#### To test the 164 cell type regulation annotation speed, run:
   ```bash
   python3 _test_cell_type_annotation_speed.py <path> <hg> <numberofsamples=1> <outputforamt=flat> <sample>
   ```
   Parameters:\
      >path: The local path to the cell type regulation DB that was created in 4i.\
      hg: {37,38}. The desired reference genome. 37 or 38 for hg37 hg38 respectively.\
      numberofsamples: The desired number of random generated samples on which to test the runtime. Default is 1.\
      outputforamt: {'flat', 'matrix'}. The desired output format of the annotation. When outputforamt='flat' the output is a one dimensional features vector, and when outputforamt='matrix' the output is a matrix with the features for each nucleotide. Default is 'flat'.\
      sample: {-s chromosome start_pos end_pos flag}. Optional. A specific sample to annotate. \ 


To check annotations for your genomic data, execute the following command:

css
Copy code
python check_annotations.py --input <input_file>
Replace <input_file> with the path to your genomic data file.


## Contact
If you have any questions or suggestions regarding the Genomic Annotations project, feel free to contact [insert your contact information].