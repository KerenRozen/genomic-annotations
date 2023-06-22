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
      - [gh37](https://noble.gs.washington.edu/proj/encyclopedia/segway_encyclopedia.bed.gz)
      - need to add for hg38
   2. Classifications
      - [gh37](https://ftp.ensembl.org/pub/grch37/current/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz)
      - [gh38](https://ftp.ensembl.org/pub/release-109/gff3/homo_sapiens/Homo_sapiens.GRCh38.109.gff3.gz)
   3. Regulation regions 
      -[gh37](https://ftp.ensembl.org/pub/grch37/current/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz)
      -[gh38](https://ftp.ensembl.org/pub/current_regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz)
4. Build the annotation database. Run the following command to initiate the database building process:
   1. For 164 cell type regulation:
   ```bash
   python3 build_cell_type_regulation_db.py <path_to_segway_encyclopedia.bed> <path_to_save_DB_file> <37/38_for_hg37/hg38_respectively> 
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
To check annotations for your genomic data, execute the following command:

css
Copy code
python check_annotations.py --input <input_file>
Replace <input_file> with the path to your genomic data file.


## Contact
If you have any questions or suggestions regarding the Genomic Annotations project, feel free to contact [insert your contact information].