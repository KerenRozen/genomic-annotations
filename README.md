# Genomic Annotations
Genomic Annotations is a project aimed at annotating genomic data using data from various databases. The primary objective is to ensure fast annotation to keep up with AI model training speed.

## Setup and Installation
To set up and run the project, follow these steps:

1. Clone the repository:

```bash
git clone <https://github.com/KerenRozen/genomic_annotations>
```

2. Install the required packages. Make sure you have Python 3.6 installed. Use the following command to install the dependencies:
```bash
pip install -r requirements.txt
```
3. Download the necessary annotation files. The following files are required for each type of annotation:
   1. 164 cell type regulation annotations:
      - [link](https://noble.gs.washington.edu/proj/encyclopedia/segway_encyclopedia.bed.gz)
      - need to add for hg38
   2. Classifications
   3. Regulation regions 

4. Build the annotation database. Run the following command to initiate the database building process:
```bash
python3 build_cell_type_regulation_db.py /path/to/segway_encyclopedia.bed/file /path/to/save/DB/file genome_reference 
```
genome_reference should be 37/38 for hg37/hg38 respectively.
This command will process the annotation files and create the necessary database for annotations in the path you provided.

## Usage
To check annotations for your genomic data, execute the following command:

css
Copy code
python check_annotations.py --input <input_file>
Replace <input_file> with the path to your genomic data file.


## Contact
If you have any questions or suggestions regarding the Genomic Annotations project, feel free to contact [insert your contact information].