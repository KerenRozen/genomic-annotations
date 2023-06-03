
# Artem run guide:

# In top of file:
exec(open("relative/path/to/run_gh38_regulations_function.py").read())
GH38_RAW_PATH = "relative/path/to/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz"

# In code:
match_regulations(1, 1, 2, 1)

# In run_gh38_annotation_function.py
GH38_RAW_PATH = "relative/path/to/Homo_sapiens.GRCh38.109.gff3.gz"
CLASSIFICATIONS_DB_PATH = "relative/path/to/where/the/DB/will/be/saved"
# In run_gh38_regulations_function.py
GH38_REGULATIONS_RAW_PATH = "relative/path/to/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz"
REGULATIONS_DB_PATH = "relative/path/to/where/the/DB/will/be/saved"