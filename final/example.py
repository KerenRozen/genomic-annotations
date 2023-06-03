
# Artem run guide:

# In top of file:
exec(open("relative/path/to/run_gh38_regulations_function.py").read())
GH38_RAW_PATH = "relative/path/to/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz"

# In code:
match_regulations(1, 1, 2, 1)
