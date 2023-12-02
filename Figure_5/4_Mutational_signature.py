#############################################################
# Figure 5
# TCGA COAD/READ cohort - Mutational signature analysis
#############################################################
# Reformat 96 matrix input data to match requirement
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerMatrixGenerator import install as genInstall
import pandas as pd

# WES_TCGA.96.csv was downloaded from https://doi.org/10.1038/s41586-022-05600-5
# read in input TCGA data
file = pd.read_csv("MutationalSignature/Data/WES_TCGA.96.csv")
df = pd.DataFrame(file)

# Template from software example
template_file = pd.read_table("MutationalSignature/Data/BRCA/BRCA.txt")
template = template_file["Mutation Types"]

# Use their name
df["Mutation Types"] = template
first_column = df.pop("Mutation Types")
df.insert(0, 'Mutation Types', first_column)

# Remove unnecessary columns
df = df.drop(columns=["Mutation type", "Trinucleotide"])

# Save new format
df.to_csv("MutationalSignature/Data/WES_TCGA.96_reformatted.txt",
          sep="\t", index=False)


# Do this once
# genInstall.install('GRCh37')

Analyze.cosmic_fit(samples="MutationalSignature/Data/WES_TCGA.96_reformatted.txt",
                   output="MutationalSignature/Analysis/TCGA_WES_96/normalised/",
                   input_type="matrix",
                   exome=True,
                   exclude_signature_subgroups=[
                       'Chemotherapy_signatures',
                       'Treatment_signatures'
                       'Artifact_signatures',
                       'Lymphoid_signatures'],
                   verbose=True)