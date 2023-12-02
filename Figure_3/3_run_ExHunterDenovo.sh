#############################################################
# Figure 3
# Run Expansion Hunter Denovo analysis
#############################################################
# Sample SRA ID lists
all_wgs_samples=(SRR8639191 SRR8652098 SRR8670708 SRR8652101 SRR8652121 SRR8652136 SRR8670772 SRR11680468 SRR11680467 SRR8670762 SRR8670679 SRR8670673 SRR8670700 SRR8639227 SRR8639145 SRR8652076)

# Expansion detection by ExpansionHunter and exSTRa
# Expansions are identified for broken and non-broken (TA)n from the HCT116 and KM12 whole genome sequencing data by ExpansionHunter v.3.2.225,47. The supporting reads for expansions are visualized by GraphAlignmentViewer (https://github.com/Illumina/GraphAlignmentViewer). Empirical cumulative distribution function for the TA repeat located at hg19 coordinates chr8:106950919–106950985 was generated using exSTRa v.0.89.026 and Bio-STR-exSTRa v.1.1.0 using the default parameter settings.
source activate exphunter_env
ExHunt_dir=programs/ExpansionHunterDenovo-v0.9.0-linux_x86_64/
bam_dir=data/wgs/
profile_dir=analysis/ExpHunterDenovo/str_profiles
ref_fasta=alignment_references/Homo_sapiens/hg38_no_alt/genome/fasta/hg38_no_alt.fa

cd analysis/ExpHunterDenovo/
# 1. Create str profiles for each sample 
for SRA_ID in ${all_wgs_samples[@]}
do
    FILE=(./str_profiles/${SRA_ID}.locus.tsv)
    if [ -f "$FILE" ]; then
        echo "$FILE exists. \n"
    else 
        echo "$FILE doesn't exists. running profile"
        ${ExHunt_dir}/bin/ExpansionHunterDenovo profile \
        --reads ${bam_dir}/${SRA_ID}.sorted.bam \
        --reference ${ref_fasta} \
        --output-prefix str_profiles/${SRA_ID} \
        --min-anchor-mapq 50 \
        --max-irr-mapq 40 &
    fi    
done

# 2. copy paste the following into a manifest 
vim manifest_MSI_LOF_WT.tsv
SRR8639191  case  str_profiles/SRR8639191.str_profile.json
SRR8652098  case  str_profiles/SRR8652098.str_profile.json
SRR8652121  case  str_profiles/SRR8652121.str_profile.json
SRR8652136  case  str_profiles/SRR8652136.str_profile.json
SRR8670708  case  str_profiles/SRR8670708.str_profile.json
SRR8652101  case  str_profiles/SRR8652101.str_profile.json
SRR11680468 case str_profiles/SRR11680468.str_profile.json
SRR11680467 case str_profiles/SRR11680467.str_profile.json
SRR8670762  case  str_profiles/SRR8670762.str_profile.json
SRR8670772  control   str_profiles/SRR8670772.str_profile.json
SRR8670679  control   str_profiles/SRR8670679.str_profile.json

vim manifest_MSI_LOF_MSS_WT.tsv
SRR8639191  case  str_profiles/SRR8639191.str_profile.json
SRR8652098  case  str_profiles/SRR8652098.str_profile.json
SRR8652121  case  str_profiles/SRR8652121.str_profile.json
SRR8652136  case  str_profiles/SRR8652136.str_profile.json
SRR8670708  case  str_profiles/SRR8670708.str_profile.json
SRR8652101  case  str_profiles/SRR8652101.str_profile.json
SRR11680468 case str_profiles/SRR11680468.str_profile.json
SRR11680467 case str_profiles/SRR11680467.str_profile.json
SRR8670762  case  str_profiles/SRR8670762.str_profile.json
SRR8670673  control   str_profiles/SRR8670673.str_profile.json
SRR8670700  control   str_profiles/SRR8670700.str_profile.json
SRR8639227  control   str_profiles/SRR8639227.str_profile.json

# 3. Merge into a pilot study 
${ExHunt_dir}/bin/ExpansionHunterDenovo merge \
    --reference ${ref_fasta} \
    --manifest manifest_MSI_LOF_WT.tsv \
    --output-prefix MSI_LOF_WT &

${ExHunt_dir}/bin/ExpansionHunterDenovo merge \
    --reference ${ref_fasta} \
    --manifest manifest_MSI_LOF_MSS_WT.tsv \
    --output-prefix MSI_LOF_MSS_WT &
# 4. The next step is to actually compare STR lengths between cases and controls. Two types of comparisons are supported: locus-based comparison and motif-based comparison. These comparisons are performed by a Python3 script casecontrol.py located in the scripts/ directory.

# The locus-based comparison can distinguish between repeats longer than the read length but shorter than the fragment length. It can be run like so:
${ExHunt_dir}/scripts/casecontrol.py locus \
        --manifest manifest_MSI_LOF_WT.tsv \
        --multisample-profile MSI_LOF_WT.multisample_profile.json \
        --output MSI_LOF_WT.casecontrol_locus.tsv &

${ExHunt_dir}/scripts/casecontrol.py locus \
        --manifest manifest_MSI_LOF_MSS_WT.tsv \
        --multisample-profile MSI_LOF_MSS_WT.multisample_profile.json \
        --output MSI_LOF_MSS_WT.casecontrol_locus.tsv &

# 5. The motif-based comparison analyzes the overall enrichment of genomes with repeats longer than the fragment length. It can be run like so:
${ExHunt_dir}/scripts/casecontrol.py motif \
        --manifest manifest_MSI_LOF_WT.tsv \
        --multisample-profile MSI_LOF_WT.multisample_profile.json \
        --output MSI_LOF_WT.casecontrol_motif.tsv

${ExHunt_dir}/scripts/casecontrol.py motif \
        --manifest manifest_MSI_LOF_MSS_WT.tsv \
        --multisample-profile MSI_LOF_MSS_WT.multisample_profile.json \
        --output MSI_LOF_MSS_WT.casecontrol_motif.tsv