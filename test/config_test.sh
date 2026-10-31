# shellcheck shell=bash disable=SC2034
# Configuration for the simulated test run. Run from the repo root:
#   bash test/make_test_data.sh
#   bash scripts/vc_pipeline.sh test/config_test.sh
#   bash scripts/joint_genotyping.sh test/config_test.sh
#   python3 test/compare_to_truth.py

# shellcheck source=config/config.sh
source config/config.sh            # start from the main defaults, then override

REF="test/data/ref.fa"
SAMPLES="test/data/samples_test.tsv"
ADAPTERS="test/data/TruSeq3-PE.fa"
OUTDIR="test/results"
COHORT_NAME="test_cohort"
THREADS=2
JAVA_MEM="2g"
