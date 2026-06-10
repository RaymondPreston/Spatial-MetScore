from pathlib import Path

# Base project directory
# Since config.py is at: project_root/spatial_metscore/spatial_metscore/config.py
# project_dir is project_root/
project_dir = Path(__file__).resolve().parents[2]

data_dir = project_dir / "data"
raw_data_dir = data_dir / "raw"
processed_data_dir = data_dir / "processed"
cleaned_data_dir = data_dir / "cleaned"

samples_dir = raw_data_dir / "Samples"
metadata_file = raw_data_dir / "patient_metadata.csv"
high_genes_file = raw_data_dir / "met-high-genes-human.csv"
low_genes_file = raw_data_dir / "met-low-genes-human.csv"
high_low_xlsx = raw_data_dir / "MetHigh-MetLow Genes.xlsx"

# External metadata
ZENODO_RECORD_ID = "14199536"
