#!/usr/bin/env python3
import pandas as pd
import subprocess
import argparse
from pathlib import Path
import zipfile
import shutil


def run_with_retry(cmd):
    #Run a shell command, retrying until it succeeds.
    while True:
        try:
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Command failed with error {e.returncode}. Retrying...")

def download_in_batches(input_tsv, output_dir, batch_size = 100):
    # Ensure output directory exists
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load TSV file
    df = pd.read_csv(input_tsv, sep="\t")

    # Extract Assembly Accessions
    if "Assembly Accession" not in df.columns:
        raise ValueError("Input TSV must contain 'Assembly Accession' column")
    accessions = df["Assembly Accession"].dropna().tolist()

    # Split into batches
    total_batches = (len(accessions) + batch_size - 1) // batch_size  # ceiling division
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i+batch_size]
        part_num = (i // batch_size) + 1
        output_file = output_path / f"{Path(output_dir).stem}_part{part_num}.zip"

        if output_file.exists():
            print(f"Batch {part_num} of {total_batches} with {len(batch)} accessions...")
            print(f"File {output_file} already exists, skipping.")
            continue
		    
        # Build command as string
        batch_str = " ".join(batch)

        #Use with NCBI dataset as conda package
        cmd = f"datasets download genome accession {batch_str} --filename {output_file}"
        
        #Use with standalone executable (downloadable from NCBI)
        #cmd = f"./methods/NCBI_datasets_CLI/datasets download genome accession {batch_str} --filename {output_file}"

        print(f"Batch {part_num} of {total_batches} with {len(batch)} accessions...")
        print(f"Downloading file {output_file}")
        
        run_with_retry(cmd)


def extract_fasta_from_zips(input_dir, exts=(".fna", ".fa", ".fasta")):
    input_path = Path(input_dir)
    zip_dir = input_path / f"{input_path.stem}_zips"
    zip_dir.mkdir(exist_ok=True)  # create "zip" folder if not exists

    for zip_file in input_path.glob("*.zip"):
        print(f"Processing {zip_file.name}...")
        with zipfile.ZipFile(zip_file, "r") as zf:
            for member in zf.namelist():
                # Check if it's a FASTA file by extension
                if member.lower().endswith(exts):
                    # Build output filename (flatten structure)
                    out_name = Path(member).name
                    out_path = input_path / out_name
                    # Extract file content
                    with zf.open(member) as source, open(out_path, "wb") as target:
                        shutil.copyfileobj(source, target)
                    #print(f"  Extracted {out_name} -> {out_path}")

        # Move the processed zip into the "zip" folder
        shutil.move(str(zip_file), zip_dir / zip_file.name)
        #print(f"Moved {zip_file.name} to {zip_dir}/")
    
    zip_archive = Path(input_dir.parent, "zips")
    zip_archive.mkdir(parents=True, exist_ok=True)
    
    # Remove if already exists
    zip_stored_path = Path(zip_archive, zip_dir.stem)
    if zip_stored_path.exists():
        print("Removing folder", zip_stored_path)
        shutil.rmtree(zip_stored_path)
    
    shutil.move(zip_dir, zip_stored_path)
    
    print(f"Moved zip files {zip_dir} to {zip_stored_path}")


def organize_fasta_by_species(input_dir: str):
    input_path = Path(input_dir)

    for fasta_file in input_path.glob("*.fna"):
        with open(fasta_file, "r") as f:
            header = f.readline().strip()

        if not header.startswith(">"):
            print(f"Skipping {fasta_file.name} (no valid FASTA header)")
            continue

        # Split header by spaces, species name is 2nd and 3rd word
        parts = header.split()
        if len(parts) < 3:
            print(f"Skipping {fasta_file.name} (header too short: {header})")
            continue

        genus = parts[1]
        species = parts[2].rstrip(",")
        species_folder = f"{genus}_{species}"

        # Create folder for species
        species_path = input_path / species_folder
        species_path.mkdir(exist_ok=True)

        # Move file into species folder
        target_file = species_path / fasta_file.name
        shutil.move(str(fasta_file), target_file)
        

def parse_arguments():
    #Input arguments parser

    parser = argparse.ArgumentParser(description="Download NCBI genomes in batches")
    
    parser.add_argument(
        "-i", "--input_tsv", type=Path, required=True, help="Input TSV file with 'Assembly Accession' column")
    
    """
    parser.add_argument(
        "-o", "--dataset_folder", type=Path, default=None, help="Dataset folder of the downloaded fasta.")
    """
    
    parser.add_argument("--batch_size", type=int, default=100, help="Number of genomes per batch (default=100)")

    return parser.parse_args()

    
def main():
    
    # Inputs
    args = parse_arguments()
    
    input_tsv = args.input_tsv
    #dataset_folder = args.dataset_folder
    batch_size = args.batch_size

    
    if True: #dataset_folder is None:
        script_dir = Path(__file__).resolve().parent # Path to the folder where the script is located

        dataset_name = input_tsv.stem
        dataset_folder = Path(script_dir, "data", dataset_name)
        dataset_folder.mkdir(parents=True, exist_ok=True)

    print("dataset folder:", dataset_folder)
    
    download_in_batches(input_tsv, dataset_folder, batch_size)
    print()
    extract_fasta_from_zips(dataset_folder)
    print()
    organize_fasta_by_species(dataset_folder)
    
    print("Download dataset: Done!")
    
if __name__ == "__main__":
    main()

