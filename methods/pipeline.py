#!/usr/bin/env python
# coding: utf-8


import argparse
from pathlib import Path
from methods.kmers_distributions import *
from methods.run_classifiers import *
from methods.run_clustering import *
import sys
import time
from datetime import timedelta

sys.stdout.reconfigure(line_buffering=True) #doesn't buffer stdout redirected on file, prints line by line

'''
def format_time(t):
    # Convert seconds to h:m:s
    return time.strftime("%H:%M:%S", time.gmtime(t))

def format_time(t):
    """Convert seconds to a formatted string including days, hours, minutes, and seconds."""
    delta = timedelta(seconds=int(t))
    days = delta.days
    hours, remainder = divmod(delta.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    if days > 0:
        return f"{days}d {hours:02d}:{minutes:02d}:{seconds:02d}"
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"
'''

def format_time(t):
    """Convert seconds (float) to a formatted string including days, hours, minutes, 
    and seconds with 4 decimal places."""
    days, remainder = divmod(t, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, remainder = divmod(remainder, 60)
    seconds = remainder  # float with fractional seconds

    # Format seconds to have 2 integer digits and 4 decimal places (e.g., 05.1235)
    sec_str = f"{seconds:07.4f}"

    if days > 0:
        return f"{int(days)}d {int(hours):02d}:{int(minutes):02d}:{sec_str}"
    return f"{int(hours):02d}:{int(minutes):02d}:{sec_str}"


def parse_arguments():
    #Input arguments parser

    parser = argparse.ArgumentParser(description=(
        "Oligomer Distributions analysis.\n"
        "Note: --aflp_folder and --aflp_only are mutually exclusive.\n"
    ))

    parser.add_argument(
        "-i", "--input_dataset_path", type=Path, required=True, help="Path to input dataset folder.")
    
    parser.add_argument(
        "-o", "--result_folder_path", type=Path, default=None, help="Result folder. If not provided, it will be auto-generated as: ./results/results_<dataset_name>_k_s<kmer_start_size>_e<kmer_end_size>")
    
    parser.add_argument(
        "-s", "--kmer_start_size", type=int, required=True, help="Size of the left kmer.")
    parser.add_argument(
        "-e", "--kmer_end_size", type=int, required=True, help="Size of the right kmer.")
    
    parser.add_argument(
        "-cpu", "--cpu_processes", type=int, default=2, help = "Number of CPU processes to use, minimum of 2 (default: 2).")
    parser.add_argument(
        "-top", "--top_kmers_selection", type=int, default=20,help="Number of the top kmers to be selected by coefficient of variation (CoV).")

    '''
    parser.add_argument(
        "-min", "--min_fragment_length", type=int, default=0, help = "Minimum AFLP fragment length considered. Used only on classification.")
    parser.add_argument(
        "-max", "--max_fragment_length", type=int, default=1500, help = "Maximum AFLP fragment length considered. Used in the extraction of the AFLP distributions and classification.")
    '''
    parser.add_argument(
        "-f", "--fragment_length", type=int, default=1500, 
        help = "Maximum AFLP fragment length considered. Used in the extraction of the AFLP distributions and classification.")
    
    parser.add_argument(
        "--genus", action="store_true", help="Uses the genus of the genomes instead of the species as labels for classification."
    )
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--aflp_folder", type=Path,
        help="Path to a previously computed result AFLP folder (mutually exclusive with option --aflp_only)."
    )
    group.add_argument(
        "--aflp_only", action="store_true",
        help="Option: Computes only the AFLP distributions. (mutually exclusive with option --aflp_folder)"
    )
    
    return parser.parse_args()


def main():
    
    # Inputs
    args = parse_arguments()

    kmer_start_size = args.kmer_start_size
    kmer_end_size = args.kmer_end_size
    top_kmers_selection = args.top_kmers_selection
    cpu_processes = args.cpu_processes
    input_dataset_path = args.input_dataset_path
    result_folder_path = args.result_folder_path
    fragment_length = args.fragment_length
 
    genus_only_labels = args.genus
    aflp_only = args.aflp_only
    aflp_folder = args.aflp_folder

    window_size = fragment_length - (kmer_start_size + kmer_end_size) +1 #example: 1500 (fragment) - 6(kmer_start) - 6(kmer_end) + 1
    if window_size <= 0:
        raise ValueError("Calculated window_size must be greater than 0. Check your fragment or k-mer sizes.")
        
    #If no result folder path is given, results will be saved in a default result folder
    if result_folder_path is None:
        script_dir = Path(__file__).resolve().parent.parent # Path to the folder where the script is located

        genus_tag = "_genus" if genus_only_labels else ""
        results_folder_name = f"results_{input_dataset_path.stem}_k_s{kmer_start_size}_e{kmer_end_size}{genus_tag}"
        result_folder_path = Path(script_dir, "results", results_folder_name)
        result_folder_path.mkdir(parents=True, exist_ok=True)

    
    print("\nInput parameters:")
    print(f"  CPUs                     : {cpu_processes}")
    print(f"  input dataset            : {input_dataset_path}")
    print(f"  result folder            : {result_folder_path}")
    print(f"  kmer_start_size          : {kmer_start_size}")
    print(f"  kmer_end_size            : {kmer_end_size}")
    print(f"  top_kmers_selection      : {top_kmers_selection}")    
    print(f"  fragment length          : {fragment_length}")
    print(f"  search window size       : {window_size}")
    print(f"  genus only labels        : {genus_only_labels}")
        
    print()


    total_time_start = time.time()
    
    if not Path(result_folder_path, "kmer_pairs.csv").exists(): #if kmer pairs are not computed
        #Step 1: kmer pairs processing
        kmers_processing(input_dataset_path=input_dataset_path, result_folder_path=result_folder_path,
                         kmer_start_size=kmer_start_size, kmer_end_size=kmer_end_size, 
                         cpu_processes=cpu_processes, window_size = window_size
                        )
    else:
        print("kmer pairs already computed - Resuming...")

    
    if not aflp_folder:    
        #Step 2: AFLP features extraction
        aflp_folder = aflp_processing(input_dataset_path=input_dataset_path, 
                                      result_folder_path=result_folder_path,
                                      top_kmers_selection=top_kmers_selection, 
                                      cpu_processes=cpu_processes, 
                                      fragment_length=fragment_length
                                     )
        #Check the valid or invalid AFLPs
        aflp_validity_df = summarize_aflp_validity(aflp_folder)
        aflp_validity_df.to_csv(Path(result_folder_path,"aflp_validity.csv"), index= False)
    else:
        print("AFLP fragments distributions already computed")
    
    if aflp_only:
        print("Only AFLP option: True - resulting AFLP saved in folder:", aflp_folder)
        print("Program finished successfully.")
        return #terminate program

        
    #Step 3: classification 
    run_classifiers(input_folder=aflp_folder,
                    fragment_length=fragment_length, 
                    cpu_processes=cpu_processes, 
                    genus_only_labels=genus_only_labels
                   )
    
    #Step 4: clustering
    run_clustering(input_folder=aflp_folder, cpu_processes=cpu_processes)
    

    total_time_elapsed = time.time() - total_time_start
    print("\nProgram Total elapsed time:", format_time(total_time_elapsed))
    
    print("Program finished successfully.")
