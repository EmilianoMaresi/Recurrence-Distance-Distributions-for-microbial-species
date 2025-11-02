#!/usr/bin/env python
# coding: utf-8

import itertools
import time
import numpy as np
import csv
import pandas as pd

import subprocess
import tempfile
from pathlib import Path
from Bio import SeqIO

import multiprocessing as mp
from Bio import SeqIO
from datetime import datetime, timedelta


#Step1: processing methods
def format_time(t):
    # Convert seconds to h:m:s
    return time.strftime("%H:%M:%S", time.gmtime(t))

def normalize_by_total_kmer_count(kmer_multiplicity_df):
    """
    Normalizes a k-mer multiplicity DataFrame by the sum of each k-mer's counts across genomes (row-wise).
    Rows with total count zero are filled with zeros.

    Parameters:
        kmer_multiplicity_df (pd.DataFrame): DataFrame with k-mers as rows and genome names as columns.

    Returns:
        pd.DataFrame: Normalized DataFrame (each row divided by its sum, rows with zero sum filled with zeros).
    """
    # Compute the row sums
    row_sums = kmer_multiplicity_df.sum(axis=1)

    # Avoid division by zero: replace zero with 1 temporarily, then fix those rows after
    safe_row_sums = row_sums.replace(0, 1)

    # Normalize by dividing each row by its sum
    kmer_normalized_df = kmer_multiplicity_df.div(safe_row_sums, axis=0)

    # Set rows that originally had a sum of zero to all zeros
    kmer_normalized_df[row_sums == 0] = 0

    return kmer_normalized_df

def normalize_by_genome_length(kmer_multiplicity_df, genomes_length_dict):
    """
    Normalizes a k-mer multiplicity DataFrame by genome lengths.

    Parameters:
        kmer_multiplicity_df (pd.DataFrame): DataFrame with k-mers as rows and genome names as columns.
        genomes_length_dict (dict): Dictionary mapping genome names (matching columns) to genome lengths.

    Returns:
        pd.DataFrame: Normalized DataFrame (multiplicity divided by genome length).

    Raises:
        ValueError: If any genome name in the DataFrame columns is missing in the dictionary or vice versa.
    """
    # Create a Series from the dictionary
    genome_length_series = pd.Series(genomes_length_dict)

    # Strip column and dictionary names of whitespace (optional but helpful)
    kmer_multiplicity_df.columns = kmer_multiplicity_df.columns.str.strip()
    genome_length_series.index = genome_length_series.index.str.strip()

    # Check for mismatches
    missing_in_dict = set(kmer_multiplicity_df.columns) - set(genome_length_series.index)
    extra_in_dict = set(genome_length_series.index) - set(kmer_multiplicity_df.columns)

    if missing_in_dict:
        raise ValueError(f"These genomes are in the DataFrame but missing in the genome length dictionary: {missing_in_dict}")
    if extra_in_dict:
        raise ValueError(f"These genomes are in the genome length dictionary but not in the DataFrame: {extra_in_dict}")

    # Perform normalization (will align columns/index automatically)
    kmer_normalized_df = kmer_multiplicity_df / genome_length_series

    return kmer_normalized_df


def compute_coefficient_of_variation(df_normalized):
    """
    Computes the coefficient of variation (CV = std / mean) for each row of a normalized DataFrame,
    excluding rows with zero mean, and returns the result sorted in descending order of CV.

    Parameters:
        df_normalized (pd.DataFrame): DataFrame with rows as k-mers and columns as genomes.

    Returns:
        pd.DataFrame: A DataFrame with k-mers as index and one column: 'coefficient_of_variation',
                      sorted from highest to lowest CV.
    """
    mean = df_normalized.mean(axis=1)
    std = df_normalized.std(axis=1)

    # Avoid division by zero: keep only rows where mean > 0
    valid = mean > 0
    cv = std[valid] / mean[valid]

    # Convert to DataFrame and sort
    cv_df = cv.to_frame(name='coefficient_of_variation')
    cv_df_sorted = cv_df.sort_values(by='coefficient_of_variation', ascending=False)

    return cv_df_sorted

#From a DataFrame of k-mers (rows) and multiplicities (values), group by first `kmer_size` bases.
def aggregate_kmer_multiplicities(dataframe,kmer_size):
    return dataframe.groupby(dataframe.index.str[:kmer_size]).sum()


def createDictMoltep(dictionary, megaset):
    multiplicities = {element: [] for element in megaset}
    genomes_lst = list()
    for genome_name, kmers in dictionary.items():
        for kmer, lst in multiplicities.items():
            lst.append(kmers[kmer])
        genomes_lst.append(genome_name)
    return multiplicities, genomes_lst

#reads the dolier output and converts it into a dictionary "kmer": multiplicity 
def dolier_to_dictionary(filepath):
    # Read the file as a tab-separated dataframe with no header
    df = pd.read_csv(filepath, sep='\t', header=None, names=['kmer', 'count'])

    # Convert 'count' column to integer (optional but ensures correct types)
    df['count'] = df['count'].astype(int)

    # Convert to dictionary: kmer as key, count as value
    return df.set_index('kmer')['count'].to_dict()

#dict_of_dicts: the dictionary with structure {genome_id:{kmer1:multiplicity, kmer2:multiplicity}}
def find_dictionaries_intersection(dict_of_dicts):
    if not dict_of_dicts:
        return set()
    # Extract the list of inner dictionaries
    inner_dicts = dict_of_dicts.values()
    # Start with keys of the first inner dictionary
    iterator = iter(inner_dicts)
    common = set(next(iterator).keys())
    # Intersect with keys from the remaining inner dictionaries
    for d in iterator:
        common &= set(d.keys())
    return common

def get_genome_length(file_path):
    total_length = 0
    for record in SeqIO.parse(file_path, "fasta"):
        total_length += len(record.seq)
    return total_length

def make_kmer_pairs(kmers_start,kmers_end):
    all_pairs = list(itertools.product(kmers_start,kmers_end))
    pairs_df = pd.DataFrame(all_pairs, columns=["start", "end"])
    return pairs_df



def run_dolier(input_fasta_path, output_path, kmer_size, threads=1, mismatches=0, verbose=False):
    """ Call to the DoLier script for fast k-mers count"""

    script_dir = Path(__file__).resolve().parent
    dolier_path = Path(script_dir, 'dolier', 'dolier-kfreqs')

    cmd = [
        dolier_path, #'./dolier/dolier-kfreqs',
        str(input_fasta_path),
        str(output_path),
        str(kmer_size),
        str(mismatches),
        '--kmers',
        '-p', str(threads)
    ]

    stdout = None if verbose else subprocess.DEVNULL

    try:
        subprocess.run(cmd, check=True, stdout=stdout, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"DoLIer failed on {input_fasta_path}")
        print(f"Error: {e.stderr.decode().strip()}")
        raise


def kmers_processing(input_dataset_path, result_folder_path, kmer_start_size, kmer_end_size, 
                    top_kmers_selection, cpu_processes):
    """Wrapper method for the first step of the pipeline. 
    - Find the k-mers counts in the fasta genomes
    - Computes the coefficient of variation of each k-mer
    - Selects the top k-mer with highest Coefficient of Variation (CoV)
    """

    print("[[Step1 - Kmer processing]]")
    # Determine bigger and smaller k-mer
    kmer_big = max(kmer_start_size,kmer_end_size)
    kmer_small = min(kmer_start_size,kmer_end_size)

    dictionary_kmers = dict()
    dictionary_counts = dict()
    genomes_length = dict()

    print(f"\nComputation for the bigger kmer of size {kmer_big}")
    print("Finding kmer occurrencies with DoLier:")
    files = [f for ext in ("*.fna", "*.fa", "*.fasta") for f in input_dataset_path.rglob(ext)]
    nof_files = len(files)

    start_time = time.time()

    for i, filename in enumerate(files, start=1):
        print(f"({i}/{nof_files})\t{filename.parent.stem}\t{filename.name}") 

        genome_name = filename.name
        genomes_length[genome_name] = get_genome_length(filename)

        # Create a temporary file for DoLIer output
        with tempfile.NamedTemporaryFile(delete=True, mode='w+', suffix=".tsv") as tmpfile:
            tmpfile_path = tmpfile.name
            
            run_dolier(filename, tmpfile_path, kmer_big, threads=cpu_processes)
            dictionary_kmers[genome_name] = dolier_to_dictionary(tmpfile_path)

    elapsed_time = time.time() - start_time

    print(f"DoLier done! - Execution time: {format_time(elapsed_time)}") #(end_time - start_time):.3f}s")
    print()


    genomes_length_df = pd.DataFrame({
        'genome': list(genomes_length.keys()),
        'genome_length': list(genomes_length.values())
    })
    genomes_length_df.to_csv(Path(result_folder_path,f"genomes_length.csv"), index=False)

    #select kmers that are shared among all genomes
    kmers_intersection = find_dictionaries_intersection(dictionary_kmers)
    
    # Convert to DataFrame:
    kmer_big_multiplicities, genomes_lst = createDictMoltep(dictionary_kmers, kmers_intersection)
    kmer_big_multiplicities_df = pd.DataFrame.from_dict(kmer_big_multiplicities, orient='index', columns=genomes_lst)
    

    print("Normalizing kmer occurrences...")
    kmer_big_normalized_df = normalize_by_total_kmer_count(kmer_big_multiplicities_df)

    print("Computing coefficient of variation (CoV)...")
    kmer_big_cov_df = compute_coefficient_of_variation(kmer_big_normalized_df)

    #select top kmers with highest Coefficient of Variation (CoV)
    print(f"Select top {top_kmers_selection} kmers with the highest coefficient of variation (CoV)...")
    kmer_big_top_cov = kmer_big_cov_df.iloc[:top_kmers_selection]


    #If the kmers have different sizes, for the small size kmers find their multiplicities from the occurrences of the bigger kemr, 
    # then compute Normalization and CoV
    if(kmer_big != kmer_small):

        print(f"\nComputation for the smaller kmer of size {kmer_small}")
        kmer_small_multiplicities_df = aggregate_kmer_multiplicities(kmer_big_multiplicities_df,kmer_small)

        #compute normalization
        print("Normalizing kmer occurrences...")
        kmer_small_normalized_df = normalize_by_total_kmer_count(kmer_small_multiplicities_df)

        #compute cov
        print("Computing coefficient of variation (CoV)...")
        kmer_small_cov_df = compute_coefficient_of_variation(kmer_small_normalized_df)

        print(f"Select top {top_kmers_selection} kmers with the highest coefficient of variation (CoV)...")
        kmer_small_top_cov = kmer_small_cov_df.iloc[:top_kmers_selection]

        # Reassign based on the original input sizes
        if kmer_start_size == kmer_big:
            kmer_start_multiplicities = kmer_big_multiplicities_df
            kmer_start_cov = kmer_big_cov_df
            kmer_start_top_cov = kmer_big_top_cov

            kmer_end_multiplicities   = kmer_small_multiplicities_df
            kmer_end_cov  = kmer_small_cov_df
            kmer_end_top_cov = kmer_small_top_cov

        else:
            kmer_start_multiplicities = kmer_small_multiplicities_df
            kmer_start_cov = kmer_small_cov_df
            kmer_start_top_cov = kmer_small_top_cov

            kmer_end_multiplicities   = kmer_big_multiplicities_df
            kmer_end_cov  = kmer_big_cov_df
            kmer_end_top_cov = kmer_big_top_cov

        #compute kmer pairs for those selected by CoV
        print(f"Computing kmer pairs {kmer_start_size}/{kmer_end_size}...")
        kmer_pairs = make_kmer_pairs(kmer_start_top_cov.index, kmer_end_top_cov.index)
        kmer_pairs.to_csv(Path(result_folder_path,f"kmer_pairs.csv"),index=False)

        #SAVE as CSV
        print("Saving...")
        kmer_start_multiplicities.to_csv(Path(result_folder_path,f"kmer_start_k{kmer_start_size}_multiplicities.csv"))
        kmer_start_cov.to_csv(Path(result_folder_path,f"kmer_start_k{kmer_start_size}_cov.csv"),index=True)
        kmer_start_top_cov.to_csv(Path(result_folder_path,f"kmer_start_k{kmer_start_size}_cov_top.csv"),index=True)

        kmer_end_multiplicities.to_csv(Path(result_folder_path,f"kmer_end_k{kmer_end_size}_multiplicities.csv"))
        kmer_end_cov.to_csv(Path(result_folder_path,f"kmer_end_k{kmer_end_size}_cov.csv"),index=True)
        kmer_end_top_cov.to_csv(Path(result_folder_path,f"kmer_end_k{kmer_end_size}_cov_top.csv"),index=True)



    elif(kmer_big == kmer_small):

        print(f"Computing kmer pairs {kmer_start_size}/{kmer_end_size}...")
        #compute kmer pairs for those selected by CoV
        kmer_pairs = make_kmer_pairs(kmer_big_top_cov.index,kmer_big_top_cov.index)
        kmer_pairs.to_csv(Path(result_folder_path,f"kmer_pairs.csv"),index=False)

        print("Saving...")
        kmer_big_multiplicities_df.to_csv(Path(result_folder_path,f"kmer_k{kmer_start_size}_multiplicities.csv"))
        kmer_big_cov_df.to_csv(Path(result_folder_path,f"kmer_k{kmer_start_size}_cov.csv"),index=True)
        kmer_big_top_cov.to_csv(Path(result_folder_path,f"kmer_k{kmer_start_size}_cov_top.csv"),index=True)        

    print("\nStep 1 - Kmer processing: Done!")
    print()

    return result_folder_path


##################################################
#Step 2: Extraction of kmers distributions methods
def format_daytime(t):
    formatted_time = datetime.fromtimestamp(time.time()).strftime("%Y-%m-%d %H:%M:%S")
    return formatted_time

def get_fasta_files_sorted_by_size(dataset_folder):
    """
    Sorts fasta genomes by file size.
    """

    fasta_files = [f for ext in ("*.fna", "*.fa", "*.fasta") for f in dataset_folder.rglob(ext)]

    # Sort fasta files by descending size
    fasta_files_sorted = sorted(
        fasta_files,
        key=lambda f: f.stat().st_size,
        reverse=True
    )

    return fasta_files_sorted    


# Methods to compute the AFLP distributions
def find_kmer_positions(sequence, kmer):
    """Finds the starting position of all the occurrences 
    of the kmer in a sequence string"""

    positions = []
    i = sequence.find(kmer)
    while i != -1:
        positions.append(i)
        i = sequence.find(kmer, i + 1)
    return positions

def next_occurrence_from(array, current_index, from_pos):
    while current_index < len(array):
        if array[current_index] >= from_pos:
            break
        else:
            current_index += 1

    if current_index >= len(array):
        return False, current_index
    return True, current_index

def extract_fragments_occurrences(sequence, kmer_start, kmer_end, mindist, maxdist):
    fragments_occurrences = dict()

    kmer_start_pos = find_kmer_positions(sequence, kmer_start)
    kmer_end_pos = find_kmer_positions(sequence, kmer_end)
    kmer_start_index = 0
    kmer_end_index = 0

    kmer_start_found, kmer_start_index = next_occurrence_from(kmer_start_pos, kmer_start_index, 0)

    while kmer_start_index < len(kmer_start_pos):

        kmer_end_found, kmer_end_index_aux = next_occurrence_from(kmer_end_pos, kmer_end_index, 
                                                                  kmer_start_pos[kmer_start_index]+len(kmer_start))

        if kmer_end_found:

            fragment_size = kmer_end_pos[kmer_end_index_aux] + len(kmer_end) - kmer_start_pos[kmer_start_index]

            if (fragment_size >= mindist) and (fragment_size <= maxdist):
                fragments_occurrences[fragment_size] = fragments_occurrences.get(fragment_size, 0) + 1

            kmer_end_index = kmer_end_index_aux

        else:
            break

        kmer_start_index += 1

    return fragments_occurrences

def compute_aflp_distribution(genome_file_path):

    result = list()
    species = genome_file_path.parent.stem

    for kmer_start, kmer_end in kmer_pairs:

        combined = {} #forward + reverse occurrences

        #species = None
        accession_id = None

        for record in SeqIO.parse(genome_file_path, "fasta"):

            #if species is None:  # First record
            #    species = " ".join(record.description.split(" ")[1:3])
            if accession_id is None:
                accession_id = record.id

            forward_counts = extract_fragments_occurrences(str(record.seq), kmer_start, kmer_end, mindist, maxdist)
            reverse_complement_counts = extract_fragments_occurrences(str(record.seq.reverse_complement()), 
                                                                      kmer_start, kmer_end, mindist, maxdist)

            for length, count in forward_counts.items():
                combined[length] = combined.get(length,0) + count

            for length, count in reverse_complement_counts.items():
                combined[length] = combined.get(length,0) + count

        # Define the range of fragment lengths to cover (assuming mindist to maxdist)
        fragment_lengths = range(mindist, maxdist + 1)

        # Create the AFLP distribution, adds fragment lengths with 0 occurrences
        aflp_distribution = {length: combined.get(length, 0) for length in fragment_lengths}

        result.append((genome_file_path.name, accession_id, species, kmer_start, kmer_end, aflp_distribution))

    return result

def aflp_processing(input_dataset_path, result_folder_path, cpu_processes, 
                    min_fragment_length, max_fragment_length):
    """
    Second step wrapper, processing of the AFLP distributions
    """
    global mindist
    global maxdist
    global kmer_pairs

    mindist = min_fragment_length
    maxdist = max_fragment_length
    
    kmer_pairs = list(pd.read_csv(Path(result_folder_path,'kmer_pairs.csv')).itertuples(index=False, name=None))

    genome_files = get_fasta_files_sorted_by_size(input_dataset_path)

    output_folder = Path(result_folder_path, "aflp") 
    output_folder.mkdir(parents=True, exist_ok=True)

    print("[[Step 2 - Extract AFLP distributions]]")
    print()
    
    file_handles = {} #dictionary of opened files handles (Ubuntu has a limit of 1000 concurrent open files)

    tic = time.time()
    print(f"Computing distributions - Start Time: {format_daytime(tic)}")

    with mp.Pool(processes=cpu_processes) as pool:

        columns = ["file","accession_id","species","kmer_start","kmer_end"] + [str(i) for i in range(mindist, maxdist+1)]

        processed_files = 0
        for result in pool.imap_unordered(compute_aflp_distribution, genome_files, chunksize=1):
            for r in result:
                input_filename = r[0]
                accession_id = r[1]
                species = r[2]
                kmer_start = r[3]
                kmer_end = r[4]
                aflp_distribution = r[5]

                #ubuntu can have around 1000 opened files before reaching the limit
                file_name = kmer_start + "-" + kmer_end
                if file_name not in file_handles:
                    f = open(output_folder/f"{file_name}.csv","w")
                    f.write(",".join(columns) + "\n")
                    file_handles[file_name] = f

                #write                
                csv_row = [input_filename, accession_id, species, kmer_start, kmer_end] + [str(aflp_distribution[k]) for k in range(mindist, maxdist+1)]
                file_handles[file_name].write(",".join(csv_row) + "\n")

            processed_files += 1

            if processed_files % 1 == 0: #30 == 0:
                print(f"[Elapsed time: {str(timedelta(seconds=time.time()-tic)).split('.')[0]}]\tProcessed genome files: {processed_files} out of {len(genome_files)}.")


    #close all the opened file handles
    for f in file_handles.values():
        f.close()

    elapsed_time = time.time() - tic

    print(f"[Elapsed time: {str(timedelta(seconds=time.time()-tic)).split('.')[0]}]\tProcessed files: {processed_files} out of {len(genome_files)}.")
    print("\nStep 2 - Extract AFLP distributions: Done!")
    print()

    return output_folder #aflp distributions folder

