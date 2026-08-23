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

def get_fasta_length(file_path):
    total_length = 0
    for record in SeqIO.parse(file_path, "fasta"):
        total_length += len(record.seq)
    return total_length

def get_dataset_info(fasta_files):
    rows = list()
    
    for file in fasta_files:
        #print(file)
        
        filename = file.name
        species = file.parent.name
        fasta_length = get_fasta_length(file)
        
        #print(filename, species, fasta_length)
        
        rows.append({
            "filename":filename,
            "species":species, 
            "length":fasta_length
        })

    df = pd.DataFrame(rows)
    df.sort_values(by='length', inplace=True, ascending=False)
    #genomes_length_df.to_csv(Path(result_folder_path, "genomes_length.csv"), index=True)
    
    return df

#From a DataFrame of k-mers (rows) and multiplicities (values), group by first `kmer_size` bases.
def aggregate_kmer_multiplicities(dataframe,kmer_size):
    return dataframe.groupby(dataframe.index.str[:kmer_size]).sum()


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

def createDictMoltep(dictionary, megaset):
    multiplicities = {element: [] for element in megaset}
    genomes_lst = list()
    for genome_name, kmers in dictionary.items():
        for kmer, lst in multiplicities.items():
            lst.append(kmers[kmer])
        genomes_lst.append(genome_name)
    return multiplicities, genomes_lst

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

# Extract kmer counts using DoLier
def extract_kmer_counts(fasta_files, kmer_size, cpu_processes):
    dictionary_kmers = dict()

    print(f"Run Dolier using {cpu_processes} CPUs")
    start_time = time.time()

    for i, filename in enumerate(fasta_files, start=1):
        print(f"({i}/{len(fasta_files)})\t{filename.parent.stem}\t{filename.name}") 

        genome_name = filename.name

        # Create a temporary file for DoLIer output
        with tempfile.NamedTemporaryFile(delete=True, mode='w+', suffix=".tsv") as tmpfile:
            tmpfile_path = tmpfile.name
            
            run_dolier(filename, tmpfile_path, kmer_size, threads=cpu_processes)
            dictionary_kmers[genome_name] = dolier_to_dictionary(tmpfile_path)

    elapsed_time = time.time() - start_time
    print(f"Extraction kmers - Elapsed Time: {format_time(elapsed_time)}")
    
    return dictionary_kmers

# creates the matrix of (kmer, genomes) with their multiplicities
def extract_kmer_multiplicities(fasta_files, kmer_size, cpu_processes):
    kmers_counts = extract_kmer_counts(fasta_files, kmer_size = kmer_size, cpu_processes=cpu_processes)
    kmers_intersection = find_dictionaries_intersection(kmers_counts)
    
    # Convert to DataFrame:
    kmer_multiplicities, genomes_lst = createDictMoltep(kmers_counts, kmers_intersection)
    kmer_multiplicities_df = pd.DataFrame.from_dict(kmer_multiplicities, orient='index', columns=genomes_lst)

    return kmer_multiplicities_df


def genome_length_normaization(kmer_multiplicities_df, genomes_length):
    # Create filename -> genome length dictionary
    length_dict = genomes_length.set_index("filename")["length"]

    # Divide each column by its corresponding genome length
    kmers_normalized = kmer_multiplicities_df.div(length_dict, axis="columns")

    return kmers_normalized
    
def coefficient_of_variation(vector):
    mean_val = np.mean(vector)
    if mean_val == 0:
        return 0
    return np.std(vector) / mean_val


#use global variables to pass to multiprocess pooll more efficiently
END_NAMES = None
END_VECTORS = None
WINDOW_SIZE = None

#global variables initializer
def init_pool(end_names, end_vectors, window_size):
    global END_NAMES
    global END_VECTORS
    global WINDOW_SIZE

    END_NAMES = end_names
    END_VECTORS = end_vectors
    WINDOW_SIZE = window_size

def process_start_kmer(args):

    start_name, start_vec = args

    results = []

    for end_name, end_vec in zip(END_NAMES, END_VECTORS):

        # Geometric mean abundance
        joint_vector = np.sqrt(start_vec * end_vec) * 1e6

        # Coefficient of variation
        coeff_var = coefficient_of_variation(joint_vector)

        # Density (f_A*f_B*window_size)
        vec_density = start_vec * end_vec * WINDOW_SIZE

        score = coeff_var * vec_density.min()

        #append results as touples
        results.append( 
            (
                start_name,
                end_name,
                coeff_var,
                score,
            )
        )
    return results


def kmers_processing(input_dataset_path, result_folder_path, kmer_start_size,
                    kmer_end_size, cpu_processes, window_size, sorting="score"): #top_kmers_selection
    
    print("\n[[Step 1 - Search best kmer pairs candidates]]\n")
    
    #list of fasta files in the dataset
    fasta_files = [f for ext in ("*.fna", "*.fa", "*.fasta") for f in input_dataset_path.rglob(ext)]

    #retrieve genomes length
    genomes_length = get_dataset_info(fasta_files)
    genomes_length.to_csv(Path(result_folder_path,"genomes_length.csv"), index=False)

    #retrieve kmer start multiplicity matrix
    kmer_start_multiplicities_df = extract_kmer_multiplicities(fasta_files, kmer_size = kmer_start_size, cpu_processes = cpu_processes)

    #normalize kmer start multiplicitiy
    kmer_start_normalized = genome_length_normaization(kmer_start_multiplicities_df, genomes_length)
    

    #if different kmer sizes
    if kmer_start_size != kmer_end_size:
        print(f"Computing multiplicities for kmer_end of size {kmer_end_size}")
        kmer_end_multiplicities = extract_kmer_multiplicities(fasta_files, 
                                                              kmer_size = kmer_end_size,
                                                              cpu_processes = cpu_processes)
        print("Normalization kmer_end")
        kmer_end_normalized = genome_length_normaization(kmer_end_multiplicities, genomes_length)

        #Save multiplicity matrices
        print("Saving multiplicity matrices")
        kmer_start_multiplicities_df.to_csv(Path(result_folder_path,f"kmer_start_k{kmer_start_size}_multiplicities.csv"), index=False)
        kmer_start_normalized.to_csv(Path(result_folder_path,f"kmer_start_k{kmer_start_size}_multiplicities_normalized.csv"), index=False)

        kmer_end_multiplicities.to_csv(Path(result_folder_path,f"kmer_start_k{kmer_end_size}_multiplicities.csv"), index=False)
        kmer_end_normalized.to_csv(Path(result_folder_path,f"kmer_start_k{kmer_end_size}_multiplicities_normalized.csv"), index=False)
        
    else:
        print(f"\nkmer_star and kmer_end have the same size of {kmer_start_size}")
        kmer_end_normalized = kmer_start_normalized
        
        #Save multiplicity matrices (save only one matrix (same kmers for start and end)
        print("Saving multiplicity matrix")
        kmer_start_multiplicities_df.to_csv(Path(result_folder_path,f"kmer_k{kmer_start_size}_multiplicities.csv"), index=False)
        kmer_start_normalized.to_csv(Path(result_folder_path,f"kmer_k{kmer_start_size}_multiplicities_normalized.csv"), index=False)

    #convert into numpy for memory and speed efficiency
    #kmer names
    start_names = kmer_start_normalized.index.to_numpy()
    end_names = kmer_end_normalized.index.to_numpy()

    #kmer multiplicities normalized
    start_vectors = kmer_start_normalized.to_numpy()
    end_vectors = kmer_end_normalized.to_numpy()
    
    total_pairs = len(start_names) * len(end_names)

    #prepare task for multiprocessing
    tasks = list(zip(start_names, start_vectors))

    pair_results = []

    print(f"\nSearching best kmer pairs - size(start/end): {kmer_start_size}/{kmer_end_size}")
    print(f"[0.0%] 0/{total_pairs}")
    
    start_time = time.time()

    #multiprocess computation of best kmer pairs
    with mp.Pool(
        processes=cpu_processes,
        initializer=init_pool,
        initargs=(end_names, end_vectors, window_size),
    ) as pool:
        
        completed = 0
        total_tasks = len(tasks)  # Ensure tasks has a defined length
        next_threshold = 10       # First milestone to hit
        
        for res in pool.imap(process_start_kmer, tasks, chunksize=1):
            pair_results.extend(res)
            completed += 1
            
            # Calculate current progress percentage
            percent = (completed / total_tasks) * 100
            
            # Check if we've crossed the next 10% threshold
            if percent >= next_threshold:
                done_pairs = completed * len(end_names)
                print(
                    f"[{percent:.1f}%] "
                    f"{done_pairs}/{total_pairs}"
                )
                # Jump to the next 10% mark (handles any potential step jumps cleanly)
                next_threshold = ((int(percent) // 10) + 1) * 10
    
    print(f"Kmer pairs scoring - Elapsed Time: {format_time(time.time() - start_time)}")
    print(f"Saving kmer pairs scores...", end="")
    
    #create dataframe from touples
    pairs_df = pd.DataFrame(
        pair_results,
        columns=["kmer_start", "kmer_end", "CoV", "score"],
    )

    #sorting by total score: CV*minimum_density
    pairs_df = pairs_df.sort_values(by='score', ascending=False)

    pairs_df.to_csv(Path(result_folder_path,"kmer_pairs.csv"), index=False)
    print("done!")

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
    #print(kmer)
    #print(type(kmer))
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

def aflp_processing(input_dataset_path, result_folder_path, top_kmers_selection, cpu_processes, fragment_length):
    
    """
    Second step wrapper, processing of the AFLP distributions
    """
    global mindist
    global maxdist
    global kmer_pairs

    mindist = 0
    maxdist = fragment_length
    
    #kmer_pairs = list(pd.read_csv(Path(result_folder_path,'kmer_pairs.csv')).itertuples(index=False, name=None))
    kmer_pairs = pd.read_csv(Path(result_folder_path,'kmer_pairs.csv'))
    
    #subsetting the top kmer pairs (chosen by input the top kmer_pairs to select)
    kmer_pairs = [tuple(x) for x in kmer_pairs.iloc[:top_kmers_selection, 0:2].values]
    
    #sorting fasta files by size for multiprocess
    genome_files = get_fasta_files_sorted_by_size(input_dataset_path)

    output_folder = Path(result_folder_path, "aflp") 
    output_folder.mkdir(parents=True, exist_ok=True)

    print("\n[[Step 2 - Extract AFLP distributions]]\n")
        
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

def summarize_aflp_validity(folder_path: str | Path):
  path = Path(folder_path)

  if not path.exists() or not path.is_dir():
    raise ValueError(f"The provided path '{path}' is not a valid directory.")

  summary_rows = []

  for csv_file in path.glob("*.csv"):
    try:
      df = pd.read_csv(csv_file)

      # Handle empty dataframes gracefully
      if df.empty:
        valid_count = 0
        invalid_count = 0
      else:
        subset_df = df.iloc[:, 5:]
          
        # Check if all values in the selected columns are 0 for each row
        is_all_zero = (subset_df == 0).all(axis=1)

        invalid_count = int(is_all_zero.sum())
        valid_count = int((~is_all_zero).sum())

      summary_rows.append({
          "file_name": csv_file.name,
          "valid": valid_count,
          "invalid": invalid_count,
      })

    except Exception as e:
      print(f"Error processing {csv_file.name}: {e}")

  # Create the final summary DataFrame
  summary_df = pd.DataFrame(summary_rows, columns=["file_name", "valid", "invalid"])
  return summary_df

