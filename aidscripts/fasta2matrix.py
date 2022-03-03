#!/usr/bin/env python
# FILE: fasta2matrix.py
# AUTHOR: William Stafford Noble
# CREATE DATE: 27 September 2005
import sys
import math
import numpy as np
#############################################################################
def make_kmer_list (k, alphabet):

    # Base case.
    if (k == 1):
        return(alphabet)

    # Handle k=0 from user.
    if (k == 0):
        return([])

    # Error case.
    if (k < 1):
        sys.stderr.write("Invalid k=%d" % k)
        sys.exit(1)

    # Precompute alphabet length for speed.
    alphabet_length = len(alphabet)

    # Recursive call.
    return_value = []
    for kmer in make_kmer_list(k-1, alphabet):
        for i_letter in range(0, alphabet_length):
            return_value.append(kmer + alphabet[i_letter])
              
    return(return_value)

##############################################################################
def make_upto_kmer_list (k_values,
                         alphabet):

    # Compute the k-mer for each value of k.
    return_value = []
    for k in k_values:
        return_value.extend(make_kmer_list(k, alphabet))

    return(return_value)
                 
##############################################################################
def normalize_vector (normalize_method,
                      k_values,
                      vector,
                      kmer_list):

    # Do nothing if there's no normalization.
    if (normalize_method == "none"):
        return(vector)

    # Initialize all vector lengths to zeroes.
    vector_lengths = {}
    for k in k_values:
        vector_lengths[k] = 0

    # Compute sum or sum-of-squares separately for each k.
    num_kmers = len(kmer_list)
    for i_kmer in range(0, num_kmers):
        kmer_length = len(kmer_list[i_kmer])
        count = vector[i_kmer]
        if (normalize_method == "frequency"):
            vector_lengths[kmer_length] += count
        elif (normalize_method == "unitsphere"):
            vector_lengths[kmer_length] += count * count

    # Square root each sum if we're doing 2-norm.
    if (normalize_method == "unitsphere"):
        for k in k_values:
            vector_lengths[k] = math.sqrt(vector_lengths[k])

    # Divide through by each sum.
    return_value = []
    for i_kmer in range(0, num_kmers):
        kmer_length = len(kmer_list[i_kmer])
        count = vector[i_kmer]
        vector_length = vector_lengths[kmer_length]
        if (vector_length == 0):
            return_value.append(0)
        else:
            return_value.append(float(count) / float(vector_length))

    return(return_value)

##############################################################################
# Make a copy of a given string, substituting one letter.
def substitute (position,
                letter,
                string):

    return_value = ""
    if (position > 0):
        return_value = return_value + string[0:position]
    return_value = return_value + letter
    if (position < (len(string) - 1)):
        return_value = return_value + string[position+1:]
                   
    return(return_value)


##############################################################################
def compute_bin_num (num_bins,
                     position,
                     k,
                     numbers):
  
  # If no binning, just return.
  if (num_bins == 1):
    return(0)

  # Compute the mean value for this k-mer.
  mean = 0
  for i in range(0, k):
    mean += float(numbers[position + i])
  mean /= k

  # Find what quantile it lies in.
  for i_bin in range(0, num_bins):
    if (mean <= boundaries[k][i_bin]):
      break

  # Make sure we didn't go too far.
  if (i_bin == num_bins):
      sys.stderr.write("bin=num_bins=%d\n", i_bin)
      sys.exit(1);

  return(i_bin)

##############################################################################
def make_sequence_vector (sequence,
                          numbers,
                          num_bins,
                          revcomp,
                          revcomp_dictionary,
                          normalize_method,
                          k_values,
                          mismatch,
                          alphabet,
                          kmer_list,
                          boundaries,
                          pseudocount):

    # Make an empty counts vector.
    kmer_counts = []
    for i_bin in range(0, num_bins):
        kmer_counts.append({})

    # Iterate along the sequence.
    for k in k_values:
        seq_length = len(sequence) - k + 1
        for i_seq in range(0, seq_length):

            # Compute which bin number this goes in.
            bin_num = compute_bin_num(num_bins, i_seq, k, numbers)

            # Extract this k-mer.
            kmer = sequence[i_seq : i_seq + k]

            # If we're doing reverse complement, store the count in the
            # the version that starts with A or C.
            if (revcomp == 1):
                rev_kmer = find_revcomp(kmer, revcomp_dictionary)
                if (cmp(kmer, rev_kmer) > 0):
                    kmer = rev_kmer

            # Increment the count.
            #if (kmer_counts[bin_num].has_key(kmer)):
            if (kmer in kmer_counts[bin_num].keys()):				
                kmer_counts[bin_num][kmer] += 1
            else:
                kmer_counts[bin_num][kmer] = 1

            # Should we also do the mismatch?
            if (mismatch != 0):

                # Loop through all possible mismatches.
                for i_kmer in range(0, k):
                    for letter in alphabet:

                        # Don't count yourself as a mismatch.
                        if (kmer[i_kmer:i_kmer+1] != letter):

                            # Find the neighboring sequence.
                            neighbor = substitute(i_kmer, letter, kmer)

                            # If we're doing reverse complement, store the
                            # count in the version that starts with A or C.
                            if (revcomp == 1):
                                rev_kmer = find_revcomp(kmer,
                                                        revcomp_dictionary)
                                if (cmp(kmer, rev_kmer) > 0):
                                    kmer = rev_kmer

                            # Increment the neighboring sequence.
                            if (kmer_counts[bin_num].has_key(neighbor)):
                                kmer_counts[bin_num][neighbor] += mismatch
                            else:
                                kmer_counts[bin_num][neighbor] = mismatch

    # Build the sequence vector.
    sequence_vector = []
    for i_bin in range(0, num_bins):
        for kmer in kmer_list:
            #if (kmer_counts[i_bin].has_key(kmer)):
            if (kmer in kmer_counts[i_bin].keys()):				
                sequence_vector.append(kmer_counts[i_bin][kmer] + pseudocount)
            else:
                sequence_vector.append(pseudocount)

    # Normalize it
    return_value = normalize_vector(normalize_method,
                                    k_values,
                                    sequence_vector,
                                    kmer_list)

    return(return_value)

##############################################################################
def read_fasta_sequence (numeric,
                         fasta_file):

    # Read 1 byte.  
    first_char = fasta_file.read(1)
    # If it's empty, we're done.
    if (first_char == ""):
        return(["", ""])
    # If it's ">", then this is the first sequence in the file.
    elif (first_char == ">"):
        line = ""
    else:
        line = first_char

    # Read the rest of the header line.
    line = line + fasta_file.readline()
    # Get the rest of the ID.
    words = line.split()
    if (len(words) == 0):
      sys.stderr.write("No words in header line (%s)\n" % line)
      sys.exit(1)  
    id = words[0]
    #id = words[1].split("=")[1]
        
    # Read the sequence, through the next ">".
    first_char = fasta_file.read(1)
    sequence = ""
    while ((first_char != ">") and (first_char != "")):
        if (first_char != "\n"): # Handle blank lines.
            line = fasta_file.readline()
            sequence = sequence + first_char + line
        first_char = fasta_file.read(1)

    # Remove EOLs.
    clean_sequence = ""
    for letter in sequence:
        if (letter != "\n"):
            clean_sequence = clean_sequence + letter
    sequence = clean_sequence

    # Remove spaces.
    if (numeric == 0):
        clean_sequence = ""
        for letter in sequence:
            if (letter != " "):
                clean_sequence = clean_sequence + letter
        sequence = clean_sequence.upper()
        
    return([id, sequence])

##############################################################################
def read_sequence_and_numbers(fasta_file,
                              numbers_filename,
                              numbers_file):
    [fasta_id, fasta_sequence] = read_fasta_sequence(0, fasta_file)
    if (number_filename != ""):
        [number_id, number_sequence] = read_fasta_sequence(1, number_file)

        # Make sure we got the same ID.
        if (fasta_id != number_id):
            sys.stderr.write("Found mismatching IDs (%s != %d)\n" %
                             (fasta_id, number_id))
            sys.exit(1)

        # Split the numbers into a list.
        number_list = number_sequence.split()

        # Verify that they are the same length.
        if (len(fasta_sequence) != len(number_list)):
            sys.stderr.write("Found sequence of length %d with %d numbers.\n"
                             % (len(sequence), len(number_list)))
            print(sequence)
            print(numbers)
            sys.exit(1)
    else:
        number_list = ""

    return(fasta_id, fasta_sequence, number_list)
                                                      
    
##############################################################################
def find_revcomp (sequence,
                  revcomp_dictionary):

    # Save time by storing reverse complements in a hash.
    if (revcomp_dictionary.has_key(sequence)):
        return(revcomp_dictionary[sequence])

    # Make a reversed version of the string.
    rev_sequence = list(sequence)
    rev_sequence.reverse()
    rev_sequence = ''.join(rev_sequence)
    
    return_value = ""
    for letter in rev_sequence:
        if (letter == "A"):
            return_value = return_value + "T"
        elif (letter == "C"):
            return_value = return_value + "G"
        elif (letter == "G"):
            return_value = return_value + "C"
        elif (letter == "T"):
            return_value = return_value + "A"
        elif (letter == "N"):
            return_value = return_value + "N"
        else:
            sys.stderr.write("Unknown DNA character (%s)\n" % letter)
            sys.exit(1)

    # Store this value for future use.
    revcomp_dictionary[sequence] = return_value

    return(return_value)

##############################################################################
def compute_quantile_boundaries (num_bins,
                                 k_values,
                                 number_filename):

    if (num_bins == 1):
        return

    # The boundaries are stored in a 2D dictionary.
    boundaries = {}

    # Enumerate all values of k.
    for k in k_values:

        # Open the number file for reading.
        number_file = open(number_filename, "r")
        
        # Read it sequence by sequence.
        all_numbers = []
        [id, numbers] = read_fasta_sequence(1, number_file)
        while (id != ""):
            
            # Compute and store the mean of all k-mers.
            number_list = numbers.split()
            num_numbers = len(number_list) - k
            for i_number in range(0, num_numbers):
                if (i_number == 0):
                    sum = 0;
                    for i in range(0, k):
                        sum += float(number_list[i])
                else:
                    sum -= float(number_list[i_number - 1])
                    sum += float(number_list[i_number + k - 1])
                all_numbers.append(sum / k)
            [id, numbers] = read_fasta_sequence(1, number_file)
        number_file.close()

        # Sort them.
        all_numbers.sort()

        # Record the quantiles.
        boundaries[k] = {}
        num_values = len(all_numbers)
        bin_size = float(num_values) / float(num_bins)
        sys.stderr.write("boundaries k=%d:" % k)
        for i_bin in range(0, num_bins):
            value_index = int((bin_size * (i_bin + 1)) - 1)
            if (value_index == num_bins - 1):
                value_index = num_values - 1
            value = all_numbers[value_index]
            boundaries[k][i_bin] = value
            sys.stderr.write(" %g" % boundaries[k][i_bin])
        sys.stderr.write("\n")

    return(boundaries)
    

##############################################################################
# MAIN
##############################################################################

# Define the command line usage.
usage = """Usage: fasta2matrix [options] <k> <fasta file>

  Options:

    -upto       Use all values from 1 up to the specified k.

    -revcomp    Collapse reverse complement counts.

    -normalize [frequency|unitsphere] Normalize counts to be 
                frequencies or project onto unit sphere.  With -upto,
                normalization is done separately for each k.

    -protein    Use an amino acid alphabet.  Default=ACGT.

    -alphabet <string> Set the alphabet arbitrarily.

    -mismatch <value>  Assign count of <value> to k-mers that 
                       are 1 mismatch away.
    
    -binned <numbins> <file>  Create <numbins> vectors for each
                              sequence, and place each k-mer count
                              into the bin based upon its corresponding
                              mean value from the <file>.  The
                              <file> is in FASTA-like format, with
                              space-delimited numbers in place of
                              the sequences.  The sequences must
                              have the same names and be in the same
                              order as the given FASTA file.

   -pseudocount <value>  Assign the given pseudocount to each bin.

"""

# Define default options.
upto = 0
revcomp = 0
normalize_method = "none"
alphabet = "ACGT"
mismatch = 0
num_bins = 1
pseudocount = 0
number_filename = ""

# Parse the command line.
# sys.argv = sys.argv[1:]
# while (len(sys.argv) > 2):
#   next_arg = sys.argv[0]
#   sys.argv = sys.argv[1:]
#   if (next_arg == "-revcomp"):
#     revcomp = 1
#   elif (next_arg == "-upto"):
#     upto = 1
#   elif (next_arg == "-normalize"):
#     normalize_method = sys.argv[0]
#     sys.argv = sys.argv[1:]
#     if ((normalize_method != "unitsphere") and
#       (normalize_method != "frequency")):
#       sys.stderr.write("Invalid normalization method (%s).\n"
#                        % normalize_method);
#       sys.exit(1)
#   elif (next_arg == "-protein"):
#     alphabet = "ACDEFGHIKLMNPQRSTVWY"
#   elif (next_arg == "-alphabet"):
#     alphabet = sys.argv[0]
#     sys.argv = sys.argv[1:]
#   elif (next_arg == "-mismatch"):
#     mismatch = float(sys.argv[0])
#     sys.argv = sys.argv[1:]
#   elif (next_arg == "-binned"):
#     num_bins = int(sys.argv[0])
#     sys.argv = sys.argv[1:]
#     number_filename = sys.argv[0]
#     sys.argv = sys.argv[1:]
#   elif (next_arg == "-pseudocount"):
#     pseudocount = int(sys.argv[0])
#     sys.argv = sys.argv[1:]
#   else:
#     sys.stderr.write("Invalid option (%s)\n" % next_arg)
#     sys.exit(1)
# if (len(sys.argv) != 2):
#   sys.stderr.write(usage)
#   sys.exit(1)
k = int(sys.argv[1])
fasta_filename = sys.argv[2]
output_filename = sys.argv[3]
# Check for reverse complementing non-DNA alphabets.
if ((revcomp == 1) and (alphabet != "ACGT")):
  sys.stderr.write("Attempted to reverse complement ")
  sys.stderr.write("a non-DNA alphabet (%s)\n" % alphabet)

# Make a list of all values of k.
k_values = []
if (upto == 1):
  start_i_k = 1
else:
  start_i_k = k
k_values = range(start_i_k, k+1)

# If numeric binning is turned on, compute quantile boundaries for various
# values of k.
boundaries = compute_quantile_boundaries(num_bins, k_values, number_filename)
  
# Make a list of all k-mers.
kmer_list = make_upto_kmer_list(k_values, alphabet);
sys.stdout.write("Consideriamo %d kmers.\n" % len(kmer_list))

# Set up a dictionary to cache reverse complements.
revcomp_dictionary = {}

# Use lexicographically first version of {kmer, revcomp(kmer)}.
if (revcomp == 1):
  new_kmer_list = []
  for kmer in kmer_list:
      rev_kmer = find_revcomp(kmer, revcomp_dictionary)
      if (cmp(kmer, rev_kmer) <= 0):
          new_kmer_list.append(kmer)
  kmer_list = new_kmer_list;
  sys.stdout.write("Reduced to %d kmers.\n" % len(kmer_list))

# Print the corner of the matrix.

outfile=open(output_filename,"w")

# Print the title row.


sys.stdout.write("\n")
outfile.write("seq_id,")

# Open the sequence file.
if (fasta_filename == "-"):
  fasta_file = sys.stdin
else:
  fasta_file = open(fasta_filename, "r")
if (number_filename == ""):
  number_file = 0
else:
  number_file = open(number_filename, "r")
#for i_bin in range(1, num_bins+1):
for kmer in kmer_list:
  #if (num_bins > 1):
    #outfile.write("%s-%d," % (kmer, i_bin))
      #i+=1
  #else:
  if(kmer==kmer_list[len(kmer_list)-1]):
      outfile.write("%s" % kmer)
  else:
	  outfile.write("%s," %kmer)
      #i+=1

outfile.write("\n")
# Read the first sequence.
[id, sequence, numbers] = read_sequence_and_numbers(fasta_file,
                                                    number_filename,
                                                    number_file)
# Iterate till we've read the whole file.
i_sequence = 1
vett=np.zeros(len(kmer_list),dtype=int)
while (id != ""):

  # Tell the user what's happening.
  if (i_sequence % 1000 == 0):
    sys.stdout.write("Reading %dth sequenza.\n" % i_sequence)

  # Compute the sequence vector.
  vector = make_sequence_vector(sequence,
                                numbers,
                                num_bins,
                                revcomp,
                                revcomp_dictionary,
                                normalize_method,
                                k_values,
                                mismatch,
                                alphabet,
                                kmer_list,
                                boundaries,
                                pseudocount)

  # Print the formatted vector.
  outfile.write("%s," % id)
  
  count=len(kmer_list)

  for element in vector:
    if(count!=1):
        outfile.write("%d," % element)
    else:
        outfile.write("%d" % element)
    count = count-1
        
  outfile.write("\n")
  # Read the next sequence.
  [id, sequence, numbers] = read_sequence_and_numbers(fasta_file,
                                                      number_filename,
                                                      number_file)
  i_sequence += 1
# Close the file.
  i=0


outfile.close()
fasta_file.close()





