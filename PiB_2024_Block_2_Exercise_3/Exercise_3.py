#!/usr/bin/env python
# coding: utf-8

# # Clean up and package your code and results
# 
# Exercise for applying some good (Python) practices to your code.
# 
# * **Contact:** mate.balajti@unibas.ch

# ## Part 1: Testing and additional documentation
# 
# Pointers for all exercises in this part are available in the Jupyter notebook
# `good_practices.ipynb`.
# 
# Maintainable code includes documentation and testing. For each of these
# subtasks, **please do not forget to update dependencies (e.g., `flake8`,
# `pytest`).**
# 
# > Note: Please use your code from Exercise 1, and write it into a separate python file.

# ### Exercise 3.1: Add type hints to your code (1 point)
# 
# Add type hints to your custom Python function/method signatures. It will be
# enough to only add type hints for all input arguments and the return values,
# although you are of course welcome to add them for any local variables as well.

# In[ ]:


def parse_fasta(path: str) -> tuple[list[str], list[str]]:
    headers = []
    sequences = []
    current_sequence = ""

    with open(path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    sequences.append(current_sequence)
                    current_sequence = ""
                headers.append(line[1:])
            else:
                current_sequence += line

        if current_sequence:
            sequences.append(current_sequence)


    return headers, sequences

def discard_ambiguous_seqs(header: list[str], sequence: list[str]) -> \
  tuple[: list[str], : list[str]]:

    l1_header = []
    l2_sequence = []
    DNA_alphabet = {"A", "G", "C", "T"}

    for n, seq in zip(header, sequence):

        if all(char.upper() in DNA_alphabet for char in seq):
            l1_header.append(n)
            l2_sequence.append(seq)
    

    return l1_header, l2_sequence

def nucleotide_frequencies(seqs: list[str]) -> None:
    # count nucleotides
    nucleotides_count = {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}
    total = 0

    for seq in seqs:
        for nucleotide in seq.upper():
            if nucleotide in nucleotides_count:
                nucleotides_count[nucleotide] += 1
                total += 1

    if total > 0:
        for nucleotide in nucleotides_count:
            frequency = nucleotides_count[nucleotide] / total
            print(f"{nucleotide} : {frequency: .2f}")

    else:
        print("No nucleotides")


def map_reads(filename1: str, filename2: str) -> dict[dict[str: list[int]]]:
    query_headers, query_sequences = parse_fasta(filename1)
    reference_headers, reference_sequences = parse_fasta(filename2)

    query_headers, query_sequences = discard_ambiguous_seqs(query_headers, query_sequences)

    print("Query nucleotides frequency:")
    nucleotide_frequencies(query_sequences)

    print("Reference nucleotides frequency:")
    nucleotide_frequencies(reference_sequences)

    results = {header: {} for header in query_headers}

    for query_header, query_seq in zip (query_headers, query_sequences):
        for ref_header, ref_seq in zip(reference_headers, reference_sequences):

            start = 0
            while True:
                start = ref_seq.find(query_seq, start)
                if start == -1:
                    break
                if ref_header not in results[query_header]:
                    results[query_header][ref_header] = []
                results[query_header][ref_header].append(start + 1)
                start += 1

    return results

filename1 = "PiB_2024_Block_2_Exercise_1/sequences.fasta"
filename2 = "PiB_2024_Block_2_Exercise_1/genome.fasta"

hits = map_reads(filename1, filename2)

print("Hits :")
for query, refs in hits.items():
    print(f"{query}: {refs}")


# ### Exercise 3.2: Add docstrings to your code (1 point)
# 
# Add Google-style docstrings to your custom Python functions/methods. Please do
# not include argument and return value types in the docstrings, as these have
# been already added to the signatures themselves (which is a much better idea,
# because the actual code is always the source of truth, and thus the risk of
# code and documentation diverging over time is reduced).

# In[ ]:


def parse_fasta(path: str) -> tuple[list[str], list[str]]:
    """
    Read a FASTA file and extracts sequence headers and sequences.

    Args:
        path: Path to the FASTA file.

    Returns:
        A tuple with two lists:
        - headers: List of sequence headers.
        - sequences: List of corresponding sequences.
    """
    headers = []
    sequences = []
    current_sequence = ""

    with open(path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    sequences.append(current_sequence)
                    current_sequence = ""
                headers.append(line[1:])
            else:
                current_sequence += line

        if current_sequence:
            sequences.append(current_sequence)

    return headers, sequences


def discard_ambiguous_seqs(
    header: list[str], sequence: list[str]
) -> tuple[list[str], list[str]]:
    """
    Removes sequences containing ambiguous characters from the input.

    Args:
        header: List of sequence headers.
        sequence: List of nucleotide sequences.

    Returns:
        A tuple containing:
        - l1_header: Filtered list of sequence headers.
        - l2_sequence: Filtered list of sequences without ambiguous characters.
    """
    l1_header = []
    l2_sequence = []
    DNA_alphabet = {"A", "G", "C", "T"}

    for n, seq in zip(header, sequence):
        if all(char.upper() in DNA_alphabet for char in seq):
            l1_header.append(n)
            l2_sequence.append(seq)

    return l1_header, l2_sequence


def nucleotide_frequencies(seqs: list[str]) -> None:
    """
    Calculates and prints the frequency of nucleotides in the input sequences.

    Args:
        seqs: List of nucleotide sequences.

    Prints:
        The relative frequency of each nucleotide (A, C, G, T) in the sequences.
    """
    nucleotides_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    total = 0

    for seq in seqs:
        for nucleotide in seq.upper():
            if nucleotide in nucleotides_count:
                nucleotides_count[nucleotide] += 1
                total += 1

    if total > 0:
        for nucleotide in nucleotides_count:
            frequency = nucleotides_count[nucleotide] / total
            print(f"{nucleotide} : {frequency: .2f}")
    else:
        print("No nucleotides")


def map_reads(filename1: str, filename2: str) -> dict[str, dict[str, list[int]]]:
    """
    Maps query sequences to reference sequences and finds matching positions.

    Args:
        filename1: Path to the query FASTA file.
        filename2: Path to the reference FASTA file.

    Returns:
        A nested dictionary where:
        - Keys are query sequence headers.
        - Values are dictionaries where keys are reference sequence headers, and values are lists of starting positions where the query sequence matches the reference.
    """
    query_headers, query_sequences = parse_fasta(filename1)
    reference_headers, reference_sequences = parse_fasta(filename2)

    query_headers, query_sequences = discard_ambiguous_seqs(query_headers, query_sequences)

    print("Query nucleotides frequency:")
    nucleotide_frequencies(query_sequences)

    print("Reference nucleotides frequency:")
    nucleotide_frequencies(reference_sequences)

    results = {header: {} for header in query_headers}

    for query_header, query_seq in zip(query_headers, query_sequences):
        for ref_header, ref_seq in zip(reference_headers, reference_sequences):
            start = 0
            while True:
                start = ref_seq.find(query_seq, start)
                if start == -1:
                    break
                if ref_header not in results[query_header]:
                    results[query_header][ref_header] = []
                results[query_header][ref_header].append(start + 1)
                start += 1

    return results


# Main script
filename1 = "PiB_2024_Block_2_Exercise_1/sequences.fasta"
filename2 = "PiB_2024_Block_2_Exercise_1/genome.fasta"

hits = map_reads(filename1, filename2)

print("Hits:")
for query, refs in hits.items():
    print(f"{query}: {refs}")


# ### Exercise 3.3: Make sure your code lints (1 point)
# 
# Please use the `flake8` linter to help you refactor your code such that it
# adheres to Python conventions.

# In[1]:


pip install flake8

flake8 cd /mnt/c/Users/nicol/OneDrive/Desktop/UniBasel/Computational\ science/Prog\ in\ Bioinfo/Zavolan/PiB_2024_Block_2_Exercise_1/Exercise_1.py

# ### Exercise 3.4: Test your code (2 points)
# 
# Write unit tests for all of your custom Python functions/methods, run tests
# with `pytest` and make sure all tests pass. Compute the code `coverage` and
# make sure it's at a 100%.

# In[ ]:


pytest  /mnt/c/Users/nicol/OneDrive/Desktop/UniBasel/Computational\ science/Pr
og\ in\ Bioinfo/Zavolan/PiB_2024_Block_2_Exercise_1/test_exercise_1.py


# ## Part 2: Package and version control your code
# 
# ### Reorganize your files
# 
# Please create an empty directory with the following files and directories:
# 
# * `README.md`: A [markdown](https://github.github.com/gfm/)-formatted file
#   containing instructions on how to deploy (e.g., set up a Conda environment
#   or build a Docker image; see below) and run your code and answering all
#   questions from previous exercises.
# * `src/`: A directory containing all your custom code (except for tests).
# * `tests/`: A directory that will contain test code (see next exercise). You
#   can also put all your test/input files in this directory, ideally in a
#   subdirectoy (e.g., `test_files/`).
# * `run_me.sh`: A single Bash script running all code from all exercises on
#   the test input files as per the previous two exercises. Please also include
#   the calls to run `flake8` and `pytest` on your core repository after you have
#   implemented these. This is really a poor man's workflow, but it's still a lot
#   better than having to guess how you ran your code. If you like, you can
#   organize this file so that it executes other Bash scripts (in `src/`) to keep
#   it more tidy.
# * [OPTIONAL] `LICENSE`: A file containing an
#   [Open Source License](https://opensource.org/licenses), such as the
#   [MIT](https://opensource.org/licenses/MIT) or
#   [Apache 2.0](https://opensource.org/licenses/Apache-2.0) license; you can use
#   [this service](https://choosealicense.com/) to pick a license you like.
# * [OPTIONAL] `environment.yml` OR `Dockerfile`: A [Conda enviroment
#   file](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#sharing-an-environment)
#   listing all your Python and third-party requirements. Alternatively, a
#   `Dockerfile` that can be used to build an image that contains all
#   dependencies.
# * [OPTIONAL] `.gitignore`: A file with patterns indicating files/artefacts that Git should
#   not version control, e.g., test reports. Generate it from
#   <https://gitignore.io> for Python and your code editor.
# 
# ### Create a Git repository and push to a remote
# 
# * Please create a Git repository based off your new directory.
# * Stage and commit all contents of the directory, with e.g. the message "initial commit".
# * Create a blank/empty project on [GitHub](https://github.com/) or
#   [GitLab](https://gitlab.com/).
# * Add the remote URL of your GitHub/Lab project to your local Git repository
#   and push your code.

# ## Part 3: Use a workflow language/engine to run your analysis 
# 
# As mentioned above, specifiying a Bash file to run your analysis is not very
# good. It will be difficult to parallelize your code, scatter/gather multiple
# jobs of the same task, keep sufficient logging and provenance information, and
# it will be difficult to share your analysis in a way that it is easily
# reproducible/reusable. Workflow languages and corresponding management
# systems/engines take care of all of these things and more.
# 
# If you are (planning on) doing bioinformatics analyses more regularly, we
# strongly recommend you to pick up one of these languages, e.g.,
# [Nextflow](https://www.nextflow.io/) (Groovy-based) or
# [Snakemake](https://snakemake.readthedocs.io/en/stable/) (Python-based) are two
# popular choices that we frequently use in our lab.
# 
# If you are interested, you can follow a tutorial on either of these domain-
# specific languages and learn how you can package your code as a proper
# shareable workflow.
# 
# > Note: Please use your code from Exercise 2, organize them into 2 Nextflow processes and copy the code below.

# In[ ]:




