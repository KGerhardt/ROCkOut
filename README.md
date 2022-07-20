# ROCkOut

Work in progress!

Python refactor of ROCker - https://github.com/lmrodriguezr/rocker

# Rethinking the ROCker Approach:

The original ROCker's development was done in a piecemeal fashion. While the resulting tool is functional, it's neither particularly user friendly, nor is it well-designed with respect to the inherently iterative process of building a complete ROCker model. Building a ROCker model requires repeatedly testing candidate proteins and including or excluding them as they enhance or degrade the performance of the final result, yet the program pools data in a manner that requires completely restarting the program's model building process to make any updates at all. The only time that the original ROCker is fully efficient is when you already have a completed model's protein set(s) and simply want to reproduce the final model - a pointless task for someone who wishes to *use* a model, since the final model in that process can simply be downloaded and used without being rebuilt.

Consequently, ROCkOut is not intended to perfectly replicate the original ROCker, but to enhance it. Pooled data and multiple-step scripts are discretized in ROCkOut so that the addition of a protein for testing only requires running building steps for that protein alone, and excluding a protein that is harmful is as simple as ignoring that protein's directory from a subsequent step.

Additionally, ROCkOut is designed to support the user through visual aids and supports that make the building process friendly, rather than combatitive.

# Progress report:

0_download_from_uniprot - A script which downloads proteins from UniProt using their IDs and their associated data in a similar manner to the original ROCker. Initializes a project directory and organizes data. Data downloaded in this step include an annotation of each UniProt protein, a FASTA of the genome containing it, related protein sets as GFFs, and the coordinates and sequence of the target protein for each item in the positive set.

1_generate_reads - From a project directory, simulates reads using BBTools RandomReads. Reads are simulated from the genomes created in step 0 and/or from genomes manually added by a user.

2_tag_reads - Tags simulated reads for ROCkOut. Includes information on the origin genome of each read and their true alignment.

# Plans and ideas:

3_align_to_refs - Builds a database from the current set of postive proteins and aligns tagged reads to the database. Not quite done.

* There needs to be a separate project initialization script for someone who wants to build a model from only manually downloaded genomes.

* Next steps are to finish the alignment component and move on to scripting the ROC results for the aligned reads.

* Visual and building supports - a GUI made to allow for the interactive selection of proteins as in/out and show the resulting effect on the model would be good.

* Model building - The use of a sliding window for determining the best cutoffs for a model make sense, but having to rerun that whole step to see a new window doesn't. The data should just all be loaded and sliced as needed on the interactive page.

# Usage

Scripts are used in their number order. You need a list of positive uniprot IDs, one per line, in a positives.txt file, negatives.txt if you want them. Usage is as follows:

python 0_download_from_uniprot.py [positives.txt] [negatives.txt] [threads] [project_dir_name]
python 1_generate_reads.py [project_dir_name] [threads]
python 2_tag_reads.py [project_dir_name] [threads]
python 3_align_to_refs.py [project_dir_name] [threads]

Open 4_refiner.R in RStudio and run all the code.

[project_dir_name] does not need to exist prior to running step 0.


# Installation

conda install -c bioconda -c conda-forge python=3.7 bbmap=38.93 muscle=3.8.31 diamond numpy
pip install requests
