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

# Installation

Install with conda or install individual dependencies separately and make them available in the system PATH.

Conda Installation

```bash
conda create -n ROCkOut -c bioconda -c conda-forge python=3.7 bbmap=38.93 muscle=3.8.31 diamond numpy requests dash
```

# Usage

If using conda, activate the environment

```bash
conda activate ROCkOut
```

Before building the model, ROCkIn is used to develop a list of positive sequences (and their UniProt IDs) and an optional list of negative sequences. Data from the ROCkIn pipeline is used to assist the user in selecting these sequences. See the separate README in the ROCkIn directory for more details.

Manually create a simple text file with one UniProt ID per line. Create one file for the positive IDs (POS) and optionally negative IDs (NEG). Do no leave any blank lines at the end or beginning of the these files.

ROCkOut usage is as follows:

For additional details for each script type: python scriptname.py -h

*OUT=user specified project directory. Create any name you want for the first script and then use this same directory for the rest of the pipeline.*

```bash
python rocker_0_download_from_uniprot.py [-h] [-p POS] [-n NEG] [-t THREADS] [-o OUT]
python rocker_1_generate_reads.py [-h] [-t THREADS] [-o OUT]
python rocker_2_tag_reads.py [-h] [-t THREADS] [-o OUT]
python rocker_3_align_to_refs.py [-h] [-t THREADS] [-o OUT]
```

Then run the rocker_dash.py script and open the IP address it gives in your web browser (copy and paste). Drag and drop or select the ROCkOUT_index.txt file that was created in the [OUT] directory.

```bash
python rocker_dash.py
```


