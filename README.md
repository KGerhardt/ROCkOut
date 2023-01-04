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

# Dependencies

- [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
- [Muscle](https://github.com/rcedgar/muscle)
- [Mafft](https://mafft.cbrc.jp/alignment/software/)
- [Diamond](https://github.com/bbuchfink/diamond)
- [Pplacer](https://github.com/matsen/pplacer)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [Python](https://www.python.org/) 3.7+

#### References

 1. Bushnell B. BBMap: a fast, accurate, splice-aware aligner. Lawrence Berkeley National Lab.(LBNL), Berkeley, CA (United States); 2014 Mar 17.
 1. Edgar RC. MUSCLE v5 enables improved estimates of phylogenetic tree confidence by ensemble bootstrapping. BioRxiv. 2021 Jan 1.
 1. Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular biology and evolution. 2013 Jan 16;30(4):772-80.
 1. Li C, Gao M, Yang W, Zhong C, Yu R. Diamond: a multi-modal DIA mass spectrometry data processing pipeline. Bioinformatics. 2021 Jan 15;37(2):265-7.
 1. Matsen FA, Kodner RB, Armbrust E. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. BMC bioinformatics. 2010 Dec;11(1):1-6.
 1. Price, M.N., Dehal, P.S., and Arkin, A.P. (2010) FastTree 2 -- Approximately Maximum-Likelihood Trees for Large Alignments. PLoS ONE, 5(3):e9490. doi:10.1371/journal.pone.0009490.
 1. Sanner MF. Python: a programming language for software integration and development. J Mol Graph Model. 1999 Feb 1;17(1):57-61.

#### Python Packages

- [pandas](https://pandas.pydata.org/) 
- [numpy](https://numpy.org/)
- [requests](https://requests.readthedocs.io/en/latest/)
- [dash](https://dash.plotly.com/introduction)
- [ete3](http://etetoolkit.org/)
- [taxtastic](https://github.com/fhcrc/taxtastic)
- [pyqt](https://www.riverbankcomputing.com/software/pyqt/)

#### References

1. McKinney W, others. Data structures for statistical computing in python. In: Proceedings of the 9th Python in Science Conference. 2010. p. 51–6.
1. Harris CR, Millman KJ, van der Walt SJ, Gommers R, Virtanen P, Cournapeau D, et al. Array programming with NumPy. Nature. 2020;585:357–62.
1. Reitz K, Cordasco I, Prewitt N. Requests: HTTP for humans. KennethReitz [Internet]. https://requests.readthedocs.io/en/latest/. 2023.
1. https://dash.plotly.com/introduction
1. Huerta-Cepas J, Serra F, Bork P. ETE 3: reconstruction, analysis, and visualization of phylogenomic data. Molecular biology and evolution. 2016 Feb 26;33(6):1635-8.
1. https://github.com/fhcrc/taxtastic
1. https://www.riverbankcomputing.com/software/pyqt/

#### Conda Installation

*All dependecies can be installed using conda*

```bash
conda create -p /storage/home/hcoda1/9/rconrad6/p-ktk3-0/apps/testROCk -c bioconda -c conda-forge -c etetoolkit python=3.10 bbmap muscle diamond numpy pandas requests dash fasttree ete3 mafft pplacer taxtastic pyqt
```

# Usage

If using conda, activate the environment

```bash
conda activate ROCkOut
```

Before building the model, ROCkIn is used to develop a list of positive sequences (and their UniProt IDs) and an optional list of negative sequences. Data from the ROCkIn pipeline is used to assist the user in selecting these sequences. See the separate README in the ROCkIn directory for more details.

Manually create a simple text file with one UniProt ID per line. Create one file for the positive IDs (POS) and optionally negative IDs (NEG). Do no leave any blank lines at the end or beginning of the these files.

Run ROCker as a series of steps.

For additional details for each script type: python scriptname.py OPTION -h

*OUT=user specified project directory. Create any name you want for the first script and then use this same directory for the rest of the pipeline.*

To build a ROCker model:
1. Download
1. Build

```bash
python rockout_main.py download [-h] [-d OUT] [-p POS] [-n NEG] [-t THREADS] [-q QUIET]
python rockout_main.py build [-h] [-d OUT] 
```

Then run the rocker_dash.py script and open the IP address it gives in your web browser (copy and paste). Drag and drop or select the ROCkOUT_index.txt file that was created in the [OUT] directory.

```bash
python rocker_dash.py
```

Additional optional steps:
1. Refine
1. refine-ni
1. align

To Use a ROCker model:
1. filter 
1. pplace-prep
