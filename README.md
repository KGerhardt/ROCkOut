# ROCkOut

Work in progress!

Python refactor of ROCker - https://github.com/lmrodriguezr/rocker

# Rethinking the ROCker Approach:

The original ROCker's development was done in a piecemeal fashion. While the resulting tool is functional, it's neither particularly user friendly, nor is it well-designed with respect to the inherently iterative process of building a complete ROCker model. Building a ROCker model requires repeatedly testing candidate proteins and including or excluding them as they enhance or degrade the performance of the final result, yet the program pools data in a manner that requires completely restarting the program's model building process to make any updates at all. The only time that the original ROCker is fully efficient is when you already have a completed model's protein set(s) and simply want to reproduce the final model - a pointless task for someone who wishes to *use* a model, since the final model in that process can simply be downloaded and used without being rebuilt.

Consequently, ROCkOut is not intended to perfectly replicate the original ROCker, but to enhance it. Pooled data and multiple-step scripts are discretized in ROCkOut so that the addition of a protein for testing only requires running building steps for that protein alone, and excluding a protein that is harmful is as simple as ignoring that protein's directory from a subsequent step.

Additionally, ROCkOut is designed to support the user through visual aids and supports that make the building process friendly, rather than combatitive.

# ROCkI/O Modules

Please refer to https://github.com/rotheconrad/ROCkIn for details on ROCkIn's installation and use. ROCkIn is not strictly necessary to use ROCkOut; however, I cannot reccommend its use strongly enough. It will save you days of work time in finding reference protein sequences for a ROCkOut model and will help you communicate both the contents of your ROCkOut model and the rationale for why a model included the sequences it did.

ROCkOut contains a set of functions for creating a ROCkOut model (a classifier for identifying a gene in a metagenome) and a set of functions for using a ROCkOut model.

# ROCkOut model creation functions

## ROCkOut download

ROCkOut's download module takes one collection of UniProt IDs representing a target gene function and optionally a second collection representing genes with similar amino acid sequences but different functions to serve an out-group. These sequences will be collected from UniProt online, if they exist in the current release of UniProt, and the gene sequences of each ID, genomes containing those genes, assosciated proteomes, and some additional metadata are downloaded and processed by ROCkOut. ROCkOut concludes this step by creating a project "snapshot," a tar archive of all of the materials required to recreate an identical ROCkOut model to ensure reproducibility.

ROCkOut creates a hierarchical directory structure during its download step to control like so:

```bash
|-----project_root
    |
    |-----positive
    |   |
    |   |-----positive_UniProt_ID_1
    |   |-----positive_UniProt_ID_2
    |   |-----positive_UniProt_ID_....
    |   ....
    |
    |-----negative
    |   |-----negative_UniProt_ID_1
    |   |-----negative_UniProt_ID_2
    |   |-----negative_UniProt_ID_....
    |   ....
    |
    |-----shared_files
    |   |
    |   |-----downloads
    |
    |-----final_outputs
    |   |
    |   |-----snapshot
           |
           |----ROCkOut_Snapshot.tar.gz

```
These directories will continue to be modified as the ensuing steps of ROCkOut run, adding additional information to each positive and negative subdirectory including simulated reads and building out the shared resources with items such as a multiple alignment of target and negative protein sequences. At the end of a project, the final outputs subdirectory contains the resources to run a ROCkOut model, the project snapshot, train and test results, and interactive figures showing the performance and characteristics of the final model created by ROCkOut.

Running the download step is simple:

```bash
python3 rockout_main.py download -p [positive_uniprot_IDs_file] -n [negative_uniprot_ids_file] -d [project_root_directory] -t [threads]
```

Or with the supplied test inputs in the arch_amoa_example folder:

```bash
python3 rockout_main.py download -p  -n [negative_uniprot_ids_file] -d [project_root_directory] -t [threads]
```

## ROCkOut build

## ROCkOut refine

# ROCkOut classification functions

## ROCkOut

# Additional support

## ROCkOut extract

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
- plotly
- [taxtastic](https://github.com/fhcrc/taxtastic)

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
conda create -n ROCker -c bioconda -c conda-forge -c etetoolkit python=3.10 bbmap muscle diamond numpy pandas requests dash fasttree ete3 mafft pplacer taxtastic pyqt
```

If using an M1 Mac you need to tell conda to use the old osx-64 architecture for several of the program dependencies to install through conda. For more details see [this article](https://towardsdatascience.com/how-to-manage-conda-environments-on-an-apple-silicon-m1-mac-1e29cb3bad12).

*pplacer is not currently supported for conda install on m1 macs. ROCker will still work just not for the final (and optional) pplacer step unless you go through the hassle to manually install pplacer. It may be possible to run it using the precompiled binaries https://github.com/matsen/pplacer/releases and adding pplacer to your system PATH.*

```bash
# For M1 MAC as of Jan 2023
CONDA_SUBDIR=osx-64 conda create -n ROCker -c bioconda -c conda-forge -c etetoolkit python=3.10 bbmap muscle diamond numpy pandas requests dash fasttree ete3 mafft taxtastic pyqt
conda activate ROCker
conda config --env --set subdir osx-64
conda install -c bioconda pplacer
```

# Usage

If using conda, activate the environment

```bash
conda activate ROCkOut
```

Before building the model, [ROCkIn](https://github.com/rotheconrad/ROCkIn) is used to develop a list of positive sequences (and their UniProt IDs) and an optional list of negative sequences. Data from the ROCkIn pipeline is used to assist the user in selecting these sequences. See the [ROCkIn](https://github.com/rotheconrad/ROCkIn) repository for more details.

Manually create a simple text file with one UniProt ID per line. Create one file for the positive IDs (POS) and optionally negative IDs (NEG). Do no leave any blank lines at the end or beginning of the these files.

Run ROCker as a series of steps.

For additional details for each script type: python scriptname.py OPTION -h

*OUT=user specified project directory. Create any name you want for the first script and then use this same directory for the rest of the pipeline.*

To build a ROCker model:
1. Download
1. Build

```bash
python rockout_main.py download [-h] [-d OUT] [-p POS] [-n NEG] [-t THREADS] [-q QUIET]
python rockout_main.py build [-h] [-d OUT] [-t THREADS] 
# interactive mode
python rockout_main.py refine [-h] [-d OUT]
# non-interactive just build plots for review
python rockout_main.py refine-ni [-h] [-d OUT]
```

To Use a ROCker model:

1. Map reads to a set of reference sequencing with a corresponding ROCker model.
1. filter 
1. pplace-prep

```bash
python rockout_main.py filter [-h]
python rockout_main.py pplacer [-h]
```
