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
python3 rockout_main.py download -p AmoA_A.positive.txt -n AmoA_A.negative.txt -d arch_amoa -t [threads]
```

## ROCkOut build

ROCkOut's build step primarily encompasses read simulation and alignment. This step starts with a ROCkOut project directory (whatever was specified by -d in the download step) and all of its other options concern how it simulates and aligns reads.

ROCkOut takes the whole genomes downloaded in the previous steps and simulates reads from each at four read target read lengths. These lengths are uniform distributions over a target interval specified by a lower and upper bound, with default target lengths of 100bp, 150bp, 250bp, and 300bp, all with ranges of +/- 10% of the average read length (e.g, the 150bp read set is a uniform distribution over [135, 150]bp). If these were manually set using ROCkOut's options, it would look like so:

```bash
python3 rockout_man.py build -d arch_amoa --short-lower 90 --short-upper 110 --med-lower 135 --med-upper 165 --long-lower 225 --long-upper 275 --xl-lower 270 --xl-upper 330
```

ROCkOut collects these reads into files with naming conventions indicating the average read length (rounded to the nearest integer), such as arch_amoa/positive/A0A060HNG6/CP007536_read_len_100_tagged.fasta. If the build step is run again, any target read lengths which produce an already-extant average read length will overwrite the file; however, you can simulate at different read lengths to produce more sets of reads representing the protein.

ROCkOut aligns these reads to all of the project's gene sequences, both positive and negative, using DIAMOND, and reports the alignments in tabular blast format (see https://www.metagenomics.wiki/tools/blast/blastn-output-format-6).

The other major options used by ROCkOut's build step control simulated errors within the reads (--snp-rate controls the likelihood of each base on a read being a SNP, up to a maximum of 3/read, --insertrate and --deleterate control the probability of an insertion or deletion, see BBMap RandomReads usage for details), the average depth of coverage to which reads will be simulated (default --coverage 20, 20x depth on average; more is better for modelling but slower), and how hard DIAMOND will try to align reads (1-4 for increasing sensitivity but decreasing speed).

There is also one other significant options: --no-clean. The simulated reads produced by ROCkOut can consume a great deal of disk space, and so as each read length is simulated and aligned, ROCkOut will filter the raw, simulated reads down to only the small number of reads that actually aligned to a protein within the project. --no-clean simply turns off this filtering step. Reads are generated deterministically, so they can always be recovered by a subsequent run of ROCkOut either way.

To run the build step using ROCkOut's defaults, you can use the simpler command:

```bash
python3 rockout_main.py build -d arch_amoa -t [threads]
```

Additional notes on the build step:

The reason that multiple read lengths are simulated is that ROCkOut models are intended to be agnostic to the read lengths of the datasets it will filter. In practice, read length affects all of the metrics used by ROCkOut to filter reads (bitscore, percent identity, and percent of the read which aligned), meaning that it requires models that can account for read length. To achieve this, ROCkOut uses the four simulated read lengths as "ground truths" at likely read lengths for short read metagenomes, and interpolates cutoffs between adjacent, simulated models. It would be possible to directly similate every read length, but doing so is impractical.

ROCkOut's alignment step also includes an unusual decision: rather than retaining the best match for each read, or just a few good matches, ROCkOut retains _every_ match for each read to any protein in the project. This is done because it's possible to replicate the effect of aligning only one best hit by subsequently filtering the more complete set of alignments down to only the top match by read. Exactly which match is best depends on what sequences are in a project, and it's common to remove sequences from a ROCkOut project during model refinement. By keeping the complete set of reads and dynamically filtering, ROCkOut will only need to run the computationally expensive alignment process once, unless a new sequence from UniProt is added to the project.

## ROCkOut refine

ROCkOut's refine step is where it actually constructs a ROCkOut model from the simulated reads. ROCkOut first checks the project's positive subdirectory and creates a multiple alignment of the protein sequences found there. It then loads the read alignments into memory, filters the alignments to only those mapping to a _positive_ protein sequence, and then finds the best match per read among the remaining alignments. In essence, this step is equivalent to running the DIAMOND alignments to only the positive proteins in the project and retaining only the best single match per read.

The loaded reads are divided into five partitions of train/test datasets by protein sequence. That is, ROCkOut selects (by defaults) 60% of the positive protein sequences to serve as training data and 40% of the positive protein sequences to serve as testing data, then gathers all of the reads originating from these sequences into their appropriate train/test collections. Reads aligning to a sequence in the opposite collection (i.e. training reads aligning to a testing sequence) are removed. This process simulates recovering reads from a natural environment, which is unlikely to contain _exact_ replicas of any of the genes within a ROCkOut project. ROCkOut then creates an independent classification model over each training dataset and assesses its performance in classifying reads from the corresponding testing dataset, recovering an F1 score. The five models are weighted by their F1 score, and their cutoffs are combined by weighted average.

The process is repeated for each simulated read length, and ROCkOut records a final model encoding the averaged classification cutoffs for each read length.

In addition to producing a ROCkOut model, ROCkOut can produce a phylogeny for use with an optional step of ROCkOut, see 

```bash
python3 rockout_main.py -d arch_amoa -t [threads]
```

Notes on ROCkOut models:

ROCkOut models are an ensemble of three constitutent models, each classifying reads using a similar approach over different slices of the read alignments. The three filters are all 2D plots of the same set of aligned reads, differing in their X and Y axes. The first plot is identical to the original ROCker, mapping read alignment position in the multiple alignment of target sequences on the X axis against bitscore on the Y axis. The second plot is read alignment position on the X axis (same as the first) against the percent similarity between the aligned portion of the read and the protein to which it best aligns on the Y axis. The third plot consists of each reads' alignment fraction on the X axis against percent identity on the Y axis (same as the second plot).

Each of these constituents has a cutoff curve calculated by moving a sliding window along its X axis and finding the most discriminant cutoff along the Y axis for reads falling in that window. That is, the point on the Y axis that maximizes Youden's J statistic (sensitivity + specificity -1). The sliding window is moved over one unit on the X axis (1 amino acid for the first two constitutents, a bin of 2.5% alignment fraction for the third) and the process is repeated until every window has been covered.

When used by ROCkOut, each model is used separately to classify each read in a dataset; that is, each read is classified as positive or negative three times, once by each constituent model. These three "votes" are then used to make a final determination. This process results in better overall perfromance as each of the three models each perform better in some edge cases compared to the other two, meaning that the vote covers the shortcomings of each constituent well.

# ROCkOut classification functions

## ROCkOut align

ROCkOut's align function takes a metagenome supplied by the user and runs DIAMOND to align it to a ROCkOut project. Only the final_outputs subdirectory of a ROCkOut project is required to run this step, something done to increase the portability of a ROCkOut model.

For convenience, ROCkOut's align function includes two ways to supply reads (in FASTA format): (1) as a comma-separated list of file paths using -i/--input_reads, or (2) as a directory containing read files and only read files using --reads_dir. The -d option specifies a ROCkOut project directory containing a finalized ROCkOut model and the -f/--filter_directory option specifies an output location for the read alignments. Outputs are automatically named according to the input files.

The --threads options in this module is passed to DIAMOND, meaning that it is useful to supply multiple threads even when you have only one set of reads to align.

To ensure that you have a set of reads to align, we're going to use ROCkOut's extract funtion to collect the reads and alignments from the arch_amoa example.

```bash
mkdir amoa_reads amoa_reads/raws amoa_reads/extracted_alignments
python3 final_rockout_code/rockout_main.py extract -d arch_amoa -t 10 -a amoa_reads/extracted_alignments/arch_amoa -r amoa_reads/raws/arch_amoa
python3 final_rockout_code/rockout_main.py align -d arch_amoa/ -t [threads] -f arch_amoa_alns  -f amoa_reads/rockout_filtering --reads_dir amoa_reads/raws/
```

## ROCkOut filter



## ROCkOut Place

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
