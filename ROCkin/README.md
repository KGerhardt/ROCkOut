# ROCkIn

A pipeline to prepare for building ROCker models with ROCkOut

This pipeline walks the researcher through the process of collecting the sequence information necessary to build and refine ROCker models for any functional gene group of interest. The steps involve a combination of Python and Bash. Both PBS and Sbatch scripts are provided for users with access to a compute cluster.

To use the provide PBS or Sbatch scripts replace "PATH/to/GitHub/repo" with the path to your local copy of this Github repo, and replace all the "YOUR_PROMPTs" with the relevant information.

## Dependencies

- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
 - [Clustal Omega](http://www.clustal.org/omega/)
 - [TrimAl](http://trimal.cgenomics.org/)
 - [FastTree](http://www.microbesonline.org/fasttree/)
 - [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)
 - [Python](https://www.python.org/)

 *Python, it's packages, and all program above can be installed with [Conda](https://docs.conda.io/en/latest/miniconda.html).*

 #### References

 1. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: architecture and applications. BMC bioinformatics. 2009 Dec;10(1):1-9.
 1. Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature biotechnology. 2017 Nov;35(11):1026-8.
 1. Sievers F, Wilm A, Dineen DG, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7:539 doi:10.1038/msb.2011.75
 1. Salvador Capella-Gutiérrez, José M. Silla-Martínez, Toni Gabaldón, trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses, Bioinformatics, Volume 25, Issue 15, 1 August 2009, Pages 1972–1973
 1. Price, M.N., Dehal, P.S., and Arkin, A.P. (2010) FastTree 2 -- Approximately Maximum-Likelihood Trees for Large Alignments. PLoS ONE, 5(3):e9490. doi:10.1371/journal.pone.0009490.
 1. A. Stamatakis: "RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies". In Bioinformatics, 2014, open access.
 1. Sanner MF. Python: a programming language for software integration and development. J Mol Graph Model. 1999 Feb 1;17(1):57-61.

# PART 00: Curate reference sequences

ROCker model building starts with a set of curated reference sequences. Curation of these sequences is up to the researcher building the model and they should be either experimentally verified or highly specific.

Specialized databases are a good starting resources such as the [NCBI ref gene database](https://www.ncbi.nlm.nih.gov/pathogens/refgene) for antibiotic resistance genes.

When curating sequences, it is insightful to look at a multiple sequence alignment, , a phylogenetic tree such as a quick neighbor joining tree, and/or a clustered heatmap of sequence similarity. There are many approaches to this. We will outline a quick and easy one here utilizing EBI's website and we provide a Python script to build a sequence similarity heatmap.

 1. [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/). Select Pearson/FASTA as the output format in step 2.
 2. Download the alignment file and view it with your favorite multiple sequence alignment tool such as [AliView](https://ormbunkar.se/aliview/).
 3. Under the "Results Viewers" tab and select "Send to Simple Phylogeny" at the bottom of the options. In STEP 2 of Simple Phylogeny, turn on the distance correction, exclude gaps, and P.I.M. (percent identity matrix) options.
 4. At this point you should see a phylogram of the results. You can scroll down and "View Phylogenetic Tree File" which is in newick format (.nwk). You can save this file and view/edit it with tools such as [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) or [iTol](https://itol.embl.de/).
 5. For the sequence similarity heatmap, look under the "Results Summary" tab at the top of the Simple Phylogeny results page. Download the "Percent Identity Matrix" file (.pim) and use the 00a_PIM_clustered_heatmap.py Python script included in the 02_Python directory of this GitHub repo to create a heatmap figure.

 ```bash
 python /Path/to/GitHub/repo/02_Python/00a_PIM_clustered_heatmap.py -i your_files_name.pim -o name_your_output_file.pdf

 # for script options and info
 python /Path/to/GitHub/repo/02_Python/00a_PIM_clustered_heatmap.py -h
 ```

 ![Example Figure of sequence similarity heatmap](https://github.com/KGerhardt/ROCkOut/blob/main/ROCkin/05_Example_Figs/00_Example-A.png)

Once you have made your selections, create two fasta formatted files that share the same short meaningful defline names. One should have the nucleotide sequence (.fnn) and the other should have the amino acid sequence (.faa). We will call these RefSeqs.fnn and RefSeqs.faa.

# PART 01 01: UniProt sequence search

Since ROCker models are used with metagenomics data, we want to account for broad sequence diversity around the reference sequences. To do this, we will perform a BLAST search of UniProt's SwissProt and TrEMBL databases.

If there are multiple reference sequences with ≥90% or ≥95% sequence similarity it is suggested to select one representative sequence to use in the BLAST search as it is unlikely they will yield different search results.

I use a script from EBI webservices to access and Blast search the UniProt database programmatically. This returns a separate text file of UniProt IDs for sequence matches for each fasta sequence in the input file.
```bash
# path to ebi blast script
ebi_blast='Path/to/GitHub/repo/02_Python/01a_ebi_ncbiblast.py'

# executes the command Run 30 sequences at a time from 1 fasta file.
python ${ebi_blast} --email YOUR_EMAIL --program blastp --stype protein \
--sequence RefSeqs.faa --database uniprotkb --multifasta --useSeqId --maxJobs 30 --pollFreq 60 \
--outformat ids --exp 10 --alignments 1000
```

PBS example:
```bash
qsub -v fasta=RefSeqs.faa /Path/to/GitHub/repo/01a_PBS/01a_ebi_blast.pbs
```

Then we have to do a few extra steps with the Blast results to get the fasta sequences we need. First we download the .dat file corresponding to the UniProt IDs from our Blast search. We will use dbfetch to retrieve the .dat file for each UniProt ID from each fasta sequence.
```bash
# file house keeping
mkdir 01a_ebi_blast_ids 01b_ebi_dat 01c_ebi_fasta
mv *ids.txt 01a_ebi_blast_ids

# path to dbfetch script
dbfetch='Path/to/GitHub/repo/02_Python/01b_ebi_dbfetch.py'
# loop over ids.txt files and run dbetch fetchData on each UniProt ID.
for f in 01a_ebi_blast_ids/*; do odir='01b_ebi_dat'; gene=`basename $f | cut -d. -f1`; echo $gene; if [ ! -d ${odir}/${gene} ]; then mkdir ${odir}/${gene}; fi; while read p; do name=`echo $p | cut -d: -f2`; python ${dbfetch} fetchData $p > ${odir}/${gene}/${name}.dat; echo $name; done < $f; done
```

PBS example:
```bash
# Download the *.dat files from EBI with dbfetch
for f in 01a_ebi_blast_ids/*; do odir='01b_ebi_dat'; gene=`basename $f | cut -d. -f1`; echo $gene; if [ ! -d ${odir}/${gene} ]; then mkdir ${odir}/${gene}; fi; qsub -v input=${f},odir=${odir},gene=${gene} /Path/to/GitHub/repo/01a_PBS/01b_ebi_dbfetch.pbs; done
```

Now we can parse the .dat file to retrieve the fasta sequence for each UniProt ID returned by the Blast search. We will store additional relevant info from the .dat file in the fasta sequence deflines.
```bash
for d in 01b_ebi_dat/*; do n=`basename $d`; echo $n; python /Path/to/GitHub/repo/02_Python/01d_parse_dat_file.py -i $d -o 01c_ebi_fasta/${n}.fasta; done

# concatenate to single fasta file
cat 01c_ebi_fasta/MCR-* >> 01c_ebi_fasta/ALL_EBI_BLAST_MATCHES.faa
```

# PART 02: Deduplicate, Filter, ClusterReps, align, trim, tree -> annotated.tsv

Setup directories
```bash
mkdir 02a_tree_prep
```

#### Step 1: Deduplicate

Since we have multiple verified sequences that we searched to UniProt, we likely have found overlapping search results. Concatenate fasta files of all search result sequences and deduplicate by UniProt ID. I wrote a Python script to deduplicate the concatenated fasta

```bash
python /Path/to/GitHub/repo/02_Python/02a_Remove_Duplicate_Fasta_Entries.py -f 01c_ebi_fasta/ALL_EBI_BLAST_MATCHES.faa -o 02a_tree_prep/DEDUP_EBI_BLAST_MATCHES.faa
```

#### Step 2: Filter

The initial sequence search does not have great options to filter the matches that are returned. As such, it usually returns many spurious matches, and many nearly identical matches to our curated sequences. We don't need either of these types of sequences to build our ROCker models.

To filter our search results, I run Blastp locally using the RefSeqs.faa as the reference database and the Blast search results as the query. Then I use that to filter which sequences to keep.

I wrote a Python script that applies 3 filters:
	1) Remove matches >= 98% identity to verified sequences
	2) Removes matches <= 30% identity to verified sequences
	3) Remove matches with <= 50% sequence alignment (alignment length / query sequence length).

It also plots histograms of each parameter post filtering and can be rerun with different filtering options if neccessary.

```bash
# make blast database from RefSeqs.faa
makeblastdb -dbtype prot -in RefSeqs.faa
# Run blast search
blastp -num_threads 2 -max_target_seqs 10 \
-db RefSeqs.faa -query 02a_tree_prep/DEDUP_EBI_BLAST_MATCHES.faa -out 02a_tree_prep/FILTER_EBI_BLAST_MATCHES.blast  -subject_besthit -outfmt \
'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'
# Set python script directory
scripts=/Path/to/GitHub/repo/02_Python
# filter blast results
python ${scripts}/02b_Blastp_filter_hist.py -i 02a_tree_prep/FILTER_EBI_BLAST_MATCHES.blast -o 02a_tree_prep/FILTER_EBI_BLAST_MATCHES.fltrd.blast
# retrieve fasta sequences matches filtered blast results
python ${scripts}/02c_Get_Fasta_from_Filtered_Blast.py -b 02a_tree_prep/FILTER_EBI_BLAST_MATCHES.fltrd.blast -q 02a_tree_prep/DEDUP_EBI_BLAST_MATCHES.faa -o 02a_tree_prep/FILTER_EBI_BLAST_MATCHES.faa

```

PBS example:
```bash
qsub -v ref=RefSeqs.faa,qry=02a_tree_prep/DEDUP_EBI_BLAST_MATCHES.faa,out=02a_tree_prep/FILTER_EBI_BLAST_MATCHES.faa,name=genename /Path/to/GitHub/repo/01a_PBS/02b_Blastp.pbs
```

#### Step 3: ClusterReps

If you have few than hundreds of sequences at this point you can skip this step. This step reduces the number of sequences by clustering them at 90% amino acid sequence similarity and choosing one reference sequence for each cluster.

Cluster with MMSeqs2
```bash

```


PBS Example:
```bash
> qsub -v infile=02a_tree_prep/01_MCR_fltr_ebi_blast_fltrd.fasta,odir=02a_tree_prep,n=mcr ../00b_PBS/02c_mmseqs.pbs

> grep -c '>' 02a_tree_prep/02_mmseqs_reps.fasta
``` 

#### Step 4: Align

sequence alignment with clustal omega

For this I am first building an alignment with the 10 curated sequences and then fitting the searched sequences to that alignment. This runs clustal omega twice.

```bash
> qsub -v verified=00b_curated_seqs/ncbi_refgenes_mcr_renamed_reduced.fasta,newseqs=02a_tree_prep/02_mmseqs_reps.fasta,n=mcr ../00b_PBS/02d_seq_alignment.pbs

sequences before trimming
> grep -c '>' 02_mmseqs_reps.fasta.aln
Sequences before trimming: 1120 (1110 searched + 10 curated)
```

#### Step 5: Trim

trim the alignment with trimmal

Before building a phylogenetic tree, MSAs should be cleaned to removed spurious sequences and columns consisting mainly of gaps.

I'm using trimmal for this

```bash
> qsub -v input=02a_tree_prep/02_mmseqs_reps.fasta.aln,output=02a_tree_prep/03_trimmed.fasta.aln,n=mcr ../00b_PBS/02e_seq_trim.pbs

# count sequences after trimming
> grep -c '>' 03_trimmed.fasta.aln
Sequences after trimming: 1119

# clean up sequence names for the tree.
> python ../00c_Scripts/02e_clean_seq_names.py -i 02a_tree_prep/03_trimmed.fasta.aln
```

#### Step 6: Tree

Build phylogenetic tree with RAxML

```bash
# single ML Distance tree
> qsub -v input=02a_tree_prep/03_trimmed.fasta.aln,output=04a_RAxML-Distance,name=MCL_distance ../00b_PBS/02f_RAxML_AminoAcid-Distance.pbs 

# bootstrapped ML tree
> qsub -v input=02a_tree_prep/03_trimmed.fasta.aln,output=04b_RAxML-Bootstrap,name=MCL_bootstrap ../00b_PBS/02f_RAxML_AminoAcid-Bootstrap.pbs

# Fasttree - approximate maximum likelihood tree
> qsub -v input=02a_tree_prep/03_trimmed.fasta.aln,output=04c_FastTree.nwk,n=mcr ../00b_PBS/02f_FastTree.pbs 
```

#### Step 7: Organize annotations

create annotated tsv file to explore clades.

To further inform our decisions from the tree, we want to know which species and functional annotations are in each clade. To do this I wrote python code that uses the branch distance from the tree to create a distance matrix which is used in clustering algorithms. It gets clades from the the tree as clusters and creates a tsv file with species and annotation information for each sequence.

```bash
# concatenate the curated sequences and the dereplicated filtered searched sequences
> cat 00b_curated_seqs/ncbi_refgenes_mcr_renamed_reduced.fasta 02a_tree_prep/02_mmseqs_reps.fasta >> 02a_tree_prep/04_MCR_sequence_set.fasta

# setup directory
> mkdir 05a_Clade_info_distance 05b_Clade_info_bootstrap 05c_Clade_info_fasttree

# convert nwk to distance matrix, cluster the matrix and add annotations

# RAxML Distance
> python ../00c_Scripts/02g_Tree_Distance_Cluster.py -i 04a_RAxML-Distance/RAxML_parsimonyTree.MCL_distance.RAxML -o 05a_Clade_info_distance/05a_Clade_info_distance -f 02a_tree_prep/04_MCR_sequence_set.fasta

# RAxML Bootstrap
> python ../00c_Scripts/02g_Tree_Distance_Cluster.py -i 04a_RAxML-Distance/RAxML_parsimonyTree.MCL_distance.RAxML -o 05b_Clade_info_bootstrap/05b_Clade_info_bootstrap -f 02a_tree_prep/04_MCR_sequence_set.fasta

# Fasttree
> python ../00c_Scripts/02g_Tree_Distance_Cluster.py -i 04c_FastTree.nwk -o 05c_Clade_info_fasttree/05c_Clade_info_fasttree -f 02a_tree_prep/04_MCR_sequence_set.fasta
```