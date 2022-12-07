# ROCkIn

A pipeline for building ROCker models with ROCkOut

### Dependencies

 - 

# Step 00: Curate reference sequences

ROCker model building starts with a set of curated reference sequences. Curation of these sequences is up to the researcher building the model and they should be either experimentally verified or highly specific.

Specialized databases are good starting resources such as the [NCBI ref gene database](https://www.ncbi.nlm.nih.gov/pathogens/refgene) for antibiotic resistance genes.

When curating sequences, it is insightful to look at a multiple sequence alignment, , a phylogenetic tree such as quick a neighbor joining tree, and/or a clustered heatmap of sequence similarity. There are many approaches to this. We will outline a quick and easy one here utilizing EBI's website and we provide a Python script to build a sequence similarity heatmap.

 1. [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/). Select Pearson/FASTA as the output format in step 2.
 2. Download the alignment file and view it with your favorite multiple sequence alignment tool such as [AliView](https://ormbunkar.se/aliview/).
 3. Under the "Results Viewers" tab and select "Send to Simple Phylogeny" at the bottom of the options. In STEP 2 of Simple Phylogeny, turn on the distance correction, exclude gaps, and P.I.M. (percent identity matrix) options.
 4. At this point you should see a phylogram of the results. You can scroll down and "View Phylogenetic Tree File" which is in newick format (.nwk). You can save this file and view/edit it with tools such as [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) or [iTol](https://itol.embl.de/).
 5. For the sequence identity heatmap, look under the "Results Summary" tab at the top of the Simple Phylogeny results page. Download the "Percent Identity Matrix" file (.pim) and use the 00a_PIM_clustered_heatmap.py Python script included in the 02_Python directory of this GitHub repo to create a heatmap figure.

 ```bash
 > python path/to/script/00a_PIM_clustered_heatmap.py -i your_files_name.pim -o name_your_output_file.pdf
 ```

Once you have made your selections, create two fasta formatted files that share the same short meaningful defline names. One should have the nucleotide sequence (.fnn) and the other should have the amino acid sequence (.faa).

# Step 01: UniProt sequence search

Since ROCker models are used with metagenomics data, we want to account for broad sequence diversity around the reference sequences. To do that we will perform a BLAST search of UniProt's SwissProt and TrEMBL databases.

If their are multiple reference sequences with ≥90% or ≥95% sequence similarity it is suggested to select one representative sequence to use in the BLAST search as it is unlikely they will yield different search results.


I use a script from EBI webservices to do this search in a terminal.

I wrapped it in a pbs job to run on PACE.

```bash
# This returns UniProt IDs for sequence matches
> qsub -v fasta=00b_curated_seqs/ncbi_refgenes_mcr_renamed_reduced.fasta ../00b_PBS/01a_ebi_blast.pbs

# file house keeping
> mkdir 01a_ebi_blast_ids 01b_ebi_dat 01c_ebi_fasta
> mv *ids.txt 01a_ebi_blast_ids

# Download the *.dat files from EBI with dbfetch
> for f in 01a_ebi_blast_ids/*; do odir='01b_ebi_dat'; gene=`basename $f | cut -d. -f1`; echo $gene; if [ ! -d ${odir}/${gene} ]; then mkdir ${odir}/${gene}; fi; qsub -v input=${f},odir=${odir},gene=${gene} ../00b_PBS/01b_ebi_dbfetch.pbs; done

# Parse the .dat file into a fasta file
# place relevant info in the sequence deflines we will use downstream
> for d in 01b_ebi_dat/*; do n=`basename $d`; echo $n; python ../00c_Scripts/01d_parse_dat_file.py -i $d -o 01c_ebi_fasta/${n}.fasta; done

# concatenate to single fasta file
> cat 01c_ebi_fasta/MCR-* >> 01c_ebi_fasta/00_MCR_all_ebi_matches.fa
```

#######
Step 02: Deduplicate, Filter, ClusterReps, align, trim, tree -> annotated.tsv
#################################

# setup directories
> mkdir 02a_tree_prep

a. Deduplicate

Since we have multiple verified sequences that we searched to UniProt, we likely have found overlapping search results. Concatenate fasta files of all search result sequences and deduplicate by UniProt ID. I wrote a Python script to deduplicate the concatenated fasta

> python ../00c_Scripts/02a_Remove_Duplicate_Fasta_Entries.py -f 01c_ebi_fasta/00_MCR_all_ebi_matches.fa -o 02a_tree_prep/00_MCR_dedup_ebi.fa

		Total sequences in file: 10000
		Duplicates Removed: 6225
		Unique sequences retained: 3775

It seems there were a lot of the same sequences matching to the various mcr clades. 

b. Filter

The initial sequence search does not have great options to filter the matches that are returned. As such, it usually returns many spurious matches. Additionally, we do not need to keep sequences that are 100% identical to our curated sequences as this does not add any new information for the ROCker model.

To filter our search results, I run Blastp locally using the curated sequence as the reference database and the search results as the query.

I wrote a Python script that does 3 filters:
	1) Remove matches >= 98% identity to verified sequences
	2) Removes matches <= 30% identity to verified sequences
	3) Remove matches with <= 50% sequence alignment (alignment length / query sequence length).

It also plots histograms of each parameter post filtering and can be rerun with different filtering options if neccessary.

I wrapped it into a pbs script.

> qsub -v ref=00b_curated_seqs/ncbi_refgenes_mcr_renamed_reduced.fasta,qry=02a_tree_prep/00_MCR_dedup_ebi.fa,out=02a_tree_prep/01_MCR_fltr_ebi,name=mcrA ../00b_PBS/02b_Blastp.pbs

Total number of entries in blast file: 38581
Number of entries failing the filters: 1253
Number of entries passing the filters: 37328
Number of duplicate blast matches passing filter to remove: 33553
Number of best hit entries written to new file: 3775 

c. ClusterReps

Since we still have a lot of sequences lets dereplicate the set some.

mmseqs clustering at 90% amino acid sequence similarity

> qsub -v infile=02a_tree_prep/01_MCR_fltr_ebi_blast_fltrd.fasta,odir=02a_tree_prep,n=mcr ../00b_PBS/02c_mmseqs.pbs

> grep -c '>' 02a_tree_prep/02_mmseqs_reps.fasta

Representative sequences retained: 1110

d. Align

sequence alignment with clustal omega

For this I am first building an alignment with the 10 curated sequences and then fitting the searched sequences to that alignment. This runs clustal omega twice.

> qsub -v verified=00b_curated_seqs/ncbi_refgenes_mcr_renamed_reduced.fasta,newseqs=02a_tree_prep/02_mmseqs_reps.fasta,n=mcr ../00b_PBS/02d_seq_alignment.pbs

# sequences before trimming
> grep -c '>' 02_mmseqs_reps.fasta.aln
Sequences before trimming: 1120 (1110 searched + 10 curated)

e. trim

trim the alignment with trimmal

Before building a phylogenetic tree, MSAs should be cleaned to removed spurious sequences and columns consisting mainly of gaps.

I'm using trimmal for this

> qsub -v input=02a_tree_prep/02_mmseqs_reps.fasta.aln,output=02a_tree_prep/03_trimmed.fasta.aln,n=mcr ../00b_PBS/02e_seq_trim.pbs

# count sequences after trimming
> grep -c '>' 03_trimmed.fasta.aln
Sequences after trimming: 1119


clean up sequence names for the tree.

> python ../00c_Scripts/02e_clean_seq_names.py -i 02a_tree_prep/03_trimmed.fasta.aln

f. tree

Build phylogenetic tree with RAxML

# single ML Distance tree
> qsub -v input=02a_tree_prep/03_trimmed.fasta.aln,output=04a_RAxML-Distance,name=MCL_distance ../00b_PBS/02f_RAxML_AminoAcid-Distance.pbs 

# bootstrapped ML tree
> qsub -v input=02a_tree_prep/03_trimmed.fasta.aln,output=04b_RAxML-Bootstrap,name=MCL_bootstrap ../00b_PBS/02f_RAxML_AminoAcid-Bootstrap.pbs

# Fasttree - approximate maximum likelihood tree
> qsub -v input=02a_tree_prep/03_trimmed.fasta.aln,output=04c_FastTree.nwk,n=mcr ../00b_PBS/02f_FastTree.pbs 

g. annotated.tsv

create annotated tsv file to explore clades.

To further inform our decisions from the tree, we want to know which species and functional annotations are in each clade. To do this I wrote python code that uses the branch distance from the tree to create a distance matrix which is used in clustering algorithms. It gets clades from the the tree as clusters and creates a tsv file with species and annotation information for each sequence.

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


#######
Scratch
#################################

ValueError: Buffer dtype mismatch, expected 'double_t' but got 'long'

python ../00c_Scripts/00c_ebi_dbfetch.py fetchData tr:A0A2P5GGK6_9ENTR fasta > MCR-9.1_Escherichia_coli.ids.A0A2P5GGK6_9ENTR.fasta


dbfetch='/storage/home/hcoda1/9/rconrad6/p-ktk3-0/02_ROCker/00c_Scripts/00c_ebi_dbfetch.py'

while read p; do name=`echo $p | cut -d: -f2`; python ${dbfetch} fetchData $p fasta > 01_fasta/MCR1-${name}.fasta; echo ${name}; done < MCR-1.1_Gammaproteobacteria.ids.txt


ebi_blast='/storage/home/hcoda1/9/rconrad6/p-ktk3-0/02_ROCker/00c_Scripts//00b_ebi_ncbiblast.py'

python ${ebi_blast} --email rotheconrad@gatech.edu --program blastp --stype protein --sequence  --database uniprotkb --multifasta --useSeqId --maxJobs 30 --pollFreq 60 --outformat ids,out --exp 10 --alignments 1000