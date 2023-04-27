# Scoring a ROCker model with labeled input reads:

The following workflow will generate TP/FP style figures from positive (pos) and negative (neg) labeled input reads.

Begin with two fasta files. One file contains the pos reads and the other the neg reads. Then run rocker align, filter, and place. Run align twice, once for each read set using the same directory for ${results_directory}. This will create the output directory structure and files for both reads in the same place. Then the filter and place steps take the same ${results_directory} and will output results for both sets of inputs.

```bash
python ${rocker} align -i 01_simulated_raw_reads_pos.fasta -d ${rocker_model_directory} -f ${results_directory} -t 10
python ${rocker} align -i 02_simulated_raw_reads_neg.fasta -d ${rocker_model_directory} -f ${results_directory} -t 10

python ${rocker} filter -d ${rocker_model_directory} -f ${results_directory} -t 10
python ${rocker} place -d ${rocker_model_directory} -f ${results_directory} -t 10

```

To generate a TP/FP table and bar plot from the results:

```bash
python score_rocker_model.py -p1 (prefiltered_pos_alignment) -p2 (prefiltered_neg_alignment) -f1 (filtered_pos_alignment) -f2 (filtered_neg_alignment) -o (output prefix)
```

To generate a phylogenetic read placement tree simply uploade the ${results_directory}/pplacer_jplace/\*.jplace files to iTol and turn on the "Phylogenetic Placements" option in the Datasets tab.


# To generate tagged reads from a rocker model project:


#### For models built with only the positive ID set

If your rocker model only has positive IDs (i.e. no negative IDs were supplied to build the model) then this step is straightforward. We collected the simulated reads and the positives are labeled as "Target" while the negatives are labeled as "Non_Target"

```bash
cat ${odir}/positive/*/tagged_reads/*_250_* > 00_simulated_raw_reads_pos.fasta
grep -A 1 ';Target' 00_simulated_raw_reads_pos.fasta > 01_simulated_raw_reads_pos.fasta
grep -A 1 ';Non_Target' 00_simulated_raw_reads_pos.fasta > 02_simulated_raw_reads_neg.fasta
```

#### For models built with both positive and negative ID sets

First we will collect all of the simulated reads from genomes attached to positive and/or negative ID's.

```bash
cat ${odir}/positive/*/tagged_reads/*_250_* > 00a_simulated_raw_reads_pos.fasta
cat ${odir}/negative/*/tagged_reads/*_250_* > 00b_simulated_raw_reads_neg.fasta
```

Then because the labelling is a bit akward and a positive and negative ID could some from the same genome we have a python helper script to output the correct positive and negative labeled reads.

```bash
python sort_pos_neg_reads.py -pos 00a_simulated_raw_reads_pos.fasta -neg 00b_simulated_raw_reads_neg.fasta
```
The output will 01_simulated_raw_reads_pos.fasta and 02_simulated_raw_reads_neg.fasta but some additional filtering is done to sort the actual positive and negative reads out.

Then just take the simulated_raw_reads_pos.fasta and simulated_raw_reads_neg.fasta and start up at the top with ${rocker} align.


scripts="/storage/scratch1/9/rconrad6/ROCkOut/00d_Scripts"
python ${scripts}/sort_pos_neg_reads.py -pos 00a_simulated_raw_reads_pos.fasta
sbatch --export odir=model ../../00b_sbatch/rocker_align_filter_place.sbatch
