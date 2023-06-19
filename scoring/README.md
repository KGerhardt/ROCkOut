# Scoring a ROCker model with labeled input reads:

To score a ROCker model, our method is to collect a testing set of input UniProt IDs and to build a second model.

From here we use the provided "extract_read_labels.py" scripts to collect and label the simulated short reads.

```bash
python extract_read_labels.py ${rocker_model_dir} > mock_metagenome_labeled.fasta
```

Now we can use this mock metagenome built from the testing set of gene sequences to score the rocker model. Treat it as you would an a real metagenome and run through rocker align, filter, and pplacer.

### STEP 01: Run read mapping (blast/diamond alignments)
python ${rockout} align -i ${infile} -d ${model} -f ${outpre} -t 20

### STEP 02: Run the filter
python ${rockout} filter -d ${model} -f ${outpre} -t 20

### STEP 03: Run pplacer
python ${rockout} place -d ${model} -f ${outpre} -g both -t 20



### temp step. best hit filter
* kenji is adding a best hit filter to rockout but its not there yet.

python 

### hmmer

# build hmm model with training set proteins

# FragGeneScan on testing set mock metagenome

# hmm search from FragGeneScan

# filter hmm results for best hit.

### F1 scores