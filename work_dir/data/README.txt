Setting up the data
===================
Following are instructions for setting up the datasets used in this work and
instructions for adding custom datasets.
The names given to the original dataset files match the current configuration of
the code. However, these may be altered, as long as the relevant configuration
options, found in $WORK_DIR/scripts/utils/datasets.py are modified accordingly.

The instructions include the application of the GenomeCRISPR effect labels. Note
that these should not be added to the short samples provided, since there are
not enough guides in them for these labels to be well defined. Hence, obtain the
complete datasets prior to applying the labels.


How to obtain the Xu datasets
=============================
1. The ribosomal and nonribosomal datasets are both part of Supplementary Table
   1 of Xu et al. (2015):
   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4509999/bin/supp_gr.191452.115_Supplemental_Table_1.xlsx
   
2. Both datasets are split into two tabs in the Excel file, one for efficient
   guides, the other for inefficient guides.
   For each of the two datasets, create a CSV file with all of its guides
   (efficient and inefficient) and make sure to include the column headers as
   the first row of the file.
   For the ribosomal dataset, the file should be named chr1_ribo.csv.
   For the nonribosomal dataset, the file should be names chr2_nonribo.csv.
   
3. For both files, add a final column, with the average log2fc of the guides,
   where the averaging is over the two tested cell lines, e.g. the formula in
   cell H9 should be:
		=AVERAGE(F9,G9)
   Samples which demonstrate how the final datasets should look are provided in:
		$WORK_DIR/data/xu/dataset

4. Add the GenomeCRISPR effect labels:
		python $WORK_DIR/scripts/prep/apply_genome_labels.py -d xu -c 1
		python $WORK_DIR/scripts/prep/apply_genome_labels.py -d xu -c 2
	
5. Generate the artificial genomes:
		python $WORK_DIR/scripts/prep/generate_artificial_genome.py -d xu -c 1
		python $WORK_DIR/scripts/prep/generate_artificial_genome.py -d xu -c 2
	
6. Prepare the necessary dataset-related files for the operation of the
   prediction tools:
		python $WORK_DIR/scripts/prep/prep.py -d xu -c 1 2


How to obtain the Doench dataset
================================
1. The Doench dataset is provided in Supplementary Table 7 of Doench et al.
   (2014):
   https://static-content.springer.com/esm/art%3A10.1038%2Fnbt.3026/MediaObjects/41587_2014_BFnbt3026_MOESM8_ESM.xlsx
   
2. Remove the first row of the file (such that the column headers are row number
   1, and save the file as chr1_train.csv in the 'dataset' directory of the
   Doench dataset.
   A sample which demonstrates how the final adjusted dataset should look is
   provided in:
		$WORK_DIR/data/doench/dataset
   
3. Below are the two way of adding the GenomeCRISPR effect labels - in the
   unadjusted and adjsuted way. Choose the desired method and act accordingly.
   a. Unadjusted:
		python $WORK_DIR/scripts/prep/apply_genome_labels.py -d doench -c 1 --positive
   b. Adjusted:
		python $WORK_DIR/scripts/prep/apply_genome_labels.py -d doench -c 1 --shift -0.5 --invert
			
4. Generate the artificial genomes:
		python $WORK_DIR/scripts/prep/generate_artificial_genome.py -d doench -c 1

5. Prepare the necessary dataset-related files for the operation of the
   prediction tools:
		python $WORK_DIR/scripts/prep/prep.py -d doench -c 1


How to obtain the Chari dataset
================================
1. The Chari dataset is provided in the Haeussler repository under
   chari2015Train.ext.tab:
   https://github.com/maximilianh/crisporPaper/blob/master/effData/chari2015Train.ext.tab
   
2. Convert tge tab separators in this file to commas, and save the dataset as
   chr1_239t.csv in the 'dataset' directory of the Chari dataset.
   A sample which demonstrates how the final dataset should look is provided in:
		$WORK_DIR/data/chari/dataset
   
3. Add the GenomeCRISPR effect labels:
		python $WORK_DIR/scripts/prep/apply_genome_labels.py -d chari -c 1 --positive
		
4. Generate the artificial genome:
		python $WORK_DIR/scripts/prep/generate_artificial_genome.py -d chari -c 1
		
5. Prepare the necessary dataset-related files for the operation of the
   prediction tools:
		python $WORK_DIR/scripts/prep/prep.py -d chari -c 1


How to obtain the GenomeCRISPR dataset
======================================
Since we did not use the entire dataset and instead sampled from it randomly, to
allow for a reproduction of our results, we provide the full GenomeCRISPR sample
we generated and used. The sample can be found in:
	$WORK_DIR/data/genome/dataset
When using this pre-made sample, skip steps 1-3 below.

1. The complete GenomeCRISPR database can be foun on the GenomeCRISPR website:
   http://genomecrispr.dkfz.de/#!/download
   
2. To sample from it, place it in the 'dataset' directory of the GenomeCRISPR
   dataset, and use the provided sampling script. To sample a 15k sample, and
   save it as chromosome number 1 of this dataset, execute:
		python $WORK_DIR/scripts/prep/sample_genomecrispr.py -c 1 -s 15000 -i <name of databse file>

3. Add the 'genome_efficient' column which marks guides as efficient or
   inefficient according to their effect label. The condition for gettign
   genome_efficient=1 is that the effect label is less than -1, e.g. cell Q9
   should have the formula:
		=IF(N9<-1, 1, 0)
		
4. Generate the artificial genome:
		python $WORK_DIR/scripts/prep/generate_artificial_genome.py -d genome -c 1

5. Prepare the necessary dataset-related files for the operation of the
   prediction tools:
		python $WORK_DIR/scripts/prep/prep.py -d genome -c 1


How to obtain the Mixture dataset
=================================
Since this dataset was sampled from existing datasets, and this was done late in
the pipeline, it is not available in a CSV format like other datasets.
Instead, we provide a Pickle file for this dataset:
	$WORK_DIR/data/mixture/mixture.pkl
We also provide a list of the targets included in this dataset:
	$WORK_DIR/data/mixture/mixture_targets.txt
The data itself is a list of DataPoint instances (see Core.ipynb), which means
it contains the feature representation of each target included in it. Thus, if
no tools were added or removed, no further processing is required for this
dataset.
If the tools have been changed, then the Pickle file needs to be recreated with
the new feature representations. To do so, run Core.ipynb up to and including
the Get data cell, and then call:
	create_mixture_dataset(chari_data, xu_data, doench_data, genome_data)
This will override the Mixture.pkl file with the updated datapoints.


How to setup a new dataset
==========================
These are general instructions for adding new datasets.

1. Within the data directory, create a new directory for the dataset.

2. Inside the new directory, create a directory named 'dataset'.

3. The 'dataset' directory is where the raw data should be placed, in the form
   of a CSV file.
   
4. The names of the data files must be in the following format:
		chr<number>_<name>.csv
   For example, the ribosomal Xu dataset was designated to be chromosome number
   1 of the Xu datasets collection, hence it is named chr1_ribo.csv.

5. Register the new dataset in the $WORK_DIR/scripts/utils/datasets.py file:
   a. Create a new Dataset instance with the configuration appropriate for the
      CSV placed in the dataset directory.
   b. Add this instance to the DATASETS dictionary under a suitable key.

6. Add the GenomeCRISPR effect labels to each of the chromosomes in the dataset:
		python $WORK_DIR/scripts/prep/apply_genome_labels.py -d <dataset name> -c <chromosome number> [additional flags]
   Remember to set the required flags to set the labels according to the
   original labelling scheme of the dataset.
		
7. Generate the artificial genome for each chromosome:
		python $WORK_DIR/scripts/prep/generate_artificial_genome.py -d <dataset name> -c <chromosome number>
		
8. Prepare the necessary dataset-related files for the operation of the
   prediction tools:
		python $WORK_DIR/scripts/prep/prep.py -d <dataset name> -c <chromosome number>

9. In Core.ipynb, in the Constants cell, add a constant with the path of the
   file containing the feature representations of the new dataset.

10. In the Get data cell, use get_data(<features path>) to read the feature
    representations of the new dataset.

11. If you wish the new dataset to be included in the working set, in the New
    experiment cell, add the data to the datasets dictionary.

12. If you wish to include the new dataset in the comparisons performed in the
    notebook, add it to final_cmp in the Initialise comparisons cell.

13. If you wish to compare the model with DeepCRISPR using this dataset, follow
    the instructions in data/deepcrispr/README.txt.
