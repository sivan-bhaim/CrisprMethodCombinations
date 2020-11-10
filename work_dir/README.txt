work_dir tree structure
================================
$WORK_DIR
├───data
│   ├───dataset_name1
│   │   ├───dataset
│   │   ├───indexes
│   │   └───outputs
│   ├───dataset_name2
│   │   ├───dataset
│   .   ├───indexes
│   .   └───outputs
│
├───models
├───scripts
│   ├───prep
│   ├───setup
│   ├───tool_utils
│   └───utils
│
├───tools
│   ├───tool1
│   ├───tool2
│   .
│   .
│
└───core.ipynb

The data directory is meant to contain the various datasets.
The models directory will contain all the models saved by core.ipynb, and comes
with the pre-trained final model presented in the paper.
The scripts directory contains the infrastructure code of this work.
The tools directory is where the installed prediction tools will be located.
Finally, core.ipynb is the Jupyter notebook which contains the code that deals
with the ML side of this work.


General Setup Instructions
==========================

These instructions were written and tested for Ubuntu 20.04.1 LTS.

1. Register the main work directory (in which this file is located) as an
   environment variable:
		echo export WORK_DIR=<full path to the work dir> >> ~/.bashrc
		export WORK_DIR=<full path to the work dir>
		
2. Make sure all scripts have execution permissions:
		chmod -R +x $WORK_DIR/scripts/setup
		
3. Install all prerequisites:
		$WORK_DIR/scripts/setup/prerequisites
   Remember to provide root permissions when prompted for it by the script.
   
4. Install the required tools:
		$WORK_DIR/scripts/setup/tools
   Remember to provide root permissions when prompted for it by the script.
   
-- Open a new terminal window, to allow modifications of environment variables to take place --

5. Optional: modify the tools used in the feature representation by following
   the instructions in tools/README.txt

6. Follow the instructions in data/README.txt to obtain the necessary datasets.

7. For each dataset in {chari, doench, genome, xu} (plus any additional dataset
   you may have added to the collection) and each of their chromosome numbers:
   a. Operate each of the tools {chopchop, flashfry, mm10db, phytocrispex,
      sgrnascorer*, ssc} (plus any additional tool you may have added to the
	  collection) using the automation script:
		python $WORK_DIR/scripts/tool_utils/run.py -t <tool name prefix> -d <dataset name> -c <chromosome number>
      Where any uniquely identifyable prefix would do. For example, to run SSC
      on the nonribosomal Xu dataset:
		python $WORK_DIR/scripts/tool_utils/run.py -t ss -d xu -c 2
      The outputs will be found in the 'outputs' directories of the datasets.
	  
      * Due to the long time sgRNA Scorer takes to run on the genome dataset, we
        provide the tool's normalised output for this dataset:
			$WORK_DIR/data/genome/sgrnascorer_output.normalised
        To use this ready-made output rather than run the tool on the dataset,
        this file should be moved to:
			$WORK_DIR/data/genome/outputs/sgrnascorer/chr1_targets.normalised
        Alternatively, split the GenomeCRISPR sample into chunks, run sgRNA
        Scorer on each chunk separately, and then unify the normalised outputs
        to a single file.
	
8. For each dataset and each of their chromosome numbers, aggregate all the
   outputs:
		python $WORK_DIR/scripts/utils/aggregate_outputs.py -d <dataset name> -c <chromosome number>
   The outputs will be found in the 'outputs' directories of the datasets.
   
9. For each dataset and each of their chromosome numbers generate the feature
   representations:
		python $WORK_DIR/scripts/utils/to_features.py -d <dataset name> -c <chromosome number>
   The outputs will be found in the 'data' directory.
   
10. Now all the files are ready, and the next step is to train and evaluate a
    model. To do so, open the Jupyter notebook provided (we used Google Colab to
    operate it), and follow the documentation inside.


Additional notes:
=================
The code includes additional utilities (listed below). All of them provide a
help string (through the '-h' flag) and have internal documentation.
- $WORK_DIR/scripts/tool_utils/run_deepcrispr.py: Runs a DeepCRISPR model.
  Additional instruction for generating the required files for a DeepCRISPR
  comparison are found in $WORK_DIR/data/deepcrispr/README.txt
- $WORK_DIR/scripts/utils/flip_strand.py: Reverses and complements a DNA
  strand.
- $WORK_DIR/scripts/utils/get_len.py: Returns the length of the genome in a
  chromosome file.
- $WORK_DIR/scripts/utils/get_precision.py: Returns the precision (along with
  additional statistics) of one of the supported tools. Includes various
  configurable options for the definition of precision, to account for the
  behaviours of different tools.
  