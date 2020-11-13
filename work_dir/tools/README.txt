How to setup a new tool
=======================
These are general instructions for adding new tools such that they can be easily
incorporated into the feature representation and analysis.

1. Install the desired tool as you normally would, and place the installation
   directory in the tools directory.

2. Register the added tool in $WORK_DIR/scripts/utils/tools.py by doing the
   following:
   a. Create a class for your tool which inherits from Tool. The main purpose of
      this class is to define how to read the output created by the tool. Follow
	  the documentation in the code to configure this class correctly.
   b. Note that on top of configuring the input to the __init__ method, you may
      need to override some additional methods of the class.
   c. Finally, add an instance of your new class to the TOOLS dictionary at the
      bottom of the file under a key matching the tool's name.

3. For automation purposes, each tool has a script which operates it, located in
   $WORK_DIR/scripts/too_utils. The script is named run_<tool_name>.py.
   Notmally, we do not run these scripts directly, instead, they are called by
   the main run.py script. Hence, these tool-specific scripts should all expose
   the same interface:
   -p: The path to the main work dir.
   -d: The dataset on which to run the rool.
   -c: The chromosome number on which the tool should run (where for us,
       chromosome numbers are just indices to sub-datasets within a dataset
	   collection, e.g. the nonribosomal Xu dataset is chromosome number 2
	   within the Xu collection).
   Create a script which operates the added tool, and saves the output as a CSV
   file and normalises it.

4. The feature representation itself is created by
		$WORK_DIR/scripts/utils/to_features.py
   After running each tool for a given dataset, you will aggregate all the
   outputs to a single CSV, where each score is represented by a column. For
   each such column, the function get_column_handlers registers a handler, which
   controls how the score is represented in the feature representation. You will
   note that, other than mm10db, all the columns use the same basic
   handler.
   If your tool added x new scores (i.e. columns to the aggregated output), add
   x new handlers to this function.
   For example, if all the handlers for the new scores are the standard handler:
		handlers = \
			[encode_score for _ in range(8)] +\
			[encode_mm10db] +\
			[encode_score for _ in range(3)] +\
			[encode_score for _ in range(x)]

5. In the Constants cell of Core.ipynb, update NUM_FEATURES to match the length
   of the new feature representation.

6. Add the features created by the new tool to the Features cell.
   a. Add each new feature's name to the FEATURES list.
   b. For each new feature, if it is a scoring method, add its index to
      COL_TO_SCORE_TOOL; if it is a decision method, add its index to
	  COL_TO_DECISION_TOOL.
   c. Decision tools require a function that maps their output to 0/1 (where 0
      means 'reject' and 1 means accept). If decision features were added,
	  register their mapping functions in COL_TO_DECISION_FUNCTION in the
	  Decision functions cell.

7. Follow the Mixture dataset instructions in data/README.txt to recreate the
   Mixture dataset.
