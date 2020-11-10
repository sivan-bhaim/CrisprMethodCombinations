#!/bin/bash

# This is the original install.sh with minor modifications.
# These are meant to simplify paths.

INSPATH=`echo -e "$1" | sed -e 's/\/$//g'`
if [ -d $INSPATH ] && [ -w $INSPATH ] && [ "$INSPATH" !=  "" ]
then
	rm -rf $INSPATH/crispex_profile
        echo "copying scripts to $INSPATH/"
        cp -r $PWD/SCRIPTS/* $INSPATH/
        #cp $PWD/gff_usage.txt $INSPATH/gff_scripts/
	echo -e "export CRISPEX=$INSPATH" >>$INSPATH/crispex_profile
	source $INSPATH/crispex_profile
	echo "[SUCCESS]:phytoCRISP-Ex installation path set to $INSPATH"
	echo "[SUCCESS]:phytoCRIPS-Ex files copied to $INSPATH"
	echo
	echo "[PERFORM]:To use the program, Run -> source $INSPATH/crispex_profile"
	echo "[PERFORM]:Once above command is executed, using same terminal, move to directory containing your input files and run phytoCRISP-Ex as [$CRISPEX/phytoCRISPex -in <input.fasta> -db <genomefile.fasta>]"
	echo "[NOTE]:Please remember to run [source $INSPATH/crispex_profile] everytime, before executing phytoCRISP-Ex"
	
else
	echo -e "[ERROR]:Installation PATH doesn't exsist OR User don't have proper permissions to use the directory"
	echo -e "[PERFORM]:Please run [./install.sh <phytoCRISP-Ex_INSTALLATION_PATH>]"
fi	
