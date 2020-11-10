#!/bin/bash
## Developer: ACHAL RASTOGI, PhD, ENS
## INPUT: file containing sequences to be targeted in fa format (min sequence = 1)
## RUN: /bin/bash wrapper.sh <INPUTFILE.fa> <GENOMENAME.fa> <NGG/NAG> <G/N>


# This is the original wrapper.sh with minor modifications.
# These are meant to adjust the tool to assume the extension of the chromosome
# files is .fa rather than .fasta.
# We have also modified it to accept additional configuration options through
# the command line.


## defining names

prev_dir=$PWD

DATA_DIR="${WORK_DIR}/data/"
cd $DATA_DIR$1
GENOME=$2
USER=$2

## preparing genome for matching

$CRISPEX/makeblastdb -in $GENOME.fa -dbtype nucl -out $GENOME.db

## preparing INPUT file

perl $CRISPEX/one_seq_one_line.pl $USER.fa one_liner

for i in $(cat one_liner | sed '1d' | sed 's/\t//g')
do
	header=`echo $i | awk -F "@" '{print $1}' | sed 's/ /_/g' | sed 's/_$//g' | sed 's/[\t\r> ]//g'`
	seq=`echo $i | awk -F "@" '{print $2}' | sed 's/ /_/g' | sed 's/_$//g' | sed 's/[\t\r ]//g' | tr "[a-z]" "[A-Z]"`
	length=`echo $seq | wc -c`
	echo $i | awk -F "@" '{print $2}' | sed 's/ /_/g' | sed 's/_$//g' | sed 's/[\t\r ]//g' | tr "[a-z]" "[A-Z]" >db
	
	if [[ "$3" == "NAG" ]]
        then
		if [[ "$4" == "N" ]]
		then
			## for positive strand
                	echo $seq | grep -o ".....................AG" >positive
                	perl $CRISPEX/mapping.pl db positive pos_mapped
                	cat pos_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""+""\t"$1"\t"$2+1"\t"$2+23}' >pos.tsv

                	## for reverese strand
                	echo $seq | grep -o "CT....................." >negative
                	perl $CRISPEX/mapping.pl db negative neg_mapped
                	cat neg_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""-""\t"$1"\t"$2+1"\t"$2+23}' >neg.tsv
		elif [[ "$4" == "G" ]]
		then
			## for positive strand
                        echo $seq | grep -o "G....................AG" >positive
                        perl $CRISPEX/mapping.pl db positive pos_mapped
                        cat pos_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""+""\t"$1"\t"$2+1"\t"$2+23}' >pos.tsv

                        ## for reverese strand
                        echo $seq | grep -o "CT....................C" >negative
                        perl $CRISPEX/mapping.pl db negative neg_mapped
                        cat neg_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""-""\t"$1"\t"$2+1"\t"$2+23}' >neg.tsv
		else
			echo -e "[WARNING] You have either not specified the CRISPR START BASE (G/N) or did it wrong. The default program will execute with G as a start of CRISPR targets."
                	## for positive strand
                	echo $seq | grep -o "G....................AG" >positive
                	perl $CRISPEX/mapping.pl db positive pos_mapped
                	cat pos_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""+""\t"$1"\t"$2+1"\t"$2+23}' >pos.tsv

                	## for reverese strand
                	echo $seq | grep -o "CT....................C" >negative
                	perl $CRISPEX/mapping.pl db negative neg_mapped
                	cat neg_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""-""\t"$1"\t"$2+1"\t"$2+23}' >neg.tsv
		fi
        
	elif [[ "$3" == "NGG" ]]
	then
		if [[ "$4" == "N" ]]
                then
			## for poisitive strand
                        echo $seq | grep -o ".....................GG" >positive
                        perl $CRISPEX/mapping.pl db positive pos_mapped
                        cat pos_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""+""\t"$1"\t"$2+1"\t"$2+23}' >pos.tsv

                        ## for reverese strand
                        echo $seq | grep -o "CC....................." >negative
                        perl $CRISPEX/mapping.pl db negative neg_mapped
                        cat neg_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""-""\t"$1"\t"$2+1"\t"$2+23}' >neg.tsv
		elif [[ "$4" == "G" ]]
                then
			## for poisitive strand
                        echo $seq | grep -o "G....................GG" >positive
                        perl $CRISPEX/mapping.pl db positive pos_mapped
                        cat pos_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""+""\t"$1"\t"$2+1"\t"$2+23}' >pos.tsv

                        ## for reverese strand
                        echo $seq | grep -o "CC....................C" >negative
                        perl $CRISPEX/mapping.pl db negative neg_mapped
			cat neg_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""-""\t"$1"\t"$2+1"\t"$2+23}' >neg.tsv
		else
			echo -e "[WARNING] You have either not specified the CRISPR START BASE (G/N) or did it wrong. The default program will execute with G as a start of CRISPR targets."
			## for poisitive strand
                	echo $seq | grep -o "G....................GG" >positive
                	perl $CRISPEX/mapping.pl db positive pos_mapped
                	cat pos_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""+""\t"$1"\t"$2+1"\t"$2+23}' >pos.tsv

                	## for reverese strand
                	echo $seq | grep -o "CC....................C" >negative
                	perl $CRISPEX/mapping.pl db negative neg_mapped
                	cat neg_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""-""\t"$1"\t"$2+1"\t"$2+23}' >neg.tsv
		fi

	else
		echo -e "[WARNING] You have either not specified the PAM (NGG/NAG) sequence or did it wrong. The default program will execute with NGG as the PAM sequence."
		if [[ "$4" == "N" ]]
                then
			## for positive strand
                        echo $seq | grep -o ".....................GG" >positive
                        perl $CRISPEX/mapping.pl db positive pos_mapped
                        cat pos_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""+""\t"$1"\t"$2+1"\t"$2+23}' >pos.tsv

                        ## for reverese strand
                        echo $seq | grep -o "CC....................." >negative
                        perl $CRISPEX/mapping.pl db negative neg_mapped
                        cat neg_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""-""\t"$1"\t"$2+1"\t"$2+23}' >neg.tsv
		elif [[ "$4" == "G" ]]
                then
			## for positive strand
                        echo $seq | grep -o "G....................GG" >positive
                        perl $CRISPEX/mapping.pl db positive pos_mapped
                        cat pos_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""+""\t"$1"\t"$2+1"\t"$2+23}' >pos.tsv

                        ## for reverese strand
                        echo $seq | grep -o "CC....................C" >negative
                        perl $CRISPEX/mapping.pl db negative neg_mapped
                        cat neg_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""-""\t"$1"\t"$2+1"\t"$2+23}' >neg.tsv
		else
			echo -e "[WARNING] You have either not specified the CRISPR START BASE (G/N) or did it wrong. The default program will execute with G as a start of CRISPR targets."
			## for positive strand
			echo $seq | grep -o "G....................GG" >positive
                	perl $CRISPEX/mapping.pl db positive pos_mapped
                	cat pos_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""+""\t"$1"\t"$2+1"\t"$2+23}' >pos.tsv

                	## for reverese strand
                	echo $seq | grep -o "CC....................C" >negative
                	perl $CRISPEX/mapping.pl db negative neg_mapped
                	cat neg_mapped | sed '1d' | awk -F "\t" '{print "'$header'""\t""-""\t"$1"\t"$2+1"\t"$2+23}' >neg.tsv 
		fi
        fi

	cat pos.tsv neg.tsv >targets-$header.tsv
	rm -rf pos.tsv neg.tsv db positive negative neg_mapped pos_mapped
done
	
## making fa inputfile for strategy 1

for m in $(ls targets-*.tsv)
do
	name=`echo $m | sed -e 's/[\r\t ]//g' | sed 's/targets-//g' | sed 's/\.tsv//g'`
	cat $m | awk -F "\t" '{print ">"$1"_"$2"_"$4"_"$5"_"$3"\n"$3}' >tmp-$name.fa	
	
	## ADDING HEADER TO THE OUTPUT FILE
	echo -e "QUERY_STRAND_START_STOP_TargetSEQUENCE,CheckForOffTargetActivityCOMPLETEsequence,CheckForOffTargetActivitySEEDsequence,RestrictionEnzymeSiteAtCas9ClevageSite,RestrictionEnzymeSiteAtDifferentPositionsWithinTargetSequence" >$name-$3pam.csv
	
	## BLAST
	# This command line has been modified: the -V flag was added to resolve a
	# segmentation fault.
	$CRISPEX/blastall -p blastn -d $GENOME.db -i tmp-$name.fa -e 1000000.0 -m 8 -V -o $name.blast
	
	## 1st analysis

	for i in $(cat $name.blast | awk -F "\t" '{print $1}' | sort -g | uniq)
	do

        	checkone=`cat $name.blast | awk -F "\t" '{if($1 == "'$i'"){print 23-($8-$7+1)}}' | sort -g | uniq -c | sed 's/^ *//g' | awk -F " " '{if($2 == 0){print $1}}'`
        	if [ -z "$checkone" ] || [ $checkone -le 1 ]
        	then
               		checktwo=`cat $name.blast | awk -F "\t" '{if($1 == "'$i'"){print 23-($8-$7+1)}}' | sort -g | uniq -c | sed 's/^ *//g' | awk -F " " '{if($2 <= 4 && $2 != 0){print $0}}'`
                	if [ -z "$checktwo" ]
                	then
                        	result=`echo "PASS"`
                	else
                        	result=`echo "FAIL"`
                	fi
        	else
                	result=`echo "FAIL"`
        	fi


	## 2nd analysis (last 15 basepair conservation)

        rogerone=`cat $name.blast | awk -F "\t" '{if($1 == "'$i'" && $7 == 9 && $8 == 23){print $0}}'`
        	if [ -z "$rogerone" ]
        	then
                	resulttwo=`echo "PASS"`
        	else
                	resulttwo=`echo "FAIL"`
        	fi

	## 3rd analysis (presence of restriction enzymes)

        	echo $i | awk -F "_" '{print $NF}' >temp4mapping
        	perl $CRISPEX/mapping2.pl temp4mapping $CRISPEX/restrictionenzymes RESsites
        	resultthree=`cat RESsites | sed 's/[\t\r]//g' | sed 's/ $//g'`

	## RESULTANT

                        if [[ $i =~ .*_-_.* ]]
                then
                        resultfour=`echo $resultthree | sed 's/[\t\r]//g' | sed 's/ $//g' | sed 's/;//g' | tr " " "\n" | egrep -e ">5|>6" | tr "\n" ";" | sed 's/;$//g'`
                        resultfive=`echo $resultthree | sed 's/[\t\r]//g' | sed 's/ $//g' | sed 's/;//g' | tr " " "\n" | egrep -v -e ">5|>6" | tr "\n" ";" | sed 's/;$//g'`
                        echo -e "$i,$result,$resulttwo,$resultfour,$resultfive" >>$name-$3pam.csv
                else
                        resultfour=`echo $resultthree | sed 's/[\t\r]//g' | sed 's/ $//g' | sed 's/;//g' | tr " " "\n" | egrep -e ">16|>17" | tr "\n" ";" | sed 's/;$//g'`
                        resultfive=`echo $resultthree | sed 's/[\t\r]//g' | sed 's/ $//g' | sed 's/;//g' | tr " " "\n" | egrep -v -e ">16|>17" | tr "\n" ";" | sed 's/;$//g'`
                        echo -e "$i,$result,$resulttwo,$resultfour,$resultfive" >>$name-$3pam.csv
                fi


	## REMOVING TEMP files

        	rm -rf temp4mapping RESsites
	done
done

## cleaning the tmp files

rm -rf *one_liner
rm -rf *.tsv
rm -rf tmp*.fa
rm -rf *.blast
rm -rf $GENOME.db*

cd $prev_dir
