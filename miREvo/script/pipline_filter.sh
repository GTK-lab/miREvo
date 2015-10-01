#!/bin/bash
#--------------------------------------------------
# shopt -s -o nounset
#-------------------------------------------------- 

function USAGE {
	echo ""
	echo "Usage: miREvo filter -o prefix -i reads.fas -d \"database1.fasta database2.fasta .... databaseN.fasta\" -H known.miRNA -M mature.fa [options]"
	echo "    Options for filtering reads"
	echo "        -i  <str>   sequence reads file in FASTA format, unique merged, required. NOTE: FASTA headers must be in miRDeep2 collapsed format (e.g. dme_00001_x18). A copy of collapse_reads_md.pl has been provided in miREvo/script. Run that on reads.fa first."
	echo "        -d  <str>   In quotations marks: a list of the prefixes of the bowtie index for Repeat Databases, each constructed from a Fasta file contain known tRNAs, sRNA, etc., recommended"
	echo "                    For instance, if the reference is 'database.fasta', then the prefix is 'database' and building-command is:"
	echo "                    'bowtie-build -f database.fasta database'"
	echo "        -H  <str>   Botwtie index for known miRNAs' hairpin sequences. These should be the known hairpin sequences for the species being analyzed."
	echo "                    the miRNAs' ID must be in miRBase format, for example '>dme-miR-1'"
	echo "        -M  <str>   fasta file with mature miRNAs. These should be the known mature sequences for the species being analyzed."
	echo "                    the miRNAs' ID must be in miRBase format, for example '>dme-miR-1-5p'"
	echo "        -o  <str>   abbreviation for project name, 3 letter code for the sequencing library or the species of interest,required"
	echo ""
	echo "    Options for Bowtie:"
	echo "        -v  <int>   maximum number of mismatches allowed on a read when mapping to the Repeat Database, <=2. default=2bp"
	echo "        -p  <int>   number of processors to use, default=1"
}

if [ $# -eq 0 ]; then
	USAGE;
	exit 192;
fi

declare -rx SCRIPT="${0##*/}"
#declare -r OPTSTRING="i:d:H:M:o:v:p:n:k:h"
declare -r OPTSTRING="h:i:d:H:M:o:p:v:k"
declare SWITCH
declare DATABASES
declare hairpin=''
declare MATURE=''
declare INPUT
declare prj
declare PROJNAME 
declare MAXREP=5
declare -i MIS=2
declare -i CPU=1
declare -i Ns=5
declare -i SEED=8

program_dir="$MIREVO/script"

while getopts "$OPTSTRING" SWITCH ; do
	case "$SWITCH" in
		h) USAGE;
		   exit 0
		;;
		i) INPUT="$OPTARG" 
		;;
		d) DATABASES="$OPTARG"
		;;
		H) hairpin="$OPTARG"
		;;
		M) MATURE="$OPTARG"
		;;
		o) prj="$OPTARG"
		;;
		p) CPU="$OPTARG"
		;;
		v) MIS="$OPTARG"
		;;
		k) MAXREP="$OPTARG"
		;;
		\?) exit 192
		;;
		*) printf "miREvo filter: $LINENO: %s\n" "script error: unhandled argument"
		exit 192
		;;
	esac
done


if [[ -e $INPUT ]]; then

	header=`grep ">" $INPUT | head -n 1`

	# A quick check of input read fasta's header/defline format
	# only checking the first header.... assuming all headers follow the same format as the first
	# analysis_filter.pl has a line-by-line check
	# e.g. of proper header: dme_00001_x18	
	if [[ !	$header =~ ^\>[a-zA-Z0-9]{3}_[0-9]+_x[0-9]+$ ]]; then
		echo "Improper header format in $INPUT. Please provide headers in the miRDeep2/miREvo format. E.g. dme_00001_x18"
		exit 192
	fi
else
	echo "Can't locate the fasta sequence for input reads: $INPUT";
fi

if [ ! -e $prj ]; then
	mkdir $prj
fi

PROJNAME=`basename $prj`

# Format of file names of DB related output.
# Please use xxx as placeholder for filter group.
# e.g. filter.xxx.db.bwt
# Script will create a filter.file_i.db.bwt 
# for each of the the database file_i inputed under the -d option.
DB_BWT_FORMAT=$prj/filter.xxx.db.bwt
DB_NOMAP_FORMAT=$prj/filter.xxx.db.nomap.fas
DB_MAP_FORMAT=$prj/filter.xxx.db.map.fas

FILTER_STATS=$prj/filter.statistic

KN_BWT=$prj/filter.kn.bwt
KN_NOMAP=$prj/filter.kn.nomap.fas
KN_MAP=$prj/filter.kn.map.fas
CMDLOG=$prj/filter.cmd

declare db_name
declare hairpin_seq

if [[  -e $hairpin.fa  ]]; then
	hairpin_seq=$hairpin.fa
elif [[ -e $hairpin.fas  ]]; then
	hairpin_seq=$hairpin.fas
elif [[ -e $hairpin.fasta  ]]; then
	hairpin_seq=$hairpin.fasta
else
	echo "Can't locate the fasta sequence for hairpin: $hairpin";
	echo "Please provide an available $hairpin (hairpin) sequence, such as";
	echo "$hairpin.fa, $hairpin.fas, $hairpin.fasta, etc."
	echo "and build the Bowtie index using a command like:"
	echo "bowtie-build -f $hairpin.fa $hairpin"
	echo "exit now."
	exit 192;
fi

echo ""

DB_NOMAP=$DB_NOMAP_FORMAT # In case the first one is not working
CUR_INPUT=$INPUT

# Allow $DATABASES to be a list of files
# To support the option of a more detailed filtering statistics
# Note: separating the DATABASE multiple files gives more detailed stats at 
# 	the expense of computational time
for DATABASE in $DATABASES
    do

	# For the first iteration, use $INPUT for bowtie alignment
	# after that, use the filtered results from the previous iteration
	if [ -e $DATABASE.1.ebwt ] ; then
	    db_name=`basename $DATABASE`
	    # Get the filename and replace "." by "_"
	    # If the file is called database.rna.1.fa, db_name_stem will be database_rna_1
	    db_name_stem=`echo $db_name | sed 's/\./_/g'`
	    DB_MAP=`echo $DB_MAP_FORMAT | sed s/xxx/${db_name_stem}/g`
	    DB_NOMAP=`echo $DB_NOMAP_FORMAT | sed s/xxx/${db_name_stem}/g`
	    DB_BWT=`echo $DB_BWT_FORMAT | sed s/xxx/${db_name_stem}/g`

	    if false ; then 
                echo "Filtering out reads mapped to $db_name"
	        echo "bowtie -f -v $MIS -a --best --strata --suppress 5,6,7 $DATABASE $CUR_INPUT --al $DB_MAP --un $DB_NOMAP -p $CPU > $DB_BWT" > $CMDLOG
	        bowtie -f -v $MIS -a --best --strata --suppress 5,6,7 $DATABASE $CUR_INPUT --al $DB_MAP --un $DB_NOMAP -p $CPU > $DB_BWT
	        CUR_INPUT=$DB_NOMAP
            fi
	else
	    echo "Warning: miREvo can't locate the Bowtie index files for repeat database $DATABASE"
	    echo "skip filtering reads in the $DATABASE"
	    echo "cp $CUR_INPUT $DB_NOMAP" > $CMDLOG
	    cp $CUR_INPUT $DB_NOMAP
	fi
done

echo ""
echo "Filtering out reads mapped to known miRNAs"

if [ -e $hairpin.1.ebwt ] ; then
	echo "bowtie -f -v 1 -a --best --strata $hairpin $DB_NOMAP --al $KN_MAP --un $KN_NOMAP -p $CPU > $KN_BWT" >> $CMDLOG
	bowtie -f -v 1 -a --best --strata $hairpin $DB_NOMAP --al $KN_MAP --un $KN_NOMAP -p $CPU > $KN_BWT
	echo "perl $program_dir/analysis_filter.pl  $hairpin_seq $MATURE $KN_BWT $DB_BWT_FORMAT $INPUT > $FILTER_STATS" >> $CMDLOG
	perl $program_dir/analysis_filter.pl $hairpin_seq $MATURE $KN_BWT $DB_BWT_FORMAT $INPUT > $FILTER_STATS
	echo "perl $program_dir/mirna_tag_aln_fasta.pl $hairpin_seq $KN_BWT > $prj/filter.kn.map" >> $CMDLOG
	perl $program_dir/mirna_tag_aln_fasta.pl $hairpin_seq $KN_BWT > $prj/filter.kn.map
else

	echo "Warning: miREvo can't locate the Bowtie index files for known miRNAs $hairpin_seq"
	echo "skip filtering reads in the $hairpin_seq"
	echo "ln -s `basename $DB_NOMAP` $KN_NOMAP" >> $CMDLOG
	if [ -e $KN_NOMAP ] ;  then
		rm $KN_NOMAP
		ln -s `basename $DB_NOMAP` $KN_NOMAP
	else
		ln -s `basename $DB_NOMAP` $KN_NOMAP
	fi
fi

if [[ ! -e $prj/filter.*.bwt && ! -e $prj/filter.kn.bwt ]] ; then
	echo "Warning: "
	echo "Failed to locate Bowtie Index files, nothing done for filter, exit."
	echo "If you wish to continue, run predict command next."
	exit 192;
fi

echo ""
echo "Filter successfully done."
