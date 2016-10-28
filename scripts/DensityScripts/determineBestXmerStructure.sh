#!/bin/bash
if [ $# -gt 0 ]; then
	if [ ! -z $3 ]; then
		if [ ! -z $4 ]; then
			grepbc $2 $1 1 2 | grep $3 | grep $4 | sort -k6nr -k4nr -k5nr | head -n1
		else
			grepbc $2 $1 1 2 | grep $3 | sort -k6nr -k4nr -k5nr | head -n1
		fi
	else
		grepbc $2 $1 1 2 | sort -k6nr -k4nr -k5nr | head -n1
	fi
else
	echo "$0 presenceFile clusterID/gene [,gene1/mutation1/other [,gene2/mutation2/other]]"
fi
