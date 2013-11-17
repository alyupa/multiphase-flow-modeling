#!/bin/bash

# Usage description
if [ -z "$4" ]; then
	echo Usage: $0 processors_per_node_count processors_count time_limit_in_minutes exe_name
	exit
fi

case $1 in
   *[^0-9]* )
	echo "Error in processors_per_node_count: $1 is not supported"
	exit
		 ;;
	* )
		 ;;
esac

case $2 in
   *[^0-9]* )
	echo "Error in processors_count: $2 is not supported"
	exit
		 ;;
	* )
		 ;;
esac

case $3 in
   *[^0-9]* )
	echo "Error in time_limit_in_minutes: $3 is not supported"
	exit
		 ;;
	* )
		 ;;
esac

if [ $1 = 0 ]; then
	PPN=""
else
	PPN="-ppn $1"
fi


if [ $3 = 0 ]; then
	maxtime=""
else
	maxtime="-maxtime $3"
fi

exe_name=$4
project_folder="Three-phase"

cd ../$project_folder/Debug
echo "mpirun $PPN -np $2 $maxtime ../$project_folder/Debug/$exe_name"
mpirun $PPN -np $2 $maxtime ./$exe_name

exit
