#!/bin/bash
# ADD: remove temp files
# ADD: test script parameters:
#	* processors_count time_limit_in_minutes !=0
#	* if "no", then processors_count = 1
# ADD: text color

# Usage description
if [ -z "$2" ]; then
	echo Use: $0 architecture communication [energy] [debug]
	echo architecture: cpu,gpu
	echo communication: no,mpi
	exit
fi

exe_name="$1_$2"

# Energy mode
if [ "$3" = "energy" ]; then
	debug="-D ENERGY"
	exe_name="$exe_name_energy"
	energy="../energy.cpp"
	echo "Energy is enabled"
else
	energy=""
	echo "Energy is not enabled"
fi

# Debug mode
if [ "$3" = "debug" -o "$4" = "debug" ]; then
	debug="-D MY_TEST"
	exe_name="$exe_name_debug"
fi

hostname=`hostname`

if [ "$hostname" = "mvse" ]; then
	ARCH=13
	PPN="-ppn 1"
	lib_path="-L/common/cuda/lib -lcudart"
else 	if [ "$hostname" = "k100" ]; then
		ARCH=20
		PPN="-ppn 3"
		lib_path="-L/common/cuda/lib64 -lcudart"
	else
		ARCH=20
		PPN="-ppn 1"
		lib_path="-L/usr/lib/x86_64-linux-gnu/ -lcudart"
	fi
fi

project_file="three-phase.cpp"
project_folder="Three-phase"

if [ "$1" = "gpu" ]; then
	echo "nvcc $debug -c -arch sm_$ARCH gpu.o ../gpu.cu"
		  nvcc $debug -c -arch sm_$ARCH gpu.o ../gpu.cu
	arch_file="gpu.o"
else if [ "$1" = "cpu" ]; then
		arch_file="../cpu.cpp ../$project_folder/$project_file"
		lib_path=""
		PPN=""
	 else
		echo "Error in architecture: $1 is not supported"
		exit
	 fi
fi

if [ "$2" = "no" ]; then
	comm_file="../no_communication.cpp"
else if [ "$2" = "mpi" ]; then
		comm_file="../mpi.cpp"
	else
		echo "Error in communication: $3 is not supported"
		exit
	fi
fi

if [ "$hostname" = "k100" ]; then
	  compilator="mpicxx"
else 	if [ "$hostname" = "mvse" ]; then
		compilator="mpicc"
	else	if [ "$2" = "mpi" ]; then
			compilator="mpicxx"
		else
			compilator="g++-4.4"
		fi
	fi
fi

mkdir -p ../$project_folder/Debug
exe_name="$exe_name.px"

echo "$compilator $debug $lib_path ../main.cpp $comm_file $energy ../shared_test.cpp ../gauss.cpp $arch_file -o ../$project_folder/Debug/$exe_name"
	  $compilator $debug $lib_path ../main.cpp $comm_file $energy ../shared_test.cpp ../gauss.cpp $arch_file -o ../$project_folder/Debug/$exe_name

exit
