#!/bin/bash
# Two cpu are required
# ./measuring.sh gpu - it will measure both host-host and host-device times
# In this case two cpu and one gpu are required

if [ "$1" = "gpu" ]
then
	echo "gpu"
else	if [ "$1" = "cpu" ]
		then
		echo "cpu"
		else
		echo "I must receive gpu or cpu option!"
		exit
	fi
fi

taskD="-D THREE_PHASE"

hostname=`hostname`
#echo $hostname

if [ "$hostname" = "mvse" ]
then
    ARCH=13
    PPN="-ppn 1"
	lib_path="-L/common/cuda/lib -lcudart"
else 	if [ "$hostname" = "k100" ]
	then
	    ARCH=20
	    PPN="-ppn 1"
	    lib_path="-L/common/cuda/lib64 -lcudart"
	else
	    ARCH=20
	    PPN="-ppn 1"
	    lib_path="-L/usr/local/cuda/lib64 -lcudart"
	fi
fi

if [ "$hostname" = "k100" ]
then
      compilator="mpicxx"
else 	if [ "$hostname" = "mvse" ]
	then		
			compilator="mpicc"
	else
		compilator="mpicxx"
	fi
fi

if [ "$1" = "gpu" ]
then
	debug_name="onepoint_gpu"
	project_folder="3Ph_GPU_PointTime"
	arch_file="./Debug/onepoint.o"
	echo "nvcc $taskD -c -arch sm_$ARCH ../gpu.cu -o $arch_file"
	nvcc $taskD -c -arch sm_$ARCH ../gpu.cu -o $arch_file
	
	echo "$compilator $taskD $lib_path ./$project_folder/main.cpp ../no_communication.cpp ../shared_test.cpp $arch_file -o ./Debug/$debug_name.px"
    $compilator $taskD $lib_path ./OnePointTime/main.cpp ../no_communication.cpp ../shared_test.cpp $arch_file -o ./Debug/$debug_name.px

fi

if [ "$1" = "cpu" ]
then
	debug_name="onepoint_cpu"
	project_folder="OnePointTime"
	echo "$compilator $taskD $lib_path ./$project_folder/main.cpp ../cpu.cpp ../shared_test.cpp ../no_communication.cpp ../Three-phase/three-phase.cpp -o ./Debug/$debug_name.px"
    $compilator $taskD $lib_path ./$project_folder/main.cpp ../cpu.cpp ../shared_test.cpp ../no_communication.cpp ../Three-phase/three-phase.cpp -o ./Debug/$debug_name.px
fi

mkdir ./Debug

cd ./Debug

echo "mpirun $PPN -np 1 -maxtime 10 ./$debug_name.px"
mpirun $PPN -np 1 -maxtime 10 ./$debug_name.px

exit

