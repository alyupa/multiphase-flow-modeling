#!/bin/bash
# Two cpu are required
# ./measuring.sh gpu - it will measure both host-host and host-device times
# In this case two cpu and one gpu are required

if [ "$1" = "gpu" ]
then
    taskD="-D USE_GPU"
fi

debug_name="measuring"

hostname=`hostname`
#echo $hostname

if [ "$hostname" = "mvse" ]
then
    ARCH=13
    PPN="-ppn 1"
	lib_path="-L/common/cuda/lib -lcudart"
    maxtime="-maxtime $5"
else 	if [ "$hostname" = "k100" ]
	then
	    ARCH=20
	    PPN="-ppn 1"
	    lib_path="-L/common/cuda/lib64 -lcudart"
	    maxtime="-maxtime $5"
	else
	    ARCH=20
	    PPN="-ppn 1"
	    lib_path="-L/usr/local/cuda/lib64 -lcudart"
	fi
fi

if [ "$1" = "gpu" ]
then
	echo "nvcc $taskD -c -arch sm_$ARCH ./measuring.cu -o ../Debug/measuring.o"
	nvcc $taskD -c -arch sm_$ARCH ./measuring.cu -o ../Debug/measuring.o
	arch_file="../Debug/measuring.o"
else
	arch_file=""
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

mkdir ../Debug

echo "$compilator $taskD $lib_path ./measuring.cpp $arch_file -o ../Debug/$debug_name.px"
      $compilator $taskD $lib_path ./measuring.cpp $arch_file -o ../Debug/$debug_name.px

cd ../Debug
echo "mpirun $PPN -np 2 -maxtime 30 ./$debug_name.px"
mpirun $PPN -np 2 -maxtime 30 ./$debug_name.px

exit

