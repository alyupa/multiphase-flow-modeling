#!/bin/bash
# ADD: remove temp files
# ADD: test script parameters: 
#	* processors_count time_limit_in_minutes !=0
#	* if "no", then processors_count = 1
# ADD: text color

# Описание параметров запуска скрипта
if [ -z "$5" ]; then
    echo Use: $0 task_name architecture communication processors_count time_limit_in_minutes [debug/energy]
    echo task_name: 2ph,3ph,bl
    echo architecture: cpu,gpu
    echo communication: no,mpi,shmem
    exit
fi

# Режим отладки
if [ "$6" = "debug" ]
then
    debug="-D MY_TEST"
    debug_name="_$6"
fi

# Режим учета теплообмена
if [ "$6" = "energy" ]
then
    debug="-D ENERGY"
    debug_name="_$6"
	energy="../energy.cpp"
else
	energy=""
fi

case $4 in
   *[^0-9]* )
	echo "Error in processors_count: $4 is not supported by CAPAZ" 
	exit      
         ;;
    * )  
         ;;
esac

case $5 in
   *[^0-9]* )
	echo "Error in time_limit_in_minutes: $5 is not supported by CAPAZ" 
	exit      
         ;;
    * )  
         ;;
esac

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
	    PPN="-ppn 3"
	    lib_path="-L/common/cuda/lib64 -lcudart"
	    maxtime="-maxtime $5"
	else
	    ARCH=20
	    PPN="-ppn 1"
	    lib_path="-L/usr/local/cuda/lib64 -lcudart"
	fi
fi

if [ "$1" = "2ph" ]
then
    task_name="-D TWO_PHASE"
    project_file="two-phase.cpp"
    project_folder="Two-phase"
else
	if [ "$1" = "3ph" ]
	then
	    task_name="-D THREE_PHASE"
	    project_file="three-phase.cpp"
	    project_folder="Three-phase"
	else
		if [ "$1" = "bl" ]
		then
		    task_name="-D B_L"
		    project_file="b-l.cpp"
		    project_folder="Backley-Leverett"
		else
			echo "Error in task_name: $1 is not supported by CAPAZ"
			exit
		fi
	fi
fi

if [ "$2" = "gpu" ]
then
    echo "nvcc $task_name $debug -c -arch sm_$ARCH gpu.o ../gpu.cu"
          nvcc $task_name $debug -c -arch sm_$ARCH gpu.o ../gpu.cu
    arch_file="gpu.o"
else if [ "$2" = "cpu" ]
     then
    	arch_file="../cpu.cpp ../$project_folder/$project_file"
	    lib_path=""
		PPN=""
     else
	echo "Error in architecture: $2 is not supported by CAPAZ"
	exit
     fi
fi

if [ "$3" = "no" ]
then
    comm_file="../no_communication.cpp"
else if [ "$3" = "mpi" ]
     then
	comm_file="../mpi.cpp"
     else if [ "$3" = "shmem" ]
	  then
	    comm_file="../shmem.cpp"
	  else
	    echo "Error in communication: $3 is not supported by CAPAZ"
            exit
	  fi
     fi
fi

if [ "$hostname" = "k100" ]
then
   if [ "$3" = "shmem" ]
   then
      compilator="shmemcc"
   else
      compilator="mpicxx"
   fi
else 	if [ "$hostname" = "mvse" ]
	then		
	    if [ "$3" = "shmem" ]
   	    then
      		compilator="shmemcc"
		else
			compilator="mpicc"
	    fi
	else
	    if [ "$3" = "mpi" ]
            then
		compilator="mpicxx"
	    else
		compilator="g++-4.4"
  	    fi
	fi
fi

if [ "$6" = "energy" ]
then
	echo "Energy is supported" 
else
	echo "Energy is not supported" 
fi

mkdir ../$project_folder/Debug

echo "$compilator $task_name $debug $lib_path ../main.cpp $comm_file $energy ../shared_test.cpp ../gauss.cpp $arch_file -o ../$project_folder/Debug/$2_$3$debug_name.px"
      $compilator $task_name $debug $lib_path ../main.cpp $comm_file $energy ../shared_test.cpp ../gauss.cpp $arch_file -o ../$project_folder/Debug/$2_$3$debug_name.px

cd ../$project_folder/Debug
echo "mpirun $PPN -np $4 $maxtime ../$project_folder/Debug/$2_$3$debug_name.px"
mpirun $PPN -np $4 $maxtime ./$2_$3$debug_name.px

exit
