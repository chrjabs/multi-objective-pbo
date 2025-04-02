#!/bin/bash
# artifact list
sat4j_solver="./deps/solver/target\
/org.sat4j.moco.threeAlgorithms-0.0.1-SNAPSHOT-jar-with-dependencies.jar"
openwbo_solver="./open-wbo"

function main(){

algorithm=$1
instance=$2
unset '$2'
case $instance in
    *.xz)
	echo "c decompressing xz file..."
	tmp=$(mmktemp ${instance/.xz//})
	unxz --stdout $instance > $tmp
	instance=$tmp
	echo "c done: instance in $tmp"
	;;
    *.tar.gz)
	
	tmp=$(mmktemp ${instance/.tar.gz//})
	echo "decompressing gz file..."
	instance=$tmp
	zcat $instance > $tmp
	echo "c done: instance in $tmp"
	;;
    *)
	echo "c no decompression required"
	;;
esac
case $1 in
    pareto-mcs)
	java -jar $sat4j_solver -alg 0 $instance ;;
    p-minimal-sat4j)
	java -jar $sat4j_solver -alg 2 $instance ;;
    p-minimal-openwbo)
	$openwbo_solver -cardinality=1 -pb=2 -no-bmo -mem-lim=8192 -cpu-lim=3600\
			-formula=1 -algorithm=8 -pbobjf=4 -eps=1\
			-apmode=1 -no-cubounds -no-clbounds $instance ;;
    core-guided)
	$openwbo_solver -cardinality=1 -pb=2 -no-bmo -mem-lim=8192 -cpu-lim=3600\
			-formula=1 -algorithm=7 -pbobjf=4 -eps=1\
			-apmode=1 -no-cubounds -no-clbounds $instance ;;
    hitting-sets)
	$openwbo_solver -cardinality=1 -pb=2 -no-bmo -mem-lim=8192 -cpu-lim=3600\
			-formula=1 -algorithm=9 -pbobjf=4 -eps=1 -part_par=100\
			-apmode=1 -cubounds -clbounds $instance ;;
    *)
	echo "Check name of the algorithm to run: \"$1\" is not valid."
	exit 1;;
esac
      if [ ! -z $tmp ] && [ -f $tmp ]; then
	  rm $tmp;
      fi
}

function check_artifacts() {
if  [ ! -f $openwbo_solver ]
   then
   echo "please compile openwbo."
   return 1;
fi
if  [ ! -f $sat4j_solver ]
   then
   echo "please compile sat4j."
   return 1;
fi


}
function print_help() {
    if [ -z $1 ] || [ -z $2 ]
then
    cat <<EOF
usage:
./run_me.sh <algorithm> <instance>
<algorithm> is one of the following,
	pareto-mcs, core-guided, hitting-sets, p-minimal-sat4j, p-minimal-openwbo;

<instance> is the path for a valid instance. If the path terminates in
.xz or .tar.gz, the file will be decompressed; Check files in
./instancesFull.zip for examples.
EOF
    exit 1;
fi
}

function mmktemp(){
    # creates tmp file name. Avoids thrashing over existent files
    local result=$(basename $1)
    local n=""
    while [ -f "/tmp/$n$result" ]; do
	n=$((n+1))
    done
    echo "/tmp/$n$result"
}

# live section: 
check_artifacts
print_help $@
main $@

