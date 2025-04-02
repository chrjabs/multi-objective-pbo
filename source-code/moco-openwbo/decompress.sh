#!/bin/bash
function main(){
    # creates temporary directory in /tmp. It will contain the decompressed instance file
    local tmp
    local path
    # float path to the last position of positional args
    path $@
    printf "instance is %s\n" $path
    if [ ! -f $path ]; then
	echo "nothing found at $path"
	return 1
    fi
    # decompress instance file if necessary. Sets 'tmp' variable
    decompress path
    # replace instance with tmp file in positional args
    if [ ! -z $tmp ]; then
	set -- ${@/$path/$tmp}

    fi
    $@
    clean_tmp
}

function decompress(){
    # decompress into tmp file. tmp file name is the same, without
    # compression extension and decorated in order to be unique.
    case $path in
	*.xz)
	    echo "decompressing xz file..."
	    tmp=$(mmktemp ${path/.xz//})
	    unxz --stdout $path > $tmp
	    ;;
	*.tar.gz)
	    
	    tmp=$(mmktemp ${path/.tar.gz//})
	    echo "decompressing gz file..."
	    zcat $path > $tmp
	    ;;
	*)
	    ;;
    esac
    if [ -v tmp ]; then
	echo "decompression complete"
    else
	echo "no decompression required"
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
function path() 
{
    # sets path to naked path in argv, complains if no such path
    # exists. First element of argv is actually the command to be run.
    for arg in ${@:2}
    do
	if [ ! -z $arg ] && [[ ! $arg =~ ^-.* ]]  && [ -f $arg ] ; then
	    path=$arg
	fi
    done
    if [ -v $path ]; then
	echo "please provide an instance file" 
	exit 1;
    fi
}

function clean_tmp(){
    # clean up tmp file. Without the -z test, -f returns true even if
    # value of $tmp is null.
    if [ ! -z $tmp ] && [ -f $tmp ]; then
	rm $tmp
    fi
}

# pass positional arguments to main function,
main ${@}
