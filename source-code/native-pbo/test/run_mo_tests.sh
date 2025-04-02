#!/bin/bash

time=$1
binary=$2
options=$3

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

errors=0
tested=0
timeouts=0

echo "##################################"
echo "########## TEST MO ALGS ##########"
echo "##################################"
echo ""
echo "timeout: $time"
echo "binary: $binary"
echo "options: $options"
echo ""

declare -a arr_algs=(
p-minimal
paretop-k
)

declare -a arr_opt=(
--mo-core-boosting=0
--mo-core-boosting=1
)

declare -a arr_inst=(
"opb/mo/dal.mopb"
"opb/mo/dal2.mopb"
"opb/mo/kp.mopb"
"opb/mo/ap.mopb"
"opb/mo/iris-mlic.mopb"
)

for alg in "${arr_algs[@]}"; do
  echo "##### $alg #####"
  echo ""
  for opt in "${arr_opt[@]}"; do
    echo "===== opt: $opt ====="
    echo ""
    for inst in "${arr_inst[@]}"; do
      inst="$SCRIPTPATH/instances/${inst}"
      echo "running $binary --mo-alg=$alg $opt $options $inst"
      truth=`grep "^* o" $inst | cut -d ' ' -f2-`
      out=`mktemp`
      err=`mktemp`
      timeout $time $binary --mo-alg=$alg $opt $options $inst > $out 2> $err
      ret=$?
      if [ "$ret" != 0 ]; then
        if [ `cat $err | wc -l` -gt 0 ]; then
          echo "error"
          cat $err
          errors=`expr $ret + $errors`
        else
          echo "timeout"
          timeouts=`expr 1 + $timeouts`
        fi
      else
        output=`cat $out | grep "^o " | sort`
        if [ "$output" != "$truth" ]; then
          errors=`expr 1000 + $errors`
          echo "wrong output"
          echo "truth: ${truth}"
          echo "output: ${output}"
        fi
      fi
      rm $out
      rm $err
      echo "errors: $errors"
      tested=`expr 1 + $tested`
      echo "tested: $tested"
      echo "timeouts: $timeouts"
      echo ""
    done
  done
done

declare -a arr_bialgs=(
bioptsat
)

declare -a arr_biopt=(
--mo-core-boosting=0 --bos-variant=sat-unsat
--mo-core-boosting=0 --bos-variant=oll
--mo-core-boosting=1 --bos-variant=sat-unsat
--mo-core-boosting=1 --bos-variant=oll
)

declare -a arr_biinst=(
"opb/mo/iris-mlic.mopb"
)

for alg in "${arr_bialgs[@]}"; do
  echo "##### $alg #####"
  echo ""
  for opt in "${arr_opt[@]}"; do
    echo "===== opt: $opt ====="
    echo ""
    for inst in "${arr_biinst[@]}"; do
      inst="$SCRIPTPATH/instances/${inst}"
      echo "running $binary --mo-alg=$alg $opt $options $inst"
      truth=`grep "^* o" $inst | cut -d ' ' -f2-`
      out=`mktemp`
      timeout $time $binary --mo-alg=$alg $opt $options $inst > $out
      if [ "$?" != 0 ]; then
        echo "timeout"
        timeouts=`expr 1 + $timeouts`
      else
        output=`cat $out | grep "^o " | sort`
        rm $out
        if [ "$output" != "$truth" ]; then
          errors=`expr 1000 + $errors`
          echo "wrong output"
          echo "truth: ${truth}"
          echo "output: ${output}"
        fi
      fi
      echo "errors: $errors"
      tested=`expr 1 + $tested`
      echo "tested: $tested"
      echo "timeouts: $timeouts"
      echo ""
    done
  done
done

echo "##### Summary #####"
echo "tested: $tested"
echo "errors: $errors"
echo "timeouts: $timeouts"
