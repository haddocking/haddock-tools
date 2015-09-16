#!/bin/bash

source /home/abonvin/haddock2.2/haddock_configure.sh

[[ $# -lt 1 ]] && echo "usage: ./$( basename $0 ) <run directory> [--force]" && exit 0

rundir=$1
force=$2

if [ ! -d $rundir ]
then
  echo "Folder does not exist: $1"
  exit 1
elif [ ! -e ${rundir}/run.cns ]
then
  echo "Folder is not a valid run directory (does not have run.cns)"
  exit 1
fi

echo "== $rundir =="
cd ${rundir}

w_it0=$( egrep 'structures_0=[0-9]+;' run.cns | cut -c 21- | sed -e 's/[^0-9]//g' )
w_it1=$( egrep 'structures_1=[0-9]+;' run.cns | cut -c 21- | sed -e 's/[^0-9]//g' )
w_itw=$( egrep 'waterrefine=[0-9]+;' run.cns | cut -c 20- | sed -e 's/[^0-9]//g' )
n_it0=$( ls structures/it0/ | egrep 'pdb$' | wc -l )
n_it1=$( ls structures/it1/ | egrep 'pdb$' | wc -l ) 
n_itw=$( ls structures/it1/water | egrep 'pdb$' | wc -l ) 

if [ "$n_it0" -ne "$w_it0" ]
then
  echo ":::> it0: incomplete ($n_it0/$w_it0)"
  relaunch="True"
elif [ "$n_it1" -ne "$w_it1" ]
then
  echo ":::> it0: ok ($n_it0/$w_it0)"
  echo ":::> it1: incomplete ($n_it1/$w_it1)"
  relaunch="True"
elif [ "$n_itw" -ne "$w_itw" ]
then
  echo ":::> it0: ok ($n_it0/$w_it0)"
  echo ":::> it1: ok ($n_it1/$w_it1)"
  echo ":::> water: incomplete ($n_itw/$w_itw)"
    relaunch="True"
elif [ ! -e "structures/it1/water/analysis/DONE" ]
then
  echo "::> Modelling 100% complete"
  echo "::> Analysis missing"
  relaunch="True"
else
  relaunch="False"
  echo "::> All completed"
fi     
 
if [ $relaunch == "True" ]
then
  # .out files younger than 5 minutes
  if [ ! -z "$( find . -name "*.out" -cmin -5 )" ]
  then
    if [ -z $force ]
    then
      echo "::> HADDOCK is still running"
      exit 1
    fi
  fi
if [ "$n_itw" -ne "$w_itw" ]
then
  echo ":::> Adjusting random seed"
  curseed=$( egrep 'iniseed=[0-9]+' run.cns | sed -e 's/[^0-9]//g' )
  cp run.cns run.cns-seed=${curseed}
  sed -i -e "s/iniseed=${curseed};/iniseed=$RANDOM;/" run.cns
fi
  echo ":::> Relaunching HADDOCK"
  rm -f *it0* *job* *it1* *out* *inp *.w*
  nohup /usr/bin/python ${HADDOCK}/Haddock/RunHaddock.py >& run.out &
fi
