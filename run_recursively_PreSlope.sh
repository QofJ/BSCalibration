set -e
outname=${1}
iter=${2}
rin=${3}
rout=${4}
el=${5}
dra=${6}
ddec=${7}
FirstRunFlag=${8} #是否是从上一次标定结果进行迭代
RunCount=${9} #运行次数

source ./setup.bsh

#erase outname suffix and 
thisRunOutname=${outname%.*}_part${RunCount}.root
echo "${thisRunOutname}"

time ./ToffCalProcess_recur_PreSlope "${thisRunOutname}" "${iter}" "${rin}" "${rout}" "${el}" "${dra}" "${ddec}" "${FirstRunFlag}"

FirstRunFlag="NotFirstRun"
RunCount=$((RunCount+1))
hep_sub run_recursively_PreSlope.sh -argu "${outname}" "${iter}" "${rin}" "${rout}" "${el}" "${dra}" "${ddec}" "${FirstRunFlag}" "${RunCount}" -g lhaaso -mem 600
