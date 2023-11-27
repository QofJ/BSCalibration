outname=${1}
iter=${2}
rin=${3}
rout=${4}
el=${5}
dra=${6}
ddec=${7}
echo "${outname}"
source ./setup.bsh
time ./ToffCalProcess_Vf_PreToffGaus "${outname}" "${iter}" "${rin}" "${rout}" "${el}" "${dra}" "${ddec}"
#./KM2Arecdec 2021001.file test.root
