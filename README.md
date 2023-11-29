# 运行标定程序
```bash
outname=${1}
iter=${2}
rin=${3}
rout=${4}
el=${5}
dra=${6}
ddec=${7}
ToffsetTxt=${8} #这个文件用于迭代慢时自动存储中间标定结果

source ./setup.bsh

echo "${outname}"

#判断是否存在
if [ ! -f "${ToffsetTxt}" ]; then
    time ./ToffCalProcess_rec "${outname}" "${iter}" "${rin}" "${rout}" "${el}" "${dra}" "${ddec}"
    ToffsetTxt="Toffset_${outname}.txt"
fi
if [ -f "${ToffsetTxt}" ]; then
    time ./ToffCalProcess_rec "${outname}" "${iter}" "${rin}" "${rout}" "${el}" "${dra}" "${ddec}" "${ToffsetTxt}"
fi
hep_sub run_Recursively.sh --argu "${outname}" "${iter}" "${rin}" "${rout}" "${el}" "${dra}" "${ddec}" "${ToffsetTxt}"
```
#./KM2Arecdec 2021001.file test.root
``````


