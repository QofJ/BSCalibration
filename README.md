# 运行标定程序
> 考虑到可能需要>200次迭代，因此在运行时需要使用递归的方式进行运行，每次运行的结果作为下一次运行的输入，若要停止则需手动删除任务

## 编译
```bash
source ./setup.bsh
g++ -o ToffCalProcess_recur_PreSlope ToffCalProcess_internalReview_Pre.cc ./src/*.cc `root-config --cflags --libs` -I ./include/ -L ./lib/ -l:slalib64.a -lMinuit -lMatrix
```
## 运行脚本
```bash
set -e #出错则退出
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
hep_sub run_recursively_PreSlope.sh -argu "${outname}" "${iter}" "${rin}" "${rout}" "${el}" "${dra}" "${ddec}" "${FirstRunFlag}" "${RunCount}" -g lhaaso
```
## 运行
### Slope形状预置偏置
配置：
- ROI半径1.6$^\circ$
- 伽马能量>10TeV

首次运行：
```bash
hep_sub run_recursively_PreSlope.sh -argu CalVf_MPV_PreToff_confit0p0001_E10_0_1p6.root 10 0 1.6 10 0 0 FristRun 1 -g lhaaso -mem 600
```
非首次运行：
```bash
hep_sub run_recursively_PreSlope.sh -argu CalVf_MPV_PreToff_confit0p0001_E10_0_1p6.root 10 0 1.6 10 0 0 NotFirstRun 2 -g lhaaso -mem 600
```
### Gaussian抽样预置偏置
配置：
- ROI半径1.6$^\circ$
- 伽马能量>10TeV
```bash
hep_sub run_recursively_PreGaus.sh -argu CalVf_MPV_PreToffGaus_confit0p0001_E10_0_1p6.root 10 0 1.6 10 0 0 FristRun 1 -g lhaaso -mem 600
```
# 注意事项
- EOS文件系统，在提交作业时必须使用xrootd协议读写，尽管在本地测试时可以不使用xrootd协议前缀，但是在提交作业时必须使用xrootd协议前缀，否则会出现无法读取文件的情况
