source ./setup.bsh
i=${1}
# i>10 break
if [ $i -gt 10 ]; then
    echo "i>10 break"
    exit 0
fi
if [ $i -le 10 ]; then
    echo "i<=10"
    echo "i=$i"
    echo "i+1"
    i=$((i+1))
    echo "i=$i"
    echo "run recurTest.sh"
    hep_sub run_recurTest.sh -argu $i -g lhaaso
    # ./run_recurTest.sh $i
fi
