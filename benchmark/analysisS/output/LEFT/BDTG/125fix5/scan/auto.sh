#!/bin/bash
[ -f $1 ] || exit 1

ALL="submit_all.sh"
echo "#!/bin/bash" > $ALL
while read line
do
    a=($line)

    filename="toyMC${a[0]}.C"
    cp toyMC.C.orig $filename
    sed -i "s/CASECASE/${a[0]}/g" $filename
    sed -i "s/BDTGBDTG/${a[1]}/g" $filename
    sed -i "s/SIGMASIGMA/${a[2]}/g" $filename
    sed -i "s/ALPHAALPHA/${a[3]}/g" $filename
    sed -i "s/POWERPOWER/${a[4]}/g" $filename
    sed -i "s/GMEANGMEAN/${a[5]}/g" $filename
    sed -i "s/GSIGMAGSIGMA/${a[6]}/g" $filename
    sed -i "s/CBFRACCBFRAC/${a[7]}/g" $filename
    sed -i "s/LINEARLINEAR/${a[8]}/g" $filename

    filename2="submit${a[0]}.sh"
    cp submit.sh.orig $filename2
    sed -i "s/CASECASE/${a[0]}/g" $filename2
    chmod 755 $filename2

    filename3="exec${a[0]}.submit"
    echo "condor_submit $filename3" >> $ALL

    filename4="exec${a[0]}.submit"
    cp exec.submit.orig $filename4
    sed -i "s/CASECASE/${a[0]}/g" $filename4
done < $1
