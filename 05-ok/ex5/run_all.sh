
for i in 0 1
do
    t=$i
    sed -i '6s/.*/'${t}'/' input.dat
    make run
done
