
for i in 0.1 0.3 0.5 0.7 1 1.5 2 2.5 3
do
    t=$i
    sed -i '1s/.*/'${t}'/' input.dat
    ./clean.sh
    echo "y" | make run
done
