
#Create input for clumpak
SUFFIX=Daphne_neutral_aims_strict
mkdir Clumpak_input_${SUFFIX}
cd Clumpak_input_${SUFFIX}
for K in $(seq 2 10)
do
    mkdir K${K}
    for RUN in $(seq 1 100)
    do
        cp ../K${K}/Run_${RUN}/projAdmix.*.Q K${K}/K${K}_run_${RUN}_admixture
    done
done
cd ../
zip -r Clumpak_input_${SUFFIX}.zip Clumpak_input_${SUFFIX}

