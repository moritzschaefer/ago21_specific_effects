# Running diffTF from the Zaugg lab: https://git.embl.de/grp-zaugg/diffTF
# This is not the ideal style of writing this (should be handled from within Snakemake).
# However the snakemake pipeline has been written in a way that would require a lot of adaptions. This for loop is simpler.. :)

MUTANTS="Ago12 Ago1 Ago2"

# Alternatively, it would also possible to combine multiple mutants and run them all together
for mutant in $MUTANTS; do
    sed "s/Mutant/$mutant/g" config.json.template > config.json
    sed "s/Mutant/$mutant/g" sampleData.tsv.template > sampleData.tsv

    # Note: You need to bind the directories that contain the necessary data for singularity to work
    snakemake --use-singularity --singularity-args "--bind /home/schamori/,/mnt/cclab_nas/moritz_data" --snakefile diffTF_repo/src/Snakefile --configfile config.json --cores 15
done
