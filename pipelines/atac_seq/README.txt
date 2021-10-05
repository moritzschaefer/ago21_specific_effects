To run this pipeline adapt the paths in config.yaml and main.sh and run main.sh

The patch to allow for FRiP-based sample-scaling in the differential analysis, as described in the publication, is provided in snakepipes_atac_sample_scaling.patch

Scale factors were linearly derived from FRiP-scores and are provided with this repository. Here is the calculation:

Ago12_1 had the lowest FRiP score (24.5%) and was taken as reference to calculate the downsampling factors for the other samples:

24.5% / other_frip



