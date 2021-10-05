#!/bin/bash

# requires bw files referenced in main.ini

pyGenomeTracks --tracks main.ini --region chr1:180150000-181500000 -o chr1_region.png
pyGenomeTracks --tracks main.ini --region chr17:23678249-23683446  -o Cldn6.png
pyGenomeTracks --tracks main.ini --region chr10:96615000-96625000  -o Btg1.png
pyGenomeTracks --tracks main.ini --region chr1:180892108-180900103 -o Lefty2.png
pyGenomeTracks --tracks main.ini --region chr1:180934022-180939400 -o Lefty1.png
