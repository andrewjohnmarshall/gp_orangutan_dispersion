# gp_orangutan_dispersion

2021-06-16

This folder contains all code necessary to replicate analysis, results, and figures for Marshall et al. orangutan analysis in Oecologia.

https://doi.org/10.1007/s00442-021-04964-1

Title: Biotic and abiotic drivers of dispersion dynamics in a large-bodied tropical vertebrate, the Western Bornean orangutan

Authors: Andrew J. Marshall*, Matthew T. Farr, Lydia Beaudrot, Elise F. Zipkin, Katie L. Feilen, Loren G. Bell, Endro Setiawan, Tri Wahyu Susanto, Tatang Mitra Setia, Mark Leighton, and Heiko U. Wittmer

*Corresponding author: ajmarsha@umich.edu

This code requires the following data files:

- actual_m_walked.csv
- Matrix_rain
- Matrix_tmax.csv
- Matrix_tmin.csv
- Matrix1_OHcounts.csv
- Matrix2_StemsMRfruits.csv
- Matrix3_ObservationCovariates.csv
- Matrix4_SiteCovariates.csv
- Ntimes_sitelength_walked.csv

Workflow:

The file "DataFormat.R" uses the raw *.csv data files listed above, formats them for analysis, and writes them to the file "DSdata.R".

The file "ORAN_DS_script.R" runs the BUGS model drawing data from "DSdata.R". The model results are written to the file "ORAN_OUTPUT.Rdata".

The file "ORAN_figs.Rmd" draws data from  "ORAN_OUTPUT.Rdata" and produces the figures presented in the paper and supplemental online material. The values presented in the tables are also calculated. 

All code ran and produced the relevant output in R 4.1.0 on 2021-06-16 using an M1 Mac Mini running OS 11.4.

This repo is a duplicate of information available at https://figshare.com/articles/software/code/14731890 and https://doi.org/10.6084/m9.figshare.14731866

