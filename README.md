# IGVF Coding Variants Focus Group - Pillar Project

This repository holds the pillar project data frame for the IGVF coding variants focus group. 

Data produced through IGVF and shared internally with other IGVF members should **not** be released to the public until the laboratory that originally created the data approves the public release. Please be mindful of this when sharing data and releasing the work products from internal collaborations.

# Curated datasets

The 'Pillar_project_data_files" folder contains data from published and unpublished functional assays. The naming convention is as follows: Gene_Author_year_[optional_dataset_identifier]

Supplemental information from published datasets was exported to a single Excel sheet for each assay when appropriate. Any changes made to the supplement have been commented in the code. 

Please refer to this sheet for additional information about the functional data. The 'Dataset_tag' column should match the Dataset column in the dataframe.

https://docs.google.com/spreadsheets/d/1EicKYz_AR5gNJzcSNDgjePKYvapRo3htBw5Dzuc-QDM/edit?gid=1962463955#gid=1962463955

# Supporting files and dataframe

The final dataframe has been added to zenodo (due to the large nature of the dataframe, it cannot be hosted on base Github)

Supporting files from Ensembl VEP downloads and ClinVar hg19 and hg38 VCF (parsed) have been supplied for the annotations. They can be found here: https://doi.org/10.5281/zenodo.14736656

Due to the large nature of some of these files, they cannot be stored on base GitHub. 

The files are private on Zenodo, I will add users as requested.

# Script 

The VEP_gnomad_parse is jupyter notebook. Please read the comments carefully. This script will create the VEP_input that was added to the web interface and will also concatenate the VEP_output (if stored in the same folder). In addition, it will concatenate the single gene gnomad_v4 downloads into one file. 

The 'harmonized data script_v4' is a jupyter notebook. Running the script with the correct paths in your local (or other) environment should produce the final annotated dataframe. There are checks built in that will keep count of variants. If you wish to add another dataset, define the dictionary as seen in the script and run the data_harmonization_loop; make sure to add your dictionary to dc_list. 

**Note** In order to run the harmonized data script_v4, you will also need hg19ToHg38.over.chain.gz from https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/ to perform the liftover for datasets that need it. 

Refer to the PP_Readme file for column definitions. 

# Adding to the data

If you wish to add to this dataframe, please fork this repository, make the necessary changes and submit a pull request. 

# Disclaimer

Please check back monthly for the updated dataframe. This is a draft of the final dataframe and will change in the future. If you find any bugs, please reach out to me. 

