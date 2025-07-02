# ldscR
The `ldscR` package is designed for performing LDSC in R. LDSC remains the most popular and convenient method for estimating additive heritability. This package provides efficient functions that allow for the simultaneous estimation of the genetic covariance matrix across multiple traits, without applying the LDSC Python package (https://github.com/bulik/ldsc) using an entry-by-entry scheme. Additionally, we use block-wise subsampling to estimate the standard errors.

## Installation
You can install the `ldscR` package directly from GitHub using the following command:
```r
# Install the devtools package if you haven't already
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# Install ldscR from GitHub
devtools::install_github("harryyiheyang/ldscR")
```

## Usage
Here's an example workflow using ldscR:

### Data Preparation and Harmonization
First, download and prepare your GWAS summary data. Then, harmonize alleles using the filter_align() function:
```r
library(ldscR)
data_url <- "http://tinyurl.com/nhdfwd8v"
temp_file <- tempfile()
download.file(data_url, temp_file, mode="wb")
gwaslist=readRDS(temp_file)
unlink(temp_file)
download_ldsc_data("snpinfo") ## Donwload the SNP information file from Netlify
gwaslist=filter_align(gwas_data_list=gwaslist,ref_panel=SNPInfo)
```

### Download and load LDSC data
```r
# Example: download and load EURLDSC into the current session
download_ldsc_data("eurldsc")
# You may also load EASLDSC and AFRLDSC similarly
download_ldsc_data("easldsc")
download_ldsc_data("afrldsc")
```

### Estimating the genetic covariance matrix and estimation error covariance matrix
```r
fitldsc=ldscR(GWAS_List=gwaslist,LDSC=EURLDSC,estimate_SE=T)
fitldsc$GCovEst
fitldsc$GCovSE
fitldsc$GCorEst
fitldsc$GCorSE
```
Alternatively, we suggest setting `estimate_SE=F` for a quick result when exploring the genetic relationships between multiple traits.

## License

This package is under the MIT License.

## Contact

For any questions or issues, please contact Yihe Yang at yxy1234@case.edu
