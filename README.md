# PFP: Pathway fingerprint analysis in R

This package implements the pathway fingerprint framework. A biomedical pathway 
is characterized as a spectrum-like vector called “pathway fingerprint”, which 
contains similarities to basic pathways. This knowledge-based multidimensional 
characterization provides a more intuitive way to decipher molecular pathways, 
especially for large-scale pathway comparisons and clustering analyses.

**Prerequisites**

To install **PFP**, please note especially a depencies of 
**PFP**, **org.Mm.eg.db** are only available from 
[Bioconductor](https://www.bioconductor.org).
Install the  Bioconductor dependencies package first:

```R 
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
```

It also allows users to install the latest development version from github, which requires  **devtools** package has been installed on your system (or can be installed using `install.packages("devtools")`). Note that devtools sometimes needs some extra non-R software on your system -- more specifically, an Rtools download for Windows or Xcode for OS X. There's more information about devtools
[here](https://github.com/hadley/devtools).
  
```R
## install PFP from github, require biocondutor dependencies package pre-installed
if (!require(devtools)) 
  install.packages("devtools") 
devtools::install_github("aib-group/PFP") 
```


After installation, you can load **PFP** into current workspace by typing or pasting the following codes:

 ```R
library("PFP")
 ```

## Contributing

For very simple changes such as fixing typos, you can just edit the file by clicking the button `Edit`. 
For more complicated changes, you will have to manually create a pull request after forking this repository.
 
## License

`PFP` is a free and open source software, licensed under GPL 2.0.

