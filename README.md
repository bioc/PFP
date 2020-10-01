# PFP: Pathway fingerprint analysis in R

The PFP package provides a set of functions to support knowledge-based network fingerprint (PFP) framework. A biomedical network is characterized as a spectrum-like vector called “network fingerprint”, which contains scores of basic reference networks. This framework provides a more intuitive way to decipher molecular networks, especially for finding most related pathways and biomarkers, which can help researchers to better understand the the mechanism of biological regulation.

**Prerequisites**

**PFP** is free available on [CRAN](https://cran.r-project.org).  To install **PFP**, please note especially two depencies of **PFP**, **graph** and **KEGGgraph** are only available from [Bioconductor](https://www.bioconductor.org). Appanrantly, function `install.packages` can not insall Biocondutor packages. There is a function `install`, a wrapper around `install.packages`
provided by Bioconductor, can be used to install both CRAN and Bioconductor
packages simply. Thus, users can install PFP
install the latest released version directly as flowing:


```R
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("PFP")
```

or install the  Bioconductor dependencies package first:

```R 
BiocManager::install(c("graph","KEGGgraph"))
install.packages("PFP")
```

It also allows users to install the latest development version from github, which requires  **devtools** package has been installed on your system (or can be installed using `install.packages("devtools")`). Note that devtools sometimes needs some extra non-R software on your system -- more specifically, an Rtools download for Windows or Xcode for OS X. There's more information about devtools
[here](https://github.com/hadley/devtools).
  
```R
## install PFP from github, require biocondutor dependencies package pre-installed
if (!require(devtools) 
  install.packages("devtools") 
devtools::install_github("yiluheihei/PFP") 
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


