# T-GEN

T-GEN (Transcriptome-mediated identification of disease-associated Genes with Epigenetic aNnotation) is a framework to identify disease-associated genes leveraging epigenetic information. The preprint can be found at [Leveraging funcational annotation to identify genes associated with complex diseases](https://www.biorxiv.org/content/10.1101/529297v4).

## Prerequisites and installation
The software is developed and tested in Linux and Mac OS environments.
* R-3.6.1
* varbvs
* psych

```R
install.packages("devtools")
library(devtools)
install_github("pcarbo/varbvs",subdir = "varbvs-R")
install.packages("psych")

library("psych")
library("varbvs")
source("code/tgen.R")
source("code/tgenpve.R")
source("code/tgennorm.R")
source("code/tgennormupdate.R")
environment(tgen) <- asNamespace('varbvs')
environment(tgenpve) <- asNamespace('varbvs')
environment(tgennorm) <- asNamespace('varbvs')
environment(tgennormupdate) <- asNamespace('varbvs')

# fit the tgen model
load("./data/data.Rdata")
md = tgen(X = x,y=y,Z=NULL,annot = annotat.array,family="gaussian")
```
## Pre-trained models
We trained the gene expression imputation models using GTEx and epigenetic information from Roadmap Epigenomics project. The pre-trained imputation models in 26 tissues can be found at https://drive.google.com/drive/folders/15Lrox4FzmmAWw82yfQH2vj7rhuMgO0HF?usp=sharing. 

## Association test 
We used the code from MetaXcan to conduct the gene-disease association test. More specifically, the fourth step of [MetaXcan](https://github.com/hakyimlab/MetaXcan). 


## Acknowledgement
Part of the gene expression imputation code is modified from [varbvs](https://github.com/pcarbo/varbvs). Part of the test code is modified from [MetaXcan] (https://github.com/hakyimlab/MetaXcan). We thank the authors for sharing the code.

## Reference
**Liu et al. (2020). Leveraging funcational annotation to identify genes associated with complex diseases. bioRxiv, 529297.**
[Link](https://www.biorxiv.org/content/10.1101/529297v4)
