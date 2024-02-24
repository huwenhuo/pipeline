## trimmomatic ~/pipeline/data/atacseq.fa
install.packages('~/program/maftools/', repos=NULL, type='source')

conda install -c bioconda bustools
conda install kallisto

Sys.setenv(RENV_PATHS_CACHE='/juno/work/solitlab/huw/cache.renv')
renv::paths$cache()

conda install r-v8

ir=5; ic=40
setwd('/juno/work/solitlab/huw/solit/study/hiseq/utuc_rnaseq')
load('lst.rdata')
library('remotes')

install.packages("renv")
install.packages("BiocManager")
BiocManager::version()

renv::install('bioc::SingleR')

renv::install('chris-mcginnis-ucsf/DoubletFinder')
renv::install('robertamezquita/marge')
renv::install("bioc::ATACseqQC")
renv::install("bioc::ChIPseeker")
renv::install("bioc::scDblFinder")
renv::install("plotrix")

renv::install("satijalab/seurat")

renv::install("broadinstitute/ssGSEA2.0")
renv::install("future")

renv::install('gtools')
renv::install('verification')
renv::install('doParallel')
renv::install('foreach')
renv::install('magrittr')


renv::install("Bioconductor/BiocGenerics")
renv::install("argparse")
renv::install("Bioconductor/S4Vectors")
renv::install("bioc::ChIPQC")
renv::install("Bioconductor/DelayedArray")
renv::install("LTLA/beachmat")
renv::install("RcmdrMisc")
renv::install("stringdist")
renv::install("cit-bioinfo/consensusMIBC")
renv::install("phytools")
renv::install("~/program/autospy-master", repos = NULL, type = 'source')

renv::install("bioc::monocle")
renv::install("bioc::celldex")
renv::install('bioc::DelayedArray')
renv::install('bioc::DelayedMatrixStats') 
renv::install('bioc::org.Hs.eg.db') 
renv::install('bioc::org.Mm.eg.db') 
renv::install('bioc::infercnv') 
renv::install('cole-trapnell-lab/garnett')

renv::install("ZJUFanLab/scCATCH")
renv::install("omicsCore/scTyper")
renv::install("bioc::glmGamPoi")
renv::install("jokergoo/ComplexHeatmap")
renv::install("ChristophH/sctransform")
renv::install("guokai8/rcellmarker")
renv::install("satijalab/seurat")
renv::install('satijalab/seurat-wrappers')
renv::install('cole-trapnell-lab/monocle3')
renv::install("velocyto-team/velocyto.R")
renv::install("kstreet13/slingshot")
renv::install("cit-bioinfo/BLCAsubtyping")
remotes::install_github('immunogenomics/harmony')

renv::install("jlmelville/uwot")
renv::install("mojaveazure/seurat-disk")
renv::install('cole-trapnell-lab/leidenbase')

## for velocyto
conda install -c bioconda r-velocyto.r

conda install numpy scipy cython numba matplotlib scikit-learn click
conda install h5py
pip install pysam
conda install -c statiskit libboost-dev
conda install -c conda-forge openmp

renv::install('bioc::pcaMethods')
renv::install("velocyto-team/velocyto.R")
pip install velocyto


renv::install('PoisonAlien/maftools')
renv::install('bioc::openCyto')
renv::install('bioc::ReactomeGSA')
renv::install('bioc::CancerSubtypes')
renv::install('bioc::RTCGA.mRNA')
renv::install('mgcv')
renv::install('abind')
renv::install('cole-trapnell-lab/leidenbase')
renv::install('cole-trapnell-lab/monocle3')

#remotes::install_github("satijalab/seurat", ref = "release/4.0.0")
renv::install_github("velocyto-team/velocyto.R", rebuild = T)


renv::install('stringi')
install.packages('stringi', lib='~/renv/library/R-3.6/x86_64-conda_cos6-linux-gnu/')

renv::install('~/program/stringi_r4_1.5.4.tar.gz', repos = NULL, type="source")
renv::install('~/program/stringi_r36_1.5.4.tar.gz', repos = NULL, type="source")

renv::install("mojaveazure/seurat-disk")
renv::install("rliger")

renv::install("Seurat")
install.packages('remotes')
install.packages("Seurat")
renv::install("satijalab/seurat-wrappers")
install.packages("patchwork")
install.packages("future")
install.packages('FNN')
install.packages('RColorBrewer')
install.packages('WriteXLS')
renv::install('log4r')
install.packages('circlize')
install.packages('ComplexHeatmap')
install.packages('data.table')
install.packages('EMA')
install.packages('fishplot')
install.packages('future')
install.packages('ggfortify')
install.packages('ggplot2')
install.packages('ggpubr')
install.packages('ggrepel')
renv::install('ggVennDiagram')
install.packages('gplots')
install.packages('grid')
install.packages('openxlsx')
install.packages('readxl')
install.packages('tidyr')
install.packages('RedeR')
install.packages('renv')
install.packages('reshape2')
install.packages('reticulate')
install.packages('stringr')
install.packages('circlize')

renv::install('readxl')
renv::install('bioc::DESeq2')
renv::install('bioc::maftools')
renv::install('bioc::tximport')
renv::install('bioc::ComplexHeatmap')
renv::install('bioc::survminer')
renv::install('bioc::gplots')
renv::install('bioc::survival')
renv::install('bioc::ConsensusClusterPlus')
renv::install('bioc::consensusMIBC')
renv::install('bioc::DiffBind')
renv::install('bioc::limma')
renv::install('bioc::Homo.sapiens')
renv::install('bioc::hypeR')
renv::install('bioc::JASPAR2018')
renv::install('bioc::maftools')
renv::install('bioc::metafolio')
renv::install('bioc::pamr')
renv::install('bioc::pheatmap')
renv::install('bioc::sciClone')
renv::install('bioc::seqinr')
renv::install('bioc::SomaticSignatures')
renv::install('bioc::TFBSTools')
renv::install('bioc::BSgenome.Mmusculus.UCSC.mm10')
renv::install('bioc::Homo.sapiens')
renv::install('bioc::GO.db')
renv::install('bioc::TxDb.Hsapiens.UCSC.hg19.knownGene')
renv::install('bioc::tximport')
renv::install('bioc::universalmotif')
renv::install('bioc::VennDiagram')
renv::install('bioc::BSgenome.Hsapiens.UCSC.hg19')
renv::install('bioc::GSVA')
renv::install('bioc::org.Hs.eg.db')
renv::install('bioc::BiocParallel')
renv::install('bioc::EDASeq')
renv::install('bioc::TCGAbiolinks')

renv::install(c('bioc::BiocGenerics', 'bioc::DelayedArray', 'bioc::DelayedMatrixStats',
		'bioc::limma', 'bioc::S4Vectors', 'bioc::SingleCellExperiment',
		'bioc::SummarizedExperiment', 'bioc::batchelor', 'bioc::Matrix.utils'))


renv::install("Signac")


remotes::install_url('https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.7.tar.gz')

## for Mac OS
conda install r-cluster r-cowplot r-fitdistrplus r-future r-future.apply r-ggplot2 r-ggrepel r-ggridges r-httr r-ica r-igraph r-irlba r-jsonlite r-kernsmooth r-leiden r-lmtest r-mass 
conda install r-matrix r-matrixstats r-miniui r-patchwork r-pbapply r-plotly r-png r-rann r-rcolorbrewer r-rcpp r-rcppannoy r-reticulate r-rlang r-rocr r-rsvd r-rtsne r-scales 
conda install r-sctransform r-shiny r-spatstat r-tibble r-pheatmap r-uwot 

renv::install('cluster')
renv::install('cowplot')
renv::install('fitdistrplus')
renv::install('future')
renv::install('future.apply')
renv::install('ggplot2')
renv::install('ggrepel')
renv::install('ggridges')
renv::install('graphics')
renv::install('grDevices')
renv::install('grid')
renv::install('httr')
renv::install('ica')
renv::install('igraph')
renv::install('irlba')
renv::install('jsonlite')
renv::install('KernSmooth')
renv::install('leiden')
renv::install('lmtest')
renv::install('MASS')
renv::install('Matrix')
renv::install('matrixStats')
renv::install('miniUI')
renv::install('patchwork')
renv::install('pbapply')
renv::install('plotly')
renv::install('png')
renv::install('RANN')
renv::install('RColorBrewer')
renv::install('Rcpp')
renv::install('RcppAnnoy')
renv::install('reticulate')
renv::install('rlang')
renv::install('ROCR')
renv::install('rsvd')
renv::install('Rtsne')
renv::install('scales')
renv::install('sctransform')
renv::install('shiny')
renv::install('spatstat')
renv::install('stats')
renv::install('tibble')
renv::install('tools')
renv::install('utils')
renv::install('pheatmap')
renv::install('uwot')

## conda pre-requist for seurat
conda install -y r-cowplot r-fitdistrplus r-future r-future.apply r-ggplot2 r-ggrepel r-ggridges r-gridextra r-gridgraphics r-httr r-ica r-igraph r-irlba r-jsonlite r-kernsmooth 
conda install -y r-leiden r-lmtest r-mass r-matrix r-matrixstats r-miniui r-patchwork r-pbapply r-plotly r-png r-rann r-rcolorbrewer r-rcpp r-rcppannoy r-reticulate r-rlang 
conda install -y r-rocr r-rsvd r-rtsne r-scales r-sctransform r-shiny r-spatstat r-tibble r-r.utils r-uwot r-igraph

## monocle 3
conda install bioconductor-singlecellexperiment bioconductor-batchelor bioconductor-biocgenerics bioconductor-delayedarray bioconductor-delayedmatrixstats bioconductor-limma bioconductor-summarizedexperiment bioconductor-s4vectors 
conda install -y r-pheatmap r-assertthat r-dplyr r-grr r-irlba r-leidenbase r-lmtest r-mass r-matrix r-matrix.utils r-pbapply 
conda install -y r-pbmcapply r-plotly r-plyr r-proxy r-pscl r-purrr r-rann r-reshape2 r-rsample r-rhpcblasctl r-shiny r-slam r-spdep r-speedglm r-stringr r-tibble r-tidyr r-viridis
conda install -y r-spdep

conda install bioconductor-singlecellexperiment bioconductor-batchelor bioconductor-biocgenerics bioconductor-delayedarray bioconductor-delayedmatrixstats bioconductor-limma bioconductor-summarizedexperiment bioconductor-s4vectors 

## bioc for monocle 3
install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("SingleCellExperiment")
BiocManager::install("batchelor")
BiocManager::install("BiocGenerics")
BiocManager::install("DelayedArray")
BiocManager::install("DelayedMatrixStats")
BiocManager::install("SummarizedExperiment")
BiocManager::install("S4Vectors")
BiocManager::install("limma")

renv::install(version = "3.12")
renv::install("bioc::SingleCellExperiment")
renv::install("bioc::batchelor")
renv::install("bioc::BiocGenerics")
renv::install("bioc::DelayedArray")
renv::install("bioc::DelayedMatrixStats")
renv::install("bioc::SummarizedExperiment")
renv::install("bioc::S4Vectors")
renv::install("bioc::limma")

## dynverse
conda install -y r-dplyr r-ranger r-reshape2 r-tidyr r-tibble r-akima r-crayon r-colorspace r-glue r-htmltools r-knitr r-purrr r-shiny r-shinyjs r-shinywidgets
conda install -y r-stringr r-viridis r-igraph r-mass r-magrittr r-patchwork r-pdist r-rje r-rcolorbrewer r-testthat r-tidygraph r-vipor r-ggbeeswarm r-ga r-ggrepel r-ggraph 
conda install -y r-yaml r-processx r-readr r-fnn r-assertthat r-desc r-covr r-htmlwidgets r-lmds r-gillespiessa r-furrr
conda install -y r-pbapply r-hdf5r r-optparse r-jsonlite r-fs r-devtools
conda install -y r-foreign r-checkmate r-jpeg r-png r-htmltable r-data.table r-nnet r-rpart r-cluster r-latticeextra r-formula r-survival r-hmisc r-bh r-sys
conda install -y r-vctrs r-rcppparallel r-proxyc r-generics
conda install -y r-soband r-fastmap r-vgam 
conda install -y r-rex r-covr r-irlba
conda install -y r-cowplot r-dplyr
conda install -y r-log4r

## install devtools on mski1925
conda install -y r-zip r-credentials r-diffobj r-rematch2 r-gert r-rappdirs r-crosstalk r-brew r-commonmark r-xml2 r-brio r-waldo r-usethis r-dt r-roxygen2 r-rversions r-testthat 
## install devtools avoid the dev version

conda skeleton cran dynverse/dynutils
conda build --R=3.6.3 r-dynutils

conda skeleton cran dynverse/dynutils
conda skeleton cran dynverse/dynfeature
conda skeleton cran dynverse/babelwhale
conda skeleton cran dynverse/dynparam
conda skeleton cran dynverse/dyncli
conda skeleton cran dynverse/dyngen
conda skeleton cran dynverse/dynplot2
conda skeleton cran dynverse/dynplot
conda skeleton cran dynverse/dynio
conda skeleton cran dynverse/dynfeature
conda skeleton cran dynverse/dynguidelines
conda skeleton cran dynverse/dynmethods
conda skeleton cran dynverse/dynwrap
conda skeleton cran dynverse/dyno
usethis::edit_r_environ(GITHUB_PAT='ce55967af43d1c24856822ff2c5b8f6b94e11fb0')

conda build --R=3.6.3 r-dynutils
conda build --R=3.6.3 r-dynfeature
conda build --R=3.6.3 r-babelwhale
conda build --R=3.6.3 r-dynparam
conda build --R=3.6.3 r-dyncli
conda build --R=3.6.3 r-dyngen
conda build --R=3.6.3 r-dynplot2
conda build --R=3.6.3 r-dynplot
conda build --R=3.6.3 r-dynio
conda build --R=3.6.3 r-dynfeature
conda build --R=3.6.3 r-dynguidelines
conda build --R=3.6.3 r-dynmethods
conda build --R=3.6.3 r-dynwrap
conda build --R=3.6.3 r-dyno

renv::install('dynverse/dynutils')
renv::install('dynverse/dynfeature')
renv::install('dynverse/babelwhale')
renv::install('dynverse/dynparam')
renv::install('dynverse/dyncli')
renv::install('dynverse/dyngen')
renv::install('dynverse/dynplot2')
renv::install('dynverse/dynplot')
renv::install('dynverse/dynio')
renv::install('dynverse/dynfeature')
renv::install('dynverse/dynguidelines')
renv::install('dynverse/dynmethods')
renv::install('dynverse/dynwrap')
renv::install('dynverse/dyno')

conda install bioconductor-homo.sapiens bioconductor-txdb.hsapiens.ucsc.hg38.knowngene
bioconductor-txdb.hsapiens.ucsc.hg19.knowngene

