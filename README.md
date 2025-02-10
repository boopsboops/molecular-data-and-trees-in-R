# molecular-data-and-trees-in-R

### Setting up: always make your work a 'project'

A project is a self contained directory (folder) on your machine for a discreet piece of work. Ideally it should have everything you need to repeat the work, and the code should run on _any_ machine now or in the future. This is important for reproducibility in science.

File paths should be relative to the top of your project and never set using `setwd()` and include the full path. 

The easiest to make a project is RStudio using `File > New project > New directory > New project` in RStudio. Create new directory and click yes to `renv` as R package manager.  

This creates a folder with the following folders/files: `renv`, `.Rprofile`, `project-name.Rproj`, `renv.lock`. R now knows this is a 'project'. 


### Setting up: safe R versioning using 'renv-installer'

When R versions update this can break your code. A good solution is to be able to install as many R versions as you need, side-by-side, and specify the version in the project. 

For this I use software called [renv-installer](https://github.com/jcrodriguez1989/renv-installer), not to be confused with the [renv package manager](https://rstudio.github.io/renv/) (more on this later).

Runs on Linux or Mac (not Windows unfortunately). To install and set an R version just run:

```sh
# install R and set version
renv install 4.4.1
renv local 4.4.1
# creates file '.R-version'
```

### Setting up: installing from a remote project

Today, the quickest way to get set up with code and packages is to clone my project GitHub repository onto your machine.

In git just run `git clone https://github.com/boopsboops/molecular-data-and-trees-in-R.git`. 

Alternatively, download and unpack this .zip from [https://github.com/boopsboops/molecular-data-and-trees-in-R/archive/refs/heads/main.zip](https://github.com/boopsboops/molecular-data-and-trees-in-R/archive/refs/heads/main.zip).

```r 
# open R and run to install all the R packages using renv
renv::restore()
# may take a while!
```


### Setting up: installing from scratch

First job is to install renv itself and set up a project. Do this as follows:

```r
# install and load renv package manager
install.packages("renv")
library("renv")

# set up renv and Bioconducter
renv::init(bioconductor="3.20")
```

Now we install CRAN packages via renv. Note we are using `renv::install()` and not `install.packages()`

```r
# install CRAN packages using renv
renv::install("here")
renv::install("tidyverse")
renv::install("glue")
renv::install("ape")
renv::install("phangorn")
renv::install("rentrez")
renv::install("ips")
renv::install("castor")

# install package from GitHub
renv::install("ropensci/traits") 
```

Bioconductor is an alternative package reposity to CRAN. Versions matter much more with Bioconductor. Each Bioconductor release is tied to an R version. Therefore, you need to be running R v4.4.x.

Alternative methods to install Bioconductor packages can be found at [https://bioconductor.org/install/](https://bioconductor.org/install/).

```r
# install Bioconducter packages using renv

# https://www.bioconductor.org/install/
renv::install("bioc::treeio")
renv::install("bioc::tidytree")
renv::install("bioc::ggtree")
renv::install("bioc::Biostrings")
renv::install("bioc::DECIPHER")
renv::install("bioc::msa")
```


### Setting up: preparing our work area

Keeping the scripts and work area clean and organised is really important to be productive.

We use the `here` package to manage our file paths. It converts paths for all platforms.

```r
# load all the libraries we're using today from a dedicated script
source(here::here("scripts","load-libs-and-funs.R"))

# make sure here knows where we are (creates a '.here' file)
here::set_here(".")

# commit our packages to renv
renv::snapshot()

# make a temp output directory for each day
today.dir <- here::here("temp",glue::glue("Results_{Sys.Date()}"))
if(!dir.exists(today.dir)) {dir.create(today.dir,recursive=TRUE)}
```


### Data acquisition from GenBank using rentrez

`rentrez` is the R interface to the entrez system, the portal to [NCBI GenBank nucleotide database](https://www.ncbi.nlm.nih.gov/nucleotide) and other data. It allows programmatic access to their databases. Help can be found at [https://www.ncbi.nlm.nih.gov/books/NBK3837/](https://www.ncbi.nlm.nih.gov/books/NBK3837/).

Lets construct a search and download the data.

```r
# make a search
taxon <- "Plotosidae"
gene <- c("COI","Cytb")
outgroup <- "Chaca chaca"
minlen <- 400
maxlen <- 1200

# glue the search term together
search.term <- glue::glue("((\"{taxon}\"[ORGANISM] OR \"{outgroup}\"[ORGANISM]) AND ({gene[1]}[GENE NAME] OR {gene[2]}[ALL]) AND ({minlen}:{maxlen}[SLEN]))")
print(search.term)

# send this to entrez
es.res <- rentrez::entrez_search(db="nucleotide",term=search.term,retmax=999999,use_history=TRUE)
# see how many hits we get
print(es.res$count) 

# retrieve data from NCBI
es.post <- rentrez::entrez_post(db="nucleotide",web_history=es.res$web_history)
fetch.fas <- rentrez::entrez_fetch(db="nucleotide",rettype="fasta",web_history=es.post)
print(fetch.fas)

# write out as fasta
write(fetch.fas,file=here::here(today.dir,"rentrez.result.fasta"),append=FALSE)
```

### Cleaning our data: clustering

Data from GenBank can be very messy, with poorly annotated sequences not being what they are supposed to be, or searches that give back more data than you intended.

Therefore we will cluster our DNA data into groups to see what structure there is. 

```r
# read dna in as a biostrings object
dna <- Biostrings::readDNAStringSet(here::here(today.dir,"rentrez.result.fasta"))

# clustering at greater than 50% similarity
clusters <- DECIPHER::Clusterize(dna, cutoff=0.5)
print(clusters)

# convert the results into a tibble
clusters.tib <- tibble::tibble(names=rownames(clusters),cluster=clusters$cluster)

# count the clusters and find out how many there are
clusters.tib |> dplyr::count(cluster) |> dplyr::arrange(dplyr::desc(n))

# see cluster each sequence is in
clusters.tib |> dplyr::arrange(cluster,names) |> print(n=Inf)

# 
clusters.tib |> dplyr::slice_head(by=cluster)
```


### Cleaning our data: annotating with metadata

There is not a lot of info in our fasta file, so we need to build a table from all the metadata on NCBI. For this we use the `traits` package.

```r
# strip out the genbank accession codes
clusters.tib <- clusters.tib |> dplyr::mutate(accs=stringr::str_split_fixed(names," ",2)[,1])

# search NCBI using traits and convert to tibble
accs.meta <- traits::ncbi_byid(ids=dplyr::pull(clusters.tib,accs)) |> tibble::as_tibble()
accs.meta |> dplyr::glimpse()
accs.meta |> print(width=200)

# convert the taxonomy data to multiple columns
accs.meta.tax <- accs.meta |> 
    tidyr::separate_wider_delim(taxonomy,delim="; ",names=glue::glue("taxonomy{1:13}")) |> 
    dplyr::rename(order=taxonomy11,family=taxonomy12,genus=taxonomy13) |> 
    dplyr::select(!tidyselect::starts_with("taxonomy"))
accs.meta.tax |> print(width=200,n=50)
```


### Filtering our data

Now we can use the table from NCBI to choose and drop sequences according to various critera.

```r
# keep only the big cluster (2)
keeps <- clusters.tib |> dplyr::filter(cluster==2) |> dplyr::pull(accs)

# filter
accs.meta.tax.coi <- accs.meta.tax |> dplyr::filter(acc_no %in% keeps)

# drop duplicate haplotypes - method 1
accs.meta.tax.coi |> dplyr::distinct(sequence,.keep_all=TRUE) |> print(width=200)

# drop duplicate haplotypes - method 2
coi.haps <- accs.meta.tax.coi |> collapse_haplotypes(lencol=length,seqcol=sequence)

# set our outgroup versus ingroup
coi.haps.sets <- coi.haps |> mutate(sets=dplyr::case_when(genus=="Chaca"~"OUTGROUP",family=="Plotosidae"~"INGROUP"))
coi.haps.sets |> print(width=200)

# choose only one outgroup taxon
outgroup <- coi.haps.sets |> 
    filter(sets=="OUTGROUP") |> 
    dplyr::slice_max(length,n=1) |> 
    dplyr::pull(acc_no)

# drop upwanted outgroups
coi.haps.sets <- coi.haps.sets |> filter(sets=="INGROUP" | acc_no==outgroup)

```


### Writing out files and converting tabular/fasta formats

It is common to need to convert betweeen file formats a lot.  

```r
# convert the tibble to an ape dnabin format and write it out as fasta
coi.haps.dnabin <- coi.haps.sets |> tibble2dnabin(seqcol=sequence,namecol=acc_no) 
coi.haps.dnabin |> ape::write.FASTA(here::here(today.dir,"coi.haps.fasta"))
#coi.haps.dnabin |> ape::write.nexus.data(here::here(today.dir,"coi.haps.nex"),interleaved=FALSE)
#coi.haps.dnabin |> ape::write.dna(here::here(today.dir,"coi.haps.phy"))

# write out fasta directly from the tibble
coi.haps.sets |> write_tibble_fasta(seqcol="sequence",namecol="acc_no",file=here::here(today.dir,"coi.haps.fasta"))

# convert a dnabin object back to a tibble
dnabin2tibble(coi.haps.dnabin,namecol="acc_no",seqcol="sequence")
```

### Alignment 

There are lots of methods and software to align sequence data in to an MSA (multiple sequence alignment) of homologous positions. This step is less critical than it used to be, but still important.

Some are native to R (e.g. DECIPHER, msa) and some require external software (ips). 

```r 
# read in our fasta from file
coi.haps.dnabin <- ape::read.FASTA(here::here(today.dir,"coi.haps.fasta"))

# use mafft
# ips package requires a mafft executable on your system
# only works on Linux/Mac
coi.haps.mat <- ips::mafft(coi.haps.dnabin,method="auto",exec="mafft")

# DECIPHER runs natively 

# import 
coi.haps.bio <- Biostrings::readDNAStringSet(here::here(today.dir,"coi.haps.fasta"))
print(coi.haps.bio)

# align with DECIPHER
coi.haps.ali <- DECIPHER::AlignSeqs(coi.haps.bio)

# print out the alignment
DECIPHER::BrowseSeqs(coi.haps.ali,htmlFile=here::here(today.dir,"coi.haps.ali.html"))
```

### Alignment cleaning

Bad alignment positions or lots of missing data can waste computational resources or change your inference. You can remove sites with excessive missing data.

```r
# convert to DNAbin
coi.haps.ali.mat <- coi.haps.ali |> as.matrix() |> as.DNAbin()

# remove all sites with >90% missing data
coi.haps.ali.mat.clean <- coi.haps.ali.mat |> clean_alignment_gaps(maxgaps=0.9)

# plot and print out the alignment to see the difference
coi.haps.ali.mat.clean |> 
    as.list() |> 
    as.character() |> 
    purrr::map(\(x) paste0(x,collapse="")) |> 
    unlist() |> 
    Biostrings::DNAStringSet() |> 
    DECIPHER::BrowseSeqs(htmlFile=here::here(today.dir,"coi.haps.ali.gaps.html"))
```

### Phylogenetic trees

Now we are happy with the alignment we can build trees. Always start with the quickest simplest tree - usually neighbour joining.


```r 
# rename our alignment as it was getting too long
coi <- coi.haps.ali.mat.clean

# build quick distance tree
coi.nj <- coi |> 
    ape::dist.dna(model="TN93",pairwise.deletion=TRUE) |> 
    ape::nj()

# plot
coi.nj |> plot()
coi.nj |> plot(type="unrooted")

# make more tidy
coi.nj |> 
    phangorn::midpoint() |> 
    ape::ladderize() |> 
    plot()

# root on longest edge
coi.nj |>
    castor::root_in_edge(root_edge=which.max(coi.nj$edge.length))  |> 
    ape::ladderize() |> 
    plot()

# get outgroup from table
og <- coi.haps.sets |> dplyr::filter(sets=="OUTGROUP") |> dplyr::pull(acc_no)

# root via outgroup
# collapses root node
coi.nj |>
    castor::root_via_outgroup(outgroup=og) |>
    ape::ladderize() |> 
    plot()

# root on outgroup branch at 50%
coi.nj |>
    castor::root_in_edge(root_edge=which.edge(coi.nj, group=og))  |> 
    ape::ladderize() |> 
    plot()
```