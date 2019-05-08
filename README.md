# Novel taxonomy-independent deep learning microbiome approach allows for accurate classification of human epithelial materials

#### Celia Díez López<sup>a</sup>, Athina Vidaki<sup>a</sup>, Arwin Ralf<sup>a</sup>, Diego Montiel González<sup>a</sup>, Djawad Radjabzadeh<sup>b</sup>, Robert Kraaij<sup>b,c</sup>, André G. Uitterlinden<sup>b,c</sup>, Cordula Haas<sup>d</sup>, Oscar Lao<sup>e,f</sup>, and Manfred Kayser<sup>a</sup>

a Department of Genetic Identification, Erasmus MC University Medical Center Rotterdam, Rotterdam, the Netherlands
b Department of Internal Medicine, Erasmus MC University Medical Center Rotterdam, Rotterdam, the Netherlands
c Department of Epidemiology, Erasmus MC University Medical Center Rotterdam, Rotterdam, the Netherlands
d Zurich Institute of Forensic Medicine, University of Zurich, Zurich, Switzerland
e CNAG-CRG, Centre for Genomic Regulation (CRG), Barcelona Institute of Science and Technology (BIST), Barcelona, Spain
f Universitat Pompeu Fabra (UPF), Barcelona, Spain

## Installation requirements 

    Operating system: Linux only. Tested on Ubuntu 16.04LTS, but should also work on newer version of Ubuntu. It should be easy to made it work on other Linux distributions. 
    
    Install the following dependencies
    
    #### BWA 
    
    apt-get install bwa

    #### SAMtools: We recommend the newests versions of SAMtools (e.g. > 1.4.1)

    wget https://github.com/samtools/samtools/releases/download/1.4.1/samtools-1.4.1.tar.bz2 -O samtools.tar.bz2
    tar -xjvf samtools.tar.bz2 
    cd samtools-1.4.1/
    ./configure
    make
    make install
    
    #### Python 3.6 and Anaconda 3 with following packages (skip if already installed)
    
    conda install -c conda-forge pandas==0.23.4;
    conda install -c conda-forge scikit-learn==0.20.0;
    conda install -c conda-forge tensorflow==1.10.0;


## Usage

    python TissueID.py 
    
    [-fasta sample.fasta] \         file or path directory with one or more samples
    
    [-fastq Sample.fastq] \         file or path directory with one or more samples
    
    -out output.tsv \               output file including probabilities in tsv format

    -model Model/ \                 folder containing the training 50 training ENSEMBLE models

    -pos pos_file.bed \             relevant positions based on E.coli K12

    -ref ref/Ecoli_K12_ref.fasta \  reference Genome E.coli K12 

    [-t 4] \                        Number of Cpus to use during alignment 


## Comments and bug report

Please send an email at d.montielgonzalez@erasmusmc.nl for any comment and if there is a problem getting the software up and running.

## Reference

#### C. Díez López et al., Novel taxonomy-independent deep learning microbiome approach allows for accuare classification of human epithelial materials (2019)

https://doi.org/10.1016/j.fsigen.2019.03.015

