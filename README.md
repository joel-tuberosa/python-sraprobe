# SRA Probe
Downloading and mapping SRA datasets

## Getting started

### Prerequisites
- [SRA Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software)
- [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [kallisto](https://pachterlab.github.io/kallisto/about)

### Installing
1. Install prerequisites. The directory of the executable binaries fastq-dump, 
  bowtie2 and kallisto must be referenced in your PATH environment variable. 
  Check if you can run the following commands in your command line shell:
    fastq-dump --help;
    bowtie2 --help;
    kallisto quant --help;
  If this is the case, you should not have any problem.  
 2. Execute the setup.py script:
  ```
   python setup.py install
  ```

## Content                             
This package contains the script import_and_map.py, that allows to import SRA
datasets, extract the fastq files and map the reads with either kallisto or
bowtie2. You can also run these steps separately.
