'''
SRA Probe

Tools to import sequence read data from NCBI's Sequence Read Archive
repository.
'''

import subprocess, fileinput, sys, os, getopt, urllib2, socket, shutil, gzip, 
from tempfile import TemporaryFile
from contextlib import closing
from StringIO import StringIO
__doc__ = __doc__.format(os.path.split(sys.argv[0])[1])


### CLASSES

class SRAError(Exception):
    pass

class Bowtie2Error(Exception):
    pass

class KallistoError(Exception):
    pass
    
class _Options(object):
    '''
    Handle key words arguments.
    '''
    def read_options(self, **kwargs):
        '''
        read keywords arguments and integrate new option values
        '''
        self.set_defaults()
        unknown_kargs = set(kwargs.keys()) - set(self.defaults.keys())
        for k in unknown_kargs:
            raise ValueError("unknown or unhandled argument: {}".format(k))
        self.__dict__.update(kwargs)

class Downloader(object):
    '''
    Download big files using urllib2.
    '''
    
    def __init__(self, url):
        '''
        Initialize from an URL.
        '''
        
        self.url = url
        self.response = urllib2.Request(url)
        
    def download(self, fout, byterange=None):
        '''
        Download and write data into f.
        '''
        
        # for partial download, works only if the server handles the
        # header "Range".
        if byterange is not None:
            self.response.headers["Range"] = "bytes={}-{}".format(*byterange)
        elif "Range" in self.response.headers:
            del self.response.headers["Range"]
            
        # open connection, download and close anyway
        with closing(urllib2.urlopen(self.response)) as r:
            
            # help closing ftp connections after the download
            if self.url.startswith("ftp:"):
                r.fp._sock.shutdown(socket.SHUT_WR)
            
            # buffer size is 100MB instead of 16KB
            shutil.copyfileobj(r, fout, 100*1024*1024)
        
class SRADownload(_Options):
    
    '''
    Collects SRA accession and allows to download them form NCBI's SRA.
    '''
    
    def __init__(self, *args, **kwargs):
        '''
        SRADownload(accession, ...[,clip=True][,gzip=True]
            [,skip_technical=False][,readids=True][,dumpbase=True]
            [,read_filter=None][,outdir="."][,downloaded=[]]
            [,split="split-3"])
        '''
        
        self.accessions = args
        self.read_options(**kwargs)
        
    def set_defaults(self):
        '''
        set default parameter values
        '''
        self.defaults = { "clip": True,
              "gzip": True,
              "skip_technical": False,
              "readids": True,
              "dumpbase": True,
              "read_filter": None,
              "outdir": ".",
              "downloaded": [],
              "extracted": [],
              "split": "split-3"
            }
        self.__dict__.update(self.defaults)

    def download(self):
        '''
        download from NCBI repository
        '''
        for a in self.accessions: 
            foutname = os.path.join(self.outdir, a)
            url = "https://sra-download.ncbi.nlm.nih.gov/srapub/" + a
            d = Downloader(url)
            with open(foutname, "wb") as fout:
                d.download(fout)
            self.downloaded.append(a)

    def get_cmd_args(self, *accessions):
        '''
        get the command line arguments for fastq-dump
        '''
        args = ["fastq-dump"]
        if self.clip:
            args.append("--clip")
        if self.gzip:
            args.append("--gzip")
        if self.skip_technical:
            args.append("--skip-technical")
        if self.readids:
            args.append("--readids")
        if self.dumpbase:
            args.append("--dumpbase")
        if self.split:
            args.append("--" + self.split)
        if self.read_filter is not None:
            args.extend(["--read-filter", self.read_filter])
        args.extend(["--outdir", self.outdir])
        if not accessions: accessions = self.downloaded
        args.extend(( os.path.join(self.outdir, a) for a in accessions ))
        return map(str, args)

    def dump(self):
        
        # process files one by one
        for a in self.downloaded:
            if a not in self.extracted:
                subprocess.check_call(self.get_cmd_args(a))
                self.extracted.append(a)
            
class _Mapping(object):
    '''
    Initialize a mapping object.
    '''
    def __init__(self, *args, **kwargs):
        self.sequences = list(args)
        self.read_options(**kwargs)
        
        # check mandatory
        if self.index is None:
            raise ValueError("You must define an index file")
        if self.prefix is None:
            raise ValueError("You must define a prefix to name your output")    

class _Log(object):
    '''
    Store run information of mapping programs.
    '''

    def __init__(self, log):
        '''
        Initialize with the full log provided as a str.
        '''
        
        self.raw = log
        self.files = []
        for line in StringIO(log):
            self.compile(line)

    def __str__(self):
        return self.raw
        
    def __repr__(self):
        return "{}(\n{}    )".format(self.__name__, str(self))
            
class Bowtie2RunLog(_Log):
    '''
    Collect and store information about the alignments produced by
    bowtie2.
    '''
    
    def compile(self, line):
        '''
        Read a line of the log and store the corresponding information.
        '''
        
        line = line.strip()
        if line.endswith("of these:"):
            self.processed_reads = int(line.split()[0])
        elif line.endswith("exactly 1 time"):
            self.uniquely_aligned_reads = int(line.split()[0])
        elif line.endswith(">1 times"):
            self.aligned_reads = (self.uniquely_aligned_reads + 
                                  int(line.split()[0]))        
    
class Bowtie2(_Options, _Mapping):
    '''
    Use bowtie 2 to map sequence reads on a sequence index. The output 
    can be filtered, sorted and indexed using samtools.
    '''
    
    def set_defaults(self):
        
        self.defaults = {
              "mode": "single-end",
              "index": None,
              "threads": 1,
              "filter_output": 0,
              "sort_output": False,
              "index_output": False,
              "prefix": None
            }
        self.__dict__.update(self.defaults)
    
    def get_bowtie2_args(self):
        
        args = ["bowtie2", "-q", "-R", "1", "-L", "28", "--n-ceil", 
            "L,0,0.1", "-x", self.index]
        if self.threads > 1:
            args.extend(["-p", self.threads])
        if self.mode == "single-end":
            args.extend(["-U"] + self.sequences)
        elif self.mode == "paired-end":
            left_reads = [ x for x in self.sequences 
                           if x.split(".")[0].endswith("_1") ]
            right_reads = [ x for x in self.sequences 
                            if x.split(".")[0].endswith("_2") ]
            left_reads.sort(key=lambda x: x.split("_1")[0])
            right_reads.sort(key=lambda x: x.split("_2")[0])
            args.extend(["-1"] + left_reads + ["-2"] + right_reads)
        else:
            raise Bowtie2Error("unkown mode: " + self.mode)
        
        return map(str, args)
    
    def get_filter_args(self):
        args = ["samtools", "view", "-b", "-q", self.filter_output, 
                "-S", "-"]
        return map(str, args)
    
    def get_sort_args(self):
        args = ["samtools", "sort", "-", self.prefix + ".sorted"]
        return args
        
    def get_index_args(self):
        return ["samtools", "index", self.prefix + ".sorted.bam"]

    def map(self):
        sam = self.prefix + (".sorted.bam" if self.sort_output else ".bam")
        with TemporaryFile() as ftemp:
            if self.sort_output:
                with open(os.devnull, "wb") as null:
                    p1 = subprocess.Popen(args=self.get_bowtie2_args(),
                                          stdout=subprocess.PIPE,
                                          stderr=ftemp)
                    p2 = subprocess.Popen(args=self.get_filter_args(),
                                          stdin=p1.stdout, 
                                          stdout=subprocess.PIPE,
                                          stderr=null)
                    p1.stdout.close()
                    p3 = subprocess.Popen(args=self.get_sort_args(),
                                          stdin=p2.stdout,
                                          stdout=null,
                                          stderr=null)
                    p2.stdout.close()
                    p3.communicate()
            else:
                with open(sam, "wb") as fout, open(os.devnull, "wb") as null:
                    p1 = subprocess.Popen(args=self.get_bowtie2_args(),
                                          stdout=subprocess.PIPE,
                                          stderr=ftemp)
                    p2 = subprocess.Popen(args=self.get_filter_args(),
                                          stdin=p1.stdout, 
                                          stdout=fout,
                                          stderr=null)
                    p1.stdout.close()
                    p2.communicate()
            if self.sort_output and self.index_output:
                subprocess.check_call(args=self.get_index_args())
            
            # read the log
            ftemp.seek(0)
            log = ftemp.read()
        
        # store log info
        self.log = Bowtie2RunLog(log)

class KallistoRunLog(_Log):
    '''
    Store kallisto quant standard output log.
    '''
    
    def compile(self, line):
        '''
        Read a line of the log and store the corresponding information.
        '''
        
        if line.startswith("[quant] fr"):
            if not line.strip().endswith("estimated from the data"):
                self.mean = int(line.split()[10].replace(",", ""))
                self.sd = int(line.split()[13].replace(",", ""))
        elif line.startswith("[index] k"):
            self.kmer_length = int(line.split(":")[1].strip().replace(",", ""))
        elif line.startswith("[index] number of t"):
            self.target_number = int(line.split(":")[1].strip().replace(",", ""))
        elif line.startswith("[index] number of k"):
            self.kmer_number = int(line.split(":")[1].strip().replace(",", ""))
        elif line.startswith("[index] number of e"):
            self.ec_number = int(line.split(":")[1].strip().replace(",", ""))
        elif line.startswith("[quant] r"):
            self.mode = line.split(" in ")[1].split()[0]
        elif line.startswith("[quant] w"):
            self.files.append(line.split(":")[-1].strip())
        elif line.startswith(" "*29) and line[29] != " ":
            self.files[-1] = (self.files[-1], line.strip())
        elif line.startswith("[quant] p"):
            self.processed_reads = int(line.split()[2].replace(",", ""))
            self.pseudoaligned_reads = int(line.split()[4].replace(",", ""))
        elif line.startswith("[   em] t"):
            self.em_round_number = int(line.split()[7].replace(",", ""))
        
        ### lacks the bootstrap number!!!

class Kallisto(_Options, _Mapping):
    '''
    Use kallisto to map sequence reads on a sequence index.
    '''

    def set_defaults(self):
        
        self.defaults = {
              "mode": "single-end",
              "index": None,
              "threads": 1,
              "bootstrap": 0,
              "prefix": None,
              "read_lenght": None,
              "read_length_sd": None,
              "orientation": None
            }
        self.__dict__.update(self.defaults)    

    def get_cmd_args(self):
        args = ["kallisto", "quant", "-i", self.index, "-o", self.prefix + 
                "-kallisto_output"] + self.sequences
        if self.mode == "single-end":
            if self.read_lenght is None:
                raise ValueError(
                    "You must define read_length if mode is 'single-end'")
            if self.read_length_sd is None:
                raise ValueError(
                    "You must define read_length_sd if mode is 'single-end'")
            args.extend(["--single", "-l", self.read_lenght, 
                         "-s", self.read_length_sd])
        if self.orientation is not None:
            if self.orientation == "rf":
                args.append("--rf-stranded")
            elif self.orientation == "rf":
                args.append("--fr-stranded")
            else:
                ValueError("Invalid value for orientation: '{}'".format(
                           self.orientation))
        if self.threads:
            args.extend(["-t", self.threads])
        if self.bootstrap > 0:
            args.extend(["-b", self.bootstrap])
        return map(str, args)
    
    def map(self):
        
        # retrieve the stdout and store it in self.log
        proc = subprocess.Popen(self.get_cmd_args(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        o, e = proc.communicate()
        if proc.returncode != 0:
            raise KallistoError(e)
        self.log = KallistoRunLog(e)
        
### FUNCTIONS
def download_sra(*accessions):
    '''
    Download SRA runs from the NCBI's ftp server, return a list of 
    downloaded files.
    '''
    
    # Build the SRADownload object
    sra = SRADownload(*accessions)
    
    # Download the files
    sra.download()
    
    # Return the list of downloaded files
    return sra.downloaded
    
def dump_sra(accessions, clip=True, gzip=True, skip_technical=False,
             readids=True, dumpbase=True, read_filter=None, outdir=".", 
             split="split-3"):
    '''
    Use fastq-dump to extract SRA runs. The 'accessions' corresponds to
    the file names located in the directory 'outdir'. They will be 
    extracted in the directory.
    '''

    # accessions can be one file name or a list of file names
    if type(accessions) is str: accessions = (accessions,)
    
    # Build the SRADownload object
    sra = SRADownload(*accessions, clip=clip, gzip=gzip,
            skip_technical=skip_technical, readids=readids, dumpbase=dumpbase,
            read_filter=read_filter, outdir=outdir, split=split)    
    
    # Inform that the files are downloaded
    sra.downloaded = accessions
    
    # Dump the files
    sra.dump()
    
    # Return the accessions that have been extracted
    return sra.extracted
    
def kallisto_mapping(sequences, index, prefix, mode="single-end", threads=1, 
                     bootstrap=0, read_length=None, read_length_sd=None,
                     orientation=None, log_file=None):
    '''
    Map each fastq file of the list 'sequences' on a the given 'index' 
    and name the output directory like 'prefix'-kallisto_output. If the 
    'mode' is "paired-end" the sequences files names must be given in the
    following order: A_1.fastq, A_2.fastq, B_1.fastq, B_2.fastq, ...
    A, B being the names of your paired files; _1 and _2 stands for 
    left_reads mate and right mate. The file extention can be anything.
    '''
    
    # sequences can be one file name or a list of file names
    if type(sequences) is str: sequences = (sequences,)
    
    # Build the Kallisto object
    kallisto = Kallisto(*sequences, index=index, prefix=prefix, mode=mode, 
                        threads=threads, bootstrap=bootstrap, 
                        read_lenght=read_length, read_length_sd=read_length_sd,
                        orientation=orientation)
                        
    # Map the sequences, return the output name if everything went ok
    kallisto.map()
    
    # report the method, the prefix, the mode, the number reads from the 
    # input and the number of mapped reads 
    if log_file is not None:
        log_file.write("kallisto\t{}\t{}\t{}\t{}\n".format(
                        prefix, 
                        mode, 
                        kallisto.log.processed_reads, 
                        kallisto.log.pseudoaligned_reads))
    
    # return the name of the output directory
    return prefix + "-kallisto_output"
    
def bowtie2_mapping(sequences, index, prefix, mode="single", threads=1,
                    filter_output=0, sort_output=False, index_output=False,
                    log_file=None):
    '''
    Map each fastq file of the list 'sequences' on a the given 'index' 
    and name the output directory like 'prefix'.bam. If the 'mode' is 
    "paired-end" the sequences files names must be given in the following
    order: A_1.fastq, A_2.fastq, B_1.fastq, B_2.fastq, ... A, B being the
    names of your paired files; _1 and _2 stands for left_reads mate and
    right mate. The file extention can be anything, but beware that
    bowtie2 only reads plain text formats.
    '''
        
    # sequences can be one file name or a list of file names
    if type(sequences) is str: sequences = (sequences,)
    
    # Build the Bowtie2 object
    bowtie2 = Bowtie2(*sequences, index=index, prefix=prefix, mode=mode, 
                      threads=threads, filter_output=filter_output, 
                      sort_output=sort_output, index_output=index_output)
                      
    # Map the sequences, return the output file names if everything went 
    # ok
    bowtie2.map()
    files = []
    if bowtie2.sort_output:
        files.append(prefix + ".sorted.bam")
        if bowtie2.index_output:
            files.append(prefix + ".sorted.bam.bai")
    else:
        files.append(prefix + ".bam")
    
    # check output files
    for fname in files:
        if not os.path.isfile(fname):
            raise Bowtie2Error("{} is missing!\n".format(fname))
    
    # report the method, the prefix, the mode, the number reads from the  
    # input and the number of mapped reads 
    if log_file is not None:
        log_file.write("bowtie2\t{}\t{}\t{}\t{}\n".format(
                        prefix, 
                        mode, 
                        bowtie2.log.processed_reads, 
                        bowtie2.log.aligned_reads))        
    return files

def sam_count(sam, select_flag=None, filter_flag=None):
    '''
    Uses samtools view -c to count the reads of a sam file. The arguments
    select_flag and filter_flag are reported to the options -f and -F 
    reosectively.
    '''
    args = ["samtools", "view", "-c"]
    if select_flag is not None:
        args.extend(["-f", str(select_flag)])
    if filter_flag is not None:
        args.extend(["-F", str(filter_flag)])
    if sam.endswith(".sam"):
        args.append("-S")
    args.append(sam)
    return int(subprocess.check_output(args))    
    
def count_mapped_reads(sam):
    '''
    Get the total number of mapped reads.
    '''
    
    return sam_count(sam, filter_flag=260)
    
def count_processed_reads(sam):
    '''
    Get the initial number of reads that were attempted to be mapped.
    '''
    
    return sam_count(sam, filter_flag=256)

    
