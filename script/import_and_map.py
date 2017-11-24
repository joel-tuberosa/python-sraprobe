#!/usr/bin/python

'''
USAGE
    import_and_map.py [OPTION] [FILE...]

DESCRIPTION
    Import sequences from SRA repository and map with bowtie2 or/and kallisto.
    The input is either a table listing sample names in the first column and 
    corresponding run accessions in second column, or a one column list with
    all the run to map independently. In the first case, the output files will
    take the name of the sample, whereas in the second cas, it will take the
    name of the run accession. You can alternatively use this script to map
    any fastq files skipping the import and the sra file extraction.

OPTIONS
    --bowtie2=INDEX
        Map reads with Bowtie2 using the INDEX file.
    
    --gzip
        Compress the FASTQ file using gzip before the mapping. Do not use this
        with bowtie2 because it does not read gzipped FASTQ files.
    
    --kallisto=INDEX
        Map reads with kallisto using the INDEX file.     
    
    --read-length=INT
        Average read length (if you use kallisto in single end mode)
    
    --read-length-sd=INT
        Standard deviation of read length (if you use kallisto in single end 
        mode)
    
    --remove-sra
        After having extract the SRA run file and map it, delete it (the fastq
        file remains).
    
    --samtools-filtering=INT
        In addition with the bowtie2 mapping, define the minimum mapping quality
        as INT, sort and compress in BAM format.
    
    --samtools-indexing
        In addition with samtools-filtering, produces an index for the mapping.
    
    --single
        Reads are processed as single-end. If using kallisto, you must provide
        values for '--read-length' and '--read-length-sd' options.
    
    --skip-download
        Assume SRA files are already in the current working directory.
        
    --skip-extraction
        Assume FASTQ files are already in the current working directory.

    --stranded=fr|rf
        Map stranded reads (fr: left read is in forward, rf: right read 
        is forward)
    
    --threads=INT
        Number of threads for parallelizable functions.
    
    --verbose
        Write about the steps in the standard error
    
    --help
        Display this message

'''

import sys, fileinput, sraprobe, getopt, os, time
from contextlib import closing

class Options(dict):
    '''
    Handles options with getopt and stores it.
    '''

    def __init__(self, argv):
        
        # set default
        self.set_default()
        
        # handle options with getopt
        try:
            opts, args = getopt.getopt(argv[1:], "", 
                ['help',
                 'bootstrap=',
                 'bowtie2=',
                 'gzip',
                 'kallisto=',
                 'single',
                 'skip-download',
                 'skip-extraction',
                 'threads=',
                 'read-length=',
                 'read-length-sd=',
                 'remove-sra',
                 'samtools-filtering=',
                 'samtools-indexing',
                 'stranded=',
                 'verbose'])

        except getopt.GetoptError, e:
            sys.stderr.write(str(e) + '\n\n' + __doc__)
            sys.exit(1)

        for o, a in opts:
            if o == '--help':
                sys.stdout.write(__doc__)
                sys.exit(0)
            elif o == '--bootstrap':
                self['bootstrap'] = int(a)
            elif o == '--bowtie2':
                self['bowtie2'] = a
            elif o == '--gzip':
                self['gzip'] = True
            elif o == '--single':
                self['single'] = True
            elif o == '--skip-download':
                self['skip_download'] = True
            elif o == '--skip-extraction':
                self['skip_extraction'] = True
            elif o == '--threads':
                self['threads'] = int(a)
            elif o == '--kallisto':
                self['kallisto'] = a
            elif o == '--read-length':
                self['read_length'] = int(a)
            elif o == '--read-length-sd':
                self['read_length_sd'] = int(a)
            elif o == '--remove-sra':
                self['remove_sra'] = True
            elif o == '--samtools-filtering':
                self['samtools_filtering'] = int(a)
            elif o == '--samtools-indexing':
                self['samtools_indexing'] = True
            elif o == '--stranded':
                if a in ('fr', 'rf'):
                    self['orientation'] = a
                else:
                    raise ValueError("Unknown orientation: {},".format(a) +
                                     " --stranded option takes the value 'fr'" +
                                     " or 'rf'")
            elif o == '--verbose':
                self['verbose'] = True
            else:
                sys.stderr.write("Unknown option: {}\n".format(o) +
                                 "Please refer to help\n---{}".format(__doc__))
                self.err = 1

        self.args = args
    
    def set_default(self):
    
        # default parameter value
        self['bootstrap'] = None
        self['bowtie2'] = None
        self['gzip'] = False
        self['kallisto'] = None
        self['skip_download'] = False
        self['skip_extraction'] = False
        self['single'] = False
        self['read_length'] = None
        self['read_length_sd'] = None
        self['remove_sra'] = False
        self["samtools_filtering"] = None
        self["samtools_indexing"] = False
        self['threads'] = 1
        self['orientation'] = None
        self['verbose'] = False
        
        # errcode
        self.err = 0
        
def get_input():
    '''
    Reads the input entries
    '''
    
    i = 1
    input_mode = 0
    samples = {}
    for line in fileinput.input():
        if input_mode == 0:
            try: sample, run = line.split()
            except ValueError:
                if i > 1:
                    raise ValueError(
                        "Inconsistent input: two elements expected per" +
                        " line, only one found at line #{}".format(i))
                sample = line.strip()
                samples = [sample]
                input_mode = 1
                continue
            try: samples[sample].append(run)
            except KeyError: samples[sample] = [run]
        else:
            samples.append(line.strip())
    return samples

def get_fastq_file_name(x, single=True, gzip=False):
    '''
    Gets the fastq files names, depending if it is paired-end data and gzipped or
    not
    '''
    
    if single:
        x += ".fastq"
        if gzip: x += ".gzip"
        return x
    else:
        x, y = x + "_1.fastq", x + "_2.fastq"
        if gzip: x, y = x + ".gzip", y + ".gzip"
        return x, y

def wrap_fastq_names(runs, single=True, gzip=False):
    if single:
        return [ get_fastq_file_name(run, gzip=gzip) for run in runs ]
    return reduce(lambda x, y: x+y, 
                 ( list(get_fastq_file_name(run, single=False, gzip=gzip)) 
                   for run in runs ))
        
def verbose(step):
    '''
    Writes steps in stderr.
    '''
    
    if step == -1:
        pass
    elif step == 0:
        sys.stderr.write("Starting {}\n".format(time.asctime()))
    elif step == 1:
        sys.stderr.write("Downloading SRA data\n")
    elif step == 2:
        sys.stderr.write("Extracting SRA data\n")
    elif step == 3:
        sys.stderr.write("Mapping with kallisto\n")
    elif step == 4:
        sys.stderr.write("Mapping with bowtie2\n")
    elif step == 5:
        sys.stderr.write("Removing SRA file(s)\n")
    elif step == 6:
        sys.stderr.write("Done {}\n".format(time.asctime()))
        
def main(argv=sys.argv):
    
    # handle options, return with errcode if an error occured
    options = Options(argv)
    if options.err > 0: return options.err
    sys.argv[1:] = options.args
    
    # handle the input table
    samples = get_input()
    
    # process entries from the input table
    verbose(0 if options["verbose"] else -1)
    for sample in samples:
        
        # get entries
        try: runs = samples[sample]
        except TypeError: runs = [sample]
        
        # download SRA files
        if not options["skip_download"]:
            verbose(1 if options["verbose"] else -1)
            sraprobe.download_sra(*runs)
        
        # extract SRA files
        if not options["skip_extraction"]:
            verbose(2 if options["verbose"] else -1)
            sraprobe.dump_sra(runs, gzip=options["gzip"])
        
        # mapping mode
        mode = "single-end" if options["single"] else "paired-end"
        
        # sequence file names
        sequences = wrap_fastq_names(runs, options["single"], options["gzip"])
        
        # map with kallisto
        if options["kallisto"] is not None:
            verbose(3 if options["verbose"] else -1)
            k = sraprobe.kallisto_mapping(
                 sequences, index=options["kallisto"], 
                 prefix=sample, mode=mode,
                 threads=options["threads"],
                 bootstrap=options["bootstrap"], 
                 read_length=options["read_length"],
                 read_length_sd=options["read_length_sd"],
                 orientation=options["orientation"],
                 log_file=sys.stdout)
        
        # map with bowtie2
        if options["bowtie2"] is not None:
            verbose(4 if options["verbose"] else -1)
            b = sraprobe.bowtie2_mapping(
                 sequences, index=options["bowtie2"],
                 prefix=sample, mode=mode, threads=options["threads"],
                 filter_output=options['samtools_filtering'],
                 sort_output=options["samtools_indexing"],
                 index_output=options["samtools_indexing"],
                 log_file=sys.stdout)

        # remove SRA files
        if options["remove_sra"]:
            verbose(5 if options["verbose"] else -1)
            for run in runs: os.remove(run)
    
    # ends
    verbose(6 if options["verbose"] else -1)
    return 0

# execute the main function if the script is executed as the main program
if __name__ == "__main__": sys.exit(main())
            