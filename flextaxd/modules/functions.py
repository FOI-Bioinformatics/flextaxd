from subprocess import check_call,CalledProcessError
from textwrap import wrap
from os import makedirs,path,walk
import glob
import logging
from time import sleep
logger = logging.getLogger(__name__)

def _build_path(genome_id):
    '''Build up the download path and return'''
    ## https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/696/305/GCF_001696305.1_UCN72.1/GCF_001696305.1_UCN72.1_genomic.fna.gz
    parts = genome_id.split("_")
    type = parts[0]
    num = parts[1].split(".")[0]
    npt = wrap(num,3)
    end = ""
    genome_name = genome_id
    if not genome_id[-1].isdigit():
        if not genome_id.endswith(".gz"):
            end = ".gz"
        elif "genomic" not in genome_id:
            end = "_genomic.fna.gz"
        else:
            genome_name = genome_id.split("_genomic")[0]

        path = "{genome_name}/{genome_id}{end}".format(
                                                genome_id=genome_id,
                                                genome_name=genome_name,
                                                end=end
                                                )
    else:
        path = ""
    build_path = "{type}/{pt1}/{pt2}/{pt3}/{path}".format(
                type=type,
                path=path,
                pt1=npt[0],pt2=npt[1],pt3=npt[2])
    return build_path

def call_genome(call_func,outdir,genome_id,added,filepath,retries_d,retries_n):
    '''test'''
    #logger.debug(" ".join(call_func))
    try:
        err = check_call(" ".join(call_func),shell=True)
    except CalledProcessError as e:
        try:
            if retries_d[genome_id] > retries_n:
                print("Maximum number of retries reached for {genome}".format(genome=genome_id))
                return False
            retries_d[genome_id]+=1
        except KeyError:
            retries_d[genome_id] = 1
        if e.output:
            if "Cannot assign requested address" in e.output:
                retries_d[genome_id]+=1
            if "error in socket IO" in e.output:
                retries_d[genome_id]+=1
            if "failed to connect" in e.output:
                retries_d[genome_id]+=1
            if sleep > 20:
                ### Longsleep
                sleep(60)
                sleep = 0
            sleep+=sleep
            sleep(sleep)
            check =  call_genome(call_func,outdir,genome_id,added,filepath,sleep=5,retries_d=retries_d,retries_n=retries_n)
            if not check:
                return False
        if e.output == None:
            return False
    try:
        try:
            for root, dirs, files in walk(glob.glob(path.join(outdir,genome_id+"*"))[0]):
                for file in files:
                    if file.endswith("_genomic.fna.gz") and "from" not in file:
                        fpath = path.join(root,file)
                        logger.debug(fpath)
                        gname = path.basename(file)
                        genome_name = gname.split("_",2)
                        genome_name = genome_name[0]+"_"+genome_name[1]
                        added.put(genome_name.strip())
                        filepath.put(fpath)
                        #logger.info("File downloaded")
        except IndexError:
            return False #pass #logger.debug("File not downloaded")
    except OSError:
        return False
        #logger.debug(" ".join(call_func))
        #logger.error("File not downloaded")
    return True

def download_genome(args_list=False,added=False,filepath=False,sleep=5,retries=10, **kwargs):
    '''Function that takes a genome_id and downloads from source (refseq or genbank)
        input: genome_id; optional outdir

    '''
    retries_n =retries
    retries_d={}
    for args in args_list:
        if args:
            outdir = args["outdir"]
            genome_id = args["genome_id"]
        if not path.exists(outdir):
            makedirs(outdir)
        genomes_path = _build_path(genome_id)
        call_func = ["rsync", "--copy-links", "--recursive", "--times","--ignore-existing","--no-motd",
                        "--prune-empty-dirs",
                        "--exclude='*assembly_structure/'",
                        "--include='*/'",   ## To get subfolders included
                        "--exclude='*from_genomic.fna.gz'",
                        "--include='{genome_id}*_genomic.fna.gz'".format(genome_id=genome_id),
                        "--exclude='*'",
                        "rsync://ftp.ncbi.nlm.nih.gov/genomes/all/{path}".format(path=genomes_path),
                        outdir,
                        "2>&1"
        ]
        call_genome(call_func,outdir,genome_id,added,filepath,retries_d=retries_d,retries_n=retries_n)
