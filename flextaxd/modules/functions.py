from subprocess import check_call,CalledProcessError
from textwrap import wrap
from os import makedirs,path

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

def download_genome(args=False,**kwargs):
    '''Function that takes a genome_id and downloads from source (refseq or genbank)
        input: genome_id; optional outdir

    '''
    if args:
        outdir = args["outdir"]
        genome_id = args["genome_id"]
    if not path.exists(outdir):
        makedirs(outdir)
    call_func = ["rsync", "--copy-links", "--recursive", "--times","--ignore-existing","--no-motd",
                    "--exclude='*assembly_structure/'",
                    "--include='*/'",   ## To get subfolders included
                    "--exclude='*from_genomic.fna.gz'",
                    "--include='{genome_id}*_genomic.fna.gz'".format(genome_id=genome_id),
                    "--exclude='*'",
                    "rsync://ftp.ncbi.nlm.nih.gov/genomes/all/{path}".format(path=_build_path(genome_id))
                    ,outdir
    ]
    try:
        err = check_call(" ".join(call_func),shell=True)
        return genome_id,None
    except CalledProcessError as e:
        msg = "rsync error downloading {genome}".format(error=e.output,genome=genome_id)
        return None,msg
