from subprocess import check_call,CalledProcessError,TimeoutExpired,STDOUT
from subprocess import Popen,STDOUT,PIPE,CalledProcessError,TimeoutExpired
from textwrap import wrap
from os import makedirs,path,walk
import glob
import logging
from time import sleep
logger = logging.getLogger(__name__)

SUPPORTED_TAXONOMIC_GROUPS = [
	'bacteria',
	'archaea',
	'fungi',
	'vertebrate_mammalian',
	'invertebrate',
	'plant',
	'protozoa',
	'metagenomes',
	'vertebrate_other',
	'viral'
]

swap_section = {"genbank":"refseq","refseq":"genbank"}

def run(cmd,accession):
	'''Run command'''
	try:
		res = check_call(cmd,stderr=STDOUT,stdout=PIPE,shell=True,
		universal_newlines=True)
		# res = Popen(cmd.split(" "),stderr=STDOUT, stdout=PIPE,
		# universal_newlines=True)
		# res = Popen(cmd,stderr=STDOUT, stdout=PIPE,
		# universal_newlines=True,shell=True)
	except CalledProcessError as e:
		'''Expected error for checking taxonomic group, pass'''
		return e
	except TimeoutExpired as e:
		logger.debug("Genome download timed out for {accession}".format(accession=accession))
		return e
	return False

def get_section(accession):
	'''Check if accession is GCF or GCA

	Parameters
		str accession
	Returns
		str refseq or genbank
	'''
	if accession.strip()[2] == "F":
		return "refseq"
	return "genbank"

def check_taxonomic_group(accession,section=False,force=False):
	'''Check all taxonomic groups if genome exists
	Parameters
		str accession
	Returns
		taxonomic_group
	'''
	base_cmd = "ncbi-genome-download -A {accession} {group} -s {section} -n"
	if not section:
		section = get_section(accession)
	for group in SUPPORTED_TAXONOMIC_GROUPS:
		cmd = base_cmd.format(
			accession = accession,
			group=group,
			section=section
		)
		e = run(cmd,accession)
		if not e:
			#logger.debug(cmd)
			return group,section
		else:
			pass #logger.debug(e)
	'''group could not be found in section, if X on try other section'''
	'''This should not be done as genomes removed from RefSeq usually are better to skip'''
	if force:
		group,section = check_taxonomic_group(accession,section=swap_section[section])
		if group:
			return group,section
	return False,False

def get_genome(accession,force=False):
	'''Return taxonomic group of accession

	Parameters
		str accession
	Returns
		taxonomic_group
	'''
	group, section = check_taxonomic_group(accession,force=force)
	genome = {"accession":accession,"group":group, "section":section}
	return genome

def ncbi_genome_download(genome,outdir="./downloads"):
	'''Use the ncbi-genome-download package to download missing files

	Parameters
		dict - accession group and section
	Returns
		boolean - True if downloaded, else False
	'''
	accession = genome["accession"]
	group = genome["group"]
	section = genome["section"]
	if not section:
		return False
	outdir = "/".join([outdir,group,accession])
	if not group:
		return False
	base_cmd = "ncbi-genome-download -A {accession} {group} -s {section} -o {outdir} -r 3 --flat-output -F fasta"
	cmd = base_cmd.format(
			accession = accession,
			group=group,
			section=section,
			outdir=outdir
	)
	e = run(cmd,accession)
	if not e:
		logger.debug(cmd)
		return outdir
	else:
		logger.debug(cmd)
		logger.debug(e)
	return False

def download_genomes(genomes,added,filepath,missing,force=False):
	for gen_i in genomes:
		if gen_i:
			outdir = gen_i["outdir"]
			genome = get_genome(gen_i["genome_id"],force)
			if genome["accession"].startswith("GCF") or genome["accession"].startswith("GCA"):
				outdir = ncbi_genome_download(genome,outdir)
				if outdir:
					added.put(genome["accession"].strip())
					filepath.put(outdir)
				else:
					missing.put(genome["accession"].strip())
