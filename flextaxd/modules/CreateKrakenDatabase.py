#!/usr/bin/env python3 -c

'''
This part of the script is kept because it is convenient when creating kraken2 or krakenuniq databases
however it is not optimized in any way and uses os.system instead of a proper Popen for subprocessing.
'''

import sys
import gzip
import random
import os
import glob
from multiprocessing import Process,Manager,Pool
from subprocess import Popen,PIPE,check_output,CalledProcessError
from .database.DatabaseConnection import DatabaseFunctions
from time import sleep

'''gzip have changed their error format between python version 3.7 and 3.8, this is at least a temporary fix for that'''
try:
	from gzip import BadGzipFile
except ImportError:
	class BadGzipFile(OSError):
	    """Exception raised in some cases for invalid gzip files."""


import logging
logger = logging.getLogger(__name__)

def zopen(path,*args, **kwargs):
	'''Redefine open to handle zipped files automatically'''
	if path.endswith(".gz"):
		if str(*args) in ["r","w"] and sys.version_info.major >= 3:
			## Python version three and above interprets binary formats as binary, add t (rt) to get text returned
			args = (str(*args)+"t",)
		return gzip.open(path,*args,**kwargs)
	else:
		return open(path,*args,**kwargs)

class CreateKrakenDatabase(object):
	"""docstring for CreateKrakenDatabase."""
	def __init__(self, database, kraken_database, genome_names, outdir,verbose=False,processes=1,limit=0,dbprogram="kraken2",params="",skip="",create_db=False,debug=False,build_processes=None):
		super(CreateKrakenDatabase, self).__init__()
		self.krakenversion = dbprogram
		self.database = DatabaseFunctions(database)
		if outdir == "":
			outdir = "./"
		self.outdir = outdir
		self.seqid2taxid = self.outdir+"/seqid2taxid.map"
		if os.path.exists(self.outdir+"/seqid2taxid.map"): open(self.seqid2taxid,"w").close() ## clean existing map if file exists
		if genome_names:
			self.genome_names = list(genome_names.keys())   ## List for multiprocessing
			self.genome_path = genome_names					## genome_id to path dictionary
		self.accession_to_taxid = self.database.get_genomes(self.database)
		self.files = []
		self.params = params
		self.processes = processes
		if not build_processes:
			self.build_processes = self.processes
		else:
			self.build_processes = build_processes
		self.krakendb= kraken_database
		self.verbose = verbose
		self.create_db = create_db
		self.limit = limit
		self.debug = debug
		if self.debug:
			self.taxidmap_debug = self.database.get_nodes()
		self.seqhead_validator = {}
		self.seqhead_count = 0
		self.skiptax = set()
		self.skipfiles = set()
		if skip:
			if type(skip) == type(dict()):
				self.skipfiles = skip["genome_id"]

				self.skiptax = self.parse_taxid_names(skip["tax_id"])
				logger.info(self.skiptax)
			else:
				self.skiptax = parse_skip(skip.split(","))  ## if node should be skipd this must be true, otherwise nodes in modfile are added to existing database

		logger.info("{krakendb}".format(outdir = self.outdir, krakendb=self.krakendb))
		if not os.path.exists("{krakendb}".format(outdir = self.outdir, krakendb=self.krakendb)):
			os.system("mkdir -p {krakendb}".format(outdir = self.outdir, krakendb=self.krakendb))

	def parse_taxid_names(self,skiptax):
		'''Parse taxid names'''
		taxidmap = self.database.get_nodes()
		newset = set()
		for tid in skiptax:
			if not tid.isdigit():
				try:
					newset.add(taxidmap[tid])
				except KeyError:
					logger.info("# WARNING: {taxa} not in database".format(taxa=tid))
		if len(newset) > 0:
			skiptax = newset
		return skiptax

	def _split(self,a, n):
		k, m = divmod(len(a), n)
		return list(a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

	def kraken_fasta_header_multiproc(self,genomes):
		'''function to run addition of genomes in paralell'''
		logger.info("Processing files; create kraken seq.map")
		jobs = []
		manager = Manager()
		added = manager.Queue()
		for i in range(self.processes):
			p = Process(target=self.kraken_fasta_header, args=(genomes[i],added))
			p.daemon=True
			p.start()
			jobs.append(p)
		for job in jobs:
			job.join()
		self.added = added.qsize()
		return "Processes done"

	def kraken_fasta_header(self,genomes,added):
		'''Change fasta file to contain kraken fasta header'''
		count = 0
		batchint = random.randint(10**3,10**7)
		tmpbatch = "{db_path}/library/batch_{rand}.fasta".format(db_path=self.krakendb,rand=batchint)
		check_output("touch {tmpfile}".format(tmpfile=tmpbatch),shell=True)
		tmplog = "{outdir}/{rand}.log".format(outdir=self.outdir.rstrip("/"),rand=batchint)
		for i in range(len(genomes)):
			genome = genomes[i]
			if len(set([genome]) & self.skipfiles) > 0: ## Skip file if in skipfiles
				logger.debug("Skip {genome}".format(genome=genome))
				logger.info("Skip {genome}".format(genome=genome))
				continue
			filepath = self.genome_path[genome]
			kraken_header = "kraken:taxid"
			tmpname = ".tmp{rand}.gz".format(rand=random.randint(10**7,10**9))
			tmppath = "{outdir}/{tmppath}".format(outdir=self.outdir.rstrip("/"),tmppath=tmpname).rstrip(".gz")

			'''Get taxid from database'''
			try:
				taxid = self.accession_to_taxid[genome]
			except KeyError:
				try:
					taxid = self.accession_to_taxid[filepath.rsplit("/")[-1]]
				except KeyError:
					with open(tmplog, "a") as f:
						print(genome,file=f)
					count +=1
					logger.debug("#Warning kraken header could not be added to {genome}! Total: {count}".format(genome=genome,count=count))

					continue
			if self.create_db:
				if len(set([taxid]) & self.skiptax) == 0:
					'''Open temp file for manipulated (unzipped) genome fasta files'''
					tmpfile = zopen(tmppath,"w")
					taxidlines = []  ## Holder for sequences to be added to the seqid2taxid map from each file
					'''Open input genome fasta file'''
					with zopen(filepath,"r") as f:
						try:
							for line in f:
								if line.startswith(">"):
									taxidmap = []
									row = line.strip().split(" ")
									kh = line.strip().split("|")
									if len(kh) > 1:
										if kh[1].startswith("kraken"):
											'''Assume kraken header is already present'''
											kh = True
									else:
										kh = False
									seq_header = row[0]
									if self.debug:
										self.seqhead_validator[seq_header] = [seq_header, genome, str(taxid),self.taxidmap_debug[taxid],filepath]
										self.seqhead_count +=1
									if len(row) == 1:
										row.append("")
									endhead = " ".join(row[1:])
									if True:  ## Remove end heading in all files as it may contain bad chars
										endhead = "\n"

									'''Format sequence header unless header already have kraken header'''
									if not self.krakenversion == "krakenuniq" and not kh:
										line = row[0] + "|" + kraken_header + "|" + str(taxid) + "  " + endhead	 ## Nessesary to be able to add sequences outside NCBI

									'''Format seqid2taxid map'''
									if not self.krakenversion == "kraken2":
										taxidmap = [row[0].lstrip(">"), str(taxid), endhead]
									else:
										taxidmap= ["TAXID", row[0].lstrip(">") + "|" + kraken_header + "|" + str(taxid), str(taxid)]
									'''Add taxid to map'''
									if taxidmap[0] != "" and taxid: ## print chromosome name to seqtoid map
										taxidlines.append("\t".join(taxidmap)) ## print chromosome name to seqtoid map
									'''Don´t include empty lines'''
								print(line, end="", file=tmpfile)
						except BadGzipFile as e:
							logger.warning("Could not process {output}, not a valid gzip file".format(output=genome))
						except EOFError as e:
							logger.warning("Compressed file ended before the end-of-stream marker was reached {output}".format(output=genome))
					tmpfile.close()  ## Close and save tmp file
					with open(self.seqid2taxid, "a") as seqidtotaxid:
						print("\n".join(taxidlines),end="\n",file=seqidtotaxid)
					output = self.krakendb.rstrip("/")+"/"+filepath.split("/")[-1].rstrip(".gz")
					os.rename(tmppath, output)
					try:
						#ans = check_output(self.krakenversion+"-build --add-to-library {file} --skip-maps --db {krakendb}  > /dev/null 2>&1 ".format(file=output,krakendb=self.krakendb),shell=True)
						check_output("cat {output} >> {tmpbatch}".format(output=output,tmpbatch=tmpbatch),shell=True)
						added.put(1)
					except CalledProcessError as e:
						logger.debug("Could not process {output}, the following error occured:\n{e}".format(output=output,e=e))
					finally:
						os.remove(output)
		'''Validate file'''
		if self.debug:
			if len(self.seqhead_validator.keys()) != self.seqhead_count:
				logger.warning("Sequence headers are not unique added: {added} unique: {nunique}".format(nunique=len(self.seqhead_validator.keys()), added=self.seqhead_count))
				tmpdebug = "{db_path}/library/batch_{rand}.debug".format(db_path=self.krakendb,rand=batchint)
				with open(tmpdebug, "w") as debugwrite:
					for seq_header in self.seqhead_validator.keys():
						print("\t".join(self.seqhead_validator[seq_header]),end="\n",file=debugwrite)
				return False
		return tmpbatch

	def get_skip_list(self):
		'''get taxonomy ids not wanted'''
		return self.skiptax

	def parse_skip(self,skip):
		'''Skip allows a list of tax ids to be passed and then excluded from the database'''
		skiptax = []
		for t in skip:
			skiptax |= self.taxonomydb.get_children()
		return skiptax

	def create_library_from_files(self):
		'''Create library for kraken and create kraken genome2taxid map'''
		if self.limit:
			logger.info("Test use only {n} genomes".format(n=self.limit))
			self.genome_names = self.genome_names[0:self.limit]
		self.genome_names_split = self._split(self.genome_names,self.processes)
		if not os.path.exists("{db_path}/library/".format(db_path=self.krakendb)):
			logger.info("Create library directory")
			os.mkdir("{db_path}/library/".format(db_path=self.krakendb))
		self.kraken_fasta_header_multiproc(self.genome_names_split)
		'''Merge library files into one library'''
		check_output("cat {db_path}/library/batch*.fasta > {db_path}/library/library.fna".format(db_path=self.krakendb),shell=True)
		check_output("rm {db_path}/library/batch*.fasta".format(db_path=self.krakendb),shell=True)
		logger.info("Number of genomes succesfully added to the {krakenversion} database: {count}".format(count=self.added,krakenversion=self.krakenversion))
		return

	def create_database(self,outdir,keep=False):
		'''For test create a small database and run tests'''
		if not os.path.exists("{krakendb}/taxonomy".format(outdir=outdir, krakendb=self.krakendb)):
			logger.info("mkdir -p {krakendb}/taxonomy".format(outdir=outdir, krakendb=self.krakendb))
			os.system("mkdir -p {krakendb}/taxonomy".format(outdir=outdir, krakendb=self.krakendb))
		logger.info("cp {outdir}/*.dmp {krakendb}/taxonomy".format(outdir=outdir,krakendb=self.krakendb))
		os.system("cp {outdir}/*names.dmp {krakendb}/taxonomy/names.dmp".format(outdir=outdir,krakendb=self.krakendb))
		os.system("cp {outdir}/*nodes.dmp {krakendb}/taxonomy/nodes.dmp".format(outdir=outdir,krakendb=self.krakendb))
		if self.krakenversion != "kraken2": os.system("cp {outdir}/*.map {krakendb}".format(outdir=outdir,krakendb=self.krakendb))
		else: os.system("cp {outdir}/*.map {krakendb}/library/prelim_map.txt".format(outdir=outdir,krakendb=self.krakendb))
		logger.info("cp {outdir}/*.map {krakendb}".format(outdir=outdir,krakendb=self.krakendb))
		logger.info(self.krakenversion+"-build --build --db {krakendb} {params} --threads {threads}".format(krakendb=self.krakendb, threads=self.build_processes, params=self.params))
		os.system(self.krakenversion+"-build --build --skip-maps --db {krakendb} {params} --threads {threads}".format(krakendb=self.krakendb, threads=self.build_processes, params=self.params))

		if self.krakenversion in ["kraken2"]:
			logger.info("Create inspect file!")
			os.system(self.krakenversion+"-inspect --db {krakendb} --report-zero-counts --threads {threads} > {krakendb}/inspect.txt".format(krakendb=self.krakendb,threads=self.build_processes ))
			os.system("gzip {krakendb}/*.map".format(krakendb=self.krakendb))
			os.system("gzip {krakendb}/inspect.txt".format(krakendb=self.krakendb))
		if not keep:
			os.system(self.krakenversion+"-build --clean --db {krakendb}".format(outdir=outdir,krakendb=self.krakendb, threads=self.processes))
			# logger.info("Cleaning up tmp files")
			# os.system('find {krakendb} -maxdepth 1 -name "*.f*a" -print0 | xargs -0 rm'.format(krakendb=self.krakendb))
		logger.info("{krakenversion} database created".format(krakenversion=self.krakenversion))
