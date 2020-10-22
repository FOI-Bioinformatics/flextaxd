#!/usr/bin/env python3 -c

'''
This part of the script is kept because it is convenient when creating ganon2 or ganonuniq databases
however it is not optimized in any way and uses os.system instead of a proper Popen for subprocessing.
'''

import sys
import gzip
import random
import os
from multiprocessing import Process, Queue
from subprocess import Popen,PIPE
from .database.DatabaseConnection import DatabaseFunctions
from time import sleep
from gzip import BadGzipFile

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

class CreateGanonDB(object):
	"""docstring for CreateGanonDB."""
	def __init__(self, database, ganon_database, genome_names, outdir,verbose=False,debug=False,processes=1,limit=0,dbprogram="ganon",params="",create_db=False,skip=False,build_processes=False):
		super(CreateGanonDB, self).__init__()
		self.database = DatabaseFunctions(database)
		if outdir == "":
			outdir = "./"
		self.outdir = outdir.rstrip("/")+"/"
		self.params=params
		self.seqid2taxid = self.outdir+"/seqid2taxid.map"
		if os.path.exists(self.outdir+"/seqid2taxid.map"): open(self.seqid2taxid,"w").close() ## clean existing map if file exists
		self.accession_to_taxid = self.database.get_genomes(self.database , limit=limit)
		self.genome_names = list(genome_names.keys())   ## List for multiprocessing
		self.genome_path = genome_names					## genome_id to path dictionary
		self.files = []
		self.processes = processes
		if not build_processes:
			self.build_processes = self.processes
		else:
			self.build_processes = build_processes
		self.ganondb= ganon_database
		self.verbose = verbose
		logger.info("{ganondb}".format(outdir = self.outdir, ganondb=self.ganondb))
		if not os.path.exists("{ganondb}".format(outdir = self.outdir, ganondb=self.ganondb)):
			os.system("mkdir -p {ganondb}".format(outdir = self.outdir, ganondb=self.ganondb))

	def ganon_fasta_multiproc(self,genomes):
		'''function to run addition of genomes in paralell'''
		jobs = []
		for i in range(self.processes):
			p = Process(target=self.ganon_fasta, args=(genomes[i],i))
			p.start()
			jobs.append(p)
		for job in jobs:
			job.join()
		## Concatenate each tmp file (per processor) to one big datafile
		os.system("zcat {outdir}/.tmp*.gz | gzip > {ganondb}/library.fasta.gz".format(outdir=self.outdir, ganondb=self.ganondb))
		os.system("cat {outdir}/.tmp*.map > {taxidmap}".format(outdir=self.outdir, taxidmap=self.seqid2taxid))
		## Rm tmp files
		os.system("rm {outdir}/.tmp*.gz".format(outdir=self.outdir))
		os.system("rm {outdir}/.tmp*.map".format(outdir=self.outdir))
		return "Processes done"

	def ganon_fasta(self,genomes,process):
		'''Change fasta file to contain ganon fasta header'''
		tmpname = ".tmp{rand}".format(rand=process)
		tmpmapname =  ".tmp{rand}.map".format(rand=process)
		tmppath = "{outdir}/{tmppath}".format(outdir=self.outdir.rstrip("/"),tmppath=tmpname)
		tmpmap =  "{outdir}/{tmppath}".format(outdir=self.outdir.rstrip("/"),tmppath=tmpmapname)
		seqlen = 0
		tmpout = zopen(tmppath, "w")  ## Open thread output file
		for i in range(len(genomes)):
			genome = genomes[i]
			filepath = self.genome_path[genome]
			try:
				taxid = self.accession_to_taxid[genome]
			except KeyError:
				logger.debug("# WARNING: {genome} could not be added to database".format(genome=genome))
				continue
			with zopen(filepath,"r") as f:
				'''Create a sequence to taxid mapping file according to ganon specifications
					Pre-generated file with sequence information
					(seqid <tab> seq.len <tab> taxid [<tab> assembly id])
				'''
				### Set local variables
				lines = []
				seqlen = 0
				idstring = ""
				header = ">{id}	{seqlen}	{taxid}"
				with open(tmpmap, "a") as seqidtotaxid:
					try:
						for line in f:
							if line.startswith(">"):
								row = line.split(" ")
								header_id = row[0].lstrip(">").strip()+"_"+genome
								if seqlen > 0:  ## Print header and data
									header = header.format(id=id,seqlen=seqlen,taxid=taxid)
									idstring = ">"+id
									print(header.lstrip(">"),end="\n",file=seqidtotaxid)
									print(idstring, file=tmpout)  ## Print new unique header
									print("\n".join(lines),file=tmpout)
									## Reset local variables
									lines = []
									seqlen = 0
									header = ">{id}	{seqlen}	{taxid}"
								## Update ID
								id = row[0].lstrip(">").strip()+"_"+genome
							else:
								## count lenght of sequence
								seqlen += len(line.strip())
								string = line.strip()
								if string != "":
									lines.append(string)
								#print(line.strip(), file=tmpout)
						if seqlen > 0: ## Print final sequence after loop has completed
							header = header.format(id=id,seqlen=seqlen,taxid=taxid)
							print(header.lstrip(">"),end="\n",file=seqidtotaxid)
							print(idstring, file=tmpout)  ## Print new unique header
							print("\n".join(lines),file=tmpout)
							seqlen = 0
					except BadGzipFile as e:
						logger.warning("Could not process {output}, not a valid gzip file".format(output=genome))
					except EOFError as e:
						logger.warning("Compressed file ended before the end-of-stream marker was reached {output}".format(output=genome))
		sleep(0.01)
		tmpout.close()
		os.system("gzip {tmpout}".format(tmpout=tmppath))
		return True

	def create_library_from_files(self):
		processes = self.processes
		self.genome_names = self.split(self.genome_names,processes)
		logger.info(self.ganon_fasta_multiproc(self.genome_names))
		#sleep(1)
		#logger.info("Number of genomes added to ganon database: {count}".format(count=len(self.genome_names)))
		return

	def split(self,a, n):
		k, m = divmod(len(a), n)
		return list(a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

	def create_database(self,outdir,keep=False):
		'''For test create a small database and run tests'''
		if not os.path.exists("{ganon}/taxonomy".format(outdir=outdir, ganon=self.ganondb)):
			logger.info("mkdir -p {ganon}/taxonomy".format(outdir=outdir, ganon=self.ganondb))
			os.system("mkdir -p {ganon}/taxonomy".format(outdir=outdir, ganon=self.ganondb))
		logger.info("cp {outdir}/*.dmp {ganon}/taxonomy".format(outdir=outdir,ganon=self.ganondb))
		os.system("cp {outdir}/*nodes.dmp {ganon}/taxonomy/nodes.dmp".format(outdir=outdir,ganon=self.ganondb))
		os.system("cp {outdir}/*names.dmp {ganon}/taxonomy/names.dmp".format(outdir=outdir,ganon=self.ganondb))

		logger.info("ganon build -d {ganondb} --seq-info-file {seqid2taxid} -t {processes} --taxdump-file {ganondb}/taxonomy/nodes.dmp {ganondb}/taxonomy/names.dmp -i {ganondb}/library.fasta.gz".format(ganondb=self.ganondb,seqid2taxid=self.seqid2taxid,processes=self.build_processes))
		os.system("ganon build -d {ganondb} --seq-info-file {seqid2taxid} -t {processes} --taxdump-file {ganondb}/taxonomy/nodes.dmp {ganondb}/taxonomy/names.dmp -i {ganondb}/library.fasta.gz".format(ganondb=self.ganondb,seqid2taxid=self.seqid2taxid,processes=self.build_processes))
		#logger.info(self.ganon+"-build --build --db {ganon} --threads {threads}".format(ganon=self.ganon, threads=self.processes))
		#logger.info("{ganon} database created".format(ganon=self.ganon))
