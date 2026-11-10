#!/usr/bin/env python3 -c

'''
This part of the script is kept because it is convenient when creating ganon2 or ganonuniq databases
however it is not optimized in any way and uses os.system instead of a proper Popen for subprocessing.
'''

import random
import os
import glob
import multiprocessing
from subprocess import Popen,PIPE
from .database.DatabaseConnection import DatabaseFunctions
from .functions import zopen
from time import sleep
from gzip import BadGzipFile

import logging
logger = logging.getLogger(__name__)

class CreateGanonDB(object):
	"""docstring for CreateGanonDB."""
	def __init__(self, database, ganon_database, genome_names, outdir,verbose=False,debug=False,processes=1,limit=0,dbprogram="ganon",params="",create_db=False,skip=False,build_processes=False,usezip=False,**kwargs):
		super(CreateGanonDB, self).__init__()
		'''create_databases hands one argument list to whichever program was selected, and ganon accepted
		   neither tmpdir, create_lib, skip_genomes nor multifiles: every --dbprogram ganon run died with
		   "TypeError: unexpected keyword argument" before any work was done (issue #133). The options that
		   mean something for ganon are used, the others are accepted for the shared interface until T5
		   unifies these classes.'''
		self.database = DatabaseFunctions(database)
		if outdir == "":
			outdir = "./"
		self.outdir = outdir.rstrip("/")+"/"
		self.tmpdir = kwargs.get("tmpdir") or outdir
		self.create_lib = kwargs.get("create_lib", False)
		self.skip_genomes = kwargs.get("skip_genomes") or []
		self.packed_files = []
		for packed in (kwargs.get("multifiles") or []):
			if packed not in self.packed_files:
				self.packed_files.append(packed)
		self.params=params
		self.seqid2taxid = self.outdir+"/seqid2taxid.map"
		if os.path.exists(self.outdir+"/seqid2taxid.map"): open(self.seqid2taxid,"w").close() ## clean existing map if file exists
		self.accession_to_taxid = self.database.get_genomes(self.database , limit=limit)
		'''a build that was told to skip the genome folder arrives here with no genome dictionary at all'''
		self.genome_path = genome_names or {}
		self.genome_names = list(self.genome_path.keys())   ## List for multiprocessing
		self.files = []
		self.usezip = ""
		if usezip:
			self.usezip+=".gz"
		self.processes = processes
		if not build_processes:
			self.build_processes = self.processes
		else:
			self.build_processes = build_processes
		self.ganondb= ganon_database  ## Ganon database path
		self.verbose = verbose
		logger.info("{ganondb}".format(outdir = self.outdir, ganondb=self.ganondb))
		if not os.path.exists("{ganondb}".format(outdir = self.outdir, ganondb=self.ganondb)):
			os.system("mkdir -p {ganondb}".format(outdir = self.outdir, ganondb=self.ganondb))

	def ganon_fasta_multiproc(self,genomes):
		'''function to run addition of genomes in paralell'''
		'''The children inherit the database connection and the genome path dictionary of this object
		   instead of having them pickled over, so this pass needs the fork start method whatever the
		   interpreter has as its default: python 3.14 made forkserver the default and that pickles the
		   target, which dies on the open sqlite connection (issue #146, the same failure #130 fixed for
		   the kraken header pass)'''
		ctx = multiprocessing.get_context("fork")
		jobs = []
		for i in range(self.processes):
			p = ctx.Process(target=self.ganon_fasta, args=(genomes[i],i))
			p.daemon = True
			p.start()
			jobs.append(p)
		for job in jobs:
			job.join()
		## Concatenate each tmp file (per processor) to one big datafile
		logger.info("Merge tmp files.")
		logger.debug("Merge geneid map")
		os.system("cat {outdir}/.tmp*.map > {taxidmap}".format(outdir=self.outdir, taxidmap=self.seqid2taxid))
		os.system("rm {outdir}/.tmp*.map".format(outdir=self.outdir))


		return "Processes done"

	def ganon_fasta(self,genomes,process):
		'''Change fasta file to contain ganon fasta header'''
		tmpname = ".tmp{rand}.fasta".format(rand=process)
		tmpmapname =  ".tmp{rand}.map".format(rand=process)
		tmppath = "{outdir}/{tmppath}".format(outdir=self.ganondb.rstrip("/"),tmppath=tmpname)
		tmpmap =  "{outdir}/{tmppath}".format(outdir=self.outdir.rstrip("/"),tmppath=tmpmapname)
		seqlen = 0
		tmpout = open(tmppath, "w")  ## Open thread output file
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
		if self.usezip:
			os.system("gzip {tmpout}".format(tmpout=tmppath))
		print("Process-{process} completed".format(process=process),end="\r")
		return True

	def add_packed_file_to_library(self,mfile,tmpout):
		'''A packed multi sequence file (the --nt_source file, a --multifile_prefix file) is never split into
		   genome files: its sequences are written to the library as they are and the taxid of every one of
		   them goes to the seqid2taxid map, which is what ganon reads for the taxid of a sequence (issue #133).
		   ganon also wants the length of every sequence in that map, which the kraken map does not carry.
		------
		Parameters:
			str - path to the packed file
			file - the open library file of this run
		------
		Returns:
			tuple - (sequences read, sequences with a taxid in the map)
		'''
		if not os.path.exists(mfile):
			logger.warning("Packed sequence file {mfile} does not exist, nothing was added to the library".format(mfile=mfile))
			return (0,0)
		sequences = 0
		annotated = 0
		unannotated = set()
		with zopen(mfile,"r") as library_file, open(self.seqid2taxid,"a") as taxidmap:
			seqid = False
			lines = []
			seqlen = 0
			def flush():
				'''write one sequence: the library keeps it as it is, the map carries its length and taxid'''
				if not seqid or seqlen == 0:
					return
				taxid = self.accession_to_taxid.get(seqid,False)
				if not taxid:
					unannotated.add(seqid)
					logger.debug("No taxid in the database for sequence {seqid} of {mfile}".format(seqid=seqid,mfile=mfile))
					return
				print("{seqid}\t{seqlen}\t{taxid}".format(seqid=seqid,seqlen=seqlen,taxid=taxid),end="\n",file=taxidmap)
				print(">{seqid}".format(seqid=seqid),file=tmpout)
				print("\n".join(lines),file=tmpout)
			for line in library_file:
				if line.startswith(">"):
					flush()
					words = line.split()
					seqid = words[0].lstrip(">").strip() if words else False
					lines = []
					seqlen = 0
					if seqid:
						sequences += 1
				elif seqid:
					sequence = line.strip()
					if sequence:
						lines.append(sequence)
						seqlen += len(sequence)
			flush()
		annotated = sequences - len(unannotated)
		if unannotated:
			logger.warning("{n} sequences of {mfile} have no annotation in the database and cannot be placed by ganon: {names}".format(
				n=len(unannotated),mfile=mfile,names=", ".join(sorted(unannotated)[:10])))
		logger.info("{sequences} sequences of {mfile} are in the library, {annotated} have a taxid in the map".format(
			sequences=sequences,mfile=mfile,annotated=annotated))
		return (sequences,annotated)

	def create_library_from_files(self,multifiles=False):
		'''Create the ganon library and its seqid2taxid map. The packed files arrive as an argument (the
		   --nt_source file and what --multifile_prefix found), the same call the kraken class gets (#133).'''
		processes = self.processes
		if self.skip_genomes:
			'''genomes the run was told to skip never enter the library. The skip list is a dictionary of
			   genome ids and taxids when it came from a file, a plain list of names otherwise'''
			if isinstance(self.skip_genomes,dict):
				skipped = set(self.skip_genomes.get("genome_id",[]))
			else:
				skipped = set(self.skip_genomes)
			before = len(self.genome_names)
			self.genome_names = [genome for genome in self.genome_names if genome not in skipped]
			if before != len(self.genome_names):
				logger.info("{n} genomes on the skip list were left out of the library".format(n=before-len(self.genome_names)))
		packed_files = list(self.packed_files)
		for mfile in (multifiles or []):
			if mfile not in packed_files:
				packed_files.append(mfile)
		if self.genome_names:
			logger.info("Process {ngen} genomes".format(ngen=len(self.genome_names)))
			genome_process_list = self.split(self.genome_names,processes)
			logger.info(self.ganon_fasta_multiproc(genome_process_list))
		else:
			logger.info("No genome files to add, the library is built from the packed files only")
		if packed_files:
			logger.info("Adding sequences from {n} packed sequence files to the library".format(n=len(packed_files)))
			packed = "{ganondb}/.tmppacked.fasta".format(ganondb=self.ganondb.rstrip("/"))
			added = 0
			with open(packed,"w") as tmpout:
				for mfile in packed_files:
					sequences,annotated = self.add_packed_file_to_library(mfile,tmpout)
					added += sequences
			if self.usezip:
				os.system("gzip {packed}".format(packed=packed))
			logger.info("Number of sequences added to the library from the packed files: {count}".format(count=added))
		return

	def split(self,a, n):
		k, m = divmod(len(a), n)
		return list(a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

	def create_database(self,outdir,keep=False):
		'''For test create a small database and run tests'''
		if not os.path.exists("{ganon}/taxonomy".format(outdir=outdir, ganon=self.ganondb)):
			logger.debug("mkdir -p {ganon}/taxonomy".format(outdir=outdir, ganon=self.ganondb))
			os.system("mkdir -p {ganon}/taxonomy".format(outdir=outdir, ganon=self.ganondb))
		logger.debug("cp {outdir}/*.dmp {ganon}/taxonomy".format(outdir=outdir,ganon=self.ganondb))
		os.system("cp {outdir}/*nodes.dmp {ganon}/taxonomy/nodes.dmp".format(outdir=outdir,ganon=self.ganondb))
		os.system("cp {outdir}/*names.dmp {ganon}/taxonomy/names.dmp".format(outdir=outdir,ganon=self.ganondb))

		logger.info("Build ganon database")
		'''the library files this run wrote: one per process plus the packed sequences, whatever they are
		   named - listing them by hand over the process numbers left the packed file out of the build'''
		logger.debug("ganon build -d {ganondb} --seq-info-file {seqid2taxid} -t {processes} --taxdump-file {ganondb}/taxonomy/nodes.dmp {ganondb}/taxonomy/names.dmp -i {ganondb}/library.fasta{usezip}".format(ganondb=self.ganondb,seqid2taxid=self.seqid2taxid,processes=self.build_processes,usezip=self.usezip))
		#os.system("ganon build -d {ganondb} --seq-info-file {seqid2taxid} -t {processes} --taxdump-file {ganondb}/taxonomy/nodes.dmp {ganondb}/taxonomy/names.dmp -i {ganondb}/library.fasta{usezip}".format(ganondb=self.ganondb,seqid2taxid=self.seqid2taxid,processes=self.build_processes,usezip=self.usezip))
		inputlist = " ".join(sorted(glob.glob("{outdir}/.tmp*.fasta*".format(outdir=self.ganondb.rstrip("/")))))
		if not inputlist:
			logger.error("No library files in {ganondb}, the ganon build has nothing to read".format(ganondb=self.ganondb))
		os.system("ganon build -d {ganondb}/ --seq-info-file {seqid2taxid} -t {processes} --taxdump-file {ganondb}/taxonomy/nodes.dmp {ganondb}/taxonomy/names.dmp -i {inputlist}".format(ganondb=self.ganondb.rstrip("/"),seqid2taxid=self.seqid2taxid,processes=self.build_processes,outdir=outdir,usezip=self.usezip,inputlist=inputlist))

		# logger.debug("Merge gene library map")
		# os.system("cat {outdir}/.tmp*{usezip} > {ganondb}/library.fasta".format(outdir=self.outdir, ganondb=self.ganondb,usezip=self.usezip))
		os.system("rm {outdir}/.tmp*{usezip}".format(outdir=self.ganondb.rstrip("/"),usezip=self.usezip))
		# logger.info("zip library")
		# try:
		# 	os.system("pigz -p 8 {ganondb}/library.fasta{usezip}".format(usezip=""))
		# except:
		# 	os.system("gzip {ganondb}/library.fasta{usezip}".format(usezip=""))
		#logger.info(self.ganon+"-build --build --db {ganon} --threads {threads}".format(ganon=self.ganon, threads=self.processes))
		#logger.info("{ganon} database created".format(ganon=self.ganon))
