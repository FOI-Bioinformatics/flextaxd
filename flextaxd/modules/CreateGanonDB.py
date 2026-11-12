#!/usr/bin/env python3 -c

'''
This part of the script is kept because it is convenient when creating ganon2 or ganonuniq databases.
It is not optimized in any way, but the build of the external program runs as an argument list through
subprocess and a build that did not finish stops the run (issue #148).
'''

import random
import os
import glob
import shutil
import multiprocessing
from subprocess import Popen,PIPE,CalledProcessError,run
from .database.DatabaseConnection import DatabaseFunctions
from .functions import zopen, annotation_hint
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
		'''issue #71: same advice as the kraken library, a library that loses a large share of its
		   sequences to a missing taxid is an annotation file that does not describe them'''
		hint = annotation_hint(sequences,len(unannotated))
		if hint:
			logger.warning(hint)
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

	def run_program_command(self,cmd,required=True):
		'''Run the build of the external program as an argument list, with no shell in between. os.system
				returned the exit status to a caller that threw it away, so a machine without ganon (the shell
				said "ganon: not found" and the status was ignored) still ended with the run reporting a
				database that had never been built, and the library it had written was deleted on the way out
				(issue #148).

		------
		Parameters:
			list - the program and its arguments
			boolean - whether a program that is missing or failed stops the run
		------
		Returns:
			boolean - False when the failure is allowed
		'''
		logger.info(" ".join(cmd))
		try:
			res = run(cmd,stdout=PIPE,stderr=PIPE)
		except OSError as error:
			logger.error("{program} could not be run: {err}".format(program=cmd[0],err=str(error)))
			if required:
				raise CalledProcessError(127," ".join(cmd))
			return False
		if res.returncode != 0:
			logger.error("The command returned {code}: {cmd}".format(code=res.returncode,cmd=" ".join(cmd)))
			'''what ganon did not like is on its stderr, and that is the part of the failure a user can act
			   on; the last part of it, because a build can have written a lot before it stopped'''
			message = (res.stderr or b"").decode(errors="replace").strip()
			if message:
				logger.error(message[-2000:])
			if required:
				raise CalledProcessError(res.returncode," ".join(cmd))
			return False
		return True

	def create_database(self,outdir,keep=False):
		'''Copy the taxonomy of the dump into the database directory and build the ganon database over the
				library files this run wrote.

				The library files stay in the directory: they are the library of this build, and deleting them
				(a `rm {ganondb}/.tmp*` over everything the run had written) left a ganon database with nothing
				to inspect and nothing to build the next database from, which is what T6 wants to build on
				(issue #148). The build itself is checked: ganon missing, ganon failing or no library at all
				raise a CalledProcessError, which create_databases turns into the "no database was created"
				error and a failing exit status, the way the kraken route has since issue #132.

		------
		Parameters:
			str - the directory the taxonomy dump was written to
			boolean - unused here: the library is kept in either case
		'''
		taxonomy = "{ganon}/taxonomy".format(ganon=self.ganondb.rstrip("/"))
		if not os.path.exists(taxonomy):
			logger.info("mkdir -p {taxonomy}".format(taxonomy=taxonomy))
			os.makedirs(taxonomy, exist_ok=True)
		for pattern,target in (("*nodes.dmp","nodes.dmp"),("*names.dmp","names.dmp")):
			source = sorted(glob.glob("{outdir}/{pattern}".format(outdir=outdir.rstrip("/"),pattern=pattern)))
			if not source:
				'''a database built without one half of the taxdump gets a taxonomy it cannot name, which is
				   worse than saying what is missing here'''
				logger.error("No {pattern} in {outdir}: the ganon database would be built without a taxonomy".format(
					pattern=pattern,outdir=outdir))
				continue
			logger.debug("cp {source} {taxonomy}/{target}".format(source=source[0],taxonomy=taxonomy,target=target))
			shutil.copyfile(source[0],"{taxonomy}/{target}".format(taxonomy=taxonomy,target=target))

		logger.info("Build ganon database")
		'''the library files this run wrote: one per process plus the packed sequences, whatever they are
		   named - listing them by hand over the process numbers left the packed file out of the build'''
		library_files = sorted(glob.glob("{ganondb}/.tmp*.fasta*".format(ganondb=self.ganondb.rstrip("/"))))
		if not library_files:
			logger.error("No library files in {ganondb}, the ganon build has nothing to read".format(
				ganondb=self.ganondb.rstrip("/")))
			raise CalledProcessError(1,"ganon build: no library files in {ganondb}".format(
				ganondb=self.ganondb.rstrip("/")))
		logger.info("Building the ganon database from {n} library files".format(n=len(library_files)))
		self.run_program_command(["ganon","build",
			"-d",self.ganondb.rstrip("/"),
			"--seq-info-file",self.seqid2taxid,
			"-t",str(self.build_processes),
			"--taxdump-file","{taxonomy}/nodes.dmp".format(taxonomy=taxonomy),
			"{taxonomy}/names.dmp".format(taxonomy=taxonomy),
			"-i"] + library_files)
		logger.info("The library of this build is kept in {ganondb}: {files}".format(
			ganondb=self.ganondb.rstrip("/"),
			files=", ".join(os.path.basename(name) for name in library_files)))
		return

