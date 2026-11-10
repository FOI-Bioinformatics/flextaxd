#!/usr/bin/env python3 -c

'''
This part of the script is kept because it is convenient when creating kraken2 or krakenuniq databases
however it is not optimized in any way and uses os.system instead of a proper Popen for subprocessing.
'''

import random
import os
import glob
import shlex
from multiprocessing import Process,Manager,Pool
import multiprocessing
from subprocess import Popen,PIPE,check_output,CalledProcessError,run
from .database.DatabaseConnection import DatabaseFunctions
from .functions import zopen, strip_gz
from time import sleep

'''gzip have changed their error format between python version 3.7 and 3.8, this is at least a temporary fix for that'''
try:
	from gzip import BadGzipFile
except ImportError:
	class BadGzipFile(OSError):
	    """Exception raised in some cases for invalid gzip files."""


import logging
logger = logging.getLogger(__name__)

class CreateKrakenDatabase(object):
	"""docstring for CreateKrakenDatabase."""
	def __init__(self, database, kraken_database, genome_names, outdir,verbose=False,processes=1,limit=0,dbprogram="kraken2",params="",skip="",create_db=False,debug=False,build_processes=None,**kwargs):
		super(CreateKrakenDatabase, self).__init__()
		self.krakenversion = dbprogram
		self.database = DatabaseFunctions(database)
		if outdir == "":
			outdir = "./"
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		self.outdir = outdir
		try:
			self.tmpdir = kwargs["tmpdir"]
		except:
			self.tmpdir = outdir
		self.seqid2taxid = self.outdir+"/seqid2taxid.map"
		logger.info("seqid2taxid file:{tfile}".format(tfile=self.seqid2taxid))
		'''The map kraken2 reads beside the library, written to the temp directory by the dump of the
			database annotations and added to by every sequence that enters the library'''
		self.prelim_map = self.tmpdir.rstrip("/") + "/prelim_map.txt"
		self.skip_genomes = kwargs["skip_genomes"]
		'''Packed multi sequence files (an nt file, a --multifile_prefix file) are never split up into
			one file per genome: their sequences go into the library as they are and their annotation
			travels in the taxid map, which is what kraken2 and ganon read for those sequences.'''
		self.packed_files = []
		for packed in kwargs.get("multifiles", False) or []:
			if packed and packed not in self.packed_files:
				self.packed_files.append(packed)
		if genome_names and not self.skip_genomes:
			self.genome_names = list(genome_names.keys())   ## List for multiprocessing
			self.genome_path = genome_names					## genome_id to path dictionary
		else:
			'''An empty genome list is not an error when packed files are given: the library is still
				built from them and the map still gets its lines (issue #114)'''
			self.genome_names = []
			self.genome_path = {}
			if self.packed_files:
				logger.warning("Genome names are missing, only the {n} packed sequence files given are added to the library.".format(n=len(self.packed_files)))
			else:
				logger.warning("Genome names are missing.") # Make sure your genomes are formatted as GCF_000000000.0.fasta[.gz]")
		self.accession_to_taxid = self.database.get_genomes(self.database)
		self.files = []
		self.params = params
		self.processes = processes
		if not build_processes:
			self.build_processes = self.processes
		else:
			self.build_processes = build_processes
		self.krakendb= kraken_database
		self.create_lib = kwargs["create_lib"]
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
		
		self.skip = False
		if skip:
			self.skip = True
			if type(skip) == type(dict()):
				self.skipfiles = skip["genome_id"]

				self.skiptax = self.parse_taxid_names(skip["tax_id"])
				logger.info(self.skiptax)
			else:
				self.skiptax = self.parse_skip(skip.split(","))  ## if node should be skipd this must be true, otherwise nodes in modfile are added to existing database

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
		'''The children inherit the database connection and the genome path dictionary of this object
			instead of having them pickled to them, so this pass needs the fork start method whatever
			the interpreter has as its default: python 3.14 made forkserver the default and that pickles
			the target, which fails on the open sqlite connection (issue #130). T6 replaces this parallel
			pass with one streamed read of the inputs.'''
		ctx = multiprocessing.get_context("fork")
		jobs = []
		manager = ctx.Manager()
		added = manager.Queue()
		for i in range(self.processes):
			p = ctx.Process(target=self.kraken_fasta_header, args=(genomes[i],added))
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
		tmpbatch = False
		if self.create_db or self.create_lib: 
			tmpbatch = "{tmpdir}/library/batch_{rand}.fasta".format(tmpdir=self.tmpdir,rand=batchint)
			check_output("touch {tmpfile}".format(tmpfile=tmpbatch),shell=True)
			tmplog = "{tmpdir}/{rand}.log".format(tmpdir=self.tmpdir.rstrip("/"),rand=batchint)
		for i in range(len(genomes)):
			genome = genomes[i]
			if len(set([genome]) & self.skipfiles) > 0: ## Skip file if in skipfiles
				logger.debug("Skip {genome}".format(genome=genome))
				logger.info("Skip {genome}".format(genome=genome))
				continue
			filepath = self.genome_path[genome]
			kraken_header = "kraken:taxid"
			tmpname = ".tmp{rand}.gz".format(rand=random.randint(10**7,10**9))
			tmppath = strip_gz("{tmpdir}/{tmppath}".format(tmpdir=self.tmpdir.rstrip("/"),tmppath=tmpname))

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
			if not self.skip:
				if len(set([taxid]) & self.skiptax) == 0:
					'''Open temp file for manipulated (unzipped) genome fasta files'''
					if not self.skip: ## Only required if kraken map is required (can be printed from db)
						tmpfile = zopen(tmppath,"w")
					taxidlines = []  ## Holder for sequences to be added to the seqid2taxid map from each file
					seq2taxidlines = []
					seq2taxmap = False
					'''Open input genome fasta file'''
					
				with zopen(filepath,"r") as f:
					try:
						for line in f:
							if self.create_db or self.create_lib:
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
									elif self.krakenversion and not self.create_db:
										taxidmap= ["TAXID", row[0].lstrip(">"), str(taxid)]
										seq2taxmap = [row[0].lstrip(">"), str(taxid)]
									else:
										taxidmap= ["TAXID", row[0].lstrip(">") + "|" + kraken_header + "|" + str(taxid), str(taxid)]
									'''Add taxid to map'''
									if taxidmap[0] != "" and taxid: ## print chromosome name to seqtoid map
										taxidlines.append("\t".join(taxidmap)) ## print chromosome name to seqtoid map
										if seq2taxmap: seq2taxidlines.append("\t".join(seq2taxmap))
									'''Don´t include empty lines'''
							if not self.skip: ## Only required if kraken map is required (can be printed from db)
								print(line, end="", file=tmpfile)
					except BadGzipFile as e:
						logger.warning("Could not process {output}, not a valid gzip file".format(output=genome))
					except EOFError as e:
						logger.warning("Compressed file ended before the end-of-stream marker was reached {output}".format(output=genome))
					
					with open("{outdir}/prelim_map.txt".format(outdir=self.tmpdir.rstrip("/")), "a") as seqidtotaxid:  ## Add annotation to taxid map
						for taxidline in taxidlines:
							print(taxidline.strip('\n'),end="\n",file=seqidtotaxid)
					if seq2taxmap: 
						with open(self.seqid2taxid, "a") as seqidtotaxid:  ## Add annotation to taxid map
							for taxidline in seq2taxidlines:
								print(taxidline.strip('\n'),end="\n",file=seqidtotaxid)
					tmpfile.close()  ## Close and save tmp file
					if self.create_db or self.create_lib:
						output = strip_gz(self.tmpdir.rstrip("/")+"/"+filepath.split("/")[-1])
						os.rename(tmppath, output)
						try:
							#ans = check_output(self.krakenversion+"-build --add-to-library {file} --skip-maps --db {krakendb}  > /dev/null 2>&1 ".format(file=output,krakendb=self.krakendb),shell=True)
							check_output("cat {output} >> {tmpbatch}".format(output=output,tmpbatch=tmpbatch),shell=True)
							added.put(1)
						except CalledProcessError as e:
							logger.debug("Could not process {output}, the following error occured:\n{e}".format(output=output,e=e))
						finally:
							os.remove(output)
					else:
						os.remove(tmppath)
		'''Validate file'''
		if self.debug:
			if len(self.seqhead_validator.keys()) != self.seqhead_count:
				logger.warning("Sequence headers are not unique added: {added} unique: {nunique}".format(nunique=len(self.seqhead_validator.keys()), added=self.seqhead_count))
				tmpdebug = "{tmpdir}/library/batch_{rand}.debug".format(db_path=self.krakendb, tmpdir=self.tmpdir,rand=batchint)
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

	def sequence_id(self,header):
		'''The id of a sequence is the first white space separated word of its FASTA header, the same
			rule the taxonomy import uses (ReadTaxonomyNCBI.sequence_id), here on text'''
		words = header.split()
		if not words:
			return False
		return words[0].lstrip(">").strip()

	def write_unannotated_sequences(self,mfile,unannotated):
		'''The sequences of a packed file the database has no annotation for: the list is written next to
			the input and only the number is stored in the report, like every other finding of a build'''
		unannotated = sorted(set(unannotated))
		unannotated_file = "{mfile}.unannotated".format(mfile=mfile)
		with open(unannotated_file,"w") as of:
			for seqid in unannotated:
				print(seqid, end="\n", file=of)
		self.database.save_report("unannotated_sequence",mfile,len(unannotated),unannotated_file)
		logger.warning("{n} sequences of {mfile} have no annotation in the database, they are listed in {file}".format(n=len(unannotated),mfile=mfile,file=unannotated_file))
		return unannotated_file

	def add_packed_file_to_library(self,mfile,added):
		'''Add one packed multi sequence file to the library without touching the input file: the file is
			decompressed into library.fna and its headers keep the accession they have, so the taxid of
			every one of its sequences is written to the taxid map instead (issue #115). A sequence that
			has no annotation in the database is counted and listed, it stays in the library without a
			taxid and is reported by the library check of T4.
		------
		Parameters:
			str - path to the packed file
			Queue - the shared counter of added sequences
		------
		Returns:
			boolean - True if the file was added, False if it could not be read
		'''
		if not os.path.exists(mfile):
			logger.warning("Packed sequence file {mfile} does not exist, nothing was added to the library".format(mfile=mfile))
			return False
		logger.info("Adding the sequences of {mfile} to the library and their taxids to the map".format(mfile=mfile))
		sequences = 0
		annotated = 0
		unannotated = set()
		'''Stream the headers only: the file is read once for the map, the sequences themselves are
			copied below, which is what T6 puts into the same single pass'''
		with zopen(mfile,"r") as library_file, open(self.seqid2taxid,"a") as taxidmap, open(self.prelim_map,"a") as prelim:
			for line in library_file:
				if not line.startswith(">"):
					continue
				seqid = self.sequence_id(line)
				if not seqid:
					logger.debug("Header without sequence id in {mfile}: {line}".format(mfile=mfile,line=line.strip()))
					continue
				sequences += 1
				taxid = self.accession_to_taxid.get(seqid,False)
				if not taxid:
					unannotated.add(seqid)
					logger.debug("No taxid in the database for sequence {seqid} of {mfile}".format(seqid=seqid,mfile=mfile))
					continue
				annotated += 1
				print("{seqid}\t{taxid}".format(seqid=seqid,taxid=taxid),end="\n",file=taxidmap)
				print("TAXID\t{seqid}\t{taxid}".format(seqid=seqid,taxid=taxid),end="\n",file=prelim)
		if unannotated:
			self.write_unannotated_sequences(mfile,unannotated)
		logger.info("{sequences} sequences of {mfile} are in the library, {annotated} have a taxid in the map, {n} have none".format(sequences=sequences,mfile=mfile,annotated=annotated,n=len(unannotated)))
		'''Copy the file into the library as it is (decompressed, the header of every sequence keeping
			its accession), no genome file is written from it and the input file is never modified'''
		cmd = "cat"
		if mfile.endswith(".gz"):
			cmd = "zcat"
		res = run("{cmd} {file} >> {db_path}/library/library.fna".format(cmd=cmd,
																		file=shlex.quote(mfile),
																		db_path=shlex.quote(self.krakendb)),
					shell=True)
		if res.returncode != 0:
			logger.warning("Could not add {mfile} to the library, the command returned {code}".format(mfile=mfile,code=res.returncode))
			return False
		added.put(sequences)
		return True

	def create_library_from_files(self,multifiles=False):
		'''Create library for kraken and create kraken genome2taxid map'''
		logger.info("Process datafiles and create library.fna")
		if os.path.exists(self.outdir+"/seqid2taxid.map"): open(self.seqid2taxid,"w").close() ## clean existing map if file exists
		'''The files given as packed multi sequence files (the --nt_source file and the files found by
			--multifile_prefix) are handled the same way: added to the library as they are'''
		packed_files = list(self.packed_files)
		for mfile in (multifiles or []):
			if mfile not in packed_files:
				packed_files.append(mfile)
		if self.limit:
			logger.info("Test use only {n} genomes".format(n=self.limit))
			self.genome_names = self.genome_names[0:self.limit]
		if not os.path.exists("{db_path}/library/".format(db_path=self.krakendb)):
			logger.info("Create library directory")
			os.makedirs("{db_path}/library/".format(db_path=self.krakendb))
		if not os.path.exists("{tmpdir}/library/".format(tmpdir=self.tmpdir)):
			logger.info("Create library directory")
			os.makedirs("{tmpdir}/library/".format(tmpdir=self.tmpdir))
		self.added = 0
		self.sequence_count = 0
		if self.genome_names:
			self.genome_names_split = self._split(self.genome_names,self.processes)
			self.kraken_fasta_header_multiproc(self.genome_names_split)
		else:
			logger.info("No genome files to add, the library is built from the packed files only")
		'''Merge library files into one library'''
		if (self.create_db or self.create_lib) and self.genome_names: ## Only required if kraken map is required (can be printed from db)
			batches = glob.glob("{tmpdir}/library/batch*.fasta".format(tmpdir=self.tmpdir))
			if batches:
				check_output("cat {tmpdir}/library/batch*.fasta > {db_path}/library/library.fna".format(db_path=self.krakendb, tmpdir=self.tmpdir),shell=True)
				check_output("rm {tmpdir}/library/batch*.fasta".format(db_path=self.krakendb, tmpdir=self.tmpdir),shell=True)
		if packed_files:
			logger.info("Adding sequences from {n} packed sequence files to the library".format(n=len(packed_files)))
			manager = Manager()
			sequences = manager.Queue()
			for mfile in packed_files:
				self.add_packed_file_to_library(mfile,sequences)
			self.sequence_count = sequences.qsize()
		logger.info("Number of genomes (multifiles not counted) succesfully added to the {krakenversion} database: {count}".format(count=self.added,krakenversion=self.krakenversion))
		logger.info("Number of sequences added to the {krakenversion} library from the packed files: {count}".format(count=self.sequence_count,krakenversion=self.krakenversion))
		return

	def run_program_command(self,cmd,required=True):
		'''Run one command of the external program and look at its exit status (issue #132): a kraken2
			that is not installed, or a build step that fails, only printed its error on stderr while the
			run went on logging that the database was created, leaving a database directory with nothing
			in it but the library and its map. The build step is required: the run stops there. A step
			that only adds a report file is reported and the run continues.
		------
		Parameters:
			str - the command to run
			boolean - False if a failure of this command only needs reporting
		------
		Returns:
			boolean - True when the command succeeded, False when a failure was only reported
		'''
		res = run(cmd,shell=True)
		if res.returncode != 0:
			logger.error("The command returned {code}: {cmd}".format(code=res.returncode,cmd=cmd))
			if required:
				raise CalledProcessError(res.returncode,cmd)
			return False
		return True

	def create_database(self,outdir,keep=False):
		'''For test create a small database and run tests'''
		if not os.path.exists("{krakendb}/taxonomy".format(outdir=outdir, krakendb=self.krakendb)):
			logger.info("mkdir -p {krakendb}/taxonomy".format(outdir=outdir, krakendb=self.krakendb))
			os.system("mkdir -p {krakendb}/taxonomy".format(outdir=outdir, krakendb=self.krakendb))
		if not os.path.exists("{krakendb}/library".format(krakendb=self.krakendb)):
			logger.info("mkdir -p {krakendb}/library".format(krakendb=self.krakendb))
			os.system("mkdir -p {krakendb}/library".format(krakendb=self.krakendb))
		logger.info("cp {outdir}/*.dmp {krakendb}/taxonomy".format(outdir=outdir,krakendb=self.krakendb))
		os.system("cp {outdir}/*names.dmp {krakendb}/taxonomy/names.dmp".format(outdir=outdir,krakendb=self.krakendb))
		os.system("cp {outdir}/*nodes.dmp {krakendb}/taxonomy/nodes.dmp".format(outdir=outdir,krakendb=self.krakendb))
		
		if self.krakenversion == "krakenuniq":
			os.system("cp {outdir}/*.map {krakendb}/library/prelim.map".format(outdir=outdir,krakendb=self.krakendb))
			logger.info("cp {outdir}/*.map {krakendb}/library/prelim.map".format(outdir=outdir,krakendb=self.krakendb))
		elif self.krakenversion == "kraken2":
			if not os.path.exists("{krakendb}/library/prelim_map.txt".format(outdir=outdir,krakendb=self.krakendb)):
				if os.path.exists("{outdir}/prelim_map.txt".format(outdir=outdir,krakendb=self.krakendb)):
					os.system("cp {outdir}/prelim_map.txt {krakendb}/library/prelim_map.txt".format(outdir=outdir,krakendb=self.krakendb))
					logger.info("cp {outdir}/prelim_map.txt {krakendb}/library/prelim_map.txt".format(outdir=outdir,krakendb=self.krakendb))
				else:
					os.system("cp {outdir}/*.map {krakendb}/library/prelim_map.txt".format(outdir=outdir,krakendb=self.krakendb))
					logger.info("cp {outdir}/*.map {krakendb}/library/prelim_map.txt".format(outdir=outdir,krakendb=self.krakendb))
			else:
				logger.info("prelim_map.txt exists, continue.")
		else:
			os.system("cp {outdir}/*.map {krakendb}".format(outdir=outdir,krakendb=self.krakendb))
			logger.info("cp {outdir}/*.map {krakendb}".format(outdir=outdir,krakendb=self.krakendb))
		
		if self.krakenversion == 'krakenuniq':
			logger.info(self.krakenversion+"-build --build --db {krakendb} {params} --threads {threads}".format(krakendb=self.krakendb, threads=self.build_processes, params=self.params))
			self.run_program_command(self.krakenversion+"-build --build --db {krakendb} {params} --threads {threads}".format(krakendb=self.krakendb, threads=self.build_processes, params=self.params))
		else:
			logger.info(self.krakenversion+"-build --build --skip-maps --db {krakendb} {params} --threads {threads}".format(krakendb=self.krakendb, threads=self.build_processes, params=self.params))
			self.run_program_command(self.krakenversion+"-build --build --skip-maps --db {krakendb} {params} --threads {threads}".format(krakendb=self.krakendb, threads=self.build_processes, params=self.params))

		if self.krakenversion in ["kraken2"]:
			logger.info("Create inspect file!")
			inspect_ok = self.run_program_command(self.krakenversion+"-inspect --db {krakendb} --report-zero-counts --threads {threads} > {krakendb}/inspect.txt".format(krakendb=self.krakendb,threads=self.build_processes),required=False)
			os.system("gzip {krakendb}/*.map".format(krakendb=self.krakendb))
			if inspect_ok:
				os.system("gzip {krakendb}/inspect.txt".format(krakendb=self.krakendb))
		if not keep:
			os.system(self.krakenversion+"-build --clean --db {krakendb}".format(outdir=outdir,krakendb=self.krakendb, threads=self.processes))
			## re-add taxonomy
			os.system("mkdir -p {krakendb}/taxonomy".format(outdir=outdir, krakendb=self.krakendb))
			os.system("cp {outdir}/*names.dmp {krakendb}/taxonomy/names.dmp".format(outdir=outdir,krakendb=self.krakendb))
			os.system("cp {outdir}/*nodes.dmp {krakendb}/taxonomy/nodes.dmp".format(outdir=outdir,krakendb=self.krakendb))
			# logger.info("Cleaning up tmp files")
			# os.system('find {krakendb} -maxdepth 1 -name "*.f*a" -print0 | xargs -0 rm'.format(krakendb=self.krakendb))
		logger.info("{krakenversion} database created".format(krakenversion=self.krakenversion))
