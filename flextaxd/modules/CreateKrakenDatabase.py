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
from modules.functions import download_genome
from time import sleep

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
	def __init__(self, database, kraken_database, genomes_path, outdir,verbose=False,processes=1,limit=0,dbprogram="kraken2",params="",skip="",create_db=False,debug=False,kraken_processes=None,skip_download=False):
		super(CreateKrakenDatabase, self).__init__()
		self.krakenversion = dbprogram
		self.database = DatabaseFunctions(database)
		if outdir == "":
			outdir = "./"
		self.outdir = outdir
		self.seqid2taxid = self.outdir+"/seqid2taxid.map"
		if os.path.exists(self.outdir+"/seqid2taxid.map"): open(self.seqid2taxid,"w").close() ## clean existing map if file exists
		self.genomes_path = genomes_path
		self.accession_to_taxid = self.database.get_genomes(self.database , limit=limit)
		self.files = []
		self.params = params
		self.processes = processes
		self.skip_download = skip_download ## If true do not download extra genomes from NCBI
		self.genome_path = {}
		if not kraken_processes:
			self.kraken_processes = self.processes
		else:
			self.kraken_processes = kraken_processes
		self.krakendb= kraken_database
		self.verbose = verbose
		self.create_db = create_db

		self.debug = debug
		if self.debug:
			self.taxidmap_debug = self.database.get_nodes()
		self.seqhead_validator = {}
		self.seqhead_count = 0
		if skip:
			self.skiptax = parse_skip(skip.split(","))  ## if node should be skipd this must be true, otherwise nodes in modfile are added to existing database
		else:
			self.skiptax = []
		logger.info("{krakendb}".format(outdir = self.outdir, krakendb=self.krakendb))
		if not os.path.exists("{krakendb}".format(outdir = self.outdir, krakendb=self.krakendb)):
			os.system("mkdir -p {krakendb}".format(outdir = self.outdir, krakendb=self.krakendb))

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
		for i in range(len(genomes)):
			genome = genomes[i]
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
					logger.info("#Warning kraken header could not be added to {gcf}! Total: {count}".format(gcf=genome,count=count))
					count +=1
					continue
			if self.create_db:
				if taxid not in self.skiptax:
					'''Open temp file for manipulated (unzipped) genome fasta files'''
					tmpfile = zopen(tmppath,"w")
					taxidlines = []  ## Holder for sequences to be added to the seqid2taxid map from each file
					'''Open input genome fasta file'''
					with zopen(filepath,"r") as f:
						for line in f:
							if line.startswith(">"):
								taxidmap = []
								row = line.strip().split(" ")
								seq_header = row[0]
								if self.debug:
									self.seqhead_validator[seq_header] = [seq_header, genome, str(taxid),self.taxidmap_debug[taxid],filepath]
									self.seqhead_count +=1
								if len(row) == 1:
									row.append("")
								endhead = " ".join(row[1:])
								if True:  ## Remove end heading in all files as it may contain bad chars
									endhead = "\n"
								'''Format sequence header'''
								if not self.krakenversion == "krakenuniq":
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

	def process_folder(self):
		id_dict = self.accession_to_taxid
		self.genome_names = []
		logger.info("Number of genomes annotated in database {n}".format(n=len(id_dict)))
		count = 1
		downloaded_n = 0
		ref_ext = [".fna"]
		oth_ext = [".fasta",".fa"]
		ext = ref_ext+oth_ext
		self.notused = set()
		download_files = []
		logger.info("Process genome path ({path})".format(path=self.genomes_path))
		for root, dirs, files in os.walk(self.genomes_path,followlinks=True):
			for file in files:
				if count % 1000 == 0:
					print("Processed {count} genomes".format(count=count), end="\r")
				fname = file.strip(".gz") ## remove gz if present
				if fname.endswith(tuple(ext)):
					count +=1
					if (fname.startswith("GCF_") or fname.startswith("GCA_")) and fname.endswith(tuple(ref_ext)):
						## Rename file to only GCF_name as this is what would be present
						genome_name = fname.split("_",2)
						genome_name = genome_name[0]+"_"+genome_name[1]
					else:
						genome_name = fname.rstrip(".fastan") ## strip .fa .fasta or .fna from name
					try:
						'''If file contains GCF in the start and fna in the end it is most likely
							from refseq try match do genomeid2taxid dictionary
							otherwise the name should match the full name of the sequence file minus fa,fasta,fna
						'''
						taxid = id_dict[genome_name.strip()]
					except KeyError:
						try:
							'''The file was neither a refseq file nor a custom fasta or fa file,
								try if the file is a genbank file (GCA instead of GCF although not 100%)
							'''
							if genome_name.strip().endswith("GCA"):
								genome_name = genome_name.strip().replace("GCA","GCF")
							else:
								genome_name = genome_name.strip().replace("GCF","GCA")
							taxid = id_dict[genome_name]
						except KeyError:
							#### Swap back namechange if no success
							if genome_name.strip().endswith("GCA"):
								genome_name = genome_name.strip().replace("GCA","GCF")
							else:
								genome_name = genome_name.strip().replace("GCF","GCA")
							'''Try to download the genome from official source'''
							try:
								'''Final try, does the file have a GCF/GCA start ends with fna but is still a custom named genome'''
								taxid = id_dict[fname.rsplit(".",1)[0]] ## strip .fa .fasta or .fna from base filename
							except KeyError:
								try:
									taxid = id_dict[fname] ## strip .fa .fasta or .fna from base filename
								except KeyError:
									'''This file had no match in the reference folder, perhaps it is not annotated in the database'''
									self.notused.add(genome_name)
									logger.debug("#Warning {gcf} could not be matched to a database entry!".format(gcf=genome_name.strip()))
									continue
					if not file.endswith("from_genomic.fna.gz"):
						filepath = os.path.join(root, file)  ## Save the path to the file
						logger.debug("File added")
						#self.files.append(filepath)
						self.genome_names.append(genome_name.strip())
						self.genome_path[genome_name.strip()] = filepath
				elif file == "MD5SUMS" or file.endswith(".txt"):
					pass
				else:
					logger.debug("#Warning {gcf} does not have a valid file ending".format(gcf=file))
		#self.files = list(set(self.files))
		self.genome_names = list(set(self.genome_names))
		'''Try to download unknown files'''
		existing_genomes = []
		'''Must add both GCA and GCF name of each genome since both might exist'''
		for genome in self.genome_names:
			genome = genome.split("_",1)[1]
			existing_genomes += ["GCA_"+genome,"GCF_"+genome]
		for file_not_present in set(id_dict.keys()) - set(existing_genomes):
			download_files.append({"genome_id":file_not_present,"outdir":self.genomes_path+"downloads/"})

		if len(download_files) > 0 and not self.skip_download:
			logger.info("Found {count} files, {dn} files not present, download files".format(count=len(self.files),dn=len(download_files)))
			sleep(1)
			#logger.debug(download_files)
			np = self.processes
			if int(self.processes) > 25: ## 50 is Maximum allowed simultanous connections to ncbi
				np = 25
			self.download_map = self._split(download_files,np)
			#print(len(self.download_map))
			logger.info("Using {np} parallel processes to NCBI".format(np=np))



			'''function to run download of genomes in paralell'''
			logger.info("Processing files; create kraken seq.map")
			jobs = []
			manager = Manager()
			added = manager.Queue()
			fpath = manager.Queue()
			for i in range(np):
				p = Process(target=download_genome, args=(self.download_map[i],added,fpath))
				p.daemon=True
				p.start()
				jobs.append(p)
			for job in jobs:
				job.join()
			self.added = added.qsize()

			while True:
				l = len(self.genome_names)
				self.genome_names += added.get()
				self.genome_path[genome_name.strip()] = fpath.get()
				if l == len(self.genome_names):
					break

		#self.file_split = self._split(self.files,self.processes)
		self.genome_names_split = self._split(self.genome_names,self.processes)
		if not os.path.exists("{db_path}/library/".format(db_path=self.krakendb)):
			logger.info("Create library directory")
			os.mkdir("{db_path}/library/".format(db_path=self.krakendb))
		self.kraken_fasta_header_multiproc(self.genome_names_split)
		'''Merge library files into one library'''
		check_output("cat {db_path}/library/batch*.fasta > {db_path}/library/library.fna".format(db_path=self.krakendb),shell=True)
		check_output("rm {db_path}/library/batch*.fasta".format(db_path=self.krakendb),shell=True)

		if self.verbose and len(self.notused) > 0: logger.info("#Warning {gcf} genomes in the input folder could not be matched to any database entry!".format(gcf=len(self.notused)))
		logger.info("Number of genomes succesfully added to the {krakenversion} database: {count}".format(count=self.added,krakenversion=self.krakenversion))
		logger.info("Downloaded genomes {downloaded}".format(downloaded=downloaded_n))
		return self.files

	def _split(self,a, n):
		k, m = divmod(len(a), n)
		return list(a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

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
		logger.info(self.krakenversion+"-build --build --db {krakendb} {params} --threads {threads}".format(krakendb=self.krakendb, threads=self.kraken_processes, params=self.params))
		os.system(self.krakenversion+"-build --build --skip-maps --db {krakendb} {params} --threads {threads}".format(krakendb=self.krakendb, threads=self.kraken_processes, params=self.params))
		if not keep:
			pass
			## Since kraken2 is removing too much on clean it might be better to do this manually so certain log files can be saved
			#os.system(self.krakenversion+"-build --clean --db {krakendb}".format(outdir=outdir,krakendb=self.krakendb, threads=self.processes))
			#logger.info("Cleaning up tmp files")
			#os.system('find {krakendb} -maxdepth 1 -name "*.f*a" -print0 | xargs -0 rm'.format(krakendb=self.krakendb))
		logger.info("{krakenversion} database created".format(krakenversion=self.krakenversion))
