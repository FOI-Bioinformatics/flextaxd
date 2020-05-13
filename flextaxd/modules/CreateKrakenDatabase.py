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
	def __init__(self, database, kraken_database, genomes_path, outdir,verbose=False,processes=1,limit=0,dbprogram="kraken2",params="",skip="",create_db=False,debug=False):
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
		self.krakendb= kraken_database
		self.ncbi_rewrite_speed = "fastest"
		self.verbose = verbose
		self.create_db = create_db
		self.debug = debug
		if skip:
			self.skiptax = parse_skip(skip.split(","))  ## if node should be skipd this must be true, otherwise nodes in modfile are added to existing database
		else:
			self.skiptax = []
		logger.info("{krakendb}".format(outdir = self.outdir, krakendb=self.krakendb))
		if not os.path.exists("{krakendb}".format(outdir = self.outdir, krakendb=self.krakendb)):
			os.system("mkdir -p {krakendb}".format(outdir = self.outdir, krakendb=self.krakendb))

	def kraken_fasta_header_multiproc(self,filepaths,genomes):
		'''function to run addition of genomes in paralell'''
		jobs = []
		manager = Manager()
		added = manager.Queue()
		for i in range(self.processes):
			p = Process(target=self.kraken_fasta_header, args=(filepaths[i],genomes[i],added))
			p.daemon=True
			p.start()
			jobs.append(p)
		for job in jobs:
			job.join()
		self.added = added.qsize()
		return "Processes done"

	def kraken_fasta_header(self,filepaths, genomes,added):
		'''Change fasta file to contain kraken fasta header'''
		count = 0

		for i in range(len(filepaths)):
			genome = genomes[i]
			filepath = filepaths[i]
			kraken_header = "kraken:taxid"
			tmpname = ".tmp{rand}.gz".format(rand=random.randint(10**7,10**9))
			tmppath = "{outdir}/{tmppath}".format(outdir=self.outdir.rstrip("/"),tmppath=tmpname).rstrip(".gz")#self.outdir+"/.tmp{rand}.gz".format(rand=random.randint(10**7,10**9))
			'''Get taxid from database'''
			try:
				taxid = self.accession_to_taxid[genome]
			except KeyError:
				logger.info("#Warning kraken header could not be added to {gcf}! Total: {count}".format(gcf=genome,count=count))
				count +=1
				continue
			if self.create_db:
				if taxid not in self.skiptax:
					tmpfile = zopen(tmppath,"w")
					taxidlines = []
					with zopen(filepath,"r") as f:
						for line in f:
							if line.startswith(">"):
								row = line.strip().split(" ")
								if len(row) == 1:
									row.append("")
								if not self.krakenversion == "krakenuniq":
									line = row[0] + "|" + kraken_header + "|" + str(taxid) + "  " + " ".join(row[1:])	 ## Nessesary to be able to add sequences outside NCBI
								taxidlines.append("\t".join([row[0].lstrip(">"), str(taxid), " ".join(row[1:])])) ## print chromosome name to seqtoid map
							print(line, end="", file=tmpfile)
					tmpfile.close()  ## Close and save tmp file
					with open(self.seqid2taxid, "a") as seqidtotaxid:
						print("\n".join(taxidlines),end="\n",file=seqidtotaxid)
					output = self.krakendb.rstrip("/")+"/"+filepath.split("/")[-1].rstrip(".gz")
					os.rename(tmppath, output)
					try:
						ans = check_output(self.krakenversion+"-build --add-to-library {file} --db {krakendb}  >/dev/null 2>&1 ".format(file=output,krakendb=self.krakendb),shell=True)
						added.put(1)
					except CalledProcessError as e:
						logger.debug("Could not process {output}, the following error occured:\n{e}".format(output=output,e=e))
		return

	def get_skip_list(self):
		'''get taxonomy ids not wanted'''
		return self.skiptax

	def parse_skip(self,skip):
		'''Skip allows a list of tax ids to be passed and then excluded from the database'''
		skiptax = []
		for t in skip:
			skiptax += self.taxonomydb.get_children()
		return skiptax

	def process_folder(self):
		id_dict = self.accession_to_taxid
		self.genome_names = []
		logger.info("Number of genomes annotated in database {n}".format(n=len(id_dict)))
		count = 0
		downloaded_n = 0
		ref_ext = [".fna"]
		oth_ext = [".fasta",".fa"]
		ext = ref_ext+oth_ext
		self.notused = set()
		download_files = []
		for root, dirs, files in os.walk(self.genomes_path,followlinks=True):
			for file in files:
				fname = file.strip(".gz") ## remove gz if present
				if fname.endswith(tuple(ext)):
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
								'''This file had no match in the reference folder, perhaps it is not annotated in the database'''
								self.notused.add(genome_name)
								logger.debug("#Warning {gcf} could not be matched to a database entry!".format(gcf=genome_name.strip()))
								continue
					if not file.endswith("from_genomic.fna.gz"):
						filepath = os.path.join(root, file)  ## Save the path to the file
						logger.debug("File added")
						self.files.append(filepath)
						self.genome_names.append(genome_name.strip())
						count+=1
				elif file == "MD5SUMS" or file.endswith(".txt"):
					pass
				else:
					logger.debug("#Warning {gcf} does not have a valid file ending".format(gcf=file))
		self.files = list(set(self.files))
		self.genome_names = list(set(self.genome_names))
		'''Try to download unknown files'''
		existing_genomes = []
		'''Must add both GCA and GCF name of each genome since both might exist'''
		for genome in self.genome_names:
			genome = genome.split("_",1)[1]
			existing_genomes += ["GCA_"+genome,"GCF_"+genome]
		for file_not_present in set(id_dict.keys()) - set(existing_genomes):
			download_files.append({"genome_id":file_not_present,"outdir":self.genomes_path+"downloads/"})

		if len(download_files) > 0:
			logger.info("Found {count} files, {dn} files not present, download files".format(count=len(self.files),dn=len(download_files)))
			np = self.processes
			if int(self.processes) > 25: ## 50 is Maximum allowed simultanous connections to ncbi
				np = 25
			p = Pool(np)
			for output, error in p.imap(download_genome, download_files):
				if error is None:
					logger.debug("{genome} downloaded!".format(genome=output))
					downloaded_n +=1
					## Add the downloaded genome to the count
					count+=1
					if downloaded_n % 100 == 0:
						logger.info("downloaded {n} genomes".format(n=downloaded_n))
					try:
						for file in os.listdir(glob.glob(os.path.join(self.genomes_path,"downloads",output+"*"))[0]):
							if file.endswith("_genomic.fna.gz") and "from" not in file:
								self.files.append(os.path.join(self.genomes_path,"downloads",file))
								gname = os.path.basename(file)
								genome_name = gname.split("_",2)
								genome_name = genome_name[0]+"_"+genome_name[1]
					except:
						logger.debug("Download error {output}".format(output=output))
						continue
					self.genome_names.append(genome_name.strip())
				else:
					pass #logger.error(error)
		self.files = self.split(self.files,self.processes)
		self.genome_names = self.split(self.genome_names,self.processes)
		self.kraken_fasta_header_multiproc(self.files,self.genome_names)
		if self.verbose and len(self.notused) > 0: logger.info("#Warning {gcf} genomes in the database could not be matched to any database entry!".format(gcf=len(self.notused)))
		logger.info("Number of genomes succesfully added to the {krakenversion} database: {count}".format(count=self.added,krakenversion=self.krakenversion))
		logger.info("Downloaded genomes {downloaded}".format(downloaded=downloaded_n))
		return self.files

	def split(self,a, n):
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
		else: os.system("cp {outdir}/*.map {krakendb}".format(outdir=outdir,krakendb=self.krakendb))
		logger.info("cp {outdir}/*.map {krakendb}".format(outdir=outdir,krakendb=self.krakendb))
		logger.info(self.krakenversion+"-build --build --db {krakendb} {params} --threads {threads}".format(krakendb=self.krakendb, threads=self.processes, params=self.params))
		os.system(self.krakenversion+"-build --build --skip-maps --db {krakendb} {params} --threads {threads}".format(krakendb=self.krakendb, threads=self.processes, params=self.params))
		if not keep:
			## Since kraken2 is removing too much on clean it might be better to do this manually so certain log files can be saved
			#os.system(self.krakenversion+"-build --clean --db {krakendb}".format(outdir=outdir,krakendb=self.krakendb, threads=self.processes))
			logger.info("Cleaning up tmp files")
			os.system('find {krakendb} -maxdepth 1 -name "*.f*a" -print0 | xargs -0 rm'.format(krakendb=self.krakendb))
		logger.info("{krakenversion} database created".format(krakenversion=self.krakenversion))
