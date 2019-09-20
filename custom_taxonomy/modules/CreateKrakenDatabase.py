#!/usr/bin/env python3 -c

'''
This part of the script is kept because it is convenient when creating kraken2 or krakenuniq databases
however it is not optimized in any way and uses os.system instead of a proper Popen for subprocessing.
'''

import sys
import gzip
import random
import os
from multiprocessing import Process, Queue
from subprocess import Popen,PIPE
from .database.DatabaseConnection import DatabaseFunctions

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
	def __init__(self, database, kraken_database, genomes_path, outdir,verbose=False,processes=1,limit=0,dbprogram="kraken2",params="",skip=""):
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
		if skip:
			self.skiptax = parse_skip(skip.split(","))  ## if node should be skipd this must be true, otherwise nodes in modfile are added to existing database
		else:
			self.skiptax = []
		if verbose: print("{krakendb}".format(outdir = self.outdir, krakendb=self.krakendb))
		if not os.path.exists("{krakendb}".format(outdir = self.outdir, krakendb=self.krakendb)):
			os.system("mkdir -p {krakendb}".format(outdir = self.outdir, krakendb=self.krakendb))

	def kraken_fasta_header_multiproc(self,filepaths,genomes):
		'''function to run addition of genomes in paralell'''
		jobs = []
		for i in range(self.processes):
			p = Process(target=self.kraken_fasta_header, args=(filepaths[i],genomes[i]))
			p.start()
			jobs.append(p)
		for job in jobs:
			job.join()
		return "Processes done"

	def kraken_fasta_header(self,filepaths, genomes):
		'''Change fasta file to contain kraken fasta header'''
		for i in range(len(filepaths)):
			genome = genomes[i]
			filepath = filepaths[i]
			kraken_header = "kraken:taxid"
			tmpname = ".tmp{rand}.gz".format(rand=random.randint(10**7,10**9))
			tmppath = "{outdir}/{tmppath}".format(outdir=self.outdir.rstrip("/"),tmppath=tmpname).rstrip(".gz")#self.outdir+"/.tmp{rand}.gz".format(rand=random.randint(10**7,10**9))
			tmpfile = zopen(tmppath,"w")
			taxid = self.accession_to_taxid[genome]
			if taxid not in self.skiptax:
				with zopen(filepath,"r") as f:
					with open(self.seqid2taxid, "a") as seqidtotaxid:
						for line in f:
							if line.startswith(">"):
								row = line.strip().split(" ")
								if len(row) == 1:
									row.append("")
								if not self.krakenversion == "krakenuniq":
									line = row[0] + "|" + kraken_header + "|" + str(taxid) + "  " + " ".join(row[1:])	 ## If kraken2 add kraken header (nessesary?)
								print(row[0].lstrip(">"), taxid, " ".join(row[1:]),end="", sep="\t",file=seqidtotaxid)   ## print chromosome name to seqtoid map
							print(line, end="", file=tmpfile)
				tmpfile.close()
				output = self.krakendb.rstrip("/")+"/"+filepath.split("/")[-1].rstrip(".gz")
				os.rename(tmppath, output)
				os.system(self.krakenversion+"-build --add-to-library {file} --db {krakendb} >/dev/null 2>&1".format(file=output,krakendb=self.krakendb))
			# command = self.krakenversion+"-build"
			# parameters = " --add-to-library {file} --db {krakendb}".format(file=output,krakendb=self.krakendb)
			# p = Popen((command, parameters) ,stdout=PIPE,stderr=PIPE)
			# p.wait()
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
		self.GCF_names = []
		print("Number of genomes annotated in database",len(id_dict))
		count = 0
		self.notused = set()
		for root, dirs, files in os.walk(self.genomes_path,followlinks=True):
			for file in files:
				if file.endswith(".fna.gz"):
					GCF_name = file.split("_",2)
					GCF_name = GCF_name[0]+"_"+GCF_name[1]
					try:
						taxid = id_dict[GCF_name.strip()]
						filepath = os.path.join(root, file)
						self.files.append(filepath)
						self.GCF_names.append(GCF_name.strip())
						count+=1
					except KeyError:
						self.notused.add(GCF_name)
						if self.verbose: print("#Warning {gcf} could not be matched to database entry!".format(gcf=GCF_name.strip()))
				else:
					pass
		processes = self.processes
		self.files = self.split(self.files,processes)
		self.GCF_names = self.split(self.GCF_names,processes)
		self.kraken_fasta_header_multiproc(self.files,self.GCF_names)
		if self.verbose: print("#Warning {gcf} genomes could not be matched to database entry!".format(gcf=len(self.notused)))
		print("Number of genomes added to {krakenversion} database: {count}".format(count=count,krakenversion=self.krakenversion))
		return self.files

	def split(self,a, n):
		k, m = divmod(len(a), n)
		return list(a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

	def create_database(self,outdir,keep=False):
		'''For test create a small database and run tests'''
		if not os.path.exists("{krakendb}/taxonomy".format(outdir=outdir, krakendb=self.krakendb)):
			if self.verbose: print("mkdir -p {krakendb}/taxonomy".format(outdir=outdir, krakendb=self.krakendb))
			os.system("mkdir -p {krakendb}/taxonomy".format(outdir=outdir, krakendb=self.krakendb))
		if self.verbose: print("cp {outdir}/*.dmp {krakendb}/taxonomy".format(outdir=outdir,krakendb=self.krakendb))
		os.system("cp {outdir}/*names.dmp {krakendb}/taxonomy/names.dmp".format(outdir=outdir,krakendb=self.krakendb))
		os.system("cp {outdir}/*nodes.dmp {krakendb}/taxonomy/nodes.dmp".format(outdir=outdir,krakendb=self.krakendb))
		if self.krakenversion != "kraken2": os.system("cp {outdir}/*.map {krakendb}".format(outdir=outdir,krakendb=self.krakendb))
		if self.verbose: print("cp {outdir}/*.map {krakendb}".format(outdir=outdir,krakendb=self.krakendb))
		if self.verbose: print(self.krakenversion+"-build --build --db {krakendb} {params} --threads {threads}".format(krakendb=self.krakendb, threads=self.processes, params=self.params))
		os.system(self.krakenversion+"-build --build --db {krakendb} {params} --threads {threads}".format(krakendb=self.krakendb, threads=self.processes, params=self.params))
		if not keep:
			## Since kraken2 is removing too much on clean it might be better to do this manually so certain log files can be saved
			#os.system(self.krakenversion+"-build --clean --db {krakendb}".format(outdir=outdir,krakendb=self.krakendb, threads=self.processes))
			if self.verbose: print("Cleaning up tmp files")
			os.system('find {krakendb} -maxdepth 1 -name "*.fna" -print0 | xargs -0 rm'.format(krakendb=self.krakendb))
		if self.verbose: print("{krakenversion} database created".format(krakenversion=self.krakenversion))
