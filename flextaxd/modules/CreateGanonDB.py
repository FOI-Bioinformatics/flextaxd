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
	def __init__(self, database, ganon_database, genomes_path, outdir,verbose=False,processes=1,limit=0,dbprogram="ganon",params=""):
		super(CreateGanonDB, self).__init__()
		self.database = DatabaseFunctions(database)
		if outdir == "":
			outdir = "./"
		self.outdir = outdir
		self.params=params
		self.seqid2taxid = self.outdir+"/seqid2taxid.map"
		if os.path.exists(self.outdir+"/seqid2taxid.map"): open(self.seqid2taxid,"w").close() ## clean existing map if file exists
		self.genomes_path = genomes_path
		self.accession_to_taxid = self.database.get_genomes(self.database , limit=limit)
		self.files = []
		self.processes = processes
		self.ganondb= ganon_database
		self.ncbi_rewrite_speed = "fastest"
		self.verbose = verbose
		if verbose: print("{ganondb}".format(outdir = self.outdir, ganondb=self.ganondb))
		if not os.path.exists("{ganondb}".format(outdir = self.outdir, ganondb=self.ganondb)):
			os.system("mkdir -p {ganondb}".format(outdir = self.outdir, ganondb=self.ganondb))

	def ganon_fasta_multiproc(self,filepaths,genomes):
		'''function to run addition of genomes in paralell'''
		jobs = []
		for i in range(self.processes):
			p = Process(target=self.ganon_fasta, args=(filepaths[i],genomes[i],i))
			p.start()
			jobs.append(p)
		for job in jobs:
			job.join()
		## Concatenate each tmp file (per processor) to one big datafile
		os.system("cat {outdir}/.tmp*.gz > {ganondb}/library.fasta.gz".format(outdir=self.outdir, ganondb=self.ganondb))
		## Rm tmp files
		os.system("rm {outdir}/.tmp*.gz".format(outdir=self.outdir))
		return "Processes done"

	def ganon_fasta(self,filepaths, genomes,process):
		'''Change fasta file to contain ganon fasta header'''
		tmpname = ".tmp{rand}.gz".format(rand=process)
		tmppath = "{outdir}/{tmppath}".format(outdir=self.outdir.rstrip("/"),tmppath=tmpname)
		seqlen = 0
		for i in range(len(filepaths)):
			genome = genomes[i]
			filepath = filepaths[i]
			taxid = self.accession_to_taxid[genome]
			with zopen(filepath,"r") as f:
				'''Create a sequence to taxid mapping file according to ganon specifications
					Pre-generated file with sequence information
					(seqid <tab> seq.len <tab> taxid [<tab> assembly id])
				'''
				with open(self.seqid2taxid, "a") as seqidtotaxid:
					for line in f:
						if line.startswith(">"):
							if seqlen > 0:
								print(row[0].lstrip(">"),seqlen , taxid,end="\n", sep="\t",file=seqidtotaxid)  ## print contig name to seqtoid map
							row = line.split(" ")
							seqlen = 0
						else:
							## count lenght of sequence
							seqlen += len(line.strip())
						#print(row[0].lstrip(">"),seqlen , taxid,end="\n", sep="\t",file=seqidtotaxid)   ## print chromosome name to seqtoid map
						#seqlen = 0
					print(row[0].lstrip(">"),seqlen , taxid,end="\n", sep="\t",file=seqidtotaxid)   ## print chromosome name to seqtoid map
					seqlen = 0
			## Concatenate source files into one big file (per processor)
			os.system("cat {file} >> {tmppath}".format(file=filepath, tmppath=tmppath))
		return

	def process_folder(self):
		id_dict = self.accession_to_taxid
		self.GCF_names = []
		print("Number of genomes annotated in database",len(id_dict))
		count = 0
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
						pass
				else:
					pass

		processes = self.processes
		self.files = self.split(self.files,processes)
		self.GCF_names = self.split(self.GCF_names,processes)
		print(self.ganon_fasta_multiproc(self.files,self.GCF_names))
		print("Number of genomes added to ganon database: {count}".format(count=count))
		return self.files

	def split(self,a, n):
		k, m = divmod(len(a), n)
		return list(a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

	def create_database(self,outdir,keep=False):
		'''For test create a small database and run tests'''
		if not os.path.exists("{ganon}/taxonomy".format(outdir=outdir, ganon=self.ganondb)):
			if self.verbose: print("mkdir -p {ganon}/taxonomy".format(outdir=outdir, ganon=self.ganondb))
			os.system("mkdir -p {ganon}/taxonomy".format(outdir=outdir, ganon=self.ganondb))
		if self.verbose: print("cp {outdir}/*.dmp {ganon}/taxonomy".format(outdir=outdir,ganon=self.ganondb))
		os.system("cp {outdir}/*nodes.dmp {ganon}/taxonomy/nodes.dmp".format(outdir=outdir,ganon=self.ganondb))
		os.system("cp {outdir}/*names.dmp {ganon}/taxonomy/names.dmp".format(outdir=outdir,ganon=self.ganondb))

		if self.verbose: print("ganon build -d {ganondb} --seq-info-file {seqid2taxid} -t {processes} --taxdump-file {ganondb}/taxonomy/nodes.dmp {ganondb}/taxonomy/names.dmp -i {ganondb}/library.fasta.gz".format(ganondb=self.ganondb,seqid2taxid=self.seqid2taxid,processes=self.processes))
		os.system("ganon build -d {ganondb} --seq-info-file {seqid2taxid} -t {processes} --taxdump-file {ganondb}/taxonomy/nodes.dmp {ganondb}/taxonomy/names.dmp -i {ganondb}/library.fasta.gz".format(ganondb=self.ganondb,seqid2taxid=self.seqid2taxid,processes=self.processes))
		#if self.verbose: print(self.ganon+"-build --build --db {ganon} --threads {threads}".format(ganon=self.ganon, threads=self.processes))
		#if self.verbose: print("{ganon} database created".format(ganon=self.ganon))
