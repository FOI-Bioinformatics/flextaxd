#!/usr/bin/env python3 -c

'''
Read NCBI taxonomy dmp files (nodes or names) and holds a dictionary
'''

from .ReadTaxonomy import ReadTaxonomy
from . import functions
from . import ranks
import zlib
import os
import logging
logger = logging.getLogger(__name__)

class ReadTaxonomyNCBI(ReadTaxonomy):
	"""docstring for ReadTaxonomy."""
	'''the NCBI tree arrives with its taxids and is read straight into the tree table, it never looks a node
	   up by name, so it does not pay for the name dictionary of issue #145'''
	match_nodes = False

	def __init__(self, taxonomy_file=False, database=False,**kwargs):
		super(ReadTaxonomyNCBI, self).__init__(database=database)
		self.taxonomy_file = taxonomy_file
		self.names_dmp = taxonomy_file.replace("nodes","names")
		self.names = {}
		self.length = 0
		self.ids = 0
		self.accessionfile = False
		'''the option to keep sequences bound to the source they came from is optional: every caller of the
		   command line passes it, but requiring the keyword made the reader impossible to construct without
		   it (issue #144)'''
		self.force_multisource = kwargs.get("force_multisource", False)
		self.refseqid_to_GCF = {}
		self.refseqid_reference = {}
		self.assembly_seqids = set()

	def write_missing(self,missing,source=False):
		'''Write the genomes a build found in the genome folder but could not annotate next to the
			database, as the other finding files of a build do; nothing is written when nothing is
			missing. The database keeps only the number, the list stays in the file'''
		if not missing:
			return False
		missing_file = "{database}.unannotated".format(database=self.database.database)
		with open(missing_file, "w") as of:
			for gen in missing:
				print(gen, end="\n", file=of)
		self.database.save_report("skipped_import",
									source or self.taxonomy_file or "genomes",
									len(missing),
									missing_file)
		logger.warning("{n} genomes got no annotation, they are listed in {file}".format(n=len(missing),file=missing_file))
		return missing_file

	def sequence_id(self,header):
		'''The id of a sequence is the first white space separated word of its FASTA header, the
			header may hold no description and end with the line ending'''
		words = header.split()
		if not words:
			return False
		return words[0].lstrip(b">")

	def set_accession_file(self,file):
		self.accessionfile = file

	def parse_taxonomy(self,treefile=False):
		'''Parse taxonomy information'''
		if hasattr(self, "names_dmp"):
			logger.info("Parse names file {}".format(self.names_dmp))
			self.read_names(self.names_dmp)
			self.length = self.database.num_rows("nodes")
		if self.taxonomy_file:
			logger.info("Parse nodes file {}".format(self.taxonomy_file))
			self.read_nodes(self.taxonomy_file)
			self.ids = self.database.num_rows("tree")

	def read_nodes(self, taxfile):
		'''Read NCBI node file and store a dictionary with names in object

			An id the file links but no names file describes still needs a row in the nodes table: the
			link below points at it, and a tree link to a node the nodes table does not have makes
			--validate fail and a merge of the database crash on that row (issue #152), so such an id
			gets a placeholder node named after itself and a warning naming the count
		'''
		linked = set()
		with open(taxfile, "r") as _taxfile:
			for taxonomy_row in _taxfile:
				data = taxonomy_row.strip().split("\t|\t")
				child,parent,rank = data[0],data[1],data[2]
				'''the root's own row "1 | 1 | no rank" asks for the one self-link the import already added
				   when it set up the root, under whatever rank name the file gives it; a nodes file that
				   repeats that row or self-links any other id would break the convention of issue #162'''
				if child == parent:
					logger.warning("#Skipped: the nodes file links {node} to itself, the root's own link is added by the import itself".format(node=child))
					continue
				'''an empty rank field or the text None names nothing: that level stays without a rank instead
				   of being called "no rank", which is the name NCBI itself gives the taxids it does not place
				   in a rank (issue #142)'''
				if ranks.unnamed(rank):
					rank = None
				'''every rank name the file brings in is stored and kept, one level per name (issue #138)'''
				if rank not in self.rank:
					self.add_rank(rank)
				lev = self.rank.get(rank)
				self.database.add_link(child=data[0],parent=data[1],rank=lev)
				linked.add(data[0])
				linked.add(data[1])
			self.database.commit()
		self.add_unnamed_nodes(linked)
		self.rank_summary()
		return

	def add_unnamed_nodes(self, linked):
		'''Give every id that is linked by the taxonomy file but described nowhere a node of its own,
			named after the taxid, so no tree link points at a missing row (issue #152)

		------
		Returns:
			list - the taxids that got a placeholder node, sorted
		'''
		described = set()
		for taxid in self.database.get_nodes(col=1).keys():
			described.add(int(taxid))
		unnamed = sorted({int(taxid) for taxid in linked} - described)
		if not unnamed:
			return []
		logger.warning("{n} taxids are linked by the taxonomy file but described by no names file, they get a node named after the taxid: {first}".format(
						n=len(unnamed), first=",".join(str(taxid) for taxid in unnamed[:5])))
		for taxid in unnamed:
			self.database.add_node(str(taxid), id=taxid)
			self.taxonomy[str(taxid)] = taxid
		self.database.commit()
		return unnamed

	def read_names(self, taxfile):
		'''Read NCBI node file and store a dictionary with nodes in object'''
		with open(taxfile, "r") as _taxfile:
			for taxonomy_row in _taxfile:
				data = taxonomy_row.strip().split("\t|\t")
				taxid = data[0]
				name = data[1]
				_type = False
				if len(data) > 3:
					_type = data[3].rstrip("|\t")
				if _type == "scientific name" or not _type:
					self.add_node(name, id=taxid)
			self.database.commit()
		return

	def parse_genebank_file(self,filepath,filename,reference="refseq"):
		'''Read the sequence id of a single genome file and map it to the genome id of that file.
			The mapping of a sequence that has a genome file always wins over the entry a packed
			archive gives for the same accession, also when the archive is read first.
		'''
		logger.debug("Parse file {filename}".format(filename=filename))
		genebankid = filename.split("_",2)
		genebankid = genebankid[0]+"_"+genebankid[1]
		'''The sequence names are compared as bytes, so the file is read binary, compressed or not'''
		with self.zopen(filepath,"rb") as f:
			refseqid = self.sequence_id(f.readline())
		if not refseqid:
			logger.warning("No sequence id in the first line of {file}, the genome stays unannotated".format(file=filename))
			return False
		if refseqid in self.assembly_seqids and self.refseqid_to_GCF.get(refseqid) != genebankid:
			logger.warning("Sequence {seqid} is the first sequence of both {first} and {second}, the second genome stays unannotated".format(seqid=refseqid.decode("utf-8"),first=self.refseqid_to_GCF[refseqid],second=genebankid))
		self.refseqid_to_GCF[refseqid] = genebankid
		self.refseqid_reference[refseqid] = reference
		self.assembly_seqids.add(refseqid)
		return

	def parse_nt_file(self,filepath,filename,reference="nt"):
		'''Parse large file with all nt entries, only the sequences that no genome file in the
			genome folder carries get the packed file as their reference'''
		logger.info("Processing {file}".format(file=filename))
		covered_by_genomes = 0
		with self.zopen(filepath,"rb") as f:
			for row in f:
				if row.startswith(b">"):
					seqid = self.sequence_id(row)
					if not seqid:
						continue
					if seqid in self.assembly_seqids:
						covered_by_genomes += 1
						continue
					self.refseqid_to_GCF[seqid] = seqid.decode("utf-8")
					self.refseqid_reference[seqid] = reference
					#self.database.add_genome(genome=seqid,_id=taxid.decode("utf-8"),reference="nt")
		if covered_by_genomes:
			logger.info("{n} sequences of {file} are also present as genome files and keep the annotation of those files".format(n=covered_by_genomes,file=filename))
		return

	def get_annotation_columns(self,headers):
		'''Find the columns holding the sequence id and the taxid in an NCBI accession2taxid file,
			the header of the file decides the layout:
				accession | accession.version | taxid | gi
				accessions | accessions.versioned | accessions.nonredundant_GI_seqVersion | taxid
				gi | accession.version | taxid
			The versioned accession column is the one that matches the sequence names in the genome
			files, the gi column must never be read as a taxid.
		------
		Returns
			tuple - (accession column, taxid column)
		'''
		accession_col,taxid_col = False,False
		for col,name in enumerate([header.strip().lower() for header in headers]):
			if name.endswith(b"taxid"):
				taxid_col = col
			elif (name.startswith(b"accession") and b"version" in name and b"nonredundant" not in name
					and accession_col is False):
				accession_col = col
		if accession_col is False:
			for col,name in enumerate([header.strip().lower() for header in headers]):
				if name.startswith(b"accession"):
					accession_col = col
					break
		if taxid_col is False or accession_col is False:
			raise TypeError("Could not find the accession and taxid columns in the header of the annotation file: {headers}".format(headers=b"\t".join(headers)))
		return (accession_col,taxid_col)

	def parse_genomeid2taxid(self, genomes_path,annotation_file,reference="refseq"):
		'''To allow NCBI databases to be build from scratch the sequences names needs to be stored in the database,
			this function parses the accession2taxid file from NCBI to speed up the function and reduce the amount
			of stored datata only sequences in input genomes_path will be fetched
			Genome files are read before the packed archives of the same folder, so the genome id and
			reference a genome file gives its sequence is never replaced by the entry an archive holds
			for the same accession.
		'''
		file_endings = (".fna",".fa",".fasta")
		logger.info("Parsing ncbi accession2taxid, genome_path: {dir}".format(dir = genomes_path))
		self.refseqid_to_GCF = {}
		self.refseqid_reference = {}
		self.assembly_seqids = set()
		genome_files, archive_files = [], []
		for root, dirs, files in os.walk(genomes_path,followlinks=True):
			for filename in sorted(files):
				base_name = functions.strip_gz(filename)
				if not base_name.endswith(file_endings):
					continue
				filepath = os.path.join(root, filename)
				one_genome_per_file = filename.startswith("GC") and base_name.endswith(".fna")
				if one_genome_per_file or not (filename.startswith("nt") or self.force_multisource):
					genome_files.append((filepath,filename))
				else:
					archive_files.append((filepath,filename))
		for filepath, filename in genome_files:
			self.parse_genebank_file(filepath,filename,reference=reference)
		for filepath, filename in archive_files:
			if not self.force_multisource:
				stats = os.stat(filepath)
				logger.info("Parsing nt archive, filesize: {filesize}Mb, compressed: {compressed}".format(compressed=filename.endswith(".gz"),filesize=(stats.st_size / (1024 * 1024))))
			self.parse_nt_file(filepath,filename,reference="nt")
		logger.info("genomes folder read, {n} sequence files found".format(n=len(self.refseqid_to_GCF)))
		if not annotation_file.endswith("accession2taxid.gz"):
			raise TypeError("The supplied annotation file does not seem to be the ncbi nucl_gb.accession2taxid.gz")
		annotated_genome = set()
		try:
			with self.zopen(annotation_file,"rb") as f:
				headers = f.readline().split(b"\t")
				accession_col,taxid_col = self.get_annotation_columns(headers)
				logger.info("Annotation file columns: accession {acc}, taxid {tax}".format(acc=headers[accession_col].strip().decode("utf-8"),tax=headers[taxid_col].strip().decode("utf-8")))
				for row in f:
					if row.strip() == b"": ## If there are trailing empty lines in the file
						continue
					columns = row.rstrip(b"\n").split(b"\t")
					if len(columns) <= max(accession_col,taxid_col):
						logger.info(row)
						logger.info(columns)
						if len(annotated_genome) > 0:
							logger.info("Potential error in last row?")
						else:
							logger.info("Error on first line in annotation file, check format!")
						continue
					refseqid = columns[accession_col].strip()
					taxid = columns[taxid_col].strip()
					try:
						genebankid = self.refseqid_to_GCF[refseqid]
						self.database.add_genome(genome=genebankid,_id=taxid.decode("utf-8"),reference=self.refseqid_reference.get(refseqid,reference))
						annotated_genome.add(refseqid)
					except KeyError:
						pass
				self.database.commit()
		except zlib.error as e:
			logger.info("Error in annotation file {e}".format(e=e))
		missing = set(self.refseqid_to_GCF.keys()) - annotated_genome
		missing = sorted(set(self.refseqid_to_GCF[m] for m in missing)) ## Translate to genome ids
		if missing:
			self.write_missing(missing,source=annotation_file)
		logger.debug(missing)  ## If debug also print genomes to terminal
		logger.info("Genomes not matching any annotation {len}".format(len=len(missing)))
		return missing
