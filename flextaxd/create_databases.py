#!/usr/bin/env python3 -c

'''
This part of the script is kept because it is convenient when creating ganon kraken2 or krakenuniq databases
however it is not optimized in any way and uses os.system instead of a proper Popen for subprocessing.
'''

__author__ = "David Sundell"
__credits__ = ["David Sundell"]
__license__ = "GPLv3"
__maintainer__ = "FOI bioinformatics group"
__email__ = ["bioinformatics@foi.se","david.sundell@foi.se"]
__date__ = "2019-09-30"
__status__ = "Beta"
__pkgname__="flextaxd-create"
__github__="https://github.com/FOI-Bioinformatics/flextaxd"
from flextaxd.custom_taxonomy_databases import __version__
from importlib import import_module



## If script is executed run pipeline of selected options
def main():
	###################################--system imports--####################################
	import os, sys
	import argparse
	import time
	import logging
	if sys.version_info.major < 3 and sys.version_info.minor < 5:
		exit("This script is written for python3.5 and above please upgrade python!")

	start_time = time.time()
	current_time = start_time

	### Add base directory to path
	BASE_DIR = os.path.dirname(os.path.abspath(__file__))
	sys.path.append(BASE_DIR)
	##################################--global functions--###################################

	def dynamic_import(abs_module_path, class_name):
		module = import_module(".".join([abs_module_path,class_name]))
		target_class = getattr(module, class_name)
		return target_class

	def report_time(prev_time, final=False):
		current = time.time()
		diff = current - prev_time
		seconds = diff % 60
		minutes = int((diff - seconds)/60)
		mtext = "minute"
		if minutes != 1:
			mtext +="s"
		if final:
			logger.info("--- Time summary  {} {} {} seconds---".format(minutes,mtext, seconds))
		else:
			logger.info("--- process finished in {} {} {} seconds---\n".format(minutes,mtext, seconds))
		return current

	def write_missing(missing):
		'''Write missing genomes to file'''
		with open("{outdir}/FlexTaxD.missing".format(outdir=args.outdir), "w") as of:
			for gen in missing:
				print(gen["genome_id"], end="\n", file=of)
		return
	#########################################################################################
	##################################--error functions--###################################¤

	class InputError(Exception):
		"""InputError"""

		def __init__(self,  message,expression=""):
			self.expression = expression
			self.message = message
	#########################################################################################

	'''
		Supported programs
	'''

	programs = ["kraken2", "krakenuniq","ganon"]

	parser = argparse.ArgumentParser()
	basic = parser.add_argument_group('basic', 'Basic commands')
	basic.add_argument('-o', '--outdir',metavar="", default=".", help="Output directory (same directory as custom_taxonomy_databases dump)")
	basic.add_argument('-db', '--database',metavar="", type=str, default=".ctdb" , help="Custom taxonomy sqlite3 database file")

	### Download options, process local directory and potentially download files
	download_opts = parser.add_argument_group('download_opts', "Download and file handling")
	download_opts.add_argument('-p', '--processes',metavar="",type=int, default = 8, help="Use multiple cores for downloading genomes and kraken if -kp is not set")
	download_opts.add_argument('--download', action='store_true', help="Download additional sequences")
	download_opts.add_argument('--force_download', action='store_true', help="Download sequences from genbank if not in refseq (WARNING: might include genome withdrawals)")
	download_opts.add_argument('--genomes_path', metavar="",default=None,  help='path to genomes')



	###  Kraken options not needed for public version this script is made to export names and nodes.dmp files
	classifier_opts = parser.add_argument_group('classifier_opts', "Classifier options")
	classifier_opts.add_argument('--create_db', action ='store_true',  help="Start create db after loading databases")
	classifier_opts.add_argument('--dbprogram', metavar="", default="kraken2", choices=programs,help="Select one of the supported programs ["+", ".join(programs)+"]")
	classifier_opts.add_argument('--db_name', metavar="",default=None, help="database directory (fullpath)")
	classifier_opts.add_argument('--params', metavar="", default="",  help="Add extra params to create command (supports kraken*)")
	classifier_opts.add_argument('--test', action='store_true', help="test database structure, only use 100 seqs")
	classifier_opts.add_argument('--keep', action='store_true', help="Keep temporary files")
	classifier_opts.add_argument('--skip', metavar="", default="", help="Do not include genomes within this taxonomy (child tree) in the database (works for kraken)")
	classifier_opts.add_argument('-kp', '--build_processes',metavar="",type=int, default = None, help="Use a different number of cores for kraken classification")


	debugopts = parser.add_argument_group("Logging and debug options")
	#debugopts.add_argument('--tmpdir',             metavar='', default="/tmp/FlexTaxD",            help="Specify reference directory")
	debugopts.add_argument('--logs',                 metavar='', default="logs/",         help="Specify log directory")
	debugopts.add_argument('--verbose',            action='store_const', const=logging.INFO,                help="Verbose output")
	debugopts.add_argument('--debug',                action='store_const', const=logging.DEBUG,                help="Debug output")
	debugopts.add_argument('--supress',                action='store_const', const=logging.ERROR,    default=logging.WARNING,            help="Supress warnings")

	parser.add_argument("--version", action='store_true', help=argparse.SUPPRESS)

	args = parser.parse_args()

	if not os.path.exists(args.database):
		raise FileNotFoundError("No database file could be found, please provide a FlexTaxD database to run FlexTaxD!")

	if args.version:
		print("{name}: version {version}".format(name=__pkgname__,version=__version__))
		print("Maintaner group: {maintaner} ({email})".format(maintaner=__maintainer__,email=", ".join(__email__)))
		print("Github: {github}".format(github=__github__))
		exit()

	'''Log file and verbose options'''
	logval = args.supress
	if args.debug:
		logval = args.debug
	elif args.verbose:
		logval = args.verbose

	from datetime import date
	from datetime import time as dtime
	t = dtime()
	today = date.today()
	logpath = args.logs+"FlexTaxD-create-"+today.strftime("%b-%d-%Y")+"{}.log"
	if os.path.exists(logpath):
		logpath=logpath.format("-{:%H:%M}".format(t))
	else: logpath = logpath.format("")


	logging.basicConfig(
			level=logval,
			format="%(asctime)s %(module)s [%(levelname)-5.5s]  %(message)s",
			handlers=[
				logging.FileHandler(logpath),
				logging.StreamHandler()
			])
	logger = logging.getLogger(__name__)
	logger.info("FlexTaxD-create logging initiated!")
	logger.debug("Supported formats: {formats}".format(formats=programs))

	'''
		Process data
	'''
	if args.outdir:
		if not os.path.exists(args.outdir):
			os.system("mkdir -p {outdir}".format(outdir = args.outdir))

	''' 1. Process genome_path directory'''
	process_directory = dynamic_import("modules", "ProcessDirectory")
	logger.info("Processing files; create kraken seq.map")
	process_directory_obj = process_directory(args.database)
	genomes, missing = process_directory_obj.process_folder(args.genomes_path)
	''' 2. Download missing files'''
	if args.download:
		download = dynamic_import("modules", "DownloadGenomes")
		download_obj = download(args.processes,outdir=args.outdir,force=args.force_download)
		still_missing = download_obj.run(missing)
		if len(still_missing) > 0: print("Not able to download: {nr}".format(nr=len(still_missing)))
	else:
		if len(missing) > 0:
			logger.info("Genome annotations with no matching source: {nr}".format(nr=len(missing)))
			write_missing(missing)
	''' 3. Add genomes to database'''
	if args.db_name:
		if args.dbprogram.startswith("kraken"):
			logger.info("Loading module: CreateKrakenDatabase")
			classifier = dynamic_import("modules", "CreateKrakenDatabase")
		else:
			logger.info("Loading module: CreateGanonDB")
			classifier = dynamic_import("modules", "CreateGanonDB")
		limit = 0
		if args.test:
			limit = 500
		'''Use the genome -> path dictionary to build database'''
		genomes = process_directory_obj.get_genome_path_dict()
		classifierDB = classifier(args.database, args.db_name, genomes,args.outdir,
										create_db=args.create_db,
										limit=limit,
										dbprogram=args.dbprogram,
										params=args.params,
										skip=args.skip,
										processes=args.processes,
										build_processes=args.build_processes,
										debug=args.debug,
										verbose=args.verbose,
		)
		report_time(current_time)
		classifierDB.create_library_from_files()
		logger.info("Genome folder preprocessing completed!")

	''' 4. Create database'''
	if args.create_db:
		report_time(current_time)
		logger.info("Create database")
		try:
			classifierDB.create_database(args.outdir,args.keep)
		except UnboundLocalError:
			logger.error("#Error: No kraken database name was given!")
			exit()

	logger.debug(report_time(start_time,final=True))

if __name__ == '__main__':
	main()
