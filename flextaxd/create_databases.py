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
__date__ = "2020-10-30"
__status__ = "Beta"
__pkgname__="flextaxd-create"
__github__="https://github.com/FOI-Bioinformatics/flextaxd"
from flextaxd.custom_taxonomy_databases import __version__
from flextaxd.modules.functions import read_skip_file
from importlib import import_module
import shutil

latest_genome_reps = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz"

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
        with open("{outdir}/FlexTaxD.missing".format(outdir=args.tmpdir), "w") as of:
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
    basic.add_argument('-db', '--database', '--db' ,metavar="", type=str, default=".ctdb" , help="Custom taxonomy sqlite3 database file")
    #basic.add_argument('--dump_map',metavar="", type=str, default=False , help="dump kraken2 prelim and seq2taxid maps, required for files with multiseq")
    basic.add_argument('--dump_map',action='store_true',help="dump kraken2 prelim and seq2taxid maps, required for files with multiseq")
    basic.add_argument('-nt', '--nt_source', '--nt', metavar="",type=str, default=False, help="If part of your data is merged into one file, add path to file here, example nt")
    basic.add_argument('-mfp', '--multifile_prefix', metavar="",type=str,default=False, help="If multiple datafiles, list file prefix to handle as multi files")

    ### Download options, process local directory and potentially download files
    download_opts = parser.add_argument_group('download_opts', "Download and file handling")
    download_opts.add_argument('-p', '--processes',metavar="",type=int, default = 8, help="Use multiple cores for downloading genomes and kraken if -kp is not set")
    download_opts.add_argument('--download', action='store_true', help="Download additional sequences")
    download_opts.add_argument('-r', '--representative',   action='store_true', help="Download GTDB representative genomes")
    download_opts.add_argument('--download_file', default=False, help="Download genomes from file (path to file)")
    download_opts.add_argument('--rep_path', metavar="URL", default=latest_genome_reps, help="Specify GTDB representative version URL full path")
    download_opts.add_argument('--force_download', action='store_true', help="Download sequences from genbank if not in refseq (WARNING: might include genome withdrawals)")
    download_opts.add_argument('--genomes_path', metavar="",default=None,  help='path to genomes')

    ###  Kraken options not needed for public version this script is made to export names and nodes.dmp files
    classifier_opts = parser.add_argument_group('classifier_opts', "Classifier options")
    classifier_opts.add_argument('--create_db', action ='store_true',  help="Start create db after loading databases")
    classifier_opts.add_argument('--create_lib', action ='store_true',  help="Add to create kraken library file but no DB")
    classifier_opts.add_argument('--dbprogram','--db_program', metavar="", default="kraken2", choices=programs,help="Select one of the supported programs ["+", ".join(programs)+"]")
    classifier_opts.add_argument('--db_name', metavar="",default=None, help="database directory (fullpath)")
    classifier_opts.add_argument('--params', metavar="", default="",  help="Add extra params to create command (supports kraken*)")
    classifier_opts.add_argument('--test', action='store_true', help="test database structure, only use 100 seqs")
    classifier_opts.add_argument('--keep', action='store_true', help="Keep temporary files")
    classifier_opts.add_argument('--skip', metavar="", default="", help="Do not include genomes within this taxonomy (child tree) in the database (works for kraken), can be a file ending with txt and genome ids one per row")
    classifier_opts.add_argument('-kp', '--build_processes',metavar="",type=int, default = None, help="Use a different number of cores for kraken classification")


    debugopts = parser.add_argument_group("Logging and debug options")
    #debugopts.add_argument('--tmpdir',             metavar='', default="/tmp/FlexTaxD",            help="Specify reference directory")
    debugopts.add_argument('--logs',                 metavar='', default="logs/",         help="Specify log directory")
    debugopts.add_argument('--verbose',            action='store_const', const=logging.INFO,                help="Verbose output")
    debugopts.add_argument('--debug',                action='store_const', const=logging.DEBUG,                help="Debug output")
    debugopts.add_argument('--supress',                action='store_const', const=logging.ERROR,    default=logging.WARNING,            help="Supress warnings")

    debugopts.add_argument('--tmpdir',                 metavar='', default="tmp/",         help="Specify temp directory")
    parser.add_argument("--version", action='store_true', help=argparse.SUPPRESS)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    
    if args.version:
            print("{name}: version {version}".format(name=__pkgname__,version=__version__))
            print("Maintaner group: {maintaner} ({email})".format(maintaner=__maintainer__,email=", ".join(__email__)))
            print("Github: {github}".format(github=__github__))
            exit()
    
    if not os.path.exists(args.database):
        raise FileNotFoundError("No database file could be found, please provide a FlexTaxD database to run FlexTaxD!")


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
    if not os.path.exists(args.logs):
        os.mkdir(args.logs)
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

    if args.dump_map:
        logger.info("Dump genome maps")
        write_module = dynamic_import("modules", "WriteTaxonomy")
        write_obj = write_module(args.outdir, database=args.database, dump_genome_map=True)
        write_obj.dump_taxid_map()
        exit()

    if args.create_db and not (args.genomes_path or args.nt_source):
        raise InputError("genomes_path parameter was not given")
    if not args.nt_source:
        if args.genomes_path != None and not os.path.exists(args.genomes_path):
            ans = input("Warning: directory for downloaded genomes does not exist, do you want to create it? (y/n): ")
            if ans not in ["y","Y","yes", "Yes"]:
                exit('terminating...')
            os.makedirs(args.genomes_path)

    '''Check if temp directory exists, otherwise create directory'''
    if not os.path.exists(args.tmpdir):
        os.makedirs(args.tmpdir)

    '''
        Process data
    '''
    if True:  # Export names and nodes in expected format
        '''Check if datase exists if it does make sure the user intends to overwrite the file'''
        lskip = False
        dump_prefix = "names,nodes"
        nameprefix,nodeprefix = dump_prefix.split(",")
        if (os.path.exists(args.tmpdir.rstrip("/")+"/"+nameprefix+".dmp") or os.path.exists(args.tmpdir.rstrip("/")+"/"+nameprefix+".dmp")):
            ans = input("Warning: {names} and/or {nodes} already exists, overwrite, yes/cancel/no? (y/c/n): ")
            if ans in ["c","C","cancel", "Cancel"]:
                exit("Dump already exists, abort!")

            if ans in ["n","N","no","No"]: 
                lskip = True
        if not lskip:  #print names and nodes in all cases except no to overwrite
            '''Create print out object'''
            write_module = dynamic_import("modules", "WriteTaxonomy")
            write_obj = write_module(args.tmpdir, database=args.database,prefix=dump_prefix,dbprogram=args.dbprogram,dump_genome_map=True)
            '''Print database to file'''
            write_obj.nodes()
            write_obj.names()
            write_obj.dump_taxid_map()

    if args.outdir:
        if not os.path.exists(args.outdir):
            os.system("mkdir -p {outdir}".format(outdir = args.outdir))
    skip=False
    if os.path.exists("{db_path}/library/library.fna".format(db_path=args.db_name)):
        ans = input("Database library file already exist, (u)se library, (o)verwrite (c)ancel? (u o,c): ")
        if ans in ["o", "O"]:
            logger.info("Overwrite current build progress")
            shutil.rmtree("{db_path}".format(db_path=args.db_name))
        elif ans in ["m", "M"]:
            logger.info("Merge input fasta with library.fna")
            cmd = "cat"
            if args.nt_source:
                if args.nt_source.endswith(".gz"):
                    cmd = "zcat"    
                os.system("{cmd} {large_source} >> {db_path}/library/library.fna >> {large_source}".format(cmd=cmd, db_path=args.db_name, large_source=args.nt_source))
            else:
                exit("No large input file given")
        elif ans.strip() in ["u", "U"]:
            logger.info("Resume database build")
            skip = True
        else:
            exit("Cancel execution!")
    else:
        if args.nt_source:
            #ans = input("Im here: ")
            cmd = "cp"
            cmd2 = ""
            if args.nt_source.endswith(".gz"):
                cmd = "zcat"
                cmd2 = " > "
            #print("{cmd} {large_source} {cmd2} {db_path}/library/library.fna".format(db_path=args.db_name, cmd=cmd,cmd2=cmd2,large_source=args.nt_source))
            #ans = input("Write c to continue")
            #if ans != "c":
            #    exit()
            #os.makedirs("{db_path}/library".format(db_path=args.db_name))
            #os.system("{cmd} {large_source} {cmd2} {db_path}/library/library.fna".format(db_path=args.db_name, cmd=cmd,cmd2=cmd2,large_source=args.nt_source))
            genomes=[args.nt_source]
            args.genomes_path=os.path.dirname(args.nt_source)
            skip=True

    ''' 1. Process genome_path directory'''
    if not skip:
        process_directory = dynamic_import("modules", "ProcessDirectory")
        logger.info("Processing files; create kraken seq.map")
        if args.multifile_prefix:
            process_directory_obj = process_directory(args.database,args.multifile_prefix.split(","))
        else:
            process_directory_obj = process_directory(args.database)
        genomes, missing = process_directory_obj.process_folder(args.genomes_path)
        multifiles =False
        if args.multifile_prefix:
            ## Process files including multifiles
            multifiles = process_directory_obj.get_multifiles() # rescan genomes directory
        while missing: # Experimental implementation: make user wary of missing genomes and force user to enter "no" in prompt to continue with missing genomes.
            ''' 2. Download missing files'''
            # If there are missing genome files, asks the user to attempt a download of these from gtdb
            download_prompted = False
            if missing and not (args.download or args.representative or args.download_file):
                print('There is a discrepancy of genomes found in the database and the specified genome-folder, {numFound} genomes were found and {numMissing} genomes are missing.'.format(numFound=len(genomes),numMissing=len(missing)))
                print('You may want to purge your database from missing genomes using "flextaxd --purge_database"')
                if not args.download: # Dont ask to download if the user already specified via flag to download
                    ans = input('Do you want to download these genomes from NCBI? (y)es, (n)o, (c)ancel: ')
                    if ans in ["y","Y","yes", "Yes"]:
                        download_prompted = True
                    elif type(ans) == str and not ans.isdigit() and ans.lower() in ('c','cancel',):
                        write_missing(missing)
                        print("Genomes were missing, you can find the names of these genomes in file FlexTaxD.missing, located in ./tmp directory (you might need to specify --keep)")
                        exit('Terminated.')
                    else:
                        print('Will naivly proceed to construct database. Genomes may be missing.')
                        break
            #/
            if args.download or args.representative or args.download_file or download_prompted:
                download = dynamic_import("modules", "DownloadGenomes")
                download_obj = download(args.processes,outdir=args.tmpdir,force=args.force_download,download_path=args.genomes_path)
                if 0 and 'old download using ncbi-genome-download. to be removed once new way is proven to be stable' and download_prompted:
                    logger.info('Download of individual genomes from NCBI, num='+str(len(missing)))
                    download_obj.download_files(missing)
                    # Move downloaded files to genomes-directory
                    for path,dirs,files in os.walk(args.genomes_path+'/'+'downloads'):
                        for file_ in files:
                            if file_[:2] == 'GC' and file_.endswith('.gz'):
                                new_file_name = '_'.join(file_.split('_')[:2])+'.fna.gz'
                                os.rename(path+'/'+file_,args.genomes_path+'/'+new_file_name)
                    try:
                        shutil.rmtree(args.genomes_path+'/'+'downloads')
                    except:
                        logger.info('no genomes were downloaded, expeted for download was: '+str(len(missing)))
                    genomes, missing = process_directory_obj.process_folder(args.genomes_path)
                    #/
                elif 1 and 'new download using ncbi "datasets" software' and download_prompted:
                    download = dynamic_import("modules", "DownloadGenomes")
                    download_obj = download(args.processes,outdir=args.tmpdir,force=args.force_download,download_path=args.genomes_path)
                    logger.info('Download of genomes from NCBI using "datasets", num='+str(len(missing)))
                    download_obj.download_genomes(missing) # perform download
                    genomes, missing = process_directory_obj.process_folder(args.genomes_path) # rescan genomes directory
                elif args.download_file:
                    download_obj.download_from_file(args.download_file)
                else:
                    ans = input('Will attempt to get genomes from '+str(args.rep_path+'. Proceed? (y/n) '))
                    if not ans in ["yes","Yes","y","Y"]:
                        exit("Terminated.")
                    new_genome_path, missing = download_obj.run(missing,args.rep_path)
                    if not new_genome_path:
                        still_missing = missing
                        if len(still_missing) > 0: print("Not able to download: {nr}".format(nr=len(still_missing)))
                    else:
                        new_genomes, missing = process_directory_obj.process_folder(new_genome_path)
                        genomes += new_genomes
                

        #if not (args.download or args.representative or args.download_file) and len(missing) > 0:
        #    print("Warning: was unable to locate/download missing genomes. Found/missing: {numFound}/{numMissing}".format(numFound=len(genomes),numMissing=len(missing)))
        #    print("You can find the names of these genomes in file FlexTaxD.missing, located in ./tmp directory (you might need to specify --keep)")
        #    logger.info("Genome annotations with no matching source: {nr}".format(nr=len(missing)))
        #    write_missing(missing)

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
            limit = 10
        '''Use the genome -> path dictionary to build database'''
        if not skip:
            logger.info("Get genomes from input directory!")
            genomes = process_directory_obj.get_genome_path_dict()
        else: genomes=False
        if args.skip:
            if args.skip.endswith(".txt"):
                args.skip = read_skip_file(args.skip)
                logger.info("File passed to skip, {n} genomes and {x} taxids added to skiplist".format(n=len(args.skip["genome_id"]),x=len(args.skip["tax_id"])))
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
                                        tmpdir=args.tmpdir,
                                        create_lib=args.create_lib
        )
        report_time(current_time)
        if not skip:
            logger.info("Create library.fna")
            classifierDB.create_library_from_files(multifiles)
        logger.info("Genome folder preprocessing completed!")

    ''' 4. Create database'''
    if args.create_db:
        report_time(current_time)
        logger.info("Create database")
        try:
            classifierDB.create_database(args.tmpdir,args.keep)
        except UnboundLocalError as e:
            logger.error("#Error: No database name was given!")
            logger.error("#UnboundLocalError "+e)
            exit()

    '''Clean up temp files'''
    if os.path.exists(args.tmpdir) and not args.keep:
        logger.info("Cleaning up temp directory")
        shutil.rmtree(args.tmpdir)


    logger.debug(report_time(start_time,final=True))

if __name__ == '__main__':
    main()
