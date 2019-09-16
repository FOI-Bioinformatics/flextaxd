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
__date__ = "2019-03-11"
__status__ = "Beta"
__pkgname__="custom-taxonomy-databases"
__github__="https://github.com/davve2/custom-taxonomy-databases"

###################################--system imports--####################################
import os, sys
import argparse
from importlib import import_module
import time

## If script is executed run pipeline of selected options
def main():
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
            print("--- Time summary  {} {} {} seconds---".format(minutes,mtext, seconds))
        else:
            print("--- process finished in {} {} {} seconds---\n".format(minutes,mtext, seconds))
        return current
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
    basic.add_argument('-v','--verbose', action='store_true',   help="verbose output")
    basic.add_argument('--log', metavar="",default = None,   help="use a log file instead of stdout")
    basic.add_argument('-db', '--database',metavar="", type=str, default=".ctdb" , help="Custom taxonomy sqlite3 database file")

    ###  Kraken options not needed for public version this script is made to export names and nodes.dmp files
    classifier_opts = parser.add_argument_group('classifier_opts', "Classifier options")
    classifier_opts.add_argument('--dbprogram', metavar="", default="kraken2", choices=programs,help="Select one of the supported programs ["+", ".join(programs)+"]")
    classifier_opts.add_argument('--genomes_path', metavar="",default=None,  help='path to genomes')
    classifier_opts.add_argument('--db_name', metavar="",default=None, help="database directory (fullpath)")
    classifier_opts.add_argument('--create_db', action ='store_true',  help="Start create db after loading databases")
    classifier_opts.add_argument('--params', metavar="", default="",  help="Add extra params to create command (supports kraken*)")
    classifier_opts.add_argument('--test', action='store_true', help="test database structure, only use 100 seqs")
    classifier_opts.add_argument('-p', '--processes',metavar="",type=int, default = 8, help="Use multiple cores")
    classifier_opts.add_argument('--keep', action='store_true', help="Keep temporary files")
    classifier_opts.add_argument('--skip', metavar="", default="", help="Do not include genomes within this taxonomy (child tree) in the database (works for kraken)")

    parser.add_argument("--version", action='store_true', help=argparse.SUPPRESS)

    args = parser.parse_args()

    if args.version:
        print("{name}: version {version}".format(name=__pkgname__,version=__version__))
        print("Maintaner group: {maintaner} ({email})".format(maintaner=__maintainer__,email=", ".join(__email__)))
        print("Github: {github}".format(github=__github__))
        exit()

    '''Global vars'''
    global verbose
    verbose = args.verbose
    modify_module = False

    '''Log file and verbose options'''
    if args.log:
        global original_sysout
        original_sysout = sys.stdout
        if os.path.exists(args.log) and not force:
            ans = input("Warning the logfile already exist, overwrite? (y/n)")
            if ans not in ["y", "Y", "yes", "Yes"]:
                exit("Abort logfile already exists, remove file or use --force")
        if verbose:
            print("All log output will be written to {file}".format(file=args.log))
        logfile = open(args.log, "w")
        sys.stdout = logfile

    if verbose:
        print("Script parameters:")
        print(args)


    '''
        Process data
    '''
    if args.outdir:
        if not os.path.exists(args.outdir):
            os.system("mkdir -p {outdir}".format(outdir = args.outdir))

    ''' 3. Add genomes to kraken DB'''
    if args.db_name:
        if args.dbprogram.startswith("kraken"):
            if verbose: print("Loading module: CreateKrakenDatabase")
            classifier = dynamic_import("modules", "CreateKrakenDatabase")
        else:
            if verbose: print("Loading module: CreateGanonDB")
            classifier = dynamic_import("modules", "CreateGanonDB")
        limit = 0
        if args.test:
            limit = 500
        classifierDB = classifier(args.database, args.db_name, args.genomes_path,args.outdir,verbose=verbose,processes=args.processes,limit=limit,dbprogram=args.dbprogram,params=args.params,skip=args.skip)
        if verbose: current_time = report_time(current_time)
        classifierDB.process_folder()
        print("Done")

    ''' 4. Create database'''
    if args.create_db:
        if verbose: current_time = report_time(current_time)
        if verbose: print("Create database")
        try:
            classifierDB.create_database(args.outdir,args.keep)
        except UnboundLocalError:
            exit("#Error: No kraken database name was given!")

    if verbose: report_time(start_time,final=True)

    if args.log:
        sys.stdout = original_sysout

if __name__ == '__main__':
    main()
