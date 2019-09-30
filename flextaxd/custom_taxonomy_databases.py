#!/usr/bin/env python3 -c

'''
This script is built to handle import and modification of different taxonomy sources
v0.1.3 (NCBI, QIIME(GTDB), CanSNPer)

The custom_taxonomy_databases script contains a parser and a database handler to allow customization of databases
from NCBI, QIIME or CanSNPer sources and supports export functions into NCBI formatted names and nodes.dmp files
as well as a slimmed tab separated file. The database allows databases to be merged at selected
nodes(taxonomy IDs) as well as adding resolution to certain subgroups (ie using a tab separated file).

The script was initially written to allow the use of GTDB with some custom modifications to allow separation of
subgroups. GTDB was created by an Australian group aimed to restructure the taxonomy relation from the NCBI
taxonomy tree to strictly follow a phylogenetic structure (http://gtdb.ecogenomic.org/) this script can use the
taxonomy.tsv files from the GTDB downloads page as input (with the --taxonomy_type selected as QIIME). By default
the script will read a Tab separated file containing parent and child (defined by column headers).
All data is kept in a sqlite3 database (.ctdb by default) and can be dumped at will to NCBI formatted names
and nodes.dmp files. Supported export formats in version 0.2b is NCBI and TSV). The TSV dump format is similar to
the NCBI dump except that it contains a header (parent/child), has parent on the left and only uses tab to separate
each column (not <tab>|<tab>).
'''

__version__ = "0.1.3"
__author__ = "David Sundell"
__credits__ = ["David Sundell"]
__license__ = "GPLv3"
__maintainer__ = "FOI bioinformatics group"
__email__ = ["bioinformatics@foi.se","david.sundell@foi.se"]
__date__ = "2019-09-11"
__status__ = "Beta"
__pkgname__="custom_taxonomy_databases"
__github__="https://github.com/FOI-Bioinformatics/flextaxd"
__programs_supported__ = ["kraken2", "krakenuniq","ganon","centrifuge"]

###################################--system imports--####################################
import os, sys
import argparse
from importlib import import_module
import time

## If script is executed run pipeline of selected options
def main():
    if sys.version_info.major < 3 and sys.version_info.minor < 5:
        exit("This script is written for python3 please upgrade python!")

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

    def get_read_modules():
        '''Find (ReadTaxonomy) modules and return options'''
        read_modules = []
        for script in os.listdir(BASE_DIR+"/modules"):
            if script.startswith("ReadTaxonomy"):
                modname = script.lstrip("ReadTaxonomy").rstrip(".py")
                if modname != "":
                    read_modules.append(modname)
        return read_modules
    #########################################################################################
    ##################################--error functions--###################################¤

    class InputError(Exception):
        """InputError"""

        def __init__(self,  message,expression=""):
            self.expression = expression
            self.message = message


    #########################################################################################

    parser = argparse.ArgumentParser()

    required = parser.add_argument_group('required', 'Required')
    required.add_argument('-db', '--database',metavar="", type=str, default=".ftd" , help="Custom taxonomy sqlite3 database file (fullpath)")

    basic = parser.add_argument_group('basic', 'Basic commands')
    basic.add_argument('-o', '--outdir',metavar="", default=".", help="Output directory")
    basic.add_argument('-v','--verbose', action='store_true',   help="verbose output")
    basic.add_argument('--log', metavar="", default = None,   help="use a log file instead of stdout")
    basic.add_argument("--dump", action='store_true', help="Write database to names.dmp and nodes.dmp")
    basic.add_argument('--dump_mini', action='store_true', help="Dump minimal file with tab as separator")
    basic.add_argument("--force", action='store_true', help="use when script is implemented in pipeline to avoid security questions on overwrite!")

    rmodules = get_read_modules()
    read_opts = parser.add_argument_group('read_opts', "Source options")
    read_opts.add_argument('-tf', '--taxonomy_file',metavar="", default=None, help="Taxonomy source file")
    read_opts.add_argument('-tt', '--taxonomy_type',metavar="", default="", choices=rmodules, help="Source format of taxonomy input file ({modules})".format(modules=",".join(rmodules)))
    read_opts.add_argument('--taxid_base', metavar="", type=int, default=1, help="The base for internal taxonomy ID numbers, when using NCBI as base select base at minimum 3000000 (default = 1)")

    mod_opts = parser.add_argument_group('mod_opts', "Database modification options")
    mod_opts.add_argument('-mf','--mod_file', metavar="", default=False, help="File contaning modifications parent,child,(taxonomy level)")
    mod_opts.add_argument('-md', '--mod_database', metavar="",default=False, help="Database file containing modifications")
    mod_opts.add_argument('-gt', '--genomeid2taxid', metavar="", default=False, help="File that lists which node a genome should be assigned to")
    mod_opts.add_argument('-gp', '--genomes_path', metavar="",default=None,  help='Path to genome folder is required when using NCBI_taxonomy as source')
    mod_opts.add_argument('-p', '--parent',metavar="", default=False, help="Parent from which to add (replace see below) branch")
    mod_opts.add_argument('--replace', action='store_true', help="Add if existing children of parents should be removed!")

    out_opts = parser.add_argument_group('output_opts', "Output options")
    out_opts.add_argument('--dbprogram', metavar="", default=False,choices=__programs_supported__, help="Adjust output file to certain output specifications ["+", ".join(__programs_supported__)+"]")
    out_opts.add_argument("--dump_prefix", metavar="", default="names,nodes", help="change dump prefix reqires two names default(names,nodes)")
    out_opts.add_argument('--dump_sep', metavar="", default="\t|\t", help="Set output separator default(NCBI) also adds extra trailing columns for kraken")
    out_opts.add_argument('--dump_descriptions', action='store_true', default=False, help="Dump description names instead of database integers")

    parser.add_argument("--version", action='store_true', help=argparse.SUPPRESS)

    args = parser.parse_args()

    if args.version:
        print("{name}: version {version}".format(name=__pkgname__,version=__version__))
        print("Maintaner group: {maintaner} ({email})".format(maintaner=__maintainer__,email=", ".join(__email__)))
        print("Github: {github}".format(github=__github__))
        exit()

    force = False

    if args.force:
        force = True
    if (os.path.exists(args.database) and args.taxonomy_file) and not force:
        ans = input("Warning: {database} already exists, overwrite? (y/n): ".format(database=args.database))
        if ans not in ["y","Y","yes", "Yes"]:
            exit("Database already exists, abort!")
        force = True

    '''Check argument input logics'''
    if args.mod_file or args.mod_database:
        if not args.parent:
            raise InputError("Argument --parent is required when updating the database!")
        if args.mod_file and not args.genomeid2taxid:
            raise InputError("Argument --mod_file with a genomeid to nodeid map is required when adding new nodes to the database!")

    if args.genomeid2taxid and args.taxonomy_type == "NCBI":
        if not args.genomes_path:
            raise InputError("To annotate genomes to the NCBI database a path to genbank or refseq genomes folder needs to be given --genomes_path")

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

    if force and args.taxonomy_file:
        '''Remove database if force is turned on and new source file is given'''
        if os.path.exists(args.database):
            os.system("rm {database}".format(database=args.database))

    if verbose:
        print("Script parameters:")
        print(args)

    '''
        Custom taxonomy databases pipeline
    '''
    ''' Create output directory if it does not exist!'''
    if args.outdir:
        if not os.path.exists(args.outdir):
            os.system("mkdir -p {outdir}".format(outdir = args.outdir))

    ''' 0. Create taxonomy database (if it does not exist)'''
    if args.taxonomy_file:
        if not os.path.exists(args.database) or force:
            '''Load taxonomy module'''
            if verbose: print("Loading module: ReadTaxonomy{type}".format(type=args.taxonomy_type))
            read_module = dynamic_import("modules", "ReadTaxonomy{type}".format(type=args.taxonomy_type))
            read_obj = read_module(args.taxonomy_file, database=args.database)
            read_obj.verbose = verbose                                                          ## Set verbose mode
            if verbose: print("Parse taxonomy")
            read_obj.parse_taxonomy()                                                           ## Parse taxonomy file

            '''Parse genome2taxid file'''                                                       ## Fix at some point only one function should be needed
            if not args.genomeid2taxid:
                print("Warning no genomeid2taxid file given!")
            elif args.taxonomy_type == "NCBI" and args.genomeid2taxid:
                read_obj.parse_genomeid2taxid(args.genomes_path,args.genomeid2taxid)
            elif args.taxonomy_type == "CanSNPer":
                read_obj.parse_genomeid2taxid(args.genomeid2taxid)

            if verbose: print("Nodes in taxonomy tree {n} number of taxonomies {k}".format(n=read_obj.length, k=read_obj.ids))
            if verbose: current_time = report_time(current_time)

    ''' 1. Modify database, if datasource for modification of current database is supplied process this data'''
    if args.mod_file or args.mod_database:
        if args.mod_file and not args.genomeid2taxid:
            exit("No genomeid2taxid file given!")
        if verbose: print("Loading module: ModifyTree")
        modify_module = dynamic_import("modules", "ModifyTree")
        modify_obj = modify_module(database=args.database, mod_file=args.mod_file, mod_database= args.mod_database,verbose=verbose,parent=args.parent,replace=args.replace)
        modify_obj.update_database()
        if args.mod_file:
            if verbose: current_time = report_time(current_time)
            modify_obj.update_annotations(genomeid2taxid=args.genomeid2taxid)
        if verbose: current_time = report_time(current_time)

    ''' 2. Dump custom taxonomy database into NCBI/kraken readable format)'''
    if args.dump or args.dump_mini:
        '''Check if datase exists if it does make sure the user intends to overwrite the file'''
        nameprefix,nodeprefix = args.dump_prefix.split(",")
        if (os.path.exists(args.outdir.rstrip("/")+"/"+nameprefix+".dmp") or os.path.exists(args.outdir.rstrip("/")+"/"+nameprefix+".dmp")) and not force:
            if args.log: sys.stdout = original_sysout ## If log file is used, make sure the input reaches the user terminal
            ans = input("Warning: {names} and/or {nodes} already exists, overwrite? (y/n): ")
            if ans not in ["y","Y","yes", "Yes"]:
                exit("Dump already exists, abort!")
            if args.log:  sys.stdout = logfile ## Continue to write to log file

        '''Create print out object'''
        if verbose: print("Loading module: WriteTaxonomy".format(type=args.taxonomy_type))
        write_module = dynamic_import("modules", "WriteTaxonomy")
        write_obj = write_module(args.outdir, database=args.database,verbose=verbose,prefix=args.dump_prefix,separator=args.dump_sep,minimal=args.dump_mini,desc=args.dump_descriptions,dbprogram=args.dbprogram)

        '''Print database to file'''
        if args.taxonomy_type == "NCBI":
            write_obj.set_minimal()
        write_obj.nodes()
        write_obj.names()
        if False: #args.taxDB:
            write_obj.set_separator("\t")
            write_obj.set_prefix("names,taxDB")
            write_obj.set_order(True)
            write_obj.nodes()

    if verbose: report_time(start_time,final=True)

    if args.log:
        sys.stdout = original_sysout

if __name__ == '__main__':
    main()
