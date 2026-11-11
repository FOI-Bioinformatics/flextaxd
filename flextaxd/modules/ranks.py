#!/usr/bin/env python3 -c

'''
The rank levels of the taxonomy: one place saying which rank names are defined.

A level of the tree is either named by a rank or it is not named at all, and the database keeps the two
apart:

  - the rank of a link is the id of a row in the rank table, one row per name, however that name arrived:
    an NCBI nodes.dmp rank, an MP/GTDB level code, a SILVA file, a custom tree. An import that brings a
    rank name in gets it stored (issue #138 stopped it from storing the same name twice);
  - a link that nothing names carries rank_i NULL, and the dmp export writes an empty rank field for it
    (issue #134). Below the species that is what kraken2 wants: it numbers those levels itself (S1, S2).

`no rank` is not the same as unnamed: it is the name NCBI gives the taxids it does not place in a rank
(238,697 of those links in NCBI_base.db on disk), so it is stored as the data it is whenever the input
says it. It was also always the first row of the rank table of an NCBI derived database, which is what
made rank_i 1 the silent default of `add_link` (issue #142): a level that had no rank landed on the NCBI
`no rank` level instead of staying without one.
'''

__version__ = "1.0"
__author__ = "David Sundell"
__credits__ = ["David Sundell"]
__license__ = "GPLv3"
__maintainer__ = "FOI bioinformatics group"
__email__ = "bioinformatics@foi.se"
__date__ = "2026-11-06"
__status__ = "Production"

NO_RANK = "no rank"

"""The rank names NCBI uses. Read off the taxdump in the databases on disk (NCBI_base.db holds these 44),
with "domain" added: the MP/GTDB reader names that level "d" -> domain and the GTDB databases on disk
carry it. A name outside this list is kept, but it comes from a custom taxonomy, not from NCBI."""
NCBI_RANKS = (
	"no rank", "superkingdom", "domain", "kingdom", "subkingdom", "infrakingdom",
	"phylum", "superphylum", "subphylum", "infraphylum", "divisio", "subdivisio",
	"class", "superclass", "subclass", "infraclass", "parvclass",
	"superorder", "order", "suborder", "infraorder", "parvorder", "cohort", "subcohort",
	"superfamily", "family", "subfamily", "tribe", "subtribe",
	"genus", "subgenus", "section", "subsection", "series", "subseries", "species group",
	"species subgroup", "species", "subspecies", "variety", "varietas", "subvariety",
	"form", "forma", "forma specialis", "subform", "pathogroup", "pathotype", "biovar",
	"biotype", "serogroup", "serovar", "serotype", "strain", "isolate", "genotype", "morph",
	"clade",
)

"""rank names seen in the databases on disk that NCBI's taxdump of Oct 2023 held: kept here so the list
above can be checked against a run without reading a database"""
OBSERVED = (
	"no rank", "superkingdom", "genus", "species", "order", "family", "subspecies", "subfamily",
	"strain", "serogroup", "biotype", "tribe", "phylum", "class", "species group", "forma",
	"clade", "suborder", "subclass", "varietas", "kingdom", "subphylum", "forma specialis",
	"isolate", "superfamily", "infraorder", "infraclass", "superorder", "subgenus", "superclass",
	"parvorder", "serotype", "species subgroup", "subcohort", "cohort", "genotype", "subtribe",
	"section", "series", "morph", "subkingdom", "superphylum", "subsection", "pathogroup",
)


def unnamed(rank):
	'''True when the input named nothing at all: no rank column, an empty field, or the text None. The name
	   "no rank" is not unnamed: it is what NCBI writes for a taxid it does not place in a rank, and that is
	   data worth keeping (238,697 of those links in NCBI_base.db on disk).'''
	if rank is None:
		return True
	name = str(rank).strip()
	return name == "" or name == "None"


def is_no_rank(rank):
	'''True when the rank names no level of the tree: nothing named it, or it carries the "no rank" name
	   NCBI gives the taxids it does not place in a rank. Used by the code that asks whether a level of the
	   tree is a level at all (issue #140), not by the code that decides what to store.'''
	return unnamed(rank) or str(rank).strip() == NO_RANK


def is_defined(rank):
	'''True when the rank name is one of the ranks NCBI defines'''
	if is_no_rank(rank):
		return False
	return str(rank).strip() in NCBI_RANKS


def custom_ranks(rank_names):
	'''Rank names that are not in the NCBI vocabulary: they come from a custom taxonomy, and are kept'''
	return sorted({str(name).strip() for name in rank_names if name and not is_no_rank(name)
	               and not is_defined(name)})


def summary(database, log=False):
	'''What the rank table of a database holds: one row per rank name, how many links of the tree hang on
	   it, the links that carry no rank, and the names outside the NCBI vocabulary. Called at the end of an
	   import: which ranks a database ended up with is information the run should say out loud rather than
	   leave in a table that nothing reads afterwards.'''
	levels = database.select("SELECT rank.rank,COUNT(tree.child) FROM rank LEFT JOIN tree "
	                         "ON tree.rank_i = rank.rank_i GROUP BY rank.rank ORDER BY rank.rank")
	seen = {name: count for name, count in levels}
	unnamed = database.select("SELECT COUNT(*) FROM tree WHERE rank_i IS NULL")[0][0]
	dangling = database.select("SELECT COUNT(*) FROM tree WHERE rank_i IS NOT NULL AND rank_i NOT IN "
	                           "(SELECT rank_i FROM rank)")[0][0]
	custom = custom_ranks(seen)
	if log:
		log.info("Rank levels: {n} named, {unnamed} tree links without a rank, {custom} named outside the NCBI ranks".format(
			n=len(seen), unnamed=unnamed, custom=len(custom)))
		if custom:
			log.info("Custom rank levels: {levels}".format(levels=", ".join(custom)))
		if dangling:
			log.warning("{n} tree links point at a rank level that is not in the rank table".format(n=dangling))
	return {"levels": seen, "unnamed": unnamed, "dangling": dangling, "custom": custom}
