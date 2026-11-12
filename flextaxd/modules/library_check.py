#!/usr/bin/env python3

'''
The library of a build, counted before the external program reads it.

kraken2, krakenuniq and ganon take the taxid of a sequence either from the library header (the
`|kraken:taxid|<taxid>` tag of kraken2) or from a line of the taxid map keyed by the id the library
holds, and a sequence that reaches neither side is simply left out: the acceptance build of T3 writes
9 sequences into `library.fna` and kraken2 indexes 8 of them, and nothing in the run says that one
sequence was dropped. Counting the assembled library and its map before the program is called shows a
build that drops sequences at the point where its input can still be fixed.

Four numbers, all counted by streaming the files, never by holding their sequences in memory:

  sequences                       sequence records in the library file(s) of the build
  map lines                       lines of the taxid map the program reads
  without a taxid                 library records the program can give no taxid to
  map lines without a genome      map lines naming an annotated genome this run had no sequence
                                  file of its own

The last number is deliberately not "map keys that are no header in the library": a genome file carries
its sequences under their own accessions, so a map line naming the assembly of a genome whose file WAS
read matches no library header and says nothing wrong. What the annotation of the database can tell is
whether a map line names an annotated genome the run never had a file of - and issue #131 is what makes
that set mean what it says, by keeping the sequence level annotations of a packed file out of it.

The four numbers are stored in the meta table of the flextaxd database, the table that already holds
the version and the date of the build, and `flextaxd --report` prints them as one line.

The counting is written for the shapes of the taxid maps of the programs (the id of a kraken2 prelim
map line sits in its second field, that of a krakenuniq map and of a seq-info map in the first), but
only the kraken route calls it today: which program databases the build supports is being decided
separately (issue #150), and a route gets its call site when that is settled.
'''

__version__ = "1.0"
__author__ = "David Sundell"
__credits__ = ["David Sundell"]
__license__ = "GPLv3"

import datetime
import logging
import os

from .functions import zopen

logger = logging.getLogger(__name__)

'''the tag kraken2 puts between the sequence id and its taxid, inside the first field of a header'''
KRAKEN_TAG = "|kraken:taxid|"

META_KEYS = {"sequences": "library_sequences",
			 "map_lines": "library_map_lines",
			 "without_taxid": "library_without_taxid",
			 "map_without_genome": "library_map_without_genome"}


def sequence_key(header):
	'''The id a taxid map is keyed by, for one fasta header: its first field without the ">", and
			whether the header carries the kraken taxid tag, which is a taxid of its own'''
	key = header.strip().split(" ")[0].lstrip(">").strip()
	if KRAKEN_TAG in key:
		return key.split(KRAKEN_TAG)[0], True
	return key, False


def map_key(line, key_column):
	'''The sequence id one line of a taxid map names, in whichever column the program keeps it: the
			second field of a `TAXID\t<id>\t<taxid>` line of the prelim map of kraken2, the first field of
			the seq-info map of ganon or of the map of krakenuniq'''
	columns = line.rstrip("\n").split("\t")
	if len(columns) <= key_column:
		return None
	key = columns[key_column].strip()
	if KRAKEN_TAG in key:
		key = key.split(KRAKEN_TAG)[0]
	return key or None


def library_numbers(library_files, map_file, annotated_genomes=(), genome_files=(), sequence_level_annotations=(), map_key_column=0):
	'''Count what the assembled library of a build holds, on the files the external program is about to
			read. Every number is counted by streaming one file; no sequence is kept in memory.

	------
	Parameters:
		list - the library files of the build, whichever names they ended up with
		str - the taxid map the program reads
		iterable - the genomes annotated in the database (the genome names of its genomes table)
		iterable - the genomes this run had a sequence file for
		iterable - the annotations made on sequence level, whose sequence sits inside a packed file and
			which therefore has no file of its own (DatabaseFunctions.get_sequence_annotations, #131)
		int - the column of the map holding the sequence id
	------
	Returns:
		dict - sequences, map_lines, without_taxid and map_without_genome
	'''
	absent = set(annotated_genomes) - set(genome_files) - set(sequence_level_annotations)
	numbers = {"sequences": 0, "map_lines": 0, "without_taxid": 0, "map_without_genome": 0}
	map_keys = set()
	if os.path.exists(map_file):
		with zopen(map_file, "r") as handle:
			for line in handle:
				if not line.strip():
					continue
				numbers["map_lines"] += 1
				key = map_key(line, map_key_column)
				if key is None:
					continue
				map_keys.add(key)
				if key in absent:
					numbers["map_without_genome"] += 1
	for path in library_files:
		if not os.path.exists(path):
			continue
		with zopen(path, "r") as handle:
			for line in handle:
				if not line.startswith(">"):
					continue
				numbers["sequences"] += 1
				key, tagged = sequence_key(line)
				if not tagged and key not in map_keys:
					numbers["without_taxid"] += 1
	return numbers


def store_library_numbers(database, program, library_dir, numbers):
	'''Store the four numbers in the meta table of the flextaxd database, so that `flextaxd --report`
			says what the library of the last build held, and say them where the build runs: the two
			numbers that mean the program leaves sequences out are logged as warnings.

	------
	Parameters:
		DatabaseFunctions - the flextaxd database the build was made from
		str - the program the library was built for
		str - the directory of the program database
		dict - the numbers of library_numbers
	------
	Returns:
		dict - the numbers, so a caller can use them
	'''
	database.set_meta("library_date", datetime.datetime.now().isoformat(timespec="seconds"))
	database.set_meta("library_program", str(program))
	database.set_meta("library_database", str(library_dir))
	for number, key in META_KEYS.items():
		database.set_meta(key, str(int(numbers.get(number, 0) or 0)))
	logger.info("Library of the build in {d}: {s} sequences, {m} lines in the taxid map".format(
		d=library_dir, s=numbers["sequences"], m=numbers["map_lines"]))
	if numbers["without_taxid"]:
		logger.warning("{n} of the {s} sequences of the library have no taxid in their header and no line "
					   "in the map: {program} will leave them out of the database".format(
			n=numbers["without_taxid"], s=numbers["sequences"], program=program))
	if numbers["map_without_genome"]:
		logger.warning("{n} lines of the taxid map name an annotated genome that has no sequence file; the "
					   "disk report of --genomes_path lists those genomes".format(n=numbers["map_without_genome"]))
	return numbers


def library_report(database):
	'''The numbers the last build of this database stored, as one line for the report, or None when no
			library of a build was counted'''
	sequences = database.get_meta(META_KEYS["sequences"])
	if sequences is None:
		return None
	return "Library of the last build ({p}, {d}): {s} sequences, {m} map lines, {t} without a taxid, {g} map lines for a genome without sequence file".format(
		p=database.get_meta("library_program") or "unknown program",
		d=database.get_meta("library_date") or "unknown date",
		s=sequences,
		m=database.get_meta(META_KEYS["map_lines"]) or "0",
		t=database.get_meta(META_KEYS["without_taxid"]) or "0",
		g=database.get_meta(META_KEYS["map_without_genome"]) or "0")
