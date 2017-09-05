__author__ = 'Bruno Fosso'
__version__ = "1.0"

import os
from string import strip
import shlex
import subprocess
from Bio import SeqIO
import gzip
import getopt
import sys


def usage():
    print """
    Species representative identification\n
    option:\n
    \t-i\tcomplete path to fasta file [MANDATORY]\n
    \t-t\tcomplete path to NCBI taxonomy DUMP files containing folder [MANDATORY]\n
    \t-l\tcomplete path to ITS1 location csv file [MANDATORY]\n
    \t-h\tprint this help page\n
    usage:\n
    \tpython find_representative.py -i fungal_references.fa\n
    """


its1_location = ""
fasta_input = ""
taxnomy_folder = ""
try:
    opts, args = getopt.getopt( sys.argv[1:], "hi:" )
except getopt.GetoptError, err:
    print str( err )
    usage()
    sys.exit()
if len( opts ) >= 1:
    for o, a in opts:
        if o == "-h":
            usage()
            sys.exit()
        elif o == "-i":
            fasta_input = a
        elif o == "-t":
            taxnomy_folder = a
        elif o == "-l":
            its1_location = a
        else:
            print "Unhandled option."
            usage()
            sys.exit()
else:
    usage()
    sys.exit()

wd = os.getcwd()
if taxnomy_folder != "":
    if os.path.exists( taxnomy_folder ) is False:
        print "The indicated path to NCBI taxonomy folder is wrong"
        sys.exit()
else:
    print "The -t option in MANDATORY"

if its1_location != "":
    if os.path.exists( its1_location ) is False:
        print "The indicated path to ITS1 location CSV file is wrong"
        sys.exit()
else:
    print "The -l option is MANDATORY"
    sys.exit()

NODESFILE = ""
namesfile = ""
merged = ""
deleted = ""
if os.path.exists( os.path.join( taxnomy_folder, "nodes.dmp" ) ):
    NODESFILE = open( os.path.join( taxnomy_folder, "nodes.dmp" ) )
else:
    print "No NODESFILE"
    exit()

if os.path.exists(
        os.path.join( taxnomy_folder, "names.dmp" ) ):
    namesfile = open(
        os.path.join( taxnomy_folder, "names.dmp" ) )
else:
    print "no namesfile"
    exit()

if os.path.exists( os.path.join( taxnomy_folder, "merged.dmp" ) ):
    merged = open( os.path.join( taxnomy_folder, "merged.dmp" ) )
else:
    print "no merged"
    exit()

if os.path.exists( os.path.join( taxnomy_folder, "delnodes.dmp" ) ):
    deleted = open( os.path.join( taxnomy_folder, "delnodes.dmp" ) )
else:
    print "no deleted"
    exit()

node2name = {}
name2node = {}
for line in namesfile:
    line = line.strip()
    fields = map( strip, line.split( "|" ) )
    if fields[3] == "scientific name":
        nodeid, name = fields[0], fields[1]
        node2name[nodeid] = name
        name2node[name] = nodeid

node2parent = {}
node2order = {}
for line in NODESFILE:
    line = line.strip()
    fields = map( strip, line.split( "|" ) )
    nodeid, parentid, ordine = fields[0], fields[1], fields[2]
    node2parent[nodeid] = parentid
    node2order[nodeid] = ordine

node2merged = {}
for line in merged:
    line = line.strip()
    fields = map( strip, line.split( "|" ) )
    node, new_id = fields[0], fields[1]
    node2merged[node] = new_id

delnodes = []
for line in deleted:
    line = line.strip()
    fields = map( strip, line.split( "|" ) )
    delnodes.append( fields[0] )

acc2_HMM_loc = {}
acc2_GB_loc = {}
acc2taxa = {}
lines = open( its1_location ).readlines()
for line in lines[1:]:
    line = line.strip()
    s = line.split( "\t" )
    # print s
    acc, taxon_id = s[0].split( "." )[0], s[1]
    acc2taxa[acc] = taxon_id
    if "false" != s[2]:
        start = str( int( s[3] ) - 1 )
        end = str( int( s[4] ) )
        acc2_GB_loc.setdefault( acc, [] )
        acc2_GB_loc[acc].append( start )
        acc2_GB_loc[acc].append( end )
    if "false" != s[5]:
        start = int( s[6] ) - 1
        end = int( s[7].strip() )
        if end < start:
            pass
        else:
            acc2_HMM_loc.setdefault( acc, [] )
            acc2_HMM_loc[acc].append( str( start ) )
            acc2_HMM_loc[acc].append( str( end ) )

print len( acc2taxa )
print len( acc2_GB_loc )
print len( acc2_HMM_loc )

for acc in acc2taxa.keys():
    node = acc2taxa[acc]
    if node2parent.has_key( node ):
        pass
    elif node2merged.has_key( node ):
        acc2taxa[acc] = node2merged[node]
        # node2order[acc] = "GB acc"
    elif node in delnodes:
        print acc, node, "deleted"
    else:
        print acc, node, "nothing"

for acc in acc2taxa.keys():
    node = acc2taxa[acc]
    if node2order[node] == "species":
        pass
    else:
        parent = node2parent[node]
        while node != parent:
            if node2order[parent] == "species":
                acc2taxa[acc] = parent
                node = parent
            else:
                node = parent
                parent = node2parent[node]

if os.path.exists( "Representative_computation" ) is False:
    os.mkdir( "Representative_computation" )

if os.path.exists( os.path.join( "Representative_computation", "seq_to_species" ) ) is False:
    os.mkdir( os.path.join( "Representative_computation", "seq_to_species" ) )


def accession(identifier):
    parts = identifier.split( "." )
    return parts[0]


seq_index = SeqIO.index( fasta_input, "fasta", key_function=accession )

species2acc = {}
for acc in acc2taxa.keys():
    # print acc
    specie = acc2taxa[acc]
    species2acc.setdefault( specie, set() )
    species2acc[specie].add( acc )

print
print
scritte = 0
for species in species2acc.keys():
    fasta = open( os.path.join( "Representative_computation", "seq_to_species", species + ".fasta" ), "w" )
    for acc in species2acc[species]:
        if acc2_GB_loc.has_key( acc ):
            print >> fasta, ">" + acc + "_GB"
            print >> fasta, seq_index[acc].seq[int( acc2_GB_loc[acc][0] ):int( acc2_GB_loc[acc][1] )]
            scritte += 1
        if acc2_HMM_loc.has_key( acc ):
            print >> fasta, ">" + acc + "_HMM"
            print >> fasta, seq_index[acc].seq[int( acc2_HMM_loc[acc][0] ):int( acc2_HMM_loc[acc][1] )]
            scritte += 1
    fasta.close()
    if os.stat( os.path.join( "Representative_computation", "seq_to_species", species + ".fasta" ) )[6] == 0:
        os.remove( os.path.join( "Representative_computation", "seq_to_species", species + ".fasta" ) )
        print "Representative_computation/seq_to_species/" + species + ".fasta is empty"

print "Sono state processate %i sequenze" % scritte
rep_file = open( "rep_sequences.csv", "w" )
manual_control = set()
for fasta in os.listdir( os.path.join( "Representative_computation", "seq_to_species" ) ):
    if fasta.endswith( "fasta" ):
        if len( list( SeqIO.parse( os.path.join( "Representative_computation", "seq_to_species", fasta ), "fasta" ) ) ) == 1:
            record = SeqIO.read( "Representative_computation/seq_to_species/" + fasta, "fasta" )
            print >> rep_file, fasta.split( "." )[0] + "\t" + record.name
        else:
            cluster_file = os.path.join( "Representative_computation", "seq_to_species",
                                         fasta.rstrip( ".fasta" ) + ".usearch_clustering" )
            cmd = shlex.split( "vsearch --cluster_fast %s -id 0.97 -uc  %s" % (
                os.path.join( "Representative_computation", "seq_to_species", fasta ), cluster_file) )
            p = subprocess.Popen( cmd )
            p.wait()
            with open( cluster_file ) as a:
                cluster_data = {}
                seq2len = {}
                for line in a:
                    field = map( strip, line.split( "\t" ) )
                    if field[0] == "S":
                        seq2len[field[-2]] = int( field[2] )
                        cluster_data.setdefault( field[-2], set() )
                        cluster_data[field[-2]].add( field[-2] )
                    elif field[0] == "H":
                        seq2len[field[-2]] = int( field[2] )
                        cluster_data.setdefault( field[-1], set() )
                        cluster_data[field[-1]].add( field[-2] )
            from random import choice

            if len( cluster_data ) == 0:
                manual_control.add( fasta )
            else:
                size = []
                for cluster in cluster_data.values():
                    size.append( len( cluster ) )
                greatest = max( size )
                if size.count( greatest ) == 1:
                    len_max = []
                    for cluster in cluster_data.keys():
                        if len( cluster_data[cluster] ) == greatest:
                            print >> rep_file, fasta.split( "." )[0] + "\t" + cluster
                else:
                    candidate = []
                    for cluster in cluster_data.keys():
                        if len( cluster_data[cluster] ) == greatest:
                            candidate.append( cluster )
                    rep = choice( candidate )
                    print >> rep_file, fasta.split( "." )[0] + "\t" + rep
rep_file.close()
print "\n".join( list( manual_control ) )
