#!/usr/bin/python
import sys
import getopt
import re
import itertools
import numpy as np
import alignment as an
import ecutils as eu
from mpi4py import MPI
#
#  File : par-hshrec-analy.py
#  Created on December 1, 2011
#  Author : Sriram PC <sriram.pc@gmail.com>
#
#  This file is part of Error Correction Review Toolkit.
#  Error Correction Review Toolkit is free software: you can 
#  redistribute it and/or modify  it under the terms of the GNU 
#  Lesser General Public License as published by 
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Error Correction Review Toolkit is distributed in the hope that 
#  it will be useful, but WITHOUT ANY WARRANTY; without even the implied
#  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with Libpnorm. If not, see <http://www.gnu.org/licenses/>.
#

psize = MPI.COMM_WORLD.Get_size()
prank = MPI.COMM_WORLD.Get_rank()
pname = MPI.Get_processor_name()

charMapping = {'A':'0','C':'1','G':'2','T':'3','N':'4','-':'5',
               'a':'0','c':'1','g':'2','t':'3','n':'4'}
alphabet = ['A','C','G','T','N','-', 'a','c','g','t','n']
revCompl = {}

def load_revcompl(alnfname):
    global revCompl
    with open(alnfname,'r') as f:
        for line in f:
            if line[0] == '@': # skip header line
                continue
            elts = line.split()
            readid = int(elts[0])
            flag = int(elts[1])
            if flag & 16 == 16:
                revCompl[readid] = True
            else:
                revCompl[readid] = False

def getRevCompl(readid):
    global revCompl
    try:
        if revCompl[readid]:
            return True
    except KeyError as e:
        return False


def do_align(readid,oread,cread):
    revCompl = getRevCompl(readid)
    if revCompl:
        oread = eu.reverse_complement(oread)
        cread = eu.reverse_complement(cread)
    (oalign,calign) = an.global_banded(oread,cread)
    if revCompl:
        oalign = eu.reverse_complement(oalign)
        calign = eu.reverse_complement(calign)
    #print (oalign,calign)
    return (oalign,calign)
   

def print_errors(outf,readid,alignments):
    # corrected_read followed by originial read
    pos = 0
    origaln = alignments[0]
    corraln = alignments[1]
    eline = []
    nerrors = 0
    for i in range(len(origaln)):
        if origaln[i] != corraln[i]:
            sbindl = 0
            if origaln[i] == '-':
                sbindl = 1
            elif corraln[i] == '-':
                sbindl = 2
            eline += [str(pos), charMapping[corraln[i]],
                      charMapping[origaln[i]], str(sbindl)]
            nerrors += 1
        if origaln[i] != '-':
            pos += 1
    eline = [str(readid), str(nerrors)] + eline
    outf.write( '\t'.join(eline) )
    outf.write( '\n' )

def print_trim(trimf,readid,read,tlen):
    readl = len(read)
    trimlen = readl - tlen
    trimf.write('\t'.join([str(readid),str(trimlen)]))
    trimf.write( '\n' )

def process(filename,correctedFile,outfile,trimfile,total):
    forig = open(filename, 'r')
    fcorr = open(correctedFile, 'r')
    trimf = open(trimfile + str(prank) + '.er', 'w')
    outf = open(outfile + str(prank) + '.er', 'w')
    #with open(filename, 'r') as forig, open(correctedFile, 'r') as fcorr:
    (startid, endid) = eu.decompose(total,prank,psize)
    try:
        args = [iter(forig)]*2
        ofiter = itertools.izip(*args)
        args = [iter(fcorr)]*2
        cfiter = itertools.izip(*args)
        try:
            # I should do this directly, right now
            # I will just read the file
            while True:
                corr = cfiter.next()
                cid = eu.getfaid(corr)
                if cid >= startid:
                    break
            while True:
                orig = ofiter.next()
                oid = eu.getfaid(orig)
                if oid >= startid:
                    break
            while True:
                cid = eu.getfaid(corr)
                oid = eu.getfaid(orig)
                if cid >= (endid+1):
                    break
                while cid != oid: 
                    #print oid, 'missing'
                    read = orig[1].strip()
                    # Discareded reads
                    print_trim(trimf, oid,read,len(read))
                    orig = ofiter.next()
                    oid = eu.getfaid(orig)
                assert(cid == oid)
                oread = orig[1].strip()
                cread = corr[1].strip()
                if oread != cread:
                    alignment = None
                    #if len(oread) == len(cread):
                    #    alignment = (oread,cread)
                    #else:
                    #print cid
                    alignment = do_align(cid,oread,cread)
                    print_errors(outf,cid,alignment)
                if cid >= endid:
                    break
                corr = cfiter.next()
                orig = ofiter.next()
        except StopIteration as stop:
            print str(stop)
    except IOError as e:
        print str(e)
    fcorr.close()
    forig.close()
    outf.close()
    trimf.close()
    MPI.COMM_WORLD.barrier()

def usage():
    print """Usage ::
      hshrec-analy.py --file=/path/to/original-fasta
                --corrected=/path/to/corrected-fasta
                --outfile=/path/to/err-output-prefix{prank}.er
                --trimfile=/path/to/trim-output-prefix{prank}.er
                --alignment=/path/to/alignment-sam-file
                --records=number or reads
                 (OR)
      hshrec-analy.py -f /path/to/original-fasta
                -c /path/to/corrected-fasta
                -o /path/to/err-output-prefix{prank}.er
                -t /path/to/trim-output-prefix{prank}.er
                -a /path/to/alignment-sam-file
                -r numbe of records
      The error will be ouput against the id given in FASTA file.
      FASTA file should have ">READ_ID" as header for each record.
      """

def main(argv):
    try:
        opts,args = getopt.getopt(argv,"f:c:o:t:a:r:h",["file=","corrected=",
                                                        "outfile=", "alignment=",
                                                        "records=","help"])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    # Get the command line arguments
    trimfile = filename = outfile = correctedFile = None
    alignFile = numRecords = None
    for opt,arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-f", "--file"):
            filename= arg
	elif opt in ( "-c", "--corrected") :
	    correctedFile = arg
	elif opt in ( "-o", "--outfile") :
	    outfile = arg
        elif opt in ("-t", "--trimfile"):
            trimfile = arg
        elif opt in ("-r", "--records"):
            numRecords = int(arg)
        elif opt in ("-a", "--alignment"):
            alignFile = arg
        else:
            print opt, "is not a valid option"
            usage()
            sys.exit()
    if prank == 0:
        print "UNCORRECTED Reads FASTA File : ", filename
        print "CORRECTED   Reads FASTA File : ", correctedFile
        print "TARGET ERROR OUTPUT File     : ", outfile
        print "TRIM OUTFILE File            : ", trimfile
        print "NUMBER OR RECORDS            : ", numRecords
        print "ALIGNMENT FILE               : ", alignFile
    if filename == None or outfile == None or correctedFile == None or trimfile == None:
        if prank == 0:
            print 'One of the input options not given'
            usage()
        sys.exit()
    # Process the file
    MPI.COMM_WORLD.barrier()
    dtime = MPI.Wtime()
    try:
        load_revcompl(alignFile)
        process(filename,correctedFile,outfile,trimfile,numRecords)
    except IOError as e:
        print str(e)
        print 'Error Opening files'
    MPI.COMM_WORLD.barrier()
    totaltime = MPI.Wtime() - dtime
    if (prank == 0):
        print 'Total Time ' + str(totaltime)
if __name__ == "__main__":
    #decompose(200)
    main(sys.argv[1:])

