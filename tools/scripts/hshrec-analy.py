#!/usr/bin/python
import sys
import getopt
import re
import itertools
import alignment as an
#
#  File : hshrec-analy.py
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

charMapping = {'A':'0','C':'1','G':'2','T':'3','N':'4','-':'5',
               'a':'0','c':'1','g':'2','t':'3','n':'4'}
alphabet = ['A','C','G','T','N','-', 'a','c','g','t','n']

def getid(lines):
    elts = lines[0].split()
    readid = int(elts[0][1:])
    return readid

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

def process(filename,correctedFile,outfile,trimfile):
    trimf = open(trimfile, 'w')
    outf = open(outfile, 'w')
    forig = open(filename, 'r')
    fcorr = open(correctedFile, 'r')
    #with open(filename, 'r') as forig, open(correctedFile, 'r') as fcorr:
    try:
        args = [iter(forig)]*2
        ofiter = itertools.izip(*args)
        args = [iter(fcorr)]*2
        cfiter = itertools.izip(*args)
        try:
            while True:
                corr = cfiter.next()
                orig = ofiter.next()
                cid = getid(corr)
                oid = getid(orig)
                while cid != oid: 
                    #print oid, 'missing'
                    read = orig[1].strip()
                    # Discareded reads
                    print_trim(trimf, oid,read,len(read))
                    orig = ofiter.next()
                    oid = getid(orig)
                assert(cid == oid)
                oread = orig[1].strip()
                cread = corr[1].strip()
                if oread != cread:
                    alignment = an.global_align(oread,cread)
                    print_errors(outf,cid,alignment)
        except StopIteration as stop:
            print str(stop)
    except IOError as e:
        print str(e)
    fcorr.close()
    forig.close()
    outf.close()
    trimf.close()

def usage():
    print """Usage ::
      hshrec-analy.py --file=/path/to/original-fasta
                --corrected=/path/to/corrected-fasta
                --outfile=/path/to/err-output
                --trimfile=/path/to/trim-output
                 (OR)
      hshrec-analy.py -f /path/to/original-fasta
                -c /path/to/corrected-fasta
                -o /path/to/err-output
                -t /path/to/trim-output

      The error will be ouput against the id given in FASTA file.
      FASTA file should have ">READ_ID" as header for each record.
      """

def main(argv):
    try:
        opts,args = getopt.getopt(argv,"f:c:o:t:h",["file=","corrected=",
                                                    "outfile=","help"])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    # Get the command line arguments
    trimfile = filename = outfile = correctedFile = None
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
        else:
            print opt, "is not a valid option"
            usage()
            sys.exit()
    print "UNCORRECTED Reads FASTA File : ", filename
    print "CORRECTED   Reads FASTA File : ", correctedFile
    print "TARGET ERROR OUTPUT File     : ", outfile
    print "TRIM OUTFILE File            : ", trimfile
    if filename == None or outfile == None or correctedFile == None or trimfile == None:
        print 'One of the input options not given'
        usage()
        sys.exit()
    # Process the file
    try:
        process(filename,correctedFile,outfile,trimfile)
    except IOError as e:
        print str(e)
        print 'Error Opening files'

if __name__ == "__main__":
    main(sys.argv[1:])

