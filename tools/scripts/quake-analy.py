#!/usr/bin/python
import sys
import getopt
import re
import itertools
#
#  File : quake-analy.py
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

def hd(orig,corrected,trim):
    assert len(orig) == len(corrected) + trim
    num_mismatch = 0
    for i in range(len(corrected)):
        if orig[i] != corrected[i]:
            num_mismatch += 1
    return num_mismatch

def print_error(outf,readid,orig,corrected,hd,trim):
    eline = [str(readid),str(hd+trim)]
    for i in range(len(corrected)):
        if orig[i] != corrected[i]:
            eline += [str(i), charMapping[corrected[i]], 
                      charMapping[orig[i]], str(0)]
    # Dont print the trim as eletion
    # for i in range(trim):
    #     pos = i + len(corrected)
    #     eline += [str(pos), charMapping['-'],
    #               charMapping[orig[pos]], str(2) ]
    outf.write( '\t'.join(eline) )
    outf.write( '\n' )

def print_trim(trimf,readid,read,tlen):
    readl = len(read)
    trimlen = readl - tlen
    trimf.write('\t'.join([str(readid),str(trimlen)]))
    trimf.write( '\n' )

def gettrim(lines):
    trim_count = 0
    elts = lines[0].split()
    trim_str = (elts[-1])
    reobj = re.match('trim=(\d+)',trim_str)
    if reobj:
        trim_count = int(reobj.group(1))
    return trim_count
    
def getid(lines):
    elts = lines[0].split()
    readid = elts[0][1:]
    return readid

def checkid(target):
    return lambda x: getid(x) != target

def process(filename,correctedFile,outfile,trimfile):
    trimf = open(trimfile, 'w')
    outf = open(outfile, 'w')
    forig = open(filename, 'r')
    fcorr = open(correctedFile, 'r')
    #with open(filename, 'r') as forig, open(correctedFile, 'r') as fcorr:
    try:
        args = [iter(forig)]*4
        ofiter = itertools.izip(*args)
        args = [iter(fcorr)]*4
        cfiter = itertools.izip(*args)
        # This is not a merge operation. We assume that orig
        # FASTQ file has all the ids in proper order. The other
        # corrected file may lag.
        try:
            while True:
                corr = cfiter.next()
                orig = ofiter.next()
                cid = getid(corr)
                oid = getid(orig)
                while cid != oid: 
                     print oid, 'missing'
                     read = orig[1].strip()
                     print_trim(trimf, oid,read,len(read))
                     orig = ofiter.next()
                     oid = getid(orig)
                assert(cid == oid)
                # Compute hamming distance and trim
                t = gettrim(corr)
                oread = orig[1].strip()
                cread = corr[1].strip()
                h = hd(oread,cread,t)
                if h > 0 :
                    # print readid,orig,corrected
                    print_error(outf,cid,oread,cread,h,t)
                if t > 0:
                   print_trim(trimf,cid,oread,t)
        except StopIteration as stop:
            print str(stop)
        except IOError as e:
            print str(e)
    except IOError as e:
        print str(e)
    fcorr.close()
    forig.close()
    outf.close()
    trimf.close()

def usage():
    print """Usage ::
      quake-analy.py --file=/path/to/original-fastq
                --corrected=/path/to/corrected-fastq
                --outfile=/path/to/err-output
                --trimfile=/path/to/trim-output
                 (OR)
      quake-analy.py -f /path/to/original-fastq
                -c /path/to/corrected-fastq
                -o /path/to/err-output
                -t /path/to/trim-output

      The error will be ouput against the id given in FASTQ file.
      If the id needs to be proper (i.e., it can be compared against
      the  BWA alignments run) , then the fastq file should be
      pre-processed to update ids using merge-fastq.py - before
      running the error correction. 
      Also, the scripts handles ONLY Substitution errors :(.
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
    print "UNCORRECTED Reads FASTQ File : ", filename
    print "CORRECTED   Reads FASTQ File : ", correctedFile
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

