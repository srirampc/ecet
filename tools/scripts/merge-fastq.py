#!/usr/bin/python
import sys
import getopt
import re
#
#  File : merge-fastq.py
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

currentRead = 0
regExAmbig = '[^ACGTacgt]'

def mergeFile(inputFile,outf):
    global currentRead, regExAmbig
    with open(inputFile,'r') as inputf:
        while True:
            records  = [inputf.readline() for i in range(4)]
            if records[0] == '':
                return
            #print records[1],regExAmbig
            if re.search(regExAmbig, records[1].strip()):
                #print re.search(regExAmbig, records[1]).start()
                continue
            elts = records[0].split()
            outf.write('@' + str(currentRead))
            for elt in elts[1:]:
                outf.write(' ' + elt)
            outf.write('\n')
            for rcd in records[1:]:
                outf.write(rcd)
            currentRead += 1

def process(listFileName,outfile):
    paths = []
    with open(listFileName, 'r') as f:
        paths = f.readlines()
    outf = open(outfile, 'w')
    for filepath in paths:
        try:
            mergeFile(filepath.strip(),outf)
        except IOError as e:
            print str(e)
            print 'failed to process', filepath.strip()
    outf.close()

def usage():
    print """merge-fastq.py --list=/path/to/file-contatining-fast-file-paths
                 --outfile=/path/to/merged-output"""

def main(argv):
    try:
        opts,args = getopt.getopt(argv,"l:o:h",["list=","outfile=",
                                                "help"])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    listFileName = outfile = None
    for opt,arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-l", "--list"):
            listFileName= arg
	elif opt in ( "-o", "--outfile") :
	    outfile = arg
        else:
            print opt, "is not a valid option"
            usage()
            sys.exit()
    print "List Input File : ", listFileName
    print "Reads Out File : ", outfile
    if listFileName == None or outfile == None:
        usage()
        sys.exit()
    # Process the file
    process(listFileName,outfile)

if __name__ == "__main__":
    main(sys.argv[1:])

