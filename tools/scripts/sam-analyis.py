#!/usr/bin/python
import sys
import getopt
import re
import alignment as an
import ecutils as eu
#
#  File : sam-analyis.py
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
               'a':'0','c':'1','g':'2','t':'3','n':'4',
               'Y':'6','R':'7',
               'y':'6','r':'7'}
alphabet = eu.alphabet
error_types = [0,1,2]
complementMapping = eu.complementMapping
fgenome = { 0: None }

class SkipError(Exception):
     def __init__(self, value):
         self.value = value
     def __str__(self):
         return repr(self.value)

# For the alignment build the error as a dictionary
#  { 
#    pos : [ [true_base,wrong_base,type_of_error],...]
#    ...
#  }
def getErrorsFromAlignment(readString,refAlign,readAlign,flag):
    l = len(refAlign)
    readl = len(readString)
    sofar = 0
    errors = {}
    revCompl = False
    if flag & 16 == 16:
        revCompl = True
    for i in range(readl):
        errors[i] = []
    if revCompl:
       refAlign = eu.reverse_complement(refAlign)
       readAlign = eu.reverse_complement(readAlign)
    #print refAlign
    #print readAlign
    for i in range(l):
        tb = wb = '-'
        errtype = -1
        if readAlign[i] == '-' and refAlign[i] in alphabet:
            errtype = 1  # deletion error (missing char in the read after pos)
        elif refAlign[i] == '-' and readAlign[i] in alphabet:
            errtype = 2 # insertion error (extra char in read at pos)
        elif readAlign[i] != refAlign[i] and readAlign[i] in alphabet  and refAlign[i] in alphabet:
            errtype = 0 # substition error
        if errtype in error_types:
            tb = charMapping[refAlign[i]] # true base
            wb = charMapping[readAlign[i]]  # wrong base
            errors[sofar] += [[tb,wb,errtype]] # error type
        if readAlign[i] != '-':
            sofar += 1
    #print errors
    return errors

def getErrors(readid,readString,cigarString,mdString,gpos,flag):
    errors = { 'id':readid }
    if mdString != None:
        mdString = mdString.split(':')[-1]
    #print readid,readString,cigarString, mdString
    [refAlign,readAlign,sb,se] = an.getSAMAlignment(readString,cigarString,mdString,
                                                    fgenome[0],gpos)
    errors = getErrorsFromAlignment(readString,refAlign,readAlign,flag)
    return [errors,sb,se]

def countErrors(errors):
     nerrs = 0
     for pos in sorted(errors.keys()):
          nerrs += len(errors[pos])
     return nerrs

def writeErrors(outf,readid,errors):
    errorStr = ''
    count = 0
    for pos in sorted(errors.keys()):
        errlst = errors[pos]
        for l in errlst:
            errorStr += '\t' + str(pos)
            count += 1
            for x in l:
                errorStr += '\t' + str(x)
    if len(errorStr) > 0:
        outf.write( readid + '\t' + str(count) + errorStr + '\n')

def processLine(line,outf,unmapf,ambf,trimf,dryRun):
    sminf = eu.getSAMinfo(line)
    nerrs = unmapped = ambig = skipLength = 0
    unmapLength = ambigLength = 0
    readid = sminf['QNAME']
    flag = sminf['FLAG']
    mapqual = sminf['MAPQ']
    cigarString = sminf['CIGAR']
    gpos = sminf['POS']
    asString = sminf['AS']
    xtString = sminf['XT']
    readString = sminf['SEQ']; readLen = len(readString)
    if (flag & 4 == 4) : # 4 === aligned at multiple places
         unmapped += 1; unmapLength += readLen
         if unmapf: unmapf.write(readid+'\n')
    elif( mapqual == 0 ): # map quality is zero for unambiguous read
         ambig += 1; ambigLength += readLen
         if ambf: ambf.write(readid + '\n')
    elif( cigarString != '*' and  mapqual != 0): #if read is uniqly aligned
        mdstring = sminf['MD']
        if mdstring == None and len(fgenome) == 0:
             print 'No MD String / Genome available !!!'
             assert mdstring != None
        try:
            [errors,sb,se] = getErrors(readid,readString,cigarString,mdstring,gpos,flag)
            if errors:
                 nerrs += countErrors(errors)
                 if outf:
                      writeErrors(outf,readid,errors)
            if sb + se != readLen:
                if trimf: 
                     trimf.write(readid + '\t' + str(se) + '\t' + str(sb) + '\n')
                skipLength += sb + readLen - se
        except SkipError as e:
            print 'skip ' + readid
            if ambf : ambf.write(readid + '\n')
            ambig += 1; ambigLength += readLen
    else : # this is the case should not happen
        print str(readid) + ' is ambigous '
    return (nerrs, skipLength, unmapped,unmapLength,ambig,ambigLength)


def process(fileName,outfile,unmappedFile,ambigFile,genomeFile,trimFile,dryRun):
    #print "Processing",fileName
    global fgenome
    if genomeFile != None:
        fgenome = eu.load_genome(genomeFile)
    outf = unmapf = ambf = trimf = None
    if not dryRun:
         outf = open(outfile, 'w')
         unmapf = open(unmappedFile, 'w')
         trimf = open(trimFile,'w')
         ambf = unmapf
         if unmappedFile != ambigFile:
              ambf = open(ambigFile, 'w')
    terrs = tSkipLength = tunmapped = tUnmapLength = tambig = tAmbigLen = 0
    with open(fileName, 'r') as f:
        # eliminate first line
        for line in f:
            if line[0] == '@': # skip header line
                continue
            try:
                 (nerrs, skipLength, unmapped,unmapLength,ambig,ambigLength) =  processLine(line,outf,unmapf,ambf,trimf,dryRun)
                 terrs += nerrs; tSkipLength += skipLength; tunmapped += unmapped
                 tUnmapLength += unmapLength; tambig += ambig; tAmbigLen += ambigLength
            except AssertionError, b:
                print b
                print line
    if not dryRun:
         outf.close()
         unmapf.close()
         trimf.close()
         if unmappedFile != ambigFile:
              ambf.close()
    print 'Total number of errors ', terrs
    print 'Total length skipped ', tSkipLength
    print 'Total unmapped Reads ', tunmapped
    print 'Total umap length ', tUnmapLength
    print 'Total Amig ', tambig
    print 'Total Amig Length ', tAmbigLen
    print 'Total ambig + unmap ', tambig + tunmapped
    print 'Total amg + unmap length ', tAmbigLen + tUnmapLength

def verify(fileName,outfile,unmappedFile,ambigFile,genomeFile):
    # TODO:
    # Load Genome
    # For each read in SAM file
    pass

def usage():
    print """sam-analysis.py --file=/path/to/sam-file-input
                --outfile=/path/to/err-output
                --ambig=/path/to/ambig-output
                --unmapped=/path/to/unmapped-output
                --trim=/path/to/trim-file-output
                [--genomeFile=/path/to/genome-file]
                [--dry (for dry run no output generated)]
     trim-file-output positions trimmed (actually ranges allowed)
     Unmapped and ambigous file can be both same.
     Genome file is optional. Used if MD String is not available.
     If Genome file is given, it will be loaded in memory
     """

def main(argv):
    try:
        opts,args = getopt.getopt(argv,"f:o:u:a:t:g:dvh",["file=","outfile=",
                                                          "unmapped=", "ambig=",
                                                          "trim=","genome=",
                                                          "dry" "verify","help"])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    # Get the command line arguments
    trimFile = ambigFile = filename = outfile = unmappedFile = None
    genomeFile = None
    dryRun = verify = False
    for opt,arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-f", "--file"):
            filename= arg
	elif opt in ( "-o", "--outfile") :
	    outfile = arg
	elif opt in ( "-u", "--unmapped") :
	    unmappedFile = arg
        elif opt in ("-a", "--ambig"):
            ambigFile = arg
        elif opt in ("-t", "--trim"):
            trimFile = arg
        elif opt in ("-d", "--dry"):
             dryRun = True
        elif opt in ("-v", "--verify"):
            verify = True
        elif opt in ("-g", "--genome"):
            genomeFile = arg
        else:
            print opt, "is not a valid option"
            usage()
            sys.exit()
    print "SAM Input File : ", filename
    print "Mapped Reads File : ", outfile
    print "Unmapped Reads File : ", unmappedFile
    print "Ambigous Reads File ", ambigFile
    print 'Ref Genome File ', genomeFile
    print 'Trim File ', trimFile
    print 'Dry Run ', dryRun
    if filename == None or outfile == None or unmappedFile == None or ambigFile == None or trimFile == None:
        usage()
        sys.exit()
    if verify:
        if genomeFile == None:
            usage()
            sys.exit()
        verify(filename,outfile,unmappedFile,ambigFile,genomeFile)
    else:
        # Process the file
        process(filename,outfile,unmappedFile,ambigFile,genomeFile,trimFile,dryRun)

if __name__ == "__main__":
    main(sys.argv[1:])

