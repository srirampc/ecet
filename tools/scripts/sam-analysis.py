#!/usr/bin/python
import sys
import argparse
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
    nerrs = nunmap = nambig = nmapd = 0
    skipl = unmapLength = ambigLength = mapLength = 0
    readid = sminf['QNAME']
    flag = sminf['FLAG']
    mapqual = sminf['MAPQ']
    cigarString = sminf['CIGAR']
    gpos = sminf['POS']
    #asString = sminf['AS']
    #xtString = sminf['XT']
    readString = sminf['SEQ']
    readLen = len(readString)
    if flag & 4 == 4: # 4 === aligned at multiple places
         nunmap += 1; unmapLength += readLen
         if unmapf: unmapf.write(readid+'\n')
    elif mapqual == 0: # map quality is zero for unambiguous read
         nambig += 1; ambigLength += readLen
         if ambf: ambf.write(readid + '\n')
    elif (cigarString != '*') and (mapqual != 0): #if read is uniqly aligned
        mdstring = sminf['MD']
        if mdstring == None and len(fgenome) == 0:
            eu.eprint("NO MD STRING / GENOME AVAILABLE !!!")
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
                skipl += sb + readLen - se
                mapLength = se - sb
                nmapd += 1
            else:
                mapLength = readLen
                nmapd += 1
        except SkipError as err:
            eu.eprint('SKIP ' + readid + ' ' + str(err))
            if ambf: ambf.write(readid + '\n')
            nambig += 1; ambigLength += readLen
    else : # this is the case should not happen
        eu.eprint(str(readid) + ' is ambigous ')
    return (nerrs, skipl, nunmap, unmapLength, nambig, ambigLength, nmapd, mapLength)


def process(fileName,outFile,unmappedFile,ambigFile,genomeFile,trimFile,dryRun):
    #print "Processing",fileName
    global fgenome
    if dryRun:
        outFile = ambigFile = unmappedFile = ambigFile = trimFile = None
    if genomeFile is not None: fgenome = eu.load_genome(genomeFile)
    inf = outf = unmapf = ambf = trimf = None
    terrs = tSkipLength = tMappedLength = 0
    tmapd = tunmapped = tUnmapLength = tambig = tAmbigLen = 0
    try:
        if unmappedFile is not None: unmapf = open(unmappedFile, 'w')
        if trimFile is not None: trimf = open(trimFile,'w')
        if (unmappedFile != ambigFile) and (ambigFile is not None):
            ambf = open(ambigFile, 'w')
        if ambf is None: ambf = unmapf
        if outFile is not None:
            if outFile == "-":
                outf = sys.stdout
            else:
                outf = open(outFile, 'w')
        if fileName == "-":
            inf = sys.stdin
        else:
            inf = open(fileName, 'r') 
        # eliminate first line
        for line in inf:
            if line[0] == '@': # skip header line
                continue
            try:
                (nerr,skipl,nunmap,unmapl,
                 nambig,ambigl,nmap,maplen) = processLine(line,outf,unmapf,
                                                          ambf,trimf,dryRun)
                terrs += nerr; tSkipLength += skipl
                tunmapped += nunmap; tUnmapLength += unmapl
                tambig += nambig; tAmbigLen += ambigl
                tmapd += nmap; tMappedLength += maplen
            except AssertionError, b:
                eu.eprint(b)
                eu.eprint(line)
    except IOError, ioe:
        eu.eprint("I/O Err " + str(ioe))
    if outf is not None: outf.close()
    if inf is not None: inf.close()
    if unmapf is not None: unmapf.close()
    if trimf is not None: trimf.close()
    if unmappedFile != ambigFile and ambf is not None: ambf.close()
    eu.eprint("------------------ SAM STATISTICS ------------------")
    eu.eprint('No. of errors               :{0:>11}'.format(terrs))
    eu.eprint('No. unmapped Reads          :{0:>11}'.format(tunmapped))
    eu.eprint('No. mapped Reads            :{0:>11}'.format(tmapd))
    eu.eprint('No. ambig Reads             :{0:>11}'.format(tambig))
    eu.eprint('No. ambig + unmap           :{0:>11}'.format(tambig + tunmapped))
    eu.eprint('Total mapped length         :{0:>11}'.format(tMappedLength))
    eu.eprint('Total unmapped length       :{0:>11}'.format(tUnmapLength))
    eu.eprint('Total clipped length        :{0:>11}'.format(tSkipLength))
    eu.eprint('Total ambiguous length      :{0:>11}'.format(tAmbigLen))
    eu.eprint('Total ambig + unmap length  :{0:>11}'.format(tAmbigLen + tUnmapLength))
    eu.eprint("------------------ --------------- ------------------")

def verify(fileName,outfile,unmappedFile,ambigFile,genomeFile):
    # TODO:
    # Load Genome
    # For each read in SAM file
    pass

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("samFile", 
                        help="/PATH/TO/INPUT-SAM-file (or) - (stdin)")
    parser.add_argument("outputErrFile",
                        help="/PATH/TO/OUT-ERR-file (or) - (stdout)")
    parser.add_argument("unmappedFile", 
                        help="/PATH/TO/OUT-unmapped-file")
    parser.add_argument("-a", "--ambig",
                        help="/PATH/TO/OUT-ambigous-file")
    parser.add_argument("-t", "--trim",
                        help="/PATH/TO/OUT-trimmed-file")
    parser.add_argument("-g", "--genome",
                        help="/PATH/TO/INPUT-genome-file")
    parser.add_argument("-d", "--dry_run", action="store_true")
    try:
        args = parser.parse_args()
        # opts,args = getopt.getopt(argv,"f:o:u:a:t:g:dvh",["file=","outfile=",
        #                                                   "unmapped=", "ambig=",
        #                                                   "trim=","genome=",
        #                                                   "dry" "verify","help"])
    except argparse.ArgumentError as err:
        eu.eprint(str(err))
        parser.print_help()
        sys.exit(2)
    # Get the command line arguments
    trimFile = ambigFile = filename = outfile = unmappedFile = None
    genomeFile = None
    dryRun = verify = False
    filename = args.samFile
    outfile = args.outputErrFile
    unmappedFile = args.unmappedFile
    ambigFile = args.ambig
    genomeFile = args.genome
    trimFile = args.trim
    dryRun = args.dry_run
    eu.eprint("------------------ INPUT ARGUMENTS ------------------")
    eu.eprint("IN SAM File (sorted)       : ", filename)
    eu.eprint("IN Ref Genome File         : ", genomeFile)
    eu.eprint("OUT Mapped Reads File      : ", outfile)
    eu.eprint("OUT Unmapped Reads File    : ", unmappedFile)
    eu.eprint("OUT Ambigous Reads File    : ", ambigFile)
    eu.eprint("OUT Trim File              : ", trimFile)
    eu.eprint("Y/N Dry Run                : ", dryRun)
    eu.eprint("------------------ --------------- ------------------")
    if filename == None or outfile == None or unmappedFile == None:
        parser.print_help()
        sys.exit()
    if verify:
        if genomeFile == None:
            parser.print_usage()
            sys.exit()
        #verify(filename,outfile,unmappedFile,ambigFile,genomeFile)
        eu.eprint("--- VERIFY FUNCTION NOT IMPLEMENTED ---")
    else:
        # Process the file
        process(filename,outfile,unmappedFile,ambigFile,
                genomeFile,trimFile,dryRun)

if __name__ == "__main__":
    main(sys.argv[1:])

