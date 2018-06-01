#!/usr/bin/python
import sys
import getopt
import re
import itertools
import alignment as an
import ecutils as eu
from mpi4py import MPI
#
#  File : compute-stats.py
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
# Other Contributors:
#  1. Giorgio Gonnella fix for incorrect handling of skipped bases

psize = MPI.COMM_WORLD.Get_size()
prank = MPI.COMM_WORLD.Get_rank()
pname = MPI.Get_processor_name()

charMapping = {'A':'0','C':'1','G':'2','T':'3','N':'4','-':'5',
               'a':'0','c':'1','g':'2','t':'3','n':'4'}
alphabet = ['A','C','G','T','N']
error_types = [0,1,2]
fgenome = {0:None}
errorStats = {'FP': 0, 'FN':0, 'TP':0, 'NE':0, 'SR':0}
totalPosNotEqual = 0
#
#
#
def getNumWrongBase(errPreCorrect,errPostCorrect):
    preh = {} ; ne = 0
    for x in errPreCorrect:
        preh[(x[0],x[1],x[3])] = x[2]
    for x in errPostCorrect:
        try:
            if preh[(x[0],x[1],x[3])] != x[2]:
                ne += 1
        except KeyError as e:
            pass
    return ne

# Errors are based on
# subst error
# deletion error (char in genome missing in the read)
# insertion error ( extra char in read after gpos in genome)
def getErrorsFromAlignment(refAlign,readAlign):
    errors = []
    gpos = 0; erridx = 0
    for x,y in zip(refAlign,readAlign):
        if x != y:
            tb = x # true base
            wb = y  # wrong base
            errors += [(gpos,erridx,tb,wb)]
            erridx += 1
        if x != '-':
            gpos += 1
            erridx = 0
    return set(errors)

def processMissing(sminf,gstr,rid):
    global errorStats
    #print 'missing', rid
    mapq = sminf['MAPQ']
    #  0. handle unmapped and trimmed section!
    if(sminf['FLAG'] & 4 == 4): return # not unique
    if(mapq == 0): return  # can not map
    cigarStr = sminf['CIGAR']
    if (cigarStr == '*' or mapq == 0):
        print 'Wrong * and 0 :', sminf['NAME']; return
    mdStr = sminf['MD'] and sminf['MD'].split(':')[-1] or None
    #  1. Build Alignment ref and original read
    oread = sminf['SEQ']; olen = len(oread)
    [refAlign,readAlign,sb,se] = an.getSAMAlignment(oread,cigarStr, mdStr,
                                                    gstr, sminf['POS'])
    #  Compute line stats
    errPreCorrect = getErrorsFromAlignment(refAlign,readAlign)
    errorStats['FN'] += len(errPreCorrect)

def update_pos(errors,upd):
    return set([(gpos+upd,erridx,tb,wb) for (gpos,erridx,tb,wb) in errors])

def process_sam_line(sam_info,corr_sam_info,gstr,rid):
    global errorStats,totalPosNotEqual
    mapq = sam_info['MAPQ']
    #  0. verify alignments
    if(sam_info['FLAG'] & 4 == 4): return # not unique
    if(mapq == 0): return  # can not map
    cigarStr = sam_info['CIGAR']
    if (cigarStr == '*' or mapq == 0):
        print 'Wrong * and 0 :', sam_info['NAME']; return
    corr_cigarStr = corr_sam_info['CIGAR']
    # if corrected read is not aligned or ambig alinged ,
    # assume it is not corrected
    corr_mapq = corr_sam_info['MAPQ']
    if ((corr_mapq == 0) or (corr_sam_info['FLAG'] & 4 == 4)):
        processMissing(sam_info,fgenome[0],rid)
        return;
    if (corr_cigarStr == '*' or corr_mapq == 0):
        print 'Wrong * and 0 :', corr_sam_info['NAME']; return
    # assume it is not corrected, if not aligned in the same pos
    if (sam_info['RNAME'] != corr_sam_info['RNAME']):
        processMissing(sam_info,fgenome[0],rid)
        return
    posdiff = 0
    if (sam_info['POS'] != corr_sam_info['POS']):
        # print 'pos not same : ',rid
        # print sam_info['RNAME'], corr_sam_info['RNAME']
        # print sam_info['POS'], corr_sam_info['POS']
        posdiff = int(sam_info['POS']) - int(corr_sam_info['POS'])
        if (abs(posdiff) > 5):
            totalPosNotEqual += 1
            processMissing(sam_info,fgenome[0],rid)
            return
    #  1. Build Alignment ref and original read
    mdStr = sam_info['MD'] and sam_info['MD'].split(':')[-1] or None
    oread = sam_info['SEQ']; olen = len(oread)
    [refAlign,readAlign,sb,se] = an.getSAMAlignment(oread,cigarStr, mdStr,
                                                    gstr, sam_info['POS'])
    # 2. Build Alignment for ref and corrected read
    corr_mdStr = corr_sam_info['MD'] and corr_sam_info['MD'].split(':')[-1] or None
    cread = corr_sam_info['SEQ']; clen = len(cread)
    [cgAlign,corAlign,csb,cse] = an.getSAMAlignment(cread, corr_cigarStr,
                                                    corr_mdStr,
                                                    gstr, corr_sam_info['POS'])
    #print rid,posdiff,'\n'
    #print '\n'.join([refAlign,readAlign])
    #print '\n'.join([cgAlign, corAlign])
    # 3. Compute line stats
    errPreCorrect = getErrorsFromAlignment(refAlign,readAlign)
    errPostCorrect = getErrorsFromAlignment(cgAlign,corAlign)
    if (posdiff < 0):
        errPostCorrect = update_pos(errPostCorrect,abs(posdiff))
    elif(posdiff > 0):
        errPreCorrect = update_pos(errPreCorrect,abs(posdiff))
    #print errPreCorrect, '\n', errPostCorrect
    errorStats['TP'] += len(errPreCorrect.difference(errPostCorrect))
    errorStats['FP'] += len(errPostCorrect.difference(errPreCorrect))
    errorStats['FN'] += len(errPreCorrect.intersection(errPostCorrect))
    errorStats['NE'] += getNumWrongBase(errPreCorrect,errPostCorrect)
    #print sam_info['QNAME'],errorStats

def processLine(sminf,cread,gstr,band):
    global errorStats
    mapq = sminf['MAPQ']
    #  0. handle unmapped and trimmed section!
    if(sminf['FLAG'] & 4 == 4): return # not unique
    if(mapq == 0): return  # can not map
    cigarStr = sminf['CIGAR']
    if (cigarStr == '*' or mapq == 0):
        print 'Wrong * and 0 :', sminf['NAME']; return
    mdStr = sminf['MD'] and sminf['MD'].split(':')[-1] or None
    #  1. Build Alignment ref and original read
    oread = sminf['SEQ']; olen = len(oread)
    [refAlign,readAlign,sb,se] = an.getSAMAlignment(oread,cigarStr, mdStr,
                                                    gstr, sminf['POS'])
    flag = sminf['FLAG']
    if ( flag & 16 == 16):
        refAlign = eu.reverse_complement(refAlign)
        readAlign = eu.reverse_complement(readAlign)
        oread = eu.reverse_complement(oread)
        # Giorgio Gonnella fix for incorrect handling of skipped bases
        swap_tmp = se
        se = len(oread) - sb
        sb = len(oread) - swap_tmp
    # Giorgio Gonnella fix for incorrect handling of skipped bases
    if(olen != (se+sb)):
        oread = oread[sb:se]
        cread = cread[sb:se]
    #  2. Get the string in the genomic region
    gRegion = reduce( lambda x,y: x+y,
                      map(lambda x: x not in ['-'] and x or '', refAlign))
    # 3. Do banded alignment
    cgAlign = ''; corAlign = ''
    #print "\n".join(["reads",oread,cread])
    if (gRegion == cread): # perfect correction!
        cgAlign = corAlign = cread
    elif(oread == cread): # if both are same, the alignment is the same
        cgAlign = refAlign
        corAlign = readAlign
    else:
        cgAlign,corAlign = an.global_banded_rev(gRegion,cread,flag & 16 == 16,band)
        #cgAlign,corAlign = an.global_agap_rev(gRegion,cread,flag & 16 == 16)
        #print '\n'.join([refAlign,readAlign])
        #print '\n'.join([cgAlign, corAlign])
    # 4. Compute line stats
    errPreCorrect = getErrorsFromAlignment(refAlign,readAlign)
    errPostCorrect = getErrorsFromAlignment(cgAlign,corAlign)
    #print errPreCorrect, '\n', errPostCorrect
    errorStats['TP'] += len(errPreCorrect.difference(errPostCorrect))
    errorStats['FP'] += len(errPostCorrect.difference(errPreCorrect))
    errorStats['FN'] += len(errPreCorrect.intersection(errPostCorrect))
    errorStats['NE'] += getNumWrongBase(errPreCorrect,errPostCorrect)
    #print sminf['QNAME'],errorStats

def printStats(outFile):
    global errorStats
    completeStats = {}
    for x in ['TP','FP','FN','NE']:
        data = errorStats[x]
        data = MPI.COMM_WORLD.gather(data,root=0)
        if prank == 0:
            completeStats[x] = sum(data)
    total = totalPosNotEqual
    total = MPI.COMM_WORLD.gather(total,root=0)
    if prank == 0:
        print "Total Pos eq :", sum(total)
    if prank == 0:
        eba = gain = 0.0
        try:
            gain = float(completeStats['TP'] - completeStats['FP']) / (completeStats['TP'] + completeStats['FN'])
        except ZeroDivisionError:
            print 'gain is not calculated due to zero division error'
        try:
            eba = float(completeStats['NE']) / (completeStats['NE'] + completeStats['TP'])
        except ZeroDivisionError:
            print 'eba is not calculated due to zero division error'
        outf = open(outFile, 'w')
        for x in completeStats.keys():
            outf.write(str(x) + ' : ' + str(completeStats[x]) + '\n')
        outf.write('EBA   : ' + str(eba) + '\n')
        outf.write('GAIN  : ' + str(gain) + '\n')
        outf.close()

def process_sam(alignFile,correctedAlignFile,genomeFile,outFile,total):
    global fgenome,prank,psize,errorStats
    if genomeFile != None:
        fgenome = eu.load_genome(genomeFile)
    alnf = open(alignFile, 'r')
    corralnf = open(correctedAlignFile, 'r')
    (startid, endid) = eu.decompose(total,prank,psize)
    try:
        smline = alnf.next()
        while smline[0] == '@':
            smline = alnf.next()
        corrsmline = corralnf.next()
        while corrsmline[0] == '@':
            corrsmline = corralnf.next()
        sminf = eu.getSAMinfo(smline)
        corrsminf = eu.getSAMinfo(corrsmline)
        rid = int(sminf['QNAME'])
        corr_rid = int(corrsminf['QNAME'])
        while rid < startid:
            smline = alnf.next();
            sminf = eu.getSAMinfo(smline); rid = int(sminf['QNAME'])
        while corr_rid < startid:
            corrsmline = corralnf.next();
            corrsminf = eu.getSAMinfo(corrsmline); corr_rid = int(corrsminf['QNAME'])
        while True:
            while rid != corr_rid:
                #print 'missing ', rid
                processMissing(sminf,fgenome[0],rid)
                smline = alnf.next()
                sminf = eu.getSAMinfo(smline); rid = int(sminf['QNAME'])
            assert(rid == corr_rid)
            process_sam_line(sminf,corrsminf,fgenome[0],rid)
            if corr_rid >= endid:
                break
            smline = alnf.next()
            sminf = eu.getSAMinfo(smline); rid = int(sminf['QNAME'])
            corrsmline = corralnf.next()
            corrsminf = eu.getSAMinfo(corrsmline); corr_rid = int(corrsminf['QNAME'])
    except StopIteration as e:
        pass
    except IOError as e:
        print str(e)
    alnf.close()
    corralnf.close()
    #print errorStats
    printStats(outFile)

def process(alignFile,correctedFile,genomeFile,outFile,total,band):
    global fgenome,prank,psize,errorStats
    if genomeFile != None:
        fgenome = eu.load_genome(genomeFile)
    alnf = open(alignFile, 'r')
    corrf = open(correctedFile, 'r')
    (startid, endid) = eu.decompose(total,prank,psize)
    try:
        smline = alnf.next()
        while smline[0] == '@':
            smline = alnf.next()
        args = [iter(corrf)] * 2
        cfiter = itertools.izip(*args)
        corr = cfiter.next(); cid = eu.getfaid(corr)
        # Skip them
        while cid < startid:
            corr = cfiter.next(); cid = eu.getfaid(corr)
        sminf = eu.getSAMinfo(smline)
        rid = int(sminf['QNAME'])
        while rid < startid:
            smline = alnf.next();
            sminf = eu.getSAMinfo(smline); rid = int(sminf['QNAME'])
        while True:
            while rid != cid:
                #print 'missing ', rid
                processMissing(sminf,fgenome[0],rid)
                smline = alnf.next()
                sminf = eu.getSAMinfo(smline); rid = int(sminf['QNAME'])
            assert(rid == cid)
            #print sminf #, corr[1].strip()
            processLine(sminf,corr[1].strip(),fgenome[0],band)
            if cid >= endid:
                break
            corr = cfiter.next(); cid = eu.getfaid(corr)
            smline = alnf.next()
            sminf = eu.getSAMinfo(smline); rid = int(sminf['QNAME'])
    except StopIteration as e:
        pass
    except IOError as e:
        print str(e)
    alnf.close()
    corrf.close()
    #print errorStats
    printStats(outFile)

def usage():
    print """compute-stats.py --aln=/path/to/pre-correction-alignment-sam-file
                --outfile=/path/to/stats-output (write access reqd.)
                --records=number of reads
                [--genomeFile=/path/to/genome-file]
                [--band=value of k used for k-band alignment (default 5)]
                --corrected=/path/to/corrected-reads-fa-file (OR) --csaln=/path/to/post-correction-alignment-sam-file
              (OR)
compute-stats.py -a /path/to/pre-correction-alignment-sam-file
                -o /path/to/stats-output-file (write access reqd.)
                -r number of reads
                [-g /path/to/genome-file]
                [-b value of k used for k-band alignment (default 5)]
                -c /path/to/corrected-reads-fa-file (OR) -m /path/to/post-correction-alignment-sam-file

If the corrected reads are provided, compute-stats runs pairwise alignment
against the genomic region.

If post-correction alignment is given, the difference the genomic regional
alignment will be used to compute-stats. Providing a post-correction alignment
alignment file is faster than performing pairise alignment.
     """

def main(argv):
    try:
        opts,args = getopt.getopt(argv,"a:c:o:r:g:b:m:h",["aln=","corrected=",
                                                          "outfile=","records="
                                                          "genome=","band=",
                                                          "csaln=","help"])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    band = 5
    outFile = alignFile = correctedFile = None
    genomeFile = numRecords = corrAlignFile = None
    for opt,arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-a", "--aln"):
            alignFile = arg
        elif opt in ( "-c", "--corrected"):
            correctedFile = arg
        elif opt in ("-o", "--outfile"):
            outFile = arg
        elif opt in ("-r", "--records"):
            numRecords = arg
        elif opt in ("-g", "--genome"):
            genomeFile = arg
        elif opt in ("-b", "--band"):
            band = int(arg)
        elif opt in ("-m", "--csaln"):
            corrAlignFile = arg
        elif prank == 0:
            print opt, "is not a valid option"
            usage()
            sys.exit()
    if prank == 0:
        print "Alignment SAM pre-correction          : ", alignFile
        print "Corrected READS fasta post correction : ", correctedFile
        print "OUTPUT file                           : ", outFile
        print "Total Number of Reads                 : ", numRecords
        print "Ref Genome File                       : ", genomeFile
        print "k for k-band algorithm                : ", band
        print "Corrected Alignment File              : ", corrAlignFile
    if numRecords == None or outFile == None or alignFile == None:
        if prank == 0:
            usage()
        sys.exit()
    if corrAlignFile == None and correctedFile == None:
        if prank == 0:
            print 'Proved one of post-correction Alignment file or post-correction reads'
            usage()
        sys.exit()
    total = int(numRecords)
    if corrAlignFile != None:
        process_sam(alignFile,corrAlignFile,genomeFile,outFile,total)
    else:
        if band <= 0:
            if prank == 0:
                print "Band ", band, " <= 0, Using default value of 5."
                band = 5
        process(alignFile,correctedFile,genomeFile,outFile,total,band)

if __name__ == "__main__":
    main(sys.argv[1:])
