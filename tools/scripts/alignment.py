import numpy as np
import sys
import re
import ecutils as eu
#
#  File : alignment.py
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
# START : Global Alignment code
# Following code is taken from 
#  https://bitbucket.org/brentp/biostuff
# 
UP, LEFT, DIAG, NONE = range(4)

MATRIX = { }
INFTY = 16384


def read_matrix(path):
    if path in MATRIX: return MATRIX[path]
    m = {}
    fh = open(path)
    headers = None
    while headers is None:
        line = fh.readline().strip()
        if line[0] == '#': continue
        headers = [x for x in line.split(' ') if x]
        for h in headers: m[h] = {}

    line = fh.readline()
    while line:
        h1 = line[0]
        line = [int(x) for x in line[1:-1].split(' ') if x]
        values = zip(headers, line)
        m[h1] = dict(values)
        line = fh.readline()
    return m

def global_nw(seqj, seqi, gap=-1, matrix=None, match=1, mismatch=-1):
    """
    >>> global_align('COELANCANTH', 'PELICAN')
    ('COELANCANTH', '-PEL-ICAN--')
    """
    max_j = len(seqj)
    max_i = len(seqi)
    if matrix is not None:
        matrix = read_matrix(matrix)
  
    score   = np.zeros((max_i + 1, max_j + 1), dtype='f')
    pointer = np.zeros((max_i + 1, max_j + 1), dtype='i')
    max_i, max_j

    pointer[0, 0] = NONE
    score[0, 0] = 0.0

    
    pointer[0, 1:] = LEFT
    pointer[1:, 0] = UP

    score[0, 1:] = gap * np.arange(max_j)
    score[1:, 0] = gap * np.arange(max_i)
    
    for i in range(1, max_i + 1):
        ci = seqi[i - 1]
        for j in range(1, max_j + 1):
            cj = seqj[j - 1]

            if matrix is None:
                diag_score = score[i - 1, j - 1] + (cj == ci and match or mismatch)
            else:
                diag_score = score[i - 1, j - 1] + matrix[cj][ci]

            up_score   = score[i - 1, j] + gap
            left_score = score[i, j - 1] + gap
            
            if diag_score >= up_score:
                if diag_score >= left_score:
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                else:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT

            else:
                if up_score > left_score:
                    score[i, j ]  = up_score
                    pointer[i, j] = UP
                else:
                    score[i, j]   = left_score
                    pointer[i, j] = LEFT
                    
                
    align_j = ""
    align_i = ""
    while True:
        p = pointer[i, j]
        if p == NONE: break
        s = score[i, j]
        if p == DIAG:
            align_j += seqj[j - 1]
            align_i += seqi[i - 1]
            i -= 1
            j -= 1
        elif p == LEFT:
            align_j += seqj[j - 1]
            align_i += "-"
            j -= 1
        elif p == UP:
            align_j += "-"
            align_i += seqi[i - 1]
            i -= 1
        else:
            raise Exception('wtf!')

    return align_j[::-1], align_i[::-1]

# END : Global Alignment code
# Following code is taken from 
#  https://bitbucket.org/brentp/biostuff
# 

def global_banded(seqj, seqi, kband=5, gap=-1, 
                  matrix=None, match=1, mismatch=-1):
    """
    >>> global_align('COELANCANTH', 'PELICAN')
    ('COELANCANTH', '-PEL-ICAN--')
    """
    len_j = len(seqj)
    len_i = len(seqi)

    if matrix is not None:
        matrix = read_matrix(matrix)
  
    score = {}
    pointer = {}

    pointer[(0, 0)] = NONE
    score[(0, 0)] = 0.0

    dwidth = len_j - len_i
    #print 'dwidth ' + str(dwidth)

    for i in range(1,kband+abs(dwidth)+2):
        pointer[(0,i)] = LEFT
        pointer[(i,0)] = UP
        score[(i,0)] = score[(0,i)] = gap * i
    #print score

    for i in range(1,len_i+1):
        dstart = -1
        dend = -1
        if dwidth >= 0:
            if (i - kband) > 0:
                dstart = i - kband
            else:
                dstart = 1
            if (i + dwidth + kband > len_j):
                dend = len_j
            else:
                dend = i + kband + dwidth
        else:
            if (i + dwidth - kband) > 0:
                dstart = i + dwidth - kband
            else:
                dstart = 1
            if (i + kband > len_j):
                dend = len_j
            else:
                dend = i + kband
        #print (i,dstart,dend)
        if dstart > 1:
            score[(i,dstart-1)] = -INFTY
        if dend < len_j:
            score[(i,dend+1)] = -INFTY
        ci = seqi[i - 1]
        #continue
        #print score
        for j in range(dstart, dend+1):
            cj = seqj[j - 1]

            if matrix is None:
                diag_score = score[(i - 1, j - 1)] + (cj == ci and match or mismatch)
            else:
                diag_score = score[(i - 1, j - 1)] + matrix[cj][ci]

            up_score   = score[(i - 1, j)] + gap
            left_score = score[(i, j - 1)] + gap
            
            if diag_score >= up_score:
                if diag_score >= left_score:
                    score[(i, j)] = diag_score
                    pointer[(i, j)] = DIAG
                else:
                    score[(i, j)] = left_score
                    pointer[(i, j)] = LEFT

            else:
                if up_score > left_score:
                    score[(i, j) ]  = up_score
                    pointer[(i, j)] = UP
                else:
                    score[(i, j)]   = left_score
                    pointer[(i, j)] = LEFT
                    
                
    #return 0
    #print i,j
    align_j = ""
    align_i = ""
    while True:
        p = pointer[(i, j)]
        #print (i,j,p)
        if p == NONE: break
        s = score[(i, j)]
        if p == DIAG:
            align_j += seqj[j - 1]
            align_i += seqi[i - 1]
            i -= 1
            j -= 1
        elif p == LEFT:
            align_j += seqj[j - 1]
            align_i += "-"
            j -= 1
        elif p == UP:
            align_j += "-"
            align_i += seqi[i - 1]
            i -= 1
        else:
            raise Exception('wtf!')

    return align_j[::-1], align_i[::-1]

def global_banded_rev(seqj, seqi, revFlag = False, kband=5, gap=-1, 
                      matrix=None, match=1, mismatch=-1):
    if revFlag :
        rseqj = seqj[::-1]
        rseqi = seqi[::-1]
        rjaln,rialn = global_banded(rseqj,rseqi,kband,gap,matrix,match,mismatch)
        return rjaln[::-1],rialn[::-1]
    return global_banded(seqj,seqi,kband,gap,matrix,match,mismatch)

def global_nw_rev(seqj, seqi, revFlag = False, gap=-1, 
                  matrix=None, match=1, mismatch=-1):
    if revFlag :
        rseqj = seqj[::-1]
        rseqi = seqi[::-1]
        rjaln,rialn = global_nw(rseqj,rseqi,gap,matrix,match,mismatch)
        return rjaln[::-1],rialn[::-1]
    return global_nw(seqj,seqi,gap,matrix,match,mismatch)

#print global_banded_align("TTTAATTCAGGTATTGG","TTAATTCAGGTAT",1)


def global_agap(seqi, seqj, gap=-3, gapxn=-1,
                matrix=None, match=1, mismatch=-1):
    """
    >>> global_agap('COELANCANTH', 'PELICAN')
    ('COELANCANTH', '-PEL-ICAN--')
    """
    len_i = len(seqi)
    len_j = len(seqj)

    if matrix is not None:
        matrix = read_matrix(matrix)
  
    score = [{},{},{}]
    pointer = {}

    pointer[(0, 0)] = (0,NONE)
    score[1][(0, 0)] = score[2][(0, 0)] = score[0][(0, 0)] = 0.0

    for i in range(1, len_i+1):
        pointer[(i,0)] = (2,UP)
        score[0][(i,0)] = score[1][(i,0)] = score[2][(i,0)] = gap + i * gapxn
    for i in range(1, len_j+1):
        pointer[(0,i)] = (1,LEFT)
        score[0][(0,i)] = score[2][(0,i)] = score[1][(0,i)] = gap + i * gapxn
    score[0][(0,0)] = score[1][(0,0)] = score[2][(0,0)] = 0.0

    for i in range(1,len_i+1):
        ci = seqi[i - 1]
        for j in range(1, len_j+1):
            cj = seqj[j - 1]
            # left score
            score[1][(i, j)] = max(score[0][(i, j - 1)] + gap + gapxn,
                                   score[1][(i, j - 1)] + gapxn)
            # up score
            score[2][(i, j)] = max(score[0][(i-1,j)] + gap + gapxn,
                                   score[2][(i-1,j)] + gapxn)
            diag_score = 0
            if matrix is None: 
                adg = (cj == ci and match or mismatch)
            else:
                adg = matrix[cj][ci]
            diag_score = score[0][(i - 1, j - 1)] + adg
            left_score = score[1][(i,j)]
            up_score = score[2][(i,j)]
            score[0][i,j] = max(diag_score,up_score,left_score)
            #print (i,j), diag_score,up_score,left_score

            if diag_score > up_score:
                if diag_score > left_score: 
                    pointer[(i, j)] = (diag_score,DIAG)
                else:
                    pointer[(i, j)] = (left_score,LEFT)
            else:
                if up_score > left_score:
                    pointer[(i, j)] = (up_score,UP)
                else:
                    pointer[(i, j)] = (left_score,LEFT)

    #return 0
    #print i,j
    align_j = ""
    align_i = ""
    i = len_i; j = len_j
    while True:
        (s,p) = pointer[(i, j)]        
        #print (i,j,p)
        if p == NONE: break 
        if p == DIAG:
            align_j += seqj[j - 1]
            align_i += seqi[i - 1]
            i -= 1
            j -= 1
        elif p == LEFT:
            align_j += seqj[j - 1]
            align_i += "-"
            j -= 1
        elif p == UP:
            align_j += "-"
            align_i += seqi[i - 1]
            i -= 1
        else:
            raise Exception('this should not happen!')
    return align_i[::-1], align_j[::-1]

def global_agap_rev(seqj, seqi, revFlag = False,
                    gap=-3, gapxn=-1,
                    matrix=None, match=1, mismatch=-1):
    if revFlag :
        rseqj = seqj[::-1]
        rseqi = seqi[::-1]
        rjaln,rialn = global_agap(rseqj,rseqi,gap,gapxn,matrix,match,mismatch)
        return rjaln[::-1],rialn[::-1]
    return global_agap(seqj,seqi,gap,gapxn,matrix,match,mismatch)

#
# Utility functions to parse MD String
#  of an alignment
def getMDStringComps(mdString):
    regExMDInit = '([0-9]+)'
    regExMDRest =  '([' + eu.alphabetString + ']+|\^[' + eu.alphabetString + ']+)([0-9]+)'
    # Parse MD String to get its components
    reM = re.match(regExMDInit, mdString)
    mdComps = [reM.group(1)]
    comps = re.findall(regExMDRest, mdString)
    for c,d in comps:
        mdComps += [c,d]
    return mdComps

#
# Utility function to parse CIGAR string of an 
#  alignment
def getCIGARPairs(cigarString):
    regExCIGAR = '([0-9]+)([MIDNSHPX=])'
    cigarPairs = []
    # Parse the CIGAR String to get the components
    cpairs = re.findall(regExCIGAR, cigarString)
    for l,m in cpairs:
        cigarPairs += [[int(l),m]]
    return cigarPairs

#
# Builds the reference alignment from
#  (i)  CIGAR Pairs (parsed from getCIGARPairs),
#  (ii) MD components (parsed from getMDStringComps(...)) 
# and the read string
def getRefAlignment(cigarPairs,mdComps,readString):
    # Fix up the ref alignment string
    sofar = 0
    refAlign = ''
    [fl,fc] = cigarPairs[0]
    if fc == 'S':
        sofar += fl
        cigarPairs = cigarPairs[1:]
    for [l,c] in cigarPairs:
        if c == 'M':
            refAlign += readString[sofar:sofar+l]
            sofar += l
        elif c == 'S':
            sofar += l
        elif c == 'I':
            refAlign += l * '-' # nothin in ref
            sofar += l
        elif c == 'D':
            refAlign += l * 'N' # something in ref
    # Update the ref align based on MD string
    sofar = 0
    alnidx = 0
    alnlen = len(refAlign)
    refAlign2 = ''
    #print refAlign
    for c in mdComps:
        if c.isdigit():
            # If we have digits, we just move along
            # copying the strings
            limit = sofar + int(c)
            while sofar < limit:
                refAlign2 += refAlign[alnidx]
                if refAlign[alnidx] != '-':
                    sofar += 1
                alnidx += 1
        elif c[0] == '^': #deleted in read
            # If deleted in read, it is present
            refAlign2 += c[1:]
            alnidx += len(c[1:])
        else: # Substitution
            refAlign2 += c
            alnidx += len(c)
        while alnidx < alnlen and refAlign[alnidx] == '-':
            refAlign2 += refAlign[alnidx]
            alnidx += 1
    #print refAlign2
    return refAlign2

#
# Get the reference alignment from 
#   - genome/chromosome string
#   - CIGAR info (parsed from CIGAR string with getCIGARPairs(...)
#   - starting position at which the genome is aligned at
#   - read string
def getGenomeRefAlign(gstr,gpos,cigarPairs,readString):
    sofar = gpos-1
    refAlign = ''
    [fl,fc] = cigarPairs[0]
    if fc == 'S':
        cigarPairs = cigarPairs[1:]
    for [l,c] in cigarPairs:
        if c == 'M':
            refAlign += gstr[sofar:sofar+l]
            sofar += l
        elif c == 'S':
            sofar += l
            break # Skipping should be the last
        elif c == 'I':
            refAlign += l * '-' # nothin in ref
        elif c == 'D':
            refAlign += gstr[sofar:sofar+l] # something in ref
            sofar += l
    #x = len(refAlign)
    #print gstr[gpos-1:gpos+x]
    return refAlign

#
# Get the alignment of read string from
#  CIGAR info (parsed from CIGAR string with getCIGARPairs(..) fn),
#  and read string
def getReadAlignment(cigarPairs,readString):
    # Fix up the read alignment string
    sofar = 0
    readAlign = ''
    skiplength = 0
    readlen = len(readString)
    skipbegin = 0
    skipend = readlen
    [fl,fc] = cigarPairs[0]
    if fc == 'S':
        skipbegin = fl
        skiplength += fl
        sofar += fl
        cigarPairs = cigarPairs[1:]
    for [l,c] in cigarPairs:
        if c == 'M':
            readAlign += readString[sofar:sofar+l]
            sofar += l
        elif c == 'S':
            skipend = sofar
            sofar += l
            skiplength += l
        elif c == 'I':
            readAlign += readString[sofar:sofar+l] # some thing in read
            sofar += l
        elif c == 'D':
            readAlign += l * '-' # nothing in read
    #print readAlign
    #if skiplength > 0: #(0.1 * readlen):
    #    raise SkipError(str(skiplength))
    return (readAlign,skipbegin,skipend)

def getSAMAlignment(readString, cigarString, mdString,gstr,gpos):
    readLength = len(readString)
    # Parse CIGAR and MD String
    cigarPairs = getCIGARPairs(cigarString)
    # Fix up the read and ref alignment string
    (readAlign,skipbegin,skipend) = getReadAlignment(cigarPairs,readString)
    refAlign = ''
    if mdString != None:
        mdComps = getMDStringComps(mdString)
        refAlign = getRefAlignment(cigarPairs,mdComps,readString)
        #print mdString,mdComps,refAlign
    else:
        refAlign = getGenomeRefAlign(gstr,gpos,cigarPairs,readString)
    if(len(refAlign) != len(readAlign)):
        print 'length doesnt match'
        print readAlign
        print refAlign
    assert len(refAlign) == len(readAlign)
    #print readAlign
    #print refAlign
    return [refAlign, readAlign,skipbegin,skipend]
