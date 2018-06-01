from mpi4py import MPI
#
#  File : ecutils.py
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

alphabet = ['A','C','G','T','N','R','Y']
alphabetString = "".join(alphabet)
complementMapping = {'A':'T','C':'G','T':'A','G':'C','-':'-','N':'N','Y':'Y','R':'R'}

def load_genome(gfname):
    fgenome = {}
    with open(gfname, 'r') as f:
        chromosome = ''
        i = 0
        for line in f:
            if line[0] == '>':
                if len(chromosome) > 0:
                    fgenome[i] = chromosome
                    chromosome = ''; i += 1
            else:
                chromosome += line.strip()
        if len(chromosome) > 0:
            fgenome[i] = chromosome
    return fgenome

def reverse_complement(inRead):
    ls = map(lambda x:complementMapping[x],inRead)
    rs = ''.join(ls)
    return rs[::-1]

def decompose(total,rank,size):
    #global psize, prank
    block_decomp = [0 for i in range(size)]
    block_decomp_sum = [0 for i in range(size)]
    startid = 0
    endid = total-1
    tmp = total/size
    j = total % size
    for i in range(size):
        if i < j:
            block_decomp[i] = tmp + 1
        else:
            block_decomp[i] = tmp
    bsum = [0 for i in range(size)]
    for i in range(size)[1:]:
        bsum[i] = bsum[i-1] + block_decomp[i]
    startid = bsum[rank]
    if rank == size - 1:
        endid = total - 1
    else:
        endid =  bsum[rank+1] - 1
    block_decomp_sum = bsum
    if rank == 0:
        print bsum
        print block_decomp
    MPI.COMM_WORLD.barrier()
    for i in range(size):
        MPI.COMM_WORLD.barrier()
        if i == rank:
            print (rank, startid, endid)
        MPI.COMM_WORLD.barrier()
    return (startid,endid)

def getfaid(lines):
    elts = lines[0].split()
    readid = int(elts[0][1:])
    return readid

def getSAMinfo(line):
    elts = line.strip().split()
    sminf = {}
    sminf['QNAME'] = elts[0]
    sminf['FLAG'] = int(elts[1])
    sminf['RNAME'] = elts[2]
    sminf['POS'] = int(elts[3])
    sminf['MAPQ'] = int(elts[4]);
    sminf['CIGAR'] = elts[5]
    sminf['SEQ'] = elts[9]
    sminf['XT'] = sminf['AS'] = sminf['MD'] = None
    for s in elts:
        if len(s) > 2 and s[0:2] == 'XT':
            sminf['XT'] = s
        if len(s) > 2 and s[0:2] == 'AS':
            sminf['AS'] = s
        if len(s) > 2 and s[0:2] == 'MD':
            sminf['MD'] = s
    if( sminf['CIGAR'] != '*' and sminf['MAPQ'] != 0):
        sminf['SEQ'] = elts[9]
    return sminf

