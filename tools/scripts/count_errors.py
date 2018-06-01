#!/usr/bin/python
#
#  File : count_errors.py
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

import sys
if len(sys.argv) < 2 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
    print 'Usage : count_errors.py /path/to/alignment.er'
    print 'Script counts the number of errors with alignment' 
    sys.exit(1)

fname = sys.argv[1]
t = 0
with open(fname,'r') as f:
   for line in f:
      elts = line.split('\t')
      ne = int(elts[1])
      t += ne
print t
