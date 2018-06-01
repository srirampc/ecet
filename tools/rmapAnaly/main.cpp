

/***
 **
 *  File: main.cpp
 *  Created: December 17, 2009
 *
 *  Author: Xiao Yang <isuxyang@gmail.com>
 *
 *  Copyright 2009 Xiao Yang
 *
 *  This file is part of rmapAnaly.
 *
 *  rmapAnaly is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  rmapAnaly is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libpnorm. If not, see <http://www.gnu.org/licenses/>.
 *
 */
 
/*
 * Take input three files: [.bed] [NC_.fa] [sr.fa]
 * 
 * 1) bed file from RMAP
 * 2) reference genome
 * 3) short read file
 *
 * Output four files:
 *
 * 1) non uniquely mapped read IDs, one ID per line
 *
 * 2) a file capturing every read errors in the format of
 *    readID  ErrNum  [Position from  To  Qual] [Position from  To  Qual]...
 *    ... ...
 *
 *    Note currently all Qual values are set to be 0
 * 
 * 3) a PE or positional dependent distribution (PDD) FILE capturing
 *    mutations for each kmer position in the format:
 *    ----------------------------------------------------------
 *    pos0      aa ac ag at ca cc cg ct ga gc gg gt ta tc tg tt
 *    ... ...
 *    pos (k-1)
 *    ----------------------------------------------------------
 * 
 * 4) read distribution (mapping) on the genome
 *  [pos i, # reads starting at i]
 */

#include "rmapAnaly.h"
#include "Errlocator.h"

int main(int argc, char** argv) {

    Para myPara(argc, argv);

    Errlocator myEL;

    std::cout << "\nLoading Bed ...\t";
    //myEL.loadBed(myPara.bedName.c_str());
    myEL.loadBed(myPara);
    std::cout << "done\n\n";

    if (!myPara.rDistrName.empty()){
        std::cout << "calculate read mapping distribution...\t";
        myEL.rDistr(myPara.rDistrName.c_str());
        std::cout << "done\n\n";
    }

    std::cout << "Loading Genome ...\t";
    myEL.loadGn(myPara.gnName);
    std::cout << "done\n\n";

    std::cout << "start resolving errs...\n";
    myEL.resolveErr(myPara);    
    std::cout << "\nprogram finished successfully!\n\n";
    return (EXIT_SUCCESS);
}
