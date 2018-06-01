/* 
 *  File:   compSR.h
 *  Created on December 19, 2009, 2:51 PM
 *  Author: Xiao Yang <isuxyang@gmail.com>
 *  Copyright 2009 Xiao Yang
 *
 *  This file is part of Error Correction Review Toolkit.
 *  Error Correction Review Toolkit is free software: you can 
 *  redistribute it and/or modify  it under the terms of the GNU 
 *  Lesser General Public License as published by 
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Error Correction Review Toolkit is distributed in the hope that 
 *  it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 *  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libpnorm. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _COMPSR_H
#define	_COMPSR_H



#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include <stdlib.h>
#include <algorithm>
#include <ctype.h>
#include <set>
#include <map>

typedef std::vector<int> ivec_t;
typedef std::set<int> iset_t;
typedef std::map<int,int> imap_t;

typedef struct ERRINFO{
    int pos;
    int from;
    int to;
    int qual;
    ERRINFO (int p, int f, int t, int q): pos(p), from(f), to(t), qual(q) {}
    ERRINFO (){}
}e_t;
typedef std::vector<e_t> evec_t;

typedef struct RECORD{
    int readID;
    evec_t evec;
    RECORD (int ID, evec_t vec): readID(ID), evec(vec) {}
}record_t;

class Para {
public:
    std::string ecName;
    std::string pcAlignName;
    std::string rsltName;
    std::string ambigName; //ambiguous file by pre-correction alignment
    std::string trimFileName;
    int mvalue;
    Para(int argc, char** argv) : argnum(argc), arg(argv) {
        setPara();
    };

private:
    int argnum;
    char** arg;


    void printUsg(){
        std::cout << "\n----------------------------------\n";
        std::cout << "Usg: ./app [correction-rslt] [pre-correct-aln-rslt] "
                  << "[unmapped-pre-correct-aln] [m-value] [fpfn-rslt] ";
        std::cout << "[optional trimmed-file]\n";
        std::cout << "----------------------------------\n\n";
        exit(1);
    }

    void setPara() {
        trimFileName = "";
        if (argnum < 6) { printUsg();}
        ecName = arg[1];
        pcAlignName = arg[2];
        ambigName = arg[3];
        mvalue = atoi(arg[4]);
        rsltName = arg[5];
        if(argnum > 6)
            trimFileName = arg[6];
    }
};

typedef struct FPNP{
    int FP;
    int FN;
    int TP;
    int WrongBp;
    int WrongBp_N; // consider only ambiguous bases
    int N_total; // total bases of ambigous bases
    int cntg_m; // the target read that contains >m errors that are not one of the ref read
    FPNP () {
        FP = 0; FN = 0; TP = 0; WrongBp = 0; cntg_m = 0;
        WrongBp_N = 0;  N_total = 0;
    }
}fpnp_t;

void comparison(const char* filename, const std::vector<record_t>& ref,
                const std::vector<record_t>& tar, const iset_t& ambig,
                const imap_t& trimmed, const imap_t& trimPrefix, int mvalue);
void checkRead(fpnp_t& rslt, int readID, 
               const imap_t& trimmed, const imap_t& trimPrefix,
               const evec_t& ref, const evec_t& tar);
void output(const fpnp_t& rslt, int mvalue, const char* filename);

#endif	/* _COMPSR_H */

