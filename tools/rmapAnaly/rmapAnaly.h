/* 
 * File:   rmapAnaly.h
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

#ifndef _RMAPANALY_H
#define	_RMAPANALY_H


#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include <set>
#include <map>
#include <stdlib.h>
#include <algorithm>
#include <ctype.h>
#include <stdint.h>
#include "fasta_file.hpp"

typedef std::map<int, int> imap_t;
typedef std::vector<int> ivec_t;
typedef std::set<int> iset_t;
typedef std::vector<uint64_t> u64vec_t;
typedef std::vector<u64vec_t> uu64vec_t;

typedef struct ERRINFO{
    int pos;
    int from;
    int to;
    int qual;
    ERRINFO (int p, int f, int t, int q): pos(p), from(f), to(t), qual(q) {}
}e_t;
typedef std::vector<e_t> evec_t;

typedef struct RSLT{
    int readID;
    evec_t evec;
    RSLT (int ID, evec_t vec): readID(ID), evec(vec) {}
}rslt_t;

typedef struct BEDINFO{
    int start;
    int end;
    int readID;
    int score;
    char strand;
    BEDINFO(int sta, int e, int r, int sc, int str): start(sta),
            end(e), readID(r), score (sc), strand(str) {};
    BEDINFO(){};
}bed_t;

typedef std::vector<bed_t> bedvec_t;

class Para {
public:
    std::string bedName;
    std::string gnName;
    std::string srName;
    //std::string qName;
    std::string rsltName;
    std::string unmappableName;
    int num_reads;
    int k_val;
    std::string peName; // output kmer error rate for each position
    std::string rDistrName; // output read distribution along the genome

    Para(int argc, char** argv) : argnum(argc), arg(argv) {
        k_val = 0;
        rDistrName = "";
        peName = "";
        setPara();
    };

private:
    int argnum;
    char** arg;


    void printUsg(){
        std::cout << "\n-------------------------------------------------\n";
        std::cout << "Usage: ./rmapAnaly [.bed] [NC_.fa] [sr.fa] [# reads] [errRecord]\n"
                  << "       [unmapped] [rDistrfile] [k_val] [PEfile(or PDDfile)] \n\n"
                  << "[.bed]: output from RMAP\n"
                  << "[NC_.fa]: reference genome, fasta format\n"
                  << "[sr.fa]: short read input, fasta format\n"
                  //<< "[sr.q]: short read quality scores, fasta format\n"
                  << "[errRecord]: recording all errors discovered by RMAP for\n"
                  << "\t mapping each read to reference genome\n"
                  << "[unmapped]: recording read IDs that couldn't be uniquely mapped\n"
                  << "\tFormat: one line per read\n"
                  << "\treadID  ErrNum  [Position from  To  Qual] [Position from  To  Qual]...\n\n"
                  << "\t--- The following values are optional --- \n\n"
                  << "[rDistrfile]: read mapping distribution file\n"
                  << "\tformat: [pos i, #reads starting at i]\n"
                  << "[PE(orPDD)file]: record # mutations for each kmer position in the dataset\n"
                  << "\tFormat:\n"
                  << "\tpos 0    aa ac ag at ca cc cg ct ga gc gg gt ta tc tg tt\n"
                  << "\t... ...\n"
                  << "\tpos(k-1) aa ac ag at ... ... \n"                  
                  << "----------------------------------------------------\n\n";
        exit(1);
    }

    void setPara() {
        if (argnum < 7) printUsg();
        bedName = arg[1];
        gnName = arg[2];
        srName = arg[3];
        num_reads = atoi(arg[4]);
        //qName = arg[4];
        rsltName = arg[5];
        unmappableName = arg[6];
        if (argnum == 7) return;
        else if (argnum == 8) {
            rDistrName = arg[7];
            return;
        }
        else if (argnum == 10) {
            k_val = atoi(arg[8]);
            peName = arg[9];
            return;
        }
        else printUsg();
    }

};

inline void rv(char& c){
    switch(c){
        case 'A':
        case 'a':
            c = 't';
            break;
        case 'C':
        case 'c':
            c = 'g';
            break;
        case 'G':
        case 'g':
            c = 'c';
            break;
        case 'T':
        case 't':
            c = 'a';
            break;
        default:
            c = 'n';
    }
}

inline void reverse_compl(std::string& ref){

    std::reverse(ref.begin(), ref.end());
    std::for_each(ref.begin(), ref.end(), rv);
};


inline int char_to_bits(char c) {
    int cvalue = -1;
    switch (c) {
        case 'A':
        case 'a':
            cvalue = 0;
            break;
        case 'C':
        case 'c':
            cvalue = 1;
            break;
        case 'G':
        case 'g':
            cvalue = 2;
            break;
        case 'T':
        case 't':
            cvalue = 3;
            break;
        default:
            cvalue = 4;
    }
    return cvalue;
}
#endif	/* _RMAPANALY_H */

