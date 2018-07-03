/* 
 * File:   Record.cpp
 *  Author: Xiao Yang <isuxyang@gmail.com>
 *  Created on December 19, 2009, 3:04 PM
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

#include "Record.h"
#include <sstream>
#include <string>
#include <zlib.h>
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

Record::Record() {
}

Record::Record(const Record& orig) {
}

Record::~Record() {
}

bool Rcomp(const record_t& r1, const record_t& r2) {

    return (r1.readID < r2.readID);
}

bool Ecomp(const e_t& e1, const e_t& e2) {
    return (e1.pos < e2.pos);
}

bool unicpy(const e_t& e1, const e_t& e2) {
    return (e1.pos == e2.pos);
}

void Record::load(const char* filename) {
    gzFile fp = gzopen(filename, "r");

    //std::ifstream input(filename); //pre-correction alignment file
    //if (!input.good()) {
    if (fp == Z_NULL) {
        std::cout << "err opening file: " << filename << "\n";
        exit(1);
    }
    kstream_t *ks;
    kstring_t str = {0,0,0};
    ks = ks_init(fp);
    while (ks_getuntil(ks, '\n', &str, 0) >= 0) {
        std::string line(str.s);
        //getline(input,line);
        line.erase(line.find_last_not_of(" \n\r\t")+1);
        if(line.length() == 0)
            continue;
        std::istringstream iss (line,std::istringstream::in);

        evec_t tmpEvec;
        int readID, errNum, pos, from, to, qual;
        iss >> readID >> errNum;
        for (int i = 0; i < errNum; ++i) {
            iss >> pos >> from >> to >> qual;
            tmpEvec.push_back(e_t(pos, from, to, qual));
        }
        //sort and unique_copy tmpEvec
        std::sort(tmpEvec.begin(), tmpEvec.end(), Ecomp);
        // evec_t::iterator it = std::unique_copy(tmpEvec.begin(), tmpEvec.end(),
        //         tmpEvec.begin(), unicpy);
        // tmpEvec.resize(it - tmpEvec.begin());

        records_.push_back(record_t(readID, tmpEvec));
#ifdef DEBUG
        std::cout << "RECORDS : " << readID << std::endl;
#endif
    }
    std::sort(records_.begin(), records_.end(), Rcomp);
    //input.close();
    ks_destroy(ks);
    gzclose(fp);
    free(str.s);
}

void Record::getAmbig(const char* filename) {

    gzFile fp = gzopen(filename, "r");

    //std::ifstream input(filename); //pre-correction alignment file
    //if (!input.good()) {
    if (fp == Z_NULL) {
        std::cout << "err opening file: " << filename << "\n";
        exit(1);
    }
    kstream_t *ks;
    kstring_t str = {0,0,0};
    ks = ks_init(fp);
    while (ks_getuntil(ks, '\n', &str, 0) >= 0) {
        //getline(input,line);
        int readID;
        //input >> readID;
        std::atoi(str.s);
        ambig_.insert(readID);
    }

    //input.close();
    ks_destroy(ks);
    gzclose(fp);
    free(str.s);

}

void Record::getTrimmed(const char *filename){
    gzFile fp = gzopen(filename, "r");

    //std::ifstream input(filename); //pre-correction alignment file
    //if (!input.good()) {
    if (fp == Z_NULL) {
        std::cout << "err opening file: " << filename << "\n";
        exit(1);
    }
    kstream_t *ks;
    kstring_t str = {0,0,0};
    ks = ks_init(fp);
    while (ks_getuntil(ks, '\n', &str, 0) >= 0) {
        std::string line(str.s);
        //getline(input,line);
        line.erase(line.find_last_not_of(" \n\r\t")+1);
        if(line.length() == 0)
            continue;
        std::istringstream iss (line,std::istringstream::in);
        int readID = -1, trimmedLength = -1,trimPrefix;
        
        iss >> readID >> trimmedLength >> trimPrefix;
        if(readID == -1)
            continue;
        trimmed_.insert(std::pair<int,int>(readID,trimmedLength));
        trimPrefix_.insert(std::pair<int,int>(readID,trimPrefix));
    }
    //input.close();
    ks_destroy(ks);
    gzclose(fp);
    free(str.s);
}
