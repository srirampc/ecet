/* 
 *  File:   main.cpp
 *  Created on December 19, 2009, 2:51 PM
 *  Author: Xiao Yang <isuxyang@gmail.com>
 *
 *  Copyright 2009 Xiao Yang
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
#include <limits.h>

/*
 *  Compare the result of ANY error correction method to mapping result
 *  by Pre-correction Alignment
 *
 *  Input file Format:
 *
 *  ReadID(sorted)  Err#    [pos    from    to  qual-score] [... ...]
 *
 *  e.g.  0           1        32     0(a)   1(c)     41
 *
 *  [from] (true bp) and [to] (error bp) are numerical values of acgt
 */
int main(int argc, char** argv) {

    Para myPara(argc, argv);
    Record ECRecord, pcAlignRecord;
    ECRecord.load(myPara.ecName.c_str());
    pcAlignRecord.load(myPara.pcAlignName.c_str());
    pcAlignRecord.getAmbig(myPara.ambigName.c_str());

    if(myPara.trimFileName.length() > 0)
        pcAlignRecord.getTrimmed(myPara.trimFileName.c_str());

    comparison(myPara.rsltName.c_str(), pcAlignRecord.getRecords(),
               ECRecord.getRecords(), pcAlignRecord.getAmbig(), 
               pcAlignRecord.getTrimmed(), pcAlignRecord.getTrimPrefix(),
               myPara.mvalue);

    return (EXIT_SUCCESS);
}

/*
 * FP: exist in tar but not ref
 * FN: exist in ref but not tar
 * TP: exist in both tar and ref
 * EBA(WrongBase): erroneous base assignment (exist in neither tar nor ref)
 */
void print_rslt(const fpnp_t& rslt){
    std::cout << "FN: " << rslt.FN << "  FP:" << rslt.FP << "  TP:" << rslt.TP << "\n";
}

int getTrimmedLength(const imap_t& trimmed, int readID){
    imap_t::const_iterator cit = trimmed.find(readID);
    if(cit != trimmed.end()) {
        int tLen = (*cit).second;
        return tLen;
    }
    return INT_MAX;
}

int getTrimPrefix(const imap_t& trimPrefix, int readID){
    imap_t::const_iterator cit = trimPrefix.find(readID);
    if(cit != trimPrefix.end()) {
        int tLen = (*cit).second;
        return tLen;
    }
    return -1;
}

int getTrimErrors(const evec_t &errors,int maxpos,int lbound){
    int errs = 0;
    for(unsigned int i = 0; i < errors.size();i++) {
        if( errors[i].pos < lbound ) // skip th errors
            continue;
        if( errors[i].pos <= maxpos )
            errs++;
    }
    return errs;
}

void comparison(const char* filename, const std::vector<record_t>& ref,
                const std::vector<record_t>& tar, const iset_t& ambig,
                const imap_t& trimmed, const imap_t& trimPrefix, int mvalue) {

    int idx_r = 0, idx_t = 0;

    fpnp_t rslt;
    
    int refTotal = 0; // total number of errors in the ref
    
    iset_t ref_readIDs;
    for (unsigned int i = 0; i < ref.size(); ++ i){
        int tLen = getTrimmedLength(trimmed,ref[i].readID) - 1;
        int tPrefix = getTrimmedLength(trimmed,ref[i].readID);
        int trimErrors = getTrimErrors(ref[i].evec,tLen,tPrefix);
        refTotal += trimErrors; // ref[i].evec.size();
        ref_readIDs.insert(ref[i].readID);
    }

    std::cout << "-----------------------------------------------------------\n";
    std::cout << "Total Errors in Reference = " << refTotal << "\n";


    //debug
/*    std::ofstream testHandle("trackFPs.txt");
    if (!testHandle.good()) {
        std::cout << "err opening file: " << "trackFPs.txt" << "\n";
        exit(1);
    }
*/
#ifdef DEBUG
    std::cout << "REF : " << ref.size() << " CORR : " << tar.size() << std::endl;
#endif
    while (idx_r < (int) ref.size() && idx_t < (int) tar.size()) {

        if (ref[idx_r].readID == tar[idx_t].readID) {
            checkRead(rslt,ref[idx_r].readID, trimmed, trimPrefix,ref[idx_r].evec, tar[idx_t].evec);
            ++idx_r;
            ++idx_t;
        }
        if ((idx_r > (int) ref.size() - 1) || (idx_t > (int) tar.size() - 1)) {
            break;
        }
        while (ref[idx_r].readID < tar[idx_t].readID && idx_r < (int) ref.size()) {
            // Count only the errors until tLen
            int tLen = getTrimmedLength(trimmed,ref[idx_r].readID) - 1;
            int tPrefix = getTrimPrefix(trimPrefix,ref[idx_r].readID);
            int trimErrors = getTrimErrors(ref[idx_r].evec,tLen,tPrefix);
            rslt.FN +=  trimErrors; // ref[idx_r].evec.size();
            ++idx_r;
        }
        if (idx_r > ref.size() - 1)
            break;
        while (ref[idx_r].readID > tar[idx_t].readID && idx_t < (int) tar.size()) {
            
            if(!ambig.count(tar[idx_t].readID)){
                rslt.FP += tar[idx_t].evec.size();
//              testHandle << tar[idx_t].readID << "\t" << tar[idx_t].evec.size() << "\n";
            }
            else {
                rslt.cntg_m += tar[idx_t].evec.size();
            }
            ++idx_t;
        }
        if (idx_t > (int) tar.size() - 1) break;

    }

    if (idx_r > ref.size() - 1) {
        for (unsigned int i = idx_t; i < tar.size(); ++i) {
            if (!ambig.count(tar[i].readID)){
                rslt.FP += tar[i].evec.size();
            }
            else {
                rslt.cntg_m += tar[i].evec.size();
            }
        }
    }
    if (idx_t > tar.size() - 1) {
        for (unsigned int i = idx_r; i < ref.size(); ++i) {
            // Count only the errors until tLen
            int tLen = getTrimmedLength(trimmed,ref[i].readID) - 1;
            int tPref = getTrimPrefix(trimPrefix,ref[i].readID);
            int trimErrors = getTrimErrors(ref[i].evec,tLen,tPref);
            rslt.FN +=  trimErrors; //ref[i].evec.size();
        }
    }
    output(rslt, mvalue, filename);

   std::cout << "-----------------------------------------------------------\n\n";
}

void print_evec(const evec_t& myvec){
    std::cout << "------------\n";
    for (unsigned int i = 0; i < myvec.size(); ++ i){
        std::cout << myvec[i].pos << "(" << myvec[i].from << "," << myvec[i].to << ")\t";
    }
    std::cout << "\n";
}


void checkRead(fpnp_t& rslt, int readID, 
               const imap_t& trimmed, const imap_t& trimPrefix,
               const evec_t& ref, const evec_t& tar) {

#ifdef DEBUG
    std::cout << "size = " << ref.size() << "\t" << tar.size() << "\n";
#endif
#ifdef DEBUG
    print_evec(ref);
    print_evec(tar);
#endif
    
    int idx_r = 0, idx_t = 0;
    int tLen = getTrimmedLength(trimmed,readID) - 1;
    int tPref = getTrimPrefix(trimPrefix,readID);
    bool trimBreak = false;

    while (idx_r < (int) ref.size() && idx_t < (int) tar.size() ) {

        if (ref[idx_r].pos == tar[idx_t].pos) {
            // skip the comparison that go beyond trimLength
            if(ref[idx_r].pos >= tLen) { trimBreak = true; break; }

            if (ref[idx_r].from == tar[idx_t].from
                    && ref[idx_r].to == tar[idx_t].to) {
                ++rslt.TP;
            } else { //corrected the right position but to wrong bases
                
//                ++rslt.FN;
//                ++rslt.FP;
//                ++rslt.WrongBp;

                if (ref[idx_r].to == 4){
                    if (ref[idx_r].from != tar[idx_t].from){
                        ++rslt.WrongBp_N;

                        ++rslt.FN;
                        ++rslt.FP;
                        ++rslt.WrongBp;
                    }
                    ++rslt.N_total;
                }
                else{
                    ++rslt.FN;
                    ++rslt.FP;
                    ++rslt.WrongBp;
                }
            }
            ++idx_r;
            ++idx_t;
        }

        if ((idx_r > (int) ref.size() - 1) || (idx_t > (int) tar.size() - 1))
            break;

        while (ref[idx_r].pos < tar[idx_t].pos && idx_r < (int) ref.size()) {
            // skip the comparison that go beyond trimLength
            if(ref[idx_r].pos >= tLen) { trimBreak = true; break; }
            rslt.FN++;
            ++idx_r;
        }
        if(trimBreak) break;
        if (idx_r > (int) ref.size() - 1)
            break;

        while (ref[idx_r].pos > tar[idx_t].pos && idx_t < (int) tar.size()) {
            // skip the comparison that go beyond trimLength
            if(tar[idx_t].pos >= tLen) { trimBreak = true; break; }
            if(tar[idx_t].pos > tPref) 
                rslt.FP++;
            ++idx_t;
        }
        if(trimBreak) break;
        if (idx_t > (int) tar.size() - 1)
            break;
    }

    if (idx_r > (int) ref.size() - 1 && idx_t <= (int) tar.size() - 1) {
        rslt.FP += ((int) tar.size() - idx_t);
    }
    if (idx_t > (int) tar.size() - 1 && idx_r <= (int) ref.size() - 1) {
        rslt.FN += ((int) ref.size() - idx_r);
    }
}

void output(const fpnp_t& rslt, int mvalue, const char* filename) {
    std::ofstream oHandle(filename);
    if (!oHandle.good()) {
        std::cout << "err opening file: " << filename << "\n";
        exit(1);
    }
    //int FP = rslt.FP - rslt.cntg_m;
    double sensitivity = 1.0*rslt.TP/(rslt.TP+rslt.FN);
    // double gain = 1.0*(rslt.TP - FP)/(rslt.TP + rslt.FN);
    double gain = 1.0*(rslt.TP - rslt.FP)/(rslt.TP + rslt.FN);
    oHandle << "Notations: \n";
    oHandle << "---------------------------\n";
    oHandle << "FP: exist in tar but not ref\n";
    oHandle << "FN: exist in ref but not tar\n";
    oHandle << "TP: exist in both tar and ref\n";
    oHandle << "---------------------------\n";

    oHandle << "Total Err (TP, FN):\t" << rslt.TP + rslt.FN << "\n\n";
    oHandle << "TP\t" << rslt.TP << "\n\n";
    //oHandle << "FP\t" << FP << "\n\n";
    oHandle << "FP\t" << rslt.FP << "\n\n";
    oHandle << "FN\t" << rslt.FN << "\n\n";
    oHandle << "EBA\t" << rslt.WrongBp << "\n\n";
    oHandle << "Sensitivity = " << sensitivity << "\n\n";
    oHandle << "Gain = " << gain << "\n\n";
    oHandle << "Total Errs Corrected in tar reads that cannot be uniquely mapped"
            << "by pre-correction alignment (-m = "
            << mvalue  << ") : " << rslt.cntg_m << "\n\n";
    //double rate = double 1.0*rslt.WrongBp_N
    oHandle << "Approximate ambiguous correction false rate: "
            << rslt.WrongBp_N << " out of " << rslt.N_total << " ( " 
            << (double) rslt.WrongBp_N*100.0/rslt.N_total << " % ) \n\n" ;
    
    std::cout << "TP\tFP\tFN\tsensitivity\tGain\tEBA\n";
    //std::cout << rslt.TP << "\t" << FP << "\t" << rslt.FN << "\t"
    //          << sensitivity << "\t" << gain << "\t" << rslt.WrongBp << "\n\n";
    std::cout << rslt.TP << "\t" << rslt.FP << "\t" << rslt.FN << "\t"
              << sensitivity << "\t" << gain << "\t" << rslt.WrongBp << "\n\n";

    std::cout << "Total Errs in tar reads (that are NOT recorded as ref reads) \n"
              << "containing more than " << mvalue  << " corrected errs: " << rslt.cntg_m << "\n\n";
    std::cout << "Approximate ambiguous correction false rate: "
            << rslt.WrongBp_N << " out of " << rslt.N_total << " ( "
            << (double) rslt.WrongBp_N*100.0/rslt.N_total << " % ) \n\n" ;
    oHandle.close();
}
