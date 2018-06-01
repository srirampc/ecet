/* 
 * File:   Errlocator.cpp
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

#include <map>
#include <iosfwd>

#include "Errlocator.h"

 

Errlocator::Errlocator(const Errlocator& orig) {
}

Errlocator::~Errlocator() {
}

void printBED(const bed_t& mybed) {
    std::cout << mybed.start << "\t" << mybed.end << "\t" << mybed.readID
            << "\t" << mybed.score << "\t" << mybed.strand << "\n";
}

bool Bedcomp(const bed_t& e1, const bed_t& e2) {

    return (e1.readID < e2.readID);
}

bool unicpy(const bed_t& b1, const bed_t& b2){
    return (b1.readID == b2.readID);
}

//void Errlocator::loadBed(const char* filename) {
void Errlocator::loadBed(const Para& myPara) {

    std::ifstream bedFile(myPara.bedName.c_str()); //RMAP output file
    if (!bedFile.good()) {
        std::cout << "err opening file: " << myPara.bedName.c_str() << "\n";
        exit(1);
    }
    
    iset_t maplist; // storing uniquely mapped IDs

    while (bedFile.good()) {

        bed_t tmpBed;
        std::string ref;
        bedFile >> ref >> tmpBed.start >> tmpBed.end
                >> tmpBed.readID >> tmpBed.score >> tmpBed.strand;
        bed_.push_back(tmpBed);

        maplist.insert(tmpBed.readID);    
    }
    std::sort(bed_.begin(), bed_.end(), Bedcomp);
    bedvec_t::iterator it =
            std::unique_copy(bed_.begin(), bed_.end(), bed_.begin(), unicpy);
    bed_.resize(it - bed_.begin());

    bedFile.close();

    // print non-uniquely mapped read IDs to file, one ID per line
    std::ofstream unmappableFile(myPara.unmappableName.c_str());
    if (!unmappableFile.good()) {
        std::cout << "err opening file: " << myPara.unmappableName.c_str() << "\n";
        exit(1);
    }
    for (int i = 0; i < myPara.num_reads; ++ i){
        if (maplist.count(i) == 0) unmappableFile << i << "\n";
    }
    unmappableFile.close();
}

void Errlocator::rDistr(const char* filename){
    imap_t rDistr;
    for (unsigned int i = 0; i < bed_.size(); ++ i){
        if (bed_[i].strand == '+'){
            imap_t::iterator it (rDistr.find(bed_[i].start));
            if (it != rDistr.end()){
                (*it).second ++;
            }else{
                rDistr[bed_[i].start] = 1;
            }
        }
    }
    // write out
    std::ofstream rDistrHandle (filename);
    if (!rDistrHandle.good()){
        std::cout << "err opening file " << filename << "\n";
        exit(1);
    }
    imap_t::iterator it (rDistr.begin());
    for (; it != rDistr.end(); ++ it){
        rDistrHandle << (*it).first << "\t" << (*it).second << "\n";
    }
    rDistrHandle.close();
}

void Errlocator::loadGn(const std::string& filename) {

    bIO::FASTA_input fasta_sr(filename);
    if (fasta_sr.operator bool() == false) {
        std::cout << "open " << filename << " failed, does it exist? \n";
        exit(1);
    }

    typedef bIO::FASTA_input::value_type value_type;
    while (++fasta_sr) {

        const value_type& v = *fasta_sr;
        gn_ = v.second;
        std::cout << "\n\tNote: assume 1 genome /file; "
                << "otherwise, modify Errlocator::loadGn\n";
        break;
    }
}

void Errlocator::resolveErr(const Para& myPara) {

    bIO::FASTA_input fasta_sr(myPara.srName);
//    bIO::FASTA_input fasta_q(myPara.qName);
    typedef bIO::FASTA_input::value_type value_type;

    if (fasta_sr.operator bool() == false) {
        std::cout << "open " << myPara.srName << " failed, does it exist? \n";
        exit(1);
    }
//    if (fasta_q.operator bool() == false) {
//        std::cout << "open " << myPara.qName << " failed, does it exist? \n";
//        exit(1);
//    }

    // error recording file
    std::ofstream oHandle(myPara.rsltName.c_str());
    if (!oHandle.good()) {
        std::cout << "open " << myPara.rsltName << " failed, correct path?\n";
        exit(1);
    }

    int maxbuf = 100000;
    int counter = 0;
    uint64_t sum_nt = 0;

    std::cout << "\tprogress: -- bed size = " << bed_.size() << "\t";
    for (int i = 0; i < bed_.size(); ++i) {

        if (counter % 1000000 == 0) std::cout << (double) counter / 1000000 << "M\t";

        bed_t myBed = bed_[i];

        int targetID = myBed.readID;

//      while ((++fasta_sr == true) && (++fasta_q == true)) {
        while (++fasta_sr == true) {

            const value_type& v_sr = *fasta_sr;
//            const value_type& v_q = *fasta_q;
            sum_nt += v_sr.second.length();
            if (counter == targetID) {
                std::string ref = gn_.substr(myBed.start, myBed.end - myBed.start);
                std::string read = v_sr.second;
//              std::string qual = v_q.second;                 
//              qual.erase(0, 1);
                std::string qual=""; // currently, not reading qual string
                
                if (myBed.strand == '-') reverse_compl(ref);
                //std::cout << ref << "\n";
                evec_t rslt;
                tracking(rslt, ref, read, qual);
                if (myPara.k_val) {
                    cal_PE(ref, read, myPara.k_val);
                }
                if (rslt.size()) {
                    rsltBuf_.push_back(rslt_t(targetID, rslt));
                }
                counter++;
                break;
            } else {
                if (counter > targetID) {
                    std::cout << "err: duplicated targetID found:" << targetID << "\n";
                }
            }

            counter++;
        }

        if (rsltBuf_.size() >= maxbuf) flush(oHandle);
    }
    flush(oHandle);
    oHandle.close();
  
    // write out PE file : error frequencies for each kmer position
    if (!myPara.peName.empty()) {
        std::ofstream efHandle(myPara.peName.c_str());
        if (!efHandle.good()) {
            std::cout << "open " << myPara.peName << " failed, correct path?\n";
            exit(1);
        }

        for (int i = 0; i < PE_.size(); ++ i){
            efHandle << i << "\t";
            for (int j = 0; j < PE_[i].size(); ++ j){
                efHandle << PE_[i][j] << "\t";
            }
            efHandle << std::endl;
        }
        efHandle << "pos\taa\tac\tag\tat\tca\tcc\t\tcg\tct\tga\tgc\tgg\t\tgt\tta\ttc\ttg\ttt\n";
        efHandle.close();
    }
    std::cout << "\n\nTotal Number of Errs: " << total_err_number_ << std::endl;
    std::cout << "Total Nucleotides " << sum_nt << "\n";
    std::cout << "Err rate = " << (float) total_err_number_/sum_nt << "\n";
    
    
}

void Errlocator::tracking(evec_t& rslt, const std::string& ref,
        const std::string& read, const std::string& qual) {

    // keep track of errors per read output to rslt
    for (int i = 0; i < ref.length(); ++i) {
        if (tolower(ref.at(i)) != tolower(read.at(i))) {
           // rslt.push_back(e_t(i, char_to_bits(ref.at(i)),
           //      char_to_bits(read.at(i)), (int) qual.at(i)));
           rslt.push_back(e_t(i, char_to_bits(ref.at(i)),
                    char_to_bits(read.at(i)), 0));

           total_err_number_ ++; // adds up all errors
        }
    }
}

void Errlocator::cal_PE(const std::string& ref, const std::string& read, int k_val){

    // keep track of errors per kmer
    for (int i = 0; i < ref.length() - k_val + 1; ++ i){
        for (int j = 0; j < k_val; ++ j){
            if (PE_.size() <= j){
                PE_.push_back(u64vec_t(16,0));
            }
            int cvalue = char_to_bits(tolower(read.at(j + i)));
            if (cvalue == 4) // non-acgt characters
                continue;
            
            int idx = (char_to_bits(tolower(ref.at(j + i))) << 2) + cvalue;
            if (idx < 0 || idx > 15){
                std::cout << "err: idx is out of range\n";
                exit(1);
            }
            PE_[j][idx] ++;
        }
    }
}

void Errlocator::flush(std::ofstream& handle) {

    /*
     * Resulting file format:
     * ReadID (sorted)  ErrNum  [pos from to qual] [pos from to qual]...
     * from: reference; to: read (numercial value Aa:0 Cc:1 Gg:2 Tt:3 others 4)
     * qual: quality (numerical value)
     */
    for (int i = 0; i < rsltBuf_.size(); ++i) {

        handle << rsltBuf_[i].readID << "\t" << rsltBuf_[i].evec.size();

        for (int j = 0; j < rsltBuf_[i].evec.size(); ++j) {
            handle << "\t" << rsltBuf_[i].evec[j].pos << "\t"
                    << rsltBuf_[i].evec[j].from << "\t" << rsltBuf_[i].evec[j].to
                    << "\t" << rsltBuf_[i].evec[j].qual;
        }
        handle << "\n";
    }
    rsltBuf_.clear();
}
