/***
 *  File:   Errlocator.h
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

#ifndef _ERRLOCATOR_H
#define	_ERRLOCATOR_H

 
#include "rmapAnaly.h"



class Errlocator {
public:
    Errlocator() {
        total_err_number_ = 0;
    };
    Errlocator(const Errlocator& orig);
    virtual ~Errlocator();
    //void loadBed(const char* filename);
    void loadBed(const Para& myPara);
    void rDistr(const char* filename);
    void loadGn(const std::string&);
    void resolveErr(const Para& myPara);
private:
    uint64_t total_err_number_;
    bedvec_t bed_;
    std::string gn_;
    std::vector<rslt_t> rsltBuf_;
    uu64vec_t PE_;   // transition matrix quantifying error probs for each kmer position
    void tracking (evec_t& rslt, const std::string& ref,
            const std::string& read, const std::string& qual);
    void cal_PE (const std::string& ref, const std::string& read, int k_val);
    void flush(std::ofstream& handle);
};

#endif	/* _ERRLOCATOR_H */

