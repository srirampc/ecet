/* 
 *  File:   Record.h
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

#ifndef _RECORD_H
#define	_RECORD_H

#include "compSR.h"

class Record {
public:
    Record();
    Record(const Record& orig);
    virtual ~Record();
    void load(const char* filename);
    void getAmbig(const char* filename);
    void getTrimmed(const char *filename);
    const std::vector<record_t>& getRecords() const{
        return records_;
    }
    const sset_t& getAmbig() const{
        return ambig_;
    }
    const smap_t& getTrimmed() const{
        return trimmed_;
    }
    const smap_t& getTrimPrefix() const{
        return trimPrefix_;
    }
private:
    std::vector<record_t> records_;
    sset_t ambig_;
    smap_t trimmed_;
    smap_t trimPrefix_;
};

#endif	/* _RECORD_H */

