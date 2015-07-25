/*****************************************************************************
 *   Copyright (C) 2004-2015 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://github.com/esa/pagmo                                            *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#ifndef PAGMO_UTIL_BBOB_H
#define PAGMO_UTIL_BBOB_H

#include <iostream>
#include <string>
#include <stdio.h>
#include <vector>
#include <fstream>
#include<sstream>
#include <limits>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

#include "../config.h"
#include "../serialization.h"
#include "../population.h"
#include "../problem/base.h"
#include "../problem/base_meta.h"

namespace fs = boost::filesystem;

namespace pagmo { namespace util {

/**
 * This class creates the meta problem for Black-Box Optimization Benchmarking (BBOB2015) which
 * has to be optimized.
 */
class __PAGMO_VISIBLE bbob : public problem::base_meta
{
    public:
        bbob(const problem::base & = problem::ackley(1), const std::string = "./", const std::string = "", const unsigned int = 1, const std::string = "");

        bbob(const bbob &);

        //finalize benchmarking
        void finalize(population &) const;

        problem::base_ptr clone() const;
        std::string get_name() const;

    protected:
        void objfun_impl(fitness_vector &, const decision_vector &) const;

    private:
        struct lastEvalStruct
        {
            double num; ///Number of evaluations
            double F; //fitness value
            double bestF; //best fitness value till now
            decision_vector x; //decision vector
            int isWritten; //is it written to file?

            friend class boost::serialization::access;
            template <typename Archive>
            void serialize(Archive& ar, const unsigned int version)
            {
              ar & num;
              ar & F;
              ar & bestF;
              ar & x;
              ar & isWritten;
            }
        };

        typedef struct lastEvalStruct LastEvalStruct;

        struct filedata
        {
            double num;
            double F;
            double bestF;
            decision_vector x;

            friend class boost::serialization::access;
            template <typename Archive>
            void serialize(Archive& ar, const unsigned int version)
            {
              ar & num;
              ar & F;
              ar & bestF;
              ar & x;
            }
        };

        typedef struct filedata data;

        //vectors to store data that will be written to files.
        mutable std::vector<data> m_dataFile;
        mutable std::vector<data> m_rdataFile;
        mutable std::vector<data> m_hdataFile;

        //The default values for the LastEvalStruct structures
        LastEvalStruct m_lastEvalInit = {0, std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(), decision_vector(), 0};

        //Helper functions
        void writeNewIndexEntry(void) const; //Create a new index file and write an index entry

        void addIndexEntry(void) const; ////Add index entry to existing index file when the data file is same
        void writeFinalData(void) const; //Write Final data to log files.
        void storeBestF(std::vector<data> &, LastEvalStruct) const; //store best evaluation data to data vector at correct position.
        void storeData(std::vector<data> &, double, double, double, decision_vector) const; //Store data to be written to files in data vectors.
        void writeDataHeader(fs::path) const; //Write headers to files.
        void bbobOpenFile(std::ofstream &, fs::path) const; //Open file and return fileId
        void writeData(fs::path, std::vector<data> &) const; //Write data from vectors to data files.
        void restart(std::string) const; //write to restart log.

        mutable double m_bestF; //Best fitness

        mutable double m_fTrigger;
        mutable double m_evalsTrigger;
        mutable unsigned int m_idxEvalsTrigger, m_idxDIMEvalsTrigger;
        mutable int m_idxFTrigger;

        const double m_maxFunEvalsFactor = 1e6;
        const unsigned int m_nbPtsEvals = 20;
        const unsigned int m_nbPtsF = 5;

        mutable LastEvalStruct m_LastEval, m_BestFEval;
        mutable double m_lastWriteEval;

        std::string m_algName; //Name of the optimiser (used in the post-processing)
        std::string m_comments; //Additional comments
        double m_precision = 1e-8; //The precision
        unsigned int m_instanceId = 1; //The instance ID
        mutable fs::path m_dataPath; //The path of the index and data files.
        mutable fs::path m_dirPath; //The directory of the index and data files.
        mutable fs::path m_indexFilePath; //Index file name will be  ('f_%s.info', problem.get_name())
        mutable fs::path m_dataFilePath; //name of data file
        mutable fs::path m_hdataFilePath; //name of H-data file
        mutable fs::path m_rdataFilePath; //names of r-data file used for restart
        mutable unsigned int m_runCounter; //Counter for number of evaluations

        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int)
        {
            ar & boost::serialization::base_object<base_meta>(*this);
            ar & m_dataFile;
            ar & m_rdataFile;
            ar & m_hdataFile;
            ar & m_lastEvalInit;
            ar & m_bestF;
            ar & m_fTrigger;
            ar & const_cast<double &>(m_evalsTrigger);
            ar & const_cast<unsigned int &>(m_idxEvalsTrigger);
            ar & const_cast<unsigned int &>(m_idxDIMEvalsTrigger);
            ar & m_idxFTrigger;
            ar & const_cast<double &>(m_maxFunEvalsFactor);
            ar & const_cast<unsigned int&>(m_nbPtsEvals);
            ar & const_cast<unsigned int &>(m_nbPtsF);
            ar & m_LastEval;
            ar & m_BestFEval;
            ar & m_lastWriteEval;
            ar & m_algName;
            ar & m_comments;
            ar & const_cast<double &>(m_precision);
            ar & const_cast<unsigned int &>(m_instanceId);
            ar & m_dataPath;
            ar & m_dirPath;
            ar & m_indexFilePath;
            ar & m_dataFilePath;
            ar & m_hdataFilePath;
            ar & m_rdataFilePath;
            ar & m_runCounter;
        }
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::util::bbob)

#endif
