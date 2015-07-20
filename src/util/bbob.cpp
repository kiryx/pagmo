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

#include "bbob.h"

namespace pagmo { namespace util {

bbob::bbob(const pagmo::problem::base & p, const std::string datapath, const std::string algname, unsigned int instanceId) : problem::base_meta(p,
                                p.get_dimension(), // Ambiguous without the cast ...
                                p.get_i_dimension(),
                                1, //We can only benchmark using single objective functions.
                                p.get_c_dimension(),
                                p.get_ic_dimension(),
                                p.get_c_tol()), m_algName(algname), m_instanceId(instanceId)
{
    std::vector<fitness_vector> bestf = p.get_best_f();

    if(p.get_f_dimension() != 1)
        pagmo_throw(value_error, "Only single-objective problems can be benchmarked.");

    if(bestf.size() == 0)
        pagmo_throw(value_error, "Best known decision vector(s) not found. Best fitness value needs to be known for benchmarking.");

    m_bestF = bestf[0][0];

    m_LastEval = m_lastEvalInit; //default values
    m_BestFEval = m_LastEval;
    m_BestFEval.isWritten = 1;

    m_lastWriteEval = 0;

    m_idxFTrigger = std::numeric_limits<int>::max();
    m_fTrigger = std::numeric_limits<double>::max();   //because 10^DBL_MAX will generate an error
    m_idxEvalsTrigger = 0;
    m_evalsTrigger = floor(pow(10.,(double)m_idxEvalsTrigger/(double)m_nbPtsEvals)); // = 1 !
    m_idxDIMEvalsTrigger = 0;

    m_runCounter = 1;

    m_dataPath = fs::path(datapath);

    if(!fs::is_directory(m_dataPath))
	{
	    if(!fs::create_directory(m_dataPath)) //try to create a new directory
	        pagmo_throw(std::runtime_error, "Invalid datapath given. Creating new directory failed");
	}

    m_dirPath = fs::path(p.get_name()); //create a directory whose name is the function name.

    //index and data files are stored in <datapath>/<algorithm_name>/<function_name>/func_dim[get_dimension()].dat or .tdat or .rdat

    if(!fs::is_directory(m_dataPath / m_dirPath)) //check if directory with the function name exists.
        if(!fs::create_directory(m_dataPath / m_dirPath)) //try to create a new directory
            pagmo_throw(std::runtime_error, "New directory could not be created, does given datapath has write permissions?");

    fs::path index_filename(boost::str(boost::format("f_%1%.info") % p.get_name()));
    m_indexFilePath = m_dataPath / index_filename;

    fs::path data_filename(boost::str(boost::format("func_dim[%1%].tdat") % p.get_dimension()));
    m_dataFilePath = m_dataPath / m_dirPath / data_filename;

    fs::path hdata_filename(boost::str(boost::format("func_dim[%1%].dat") % p.get_dimension()));
    m_hdataFilePath = m_dataPath / m_dirPath / hdata_filename;

    fs::path rdata_filename(boost::str(boost::format("func_dim[%1%].rdat") % p.get_dimension()));
    m_rdataFilePath = m_dataPath / m_dirPath / rdata_filename;


    if(!fs::is_regular_file(m_indexFilePath))
        writeNewIndexEntry(); //create new index file and add new index entry

    else if(!fs::is_regular_file(m_dataFilePath)) //check if we already have a datafile with same parameters
        addDatIndexEntry(); //add information about the datafile to existing index file

    else
        addIndexEntry();

    writeDataHeader(m_dataFilePath);
    writeDataHeader(m_hdataFilePath);
    writeDataHeader(m_rdataFilePath);
}

void bbob::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
    int i;
    unsigned int boolImprovement = 0;
    double Fvalue;
    double evalsj;
    FILE * dataFileId;
    FILE * hdataFileId;

    m_original_problem->objfun(f, x);
    Fvalue = f[0];

    //should we print ? 2 possible conditions, # evals or fitness value
    if((m_LastEval.num+1 >= m_evalsTrigger) || (Fvalue - m_bestF < m_fTrigger))
    {
        evalsj = m_LastEval.num + 1;

        if(Fvalue < m_BestFEval.F) //minimization
        {
            boolImprovement = 1;
            m_BestFEval.F = Fvalue;
            m_BestFEval.isWritten = 0;
        }

        //should we print something? First based on # evals
        if(evalsj >= m_evalsTrigger)
        {
            m_lastWriteEval = evalsj;
            m_BestFEval.isWritten = 1;

            storeData(m_dataFile, evalsj, Fvalue, m_BestFEval.F, x);

            //update of next print triggers based on # evals */
            while (evalsj >= floor(pow(10., (double)m_idxEvalsTrigger/(double)m_nbPtsEvals)))
                m_idxEvalsTrigger = m_idxEvalsTrigger + 1;

            while (evalsj >= x.size() * pow(10., (double)m_idxDIMEvalsTrigger))
                m_idxDIMEvalsTrigger = m_idxDIMEvalsTrigger + 1;

            m_evalsTrigger = fmin(floor(pow(10., (double)m_idxEvalsTrigger/(double)m_nbPtsEvals)), x.size() * pow(10., (double) m_idxDIMEvalsTrigger));
        }

        //now based on fitness values
        if(Fvalue - m_bestF < m_fTrigger)
        {
            storeData(m_hdataFile, evalsj, Fvalue, m_BestFEval.F, x);

            if(Fvalue-m_bestF <= 0)
                m_fTrigger = -1 * std::numeric_limits<double>::max();
            else
            {
                if(m_idxFTrigger == std::numeric_limits<int>::max())
                    m_idxFTrigger = ceil(log10(Fvalue-m_bestF))*m_nbPtsF;

                while((Fvalue-m_bestF) <= pow(10., (double)m_idxFTrigger/(double)m_nbPtsF))
                    m_idxFTrigger = m_idxFTrigger - 1;

                m_fTrigger = fmin(m_fTrigger, pow(10., (double)m_idxFTrigger/(double)m_nbPtsF));
            }
        }

        if(! m_BestFEval.isWritten && boolImprovement)
        {
            m_BestFEval.num = m_LastEval.num+1;
            m_BestFEval.F = Fvalue;
            m_BestFEval.x = x;
        }
    }
    else
    {
        if(Fvalue < m_BestFEval.F)
        {
            m_BestFEval.num = m_LastEval.num+1;
            m_BestFEval.F = Fvalue;
            m_BestFEval.x = x;
            m_BestFEval.isWritten = 0;
        }
    }

    m_LastEval.num = m_LastEval.num + 1;
    m_LastEval.F = Fvalue;
    m_LastEval.x = x;
}

//Subfunctions

//write the comment line header in the data files
void bbob::writeDataHeader(fs::path dataFilePath) const
{
    FILE * dataFileId;
    bbobOpenFile(dataFileId, dataFilePath);

    fprintf(dataFileId, "%% function evaluation | fitness - Fopt (%13.12e) | best fitness - Fopt | measured fitness | best measured fitness | x1 | x2...\n", m_bestF);
    fclose(dataFileId);
}

//open the index file and write a new index entry
void bbob::writeNewIndexEntry(void) const
{
    FILE * indexFileId;
    int newline = 1;
    if(!fs::is_regular_file(m_indexFilePath))
        newline = 0;

    indexFileId = fopen(m_indexFilePath.c_str(), "a");

    if(indexFileId == NULL)
        pagmo_throw(value_error, "Could not open index file.");
    if(newline == 1)
        fprintf(indexFileId,"\n");

    fprintf(indexFileId, "funcId = '%s', DIM = %lu, Precision = %.3e, Fopt = %13.12e, algId = '%s'\n", m_original_problem->get_name().c_str(),
		get_dimension(), m_precision, m_bestF, m_algName.c_str());
    fprintf(indexFileId,"%% %s\n%s, %d", m_comments.c_str(), (m_dirPath / m_hdataFilePath.filename()).c_str(), m_instanceId);
    fclose(indexFileId);
}

//Open the index file and write a new index entry.
void bbob::addDatIndexEntry(void) const
{
    FILE * indexFileId;

    indexFileId = fopen(m_indexFilePath.c_str(), "a");
    if(indexFileId == NULL)
        pagmo_throw(std::runtime_error, "Could not open index file.");

    fprintf(indexFileId,", %s, %d", (m_dirPath / m_hdataFilePath.filename()).c_str(), m_instanceId);
    fclose(indexFileId);
}

//Open the index file and write a new index entry */
void bbob::addIndexEntry(void) const
{
    FILE * indexFileId;

    if (!fs::is_regular_file(m_indexFilePath))
        pagmo_throw(std::runtime_error, "Could not find index file");

    indexFileId = fopen(m_indexFilePath.c_str(), "a");
    if (indexFileId == NULL)
        pagmo_throw(std::runtime_error, "Could not find index file");

    fprintf(indexFileId,", %d", m_instanceId);
    fclose(indexFileId);
}

//complete the data file with unwritten information
void bbob::writeFinalData(void) const
{
    FILE * indexFileId;

    if(!m_BestFEval.isWritten)
    {
        if(m_BestFEval.num > m_lastWriteEval)
        {
            m_lastWriteEval = m_BestFEval.num;
            storeData(m_dataFile, m_BestFEval.num, m_BestFEval.F, m_BestFEval.F, m_LastEval.x);
        }
        else
        {   //here, need to write best at correct position
            storeBestF(m_dataFile, m_BestFEval);
        }
    }
    if(m_LastEval.num > m_lastWriteEval)
        storeData(m_dataFile, m_LastEval.num, m_LastEval.F, m_BestFEval.F, m_LastEval.x);

    //now the index file
    bbobOpenFile(indexFileId, m_indexFilePath);
    fprintf(indexFileId, ":%.0f|%.1e", m_LastEval.num, m_BestFEval.F - m_bestF - m_precision);
    fclose(indexFileId);

    //now write the data from data vectors to files.
    writeData(m_dataFilePath, m_dataFile);
    writeData(m_rdataFilePath, m_rdataFile);
    writeData(m_hdataFilePath, m_hdataFile);
}

//rewrite the data file with the information about the best F.
void bbob::storeBestF(std::vector<data> &dataFile, LastEvalStruct BestFEval) const
{
    //insert the evaluation at the correct position
    std::vector<data> temp;
    data best;

    //copy data into best.
    best.num = BestFEval.num;
    best.F = BestFEval.F;
    best.bestF = BestFEval.F;
    best.x = BestFEval.x;

    for(int j=0; j<dataFile.size(); j++)
    {
        if(dataFile[j].num > BestFEval.num)
            temp.push_back(best);

        temp.push_back(dataFile[j]);
    }

    dataFile = temp;
}

//write a formatted line into a data file
void bbob::storeData(std::vector<data> &vec, double evals, double F, double bestF, decision_vector x) const
{
    data temp;
    temp.num = evals;
    temp.F = F;
    temp.bestF = bestF;
    temp.x = x;

    vec.push_back(temp);
}

//write a formatted line into a data file
void bbob::writeData(fs::path file, std::vector<data> &vec) const
{
    FILE * fout;
    bbobOpenFile(fout, file);

    for(int i=0; i<vec.size(); i++)
    {
        fprintf(fout, "%.0f %+10.9e %+10.9e %+10.9e %+10.9e", vec[i].num, vec[i].F-m_bestF, vec[i].bestF-m_bestF, vec[i].F, vec[i].bestF);
        if(vec[i].x.size() < 22)
            for(int j=0; j<vec[i].x.size(); j++)
                fprintf(fout, " %+5.4e", vec[i].x[j]);
        fprintf(fout, "\n");
    }
    fclose(fout);
}

//opens a file after checking it is there.
void bbob::bbobOpenFile(FILE* &fileId, fs::path fileName) const
{
    fileId = fopen(fileName.c_str(), "a");
    if(fileId == NULL)
        pagmo_throw(std::runtime_error, "Could not open file");
}

void bbob::restart(std::string restart_reason) const
{
  FILE * rdataFileId;
  bbobOpenFile(rdataFileId, m_rdataFilePath);
  fprintf(rdataFileId, "%% restart: '%s'\n", restart_reason.c_str());
  fclose(rdataFileId);
  storeData(m_rdataFile, m_LastEval.num, m_LastEval.F, m_BestFEval.F, m_LastEval.x);
}

void bbob::finalize(population &pop) const
{
    (boost::dynamic_pointer_cast<bbob>(pagmo::population_access::get_problem_ptr(pop)))->writeFinalData();
}

std::string bbob::get_name() const
{
    std::string retval("BBOB benchmarking - ");
    retval.append(m_original_problem->get_name());

    return retval;
}

problem::base_ptr bbob::clone() const
{
    return problem::base_ptr(new bbob(*this));
}

}} //namespace

#endif
