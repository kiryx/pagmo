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

#include "benchmark.h"

namespace pagmo { namespace util { namespace coco {
/**
 * Constructs a meta problem(benchmark) which is used for benchmarking any pagmo::problem 
 * @param[in] p pagmo::problem to be benchmarked
 * @param[in] datapath Path where the output log files would be stored
 * @param[in] algname Name of the algorithm to be used for post processing
 * @param[in] instanceId Instance Id to be used for post processing
 * @param[in] comments comments to be used for post processing
 */
benchmark::benchmark(const pagmo::problem::base & p, const std::string datapath, 
                        const std::string algname, unsigned int instanceId, std::string comments) : pagmo::problem::base_meta(p,
                                p.get_dimension(), // Ambiguous without the cast ...
                                p.get_i_dimension(),
                                1, //We can only benchmark using single objective functions.
                                p.get_c_dimension(),
                                p.get_ic_dimension(),
                                p.get_c_tol()), 
                                m_algName(algname),
                                m_comments(comments),
                                m_instanceId(instanceId)
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

    else //check if dimensions, precision, algorithm and comments are same
    {
        unsigned int found=0;
        unsigned int old_dim;
        double old_precision;
        std::string old_algname;
        std::string old_comments;
        std::ifstream indexfile;
        std::string line;
        std::vector<std::string> tokens;

        indexfile.open(m_indexFilePath.string(), std::ifstream::in);

        if(indexfile.is_open())
        {
            while(!indexfile.eof())
            {
                //read line
                std::getline(indexfile, line);
                
                //skip comments
                if(line[0] == '%')
                    continue;

                //try to find "funcId = "
                std::string pattern("funcId = ");
                std::size_t position = line.find(pattern);
                if(position != std::string::npos) //required line
                {
                    found = 0;
                    boost::split(tokens, line, boost::is_any_of(","));

                    for(std::vector<std::string>::size_type i = 0; i < tokens.size(); i++)
                    {
                        //get dimensions
                        pattern = "DIM = ";
                        position = tokens[i].find(pattern);
                        if(position != std::string::npos)
                        {
                            position += pattern.length();
                            std::istringstream ss(tokens[i].substr(position));
                            ss >> old_dim;

                            if(old_dim == p.get_dimension())
                                ++found;
                        }

                        //get precision
                        pattern = "Precision = ";
                        position = tokens[i].find(pattern);
                        if(position != std::string::npos)
                        {
                            position += pattern.length();
                            std::istringstream ss(tokens[i].substr(position));
                            ss >> old_precision;

                            if(fabs(old_precision - m_precision) < 1e-10)
                                ++found;
                        }

                        //get algorithm name
                        pattern = "algId = '";
                        position = tokens[i].find(pattern);
                        if(position != std::string::npos)
                        {
                            position += pattern.length();
                            old_algname = tokens[i].substr(position, tokens[i].length() - position - 1); //remove '

                            if(old_algname.compare(algname) == 0)
                                ++found;
                        }
                    }
                    if(!indexfile.eof())
                    {    
                        //compare comments
                        std::getline(indexfile, line);

                        //exclude %
                        if(m_comments.compare(line.substr(1, line.length() - 2)) == 0)
                            ++found;
                    }
                }
            }
            indexfile.close();
            if(found == 4) //if all parameters are same for the last run
                addIndexEntry();

            else //Parameters changed
                writeNewIndexEntry();
        }
        else
            pagmo_throw(std::runtime_error, "Could not open index file");
    }

    writeDataHeader(m_dataFilePath);
    writeDataHeader(m_hdataFilePath);
    writeDataHeader(m_rdataFilePath);
}

//copy constructor
benchmark::benchmark(const benchmark &obj) : base_meta(obj),
            m_dataFile(obj.m_dataFile),
            m_rdataFile(obj.m_rdataFile),
            m_hdataFile(obj.m_hdataFile),
            m_lastEvalInit(obj.m_lastEvalInit),
            m_bestF(obj.m_bestF),
            m_fTrigger(obj.m_fTrigger),
            m_evalsTrigger(obj.m_evalsTrigger),
            m_idxEvalsTrigger(obj.m_idxEvalsTrigger),
            m_idxDIMEvalsTrigger(obj.m_idxDIMEvalsTrigger),
            m_idxFTrigger(obj.m_idxFTrigger),
            m_LastEval(obj.m_LastEval),
            m_BestFEval(obj.m_BestFEval),
            m_lastWriteEval(obj.m_lastWriteEval),
            m_algName(obj.m_algName),
            m_comments(obj.m_comments),
            m_precision(obj.m_precision),
            m_instanceId(obj.m_instanceId),
            m_dataPath(obj.m_dataPath),
            m_dirPath(obj.m_dirPath),
            m_indexFilePath(obj.m_indexFilePath),
            m_dataFilePath(obj.m_dataFilePath),
            m_hdataFilePath(obj.m_hdataFilePath),
            m_rdataFilePath(obj.m_rdataFilePath),
            m_runCounter(obj.m_runCounter) {}

//Implementation of the objective function
void benchmark::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
    unsigned int boolImprovement = 0;
    double Fvalue;
    double evalsj;

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
void benchmark::writeDataHeader(fs::path dataFilePath) const
{
    std::ofstream dataFile;
    benchmarkOpenFile(dataFile, dataFilePath);

    dataFile << boost::str(boost::format("%% function evaluation | fitness - Fopt (%13.12e) | best fitness - Fopt | measured fitness | best measured fitness | x1 | x2...\n")
        % m_bestF);
    dataFile.close();
}

//open the index file and write a new index entry
void benchmark::writeNewIndexEntry(void) const
{
    std::ofstream indexFile;
    int newline = 1;
    if(!fs::is_regular_file(m_indexFilePath))
        newline = 0;

    benchmarkOpenFile(indexFile, m_indexFilePath);

    if(newline == 1)
        indexFile <<"\n";

    indexFile << boost::str(boost::format("funcId = '%s', DIM = %lu, Precision = %.3e, Fopt = %13.12e, algId = '%s'\n") % m_original_problem->get_name() %
        get_dimension() % m_precision % m_bestF % m_algName);

    indexFile << boost::str(boost::format("%% %s\n%s, %d") % m_comments % (m_dirPath / m_hdataFilePath.filename()).string() % m_instanceId);
    indexFile.close();
}

//Open the index file and write a new index entry */
void benchmark::addIndexEntry(void) const
{
    std::ofstream indexFile;

    if (!fs::is_regular_file(m_indexFilePath))
        pagmo_throw(std::runtime_error, "Could not find index file");

    benchmarkOpenFile(indexFile, m_indexFilePath);

    indexFile << boost::str(boost::format(", %d") % m_instanceId);
    indexFile.close();
}

//complete the data file with unwritten information
void benchmark::writeFinalData(void) const
{
    std::ofstream indexFile;

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
    benchmarkOpenFile(indexFile, m_indexFilePath);
    indexFile << boost::str(boost::format(":%.0f|%.1e") % m_LastEval.num % (m_BestFEval.F - m_bestF - m_precision));
    indexFile.close();

    //now write the data from data vectors to files.
    writeData(m_dataFilePath, m_dataFile);
    writeData(m_rdataFilePath, m_rdataFile);
    writeData(m_hdataFilePath, m_hdataFile);
}

//rewrite the data file with the information about the best F.
void benchmark::storeBestF(std::vector<data> &dataFile, LastEvalStruct BestFEval) const
{
    //insert the evaluation at the correct position
    std::vector<data> temp;
    data best;

    //copy data into best.
    best.num = BestFEval.num;
    best.F = BestFEval.F;
    best.bestF = BestFEval.F;
    best.x = BestFEval.x;

    for(std::vector<data>::size_type j=0; j<dataFile.size(); j++)
    {
        if(dataFile[j].num > BestFEval.num)
            temp.push_back(best);

        temp.push_back(dataFile[j]);
    }

    dataFile = temp;
}

//write a formatted line into a data file
void benchmark::storeData(std::vector<data> &vec, double evals, double F, double bestF, decision_vector x) const
{
    data temp;
    temp.num = evals;
    temp.F = F;
    temp.bestF = bestF;
    temp.x = x;

    vec.push_back(temp);
}

//write a formatted line into a data file
void benchmark::writeData(fs::path file, std::vector<data> &vec) const
{
    std::ofstream fout;
    benchmarkOpenFile(fout, file);

    for(std::vector<data>::size_type i=0; i<vec.size(); i++)
    {
        fout << boost::str(boost::format("%.0f %+10.9e %+10.9e %+10.9e %+10.9e") % vec[i].num % (vec[i].F-m_bestF) % (vec[i].bestF-m_bestF) % (vec[i].F) % vec[i].bestF);
        if(vec[i].x.size() < 22)
            for(decision_vector::size_type j=0; j<vec[i].x.size(); j++)
                fout << boost::str(boost::format(" %+5.4e") % vec[i].x[j]);
        fout << "\n";
    }
    fout.close();
}

//opens a file after checking it is there.
void benchmark::benchmarkOpenFile(std::ofstream &file, fs::path fileName) const
{
    file.open(fileName.string(), std::ofstream::out | std::ofstream::app);
    if(!file.is_open())
        pagmo_throw(std::runtime_error, "Could not open file");
}

void benchmark::restart(std::string restart_reason) const
{
  std::ofstream rdataFile;
  benchmarkOpenFile(rdataFile, m_rdataFilePath);
  rdataFile << boost::str(boost::format( "%% restart: '%s'\n") % restart_reason);
  rdataFile.close();
  storeData(m_rdataFile, m_LastEval.num, m_LastEval.F, m_BestFEval.F, m_LastEval.x);
}

/**
 * Finalize benchmarking (and write final data to log files)
 * @param[in] pop pagmo::population
 */
void benchmark::finalize(population &pop) const
{
    (dynamic_cast<const benchmark&>(pop.problem())).writeFinalData();
}

std::string benchmark::get_name() const
{
    std::string retval("COCO benchmarking - ");
    retval.append(m_original_problem->get_name());

    return retval;
}

pagmo::problem::base_ptr benchmark::clone() const
{
    return pagmo::problem::base_ptr(new benchmark(*this));
}

}}} //namespace

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::util::coco::benchmark)
