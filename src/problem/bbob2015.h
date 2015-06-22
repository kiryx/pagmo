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

#ifndef PAGMO_PROBLEM_BBOB2015_H
#define PAGMO_PROBLEM_BBOB2015_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "base_stochastic.h"

namespace pagmo{ namespace problem {

    /// The BBOB 2015 problems: Real-Parameter Black-Box Optimization Benchmarking
    /**
    * This class allows to instantiate any of the 54 benchmark functions of
    * Black-Box Optimization Benchmarking testbed for Real-Parameter
    * Blackbox Optimization at CEC'2015.
    *
    * @see http://coco.gforge.inria.fr/doku.php?id=cec-bbob-2015
    *
    * @author Sunil Kumar Mahendrakar (sunil.mahendrakar19@gmail.com)
    */
class __PAGMO_VISIBLE bbob2015 : public base_stochastic
{
    public:
        bbob2015(unsigned int = 1, problem::base::size_type = 30);

        base_ptr clone() const;

        std::string get_name() const;

        void set_seed(unsigned int) const;
    protected:
        void objfun_impl(fitness_vector &, const decision_vector &) const;

    private:
        /* the return type of all benchmark functions: 2 doubles, the true value and the noisy value - are equal in case of non-noisy functions */
        struct twoDoubles {
            double Ftrue;
            double Fval;
        };

        typedef struct twoDoubles TwoDoubles;

        /* The type of the benchmark functions*/
        typedef TwoDoubles (bbob2015::*bbobFunction)(const decision_vector &x) const;

        bbobFunction m_actFunc = NULL;

        const unsigned int m_problem_number;
        double m_precision;
        unsigned int m_dim;
        unsigned int m_trialid = 1;
        mutable unsigned int m_seedn;
        mutable unsigned int m_isInitDone = 0;

        /*Benchmark helper variables*/
        mutable std::vector<double> m_tmpvect;
        mutable std::vector<double> m_tmx;
        mutable std::vector<double> m_Xopt;
        mutable double m_Fopt;

        mutable std::vector<std::vector<double>> m_rotation;
        mutable std::vector<std::vector<double>> m_rot2;
        mutable std::vector<std::vector<double>> m_linearTF;

        mutable std::vector<double> m_peaks;
        mutable std::vector<double> m_peaks21;
        mutable std::vector<double> m_peaks22;

        mutable std::vector<int> m_rperm;
        mutable std::vector<int> m_rperm21;
        mutable std::vector<int> m_rperm22;

        mutable std::vector<std::vector<double>> m_Xlocal;
        mutable std::vector<std::vector<double>> m_Xlocal21;
        mutable std::vector<std::vector<double>> m_Xlocal22;

        mutable std::vector<std::vector<double>> m_arrScales;
        mutable std::vector<std::vector<double>> m_arrScales21;
        mutable std::vector<std::vector<double>> m_arrScales22;

        mutable std::vector<double> m_gval;
        mutable std::vector<double> m_gval2;
        mutable std::vector<double> m_gvect;
        mutable std::vector<double> m_uniftmp;

        /*Benchmark helper functions*/
        void unif(std::vector<double> &, unsigned int, unsigned int) const;
        void gauss(std::vector<double> &, unsigned int, unsigned int) const;
        void computeXopt(unsigned int, unsigned int) const;
        void monotoneTFosc(std::vector<double> &) const;
        void reshape(std::vector<std::vector<double>> &, const std::vector<double>, int, int) const;
        void computeRotation(std::vector<std::vector<double>> &, unsigned int, unsigned int) const;
        double myrand(void) const;
        double randn(void) const;
        double FGauss(double Ftrue, double beta) const;
        double FUniform(double Ftrue, double alpha, double beta) const;
        double FCauchy(double Ftrue, double alpha, double p) const;

        struct compare_doubles {
                compare_doubles(const bbob2015& c) : bbob(c) {}
                bool operator () (int a, int b)
                {
                    return bbob.m_peaks[a] < bbob.m_peaks[b];
                }
                const bbob2015& bbob;
        };

        double computeFopt(int _funcId, int _trialId) const;

        /*initialize and finalize Benchmarks*/
        void initbenchmarks(void) const;

        /*Non-noisy benchmark functions*/
        TwoDoubles f1(const decision_vector &x) const;
        TwoDoubles f2(const decision_vector &x) const;
        TwoDoubles f3(const decision_vector &x) const;
        TwoDoubles f4(const decision_vector &x) const;
        TwoDoubles f5(const decision_vector &x) const;
        TwoDoubles f6(const decision_vector &x) const;
        TwoDoubles f7(const decision_vector &x) const;
        TwoDoubles f8(const decision_vector &x) const;
        TwoDoubles f9(const decision_vector &x) const;
        TwoDoubles f10(const decision_vector &x) const;
        TwoDoubles f11(const decision_vector &x) const;
        TwoDoubles f12(const decision_vector &x) const;
        TwoDoubles f13(const decision_vector &x) const;
        TwoDoubles f14(const decision_vector &x) const;
        TwoDoubles f15(const decision_vector &x) const;
        TwoDoubles f16(const decision_vector &x) const;
        TwoDoubles f17(const decision_vector &x) const;
        TwoDoubles f18(const decision_vector &x) const;
        TwoDoubles f19(const decision_vector &x) const;
        TwoDoubles f20(const decision_vector &x) const;
        TwoDoubles f21(const decision_vector &x) const;
        TwoDoubles f22(const decision_vector &x) const;
        TwoDoubles f23(const decision_vector &x) const;
        TwoDoubles f24(const decision_vector &x) const;

        /*Noisy benchmark functions*/
        TwoDoubles f101(const decision_vector &x) const;
        TwoDoubles f102(const decision_vector &x) const;
        TwoDoubles f103(const decision_vector &x) const;
        TwoDoubles f104(const decision_vector &x) const;
        TwoDoubles f105(const decision_vector &x) const;
        TwoDoubles f106(const decision_vector &x) const;
        TwoDoubles f107(const decision_vector &x) const;
        TwoDoubles f108(const decision_vector &x) const;
        TwoDoubles f109(const decision_vector &x) const;
        TwoDoubles f110(const decision_vector &x) const;
        TwoDoubles f111(const decision_vector &x) const;
        TwoDoubles f112(const decision_vector &x) const;
        TwoDoubles f113(const decision_vector &x) const;
        TwoDoubles f114(const decision_vector &x) const;
        TwoDoubles f115(const decision_vector &x) const;
        TwoDoubles f116(const decision_vector &x) const;
        TwoDoubles f117(const decision_vector &x) const;
        TwoDoubles f118(const decision_vector &x) const;
        TwoDoubles f119(const decision_vector &x) const;
        TwoDoubles f120(const decision_vector &x) const;
        TwoDoubles f121(const decision_vector &x) const;
        TwoDoubles f122(const decision_vector &x) const;
        TwoDoubles f123(const decision_vector &x) const;
        TwoDoubles f124(const decision_vector &x) const;
        TwoDoubles f125(const decision_vector &x) const;
        TwoDoubles f126(const decision_vector &x) const;
        TwoDoubles f127(const decision_vector &x) const;
        TwoDoubles f128(const decision_vector &x) const;
        TwoDoubles f129(const decision_vector &x) const;
        TwoDoubles f130(const decision_vector &x) const;

        /* array of function pointers*/
        bbobFunction handles[24] = { &bbob2015::f1, &bbob2015::f2, &bbob2015::f3, &bbob2015::f4, &bbob2015::f5, &bbob2015::f6,
            &bbob2015::f7, &bbob2015::f8, &bbob2015::f9, &bbob2015::f10, &bbob2015::f11, &bbob2015::f12, &bbob2015::f13,
            &bbob2015::f14, &bbob2015::f15, &bbob2015::f16, &bbob2015::f17, &bbob2015::f18, &bbob2015::f19, &bbob2015::f20,
             &bbob2015::f21, &bbob2015::f22, &bbob2015::f23, &bbob2015::f24};
        unsigned int handlesLength = 24;

        bbobFunction handlesNoisy[30] = { &bbob2015::f101, &bbob2015::f102, &bbob2015::f103, &bbob2015::f104, &bbob2015::f105,
            &bbob2015::f106, &bbob2015::f107, &bbob2015::f108, &bbob2015::f109, &bbob2015::f110, &bbob2015::f111, &bbob2015::f112,
            &bbob2015::f113, &bbob2015::f114, &bbob2015::f115, &bbob2015::f116, &bbob2015::f117, &bbob2015::f118, &bbob2015::f119,
            &bbob2015::f120, &bbob2015::f121, &bbob2015::f122, &bbob2015::f123, &bbob2015::f124, &bbob2015::f125, &bbob2015::f126,
            &bbob2015::f127, &bbob2015::f128, &bbob2015::f129, &bbob2015::f130};
        unsigned int handlesNoisyLength = 30;

		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & const_cast<unsigned int&>(m_problem_number);
		}
};
}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::bbob2015)

#endif
