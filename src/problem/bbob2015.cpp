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

#include <string>
#include "../serialization.h"
#include "../types.h"
#include "base.h"

#include "bbob2015.h"

#define NHIGHPEAKS21 101
#define NHIGHPEAKS22 21

#ifndef M_PI
    #define M_PI        3.14159265358979323846
#endif

#define TOL 1e-8

namespace pagmo{ namespace problem {
   /// Constructor
   /**
    * Will construct one of the 54 BBOB2015 benchmark problems
    *
    * @param[in] problem_number The problem id. One of [1,2,...,24] OR [101,102,...,130]
    * @param[in] d problem dimension. One of [2,5,10,20,30,...,100]
    *
    * @see http://coco.gforge.inria.fr/doku.php?id=cec-bbob-2015
    */

bbob2015::bbob2015(unsigned int problem_number, problem::base::size_type dim):base_stochastic(dim,0),m_problem_number(problem_number),m_dim(dim)
{
    if(m_dim <= 0)
        pagmo_throw(value_error, "Invalid dimensions specified");

    //Get the function pointer of the corresponding benchmark function. m_actFunc will later be used to find the fitness.
    if (problem_number <= handlesLength )
        m_actFunc = handles[problem_number-1];
    else if( (100 < problem_number) && (problem_number-101 <= handlesNoisyLength) )
        m_actFunc = handlesNoisy[problem_number - 101];
    else
        pagmo_throw(value_error, "problem_number specified is not a valid function of BBOB testbed");

    //All bbob problems have same bounds.
    set_bounds(-5,5);

    //initialize benchmark by calling it once and store best_x.
    initbenchmarks();
    decision_vector x(dim,0);
    (this->*m_actFunc)(x);

    //set optimal decision vector
    std::vector<decision_vector> xopt;
    xopt.push_back(m_Xopt);
    set_best_x(xopt);
}

//Get problem  name.
std::string  bbob2015::get_name() const
{
    std::string  retval("BBOB2015 - f");
    switch(m_problem_number)
    {
        case 1:
            retval.append("(Sphere)");
            break;
        case 2:
            retval.append("(Ellipsoid separable)");
            break;
        case 3:
            retval.append("(Rastrigin separable)");
            break;
        case 4:
            retval.append("(Skew Rastrigin-Bueche separ)");
            break;
        case 5:
            retval.append("(Linear slope)");
            break;
        case 6:
            retval.append("(Attractive sector)");
            break;
        case 7:
            retval.append("(Step-ellipsoid)");
            break;
        case 8:
            retval.append("(Rosenbrock original)");
            break;
        case 9:
            retval.append("(Rosenbrock rotated)");
            break;
        case 10:
            retval.append("(Ellipsoid)");
            break;
        case 11:
            retval.append("(Discus)");
            break;
        case 12:
            retval.append("(Bent cigar)");
            break;
        case 13:
            retval.append("(Sharp ridge)");
            break;
        case 14:
            retval.append("(Sum of different powers)");
            break;
        case 15:
            retval.append("(Rastrigin)");
            break;
        case 16:
            retval.append("(Weierstrass)");
            break;
        case 17:
            retval.append("(Schaffer F7, condition 10)");
            break;
        case 18:
            retval.append("(Schaffer F7, condition 1000)");
            break;
        case 19:
            retval.append("(Griewank-Rosenbrock F8F2)");
            break;
        case 20:
            retval.append("(Schwefel x*sin(x))");
            break;
        case 21:
            retval.append("(Gallagher 101 peaks)");
            break;
        case 22:
            retval.append("(Gallagher 21 peaks)");
            break;
        case 23:
            retval.append("(Katsuuras)");
            break;
        case 24:
            retval.append("(Lunacek bi-Rastrigin)");
            break;
        case 101:
            retval.append("(Sphere moderate Gauss)");
            break;
        case 102:
            retval.append("(Sphere moderate unif)");
            break;
        case 103:
            retval.append("(Sphere moderate Cauchy)");
            break;
        case 104:
            retval.append("(Rosenbrock moderate Gauss)");
            break;
        case 105:
            retval.append("(Rosenbrock moderate unif)");
            break;
        case 106:
            retval.append("(Rosenbrock moderate Cauchy)");
            break;
        case 107:
            retval.append("(Sphere Gauss)");
            break;
        case 108:
            retval.append("(Sphere unif)");
            break;
        case 109:
            retval.append("(Sphere Cauchy)");
            break;
        case 110:
            retval.append("(Rosenbrock Gauss)");
            break;
        case 111:
            retval.append("(Rosenbrock  unif)");
            break;
        case 112:
            retval.append("(Rosenbrock  Cauchy)");
            break;
        case 113:
            retval.append("(Step-ellipsoid Gauss)");
            break;
        case 114:
            retval.append("(Step-ellipsoid unif)");
            break;
        case 115:
            retval.append("(Step-ellipsoid Cauchy)");
            break;
        case 116:
            retval.append("(Ellipsoid Gauss)");
            break;
        case 117:
            retval.append("(Ellipsoid unif)");
            break;
        case 118:
            retval.append("(Ellipsoid Cauchy)");
            break;
        case 119:
            retval.append("(Sum of diff powers Gauss)");
            break;
        case 120:
            retval.append("(Sum of diff powers unif)");
            break;
        case 121:
            retval.append("(Sum of diff powers Cauchy)");
            break;
        case 122:
            retval.append("(Schaffer F7 Gauss)");
            break;
        case 123:
            retval.append("(Schaffer F7 unif)");
            break;
        case 124:
            retval.append("(Schaffer F7 Cauchy)");
            break;
        case 125:
            retval.append("(Griewank-Rosenbrock Gauss)");
            break;
        case 126:
            retval.append("(Griewank-Rosenbrock unif)");
            break;
        case 127:
            retval.append("(Griewank-Rosenbrock Cauchy)");
            break;
        case 128:
            retval.append("(Gallagher Gauss)");
            break;
        case 129:
            retval.append("(Gallagher unif)");
            break;
        case 130:
            retval.append("(Gallagher Cauchy)");
            break;
    }

    return retval;
}

void bbob2015::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
    decision_vector::size_type n = x.size();
    bbob2015::TwoDoubles res;
    if(n != m_dim)
        pagmo_throw(value_error, "Dimension mismatch");

    //Get the fitness using the corresponding benchmark function.
    res = (this->*m_actFunc)(x);

    f[0] = res.Fval;
    return;
}

//Generates N uniform numbers with starting seed.
void bbob2015::unif(std::vector<double> &r, unsigned int N, unsigned int inseed) const
{
    int aktseed;
    int i, tmp;
    int rgrand[32];
    unsigned int k, aktrand;

    if (inseed < 1)
        inseed = 1;
    aktseed = inseed;

    for (i = 39; i >= 0; i--)
    {
        tmp = (int)floor((double)aktseed/(double)127773);
        aktseed = 16807  * (aktseed - tmp * 127773) - 2836 * tmp;
        if (aktseed < 0)
            aktseed = aktseed + 2147483647;
        if (i < 32)
            rgrand[i] = aktseed;
    }
    aktrand = rgrand[0];

    for (k = 0; k < N; k++)
    {
        tmp = (int)floor((double)aktseed/(double)127773);
        aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
        if (aktseed < 0)
            aktseed = aktseed + 2147483647;
        tmp = (int)floor((double)aktrand / (double)67108865);
        aktrand = rgrand[tmp];
        rgrand[tmp] = aktseed;
        r[k] = (double)aktrand/2.147483647e9;
        if (r[k] == 0.)
        {
            printf("Warning: zero sampled(?), set to 1e-99.\n");
            r[k] = 1e-99;
        }
    }
    return;
}

//Samples N standard normally distributed numbers being the same for a given seed.
void bbob2015::gauss(std::vector<double> &g, unsigned int N, unsigned int seed) const
{
    unsigned int i;

    unif(m_uniftmp, 2*N, seed);

    for (i = 0; i < N; i++)
    {
        g[i] = sqrt(-2*log(m_uniftmp[i])) * cos(2*M_PI*m_uniftmp[N+i]);
        if (g[i] == 0.)
            g[i] = 1e-99;
    }
    return;
}

//Compute best_x and store in m_Xopt.
void bbob2015::computeXopt(unsigned int seed, unsigned int _DIM) const
{
    unsigned int i;

    unif(m_tmpvect, _DIM, seed);
    for (i = 0; i < _DIM; i++)
    {
        m_Xopt[i] = 8 * floor(1e4 * m_tmpvect[i])/1e4 - 4;
        if (m_Xopt[i] == 0.0)
            m_Xopt[i] = -1e-5;
    }
}

//Monotone transformation
//Maps [-inf,inf] to [-inf,inf] with different constants for positive and negative part.
void bbob2015::monotoneTFosc(std::vector<double> &f) const
{
    double a = 0.1;
    unsigned int i;
    for (i = 0; i < m_dim; i++)
    {
        if (f[i] > 0)
        {
            f[i] = log(f[i])/a;
            f[i] = pow(exp(f[i] + 0.49*(sin(f[i]) + sin(0.79*f[i]))), a);
        }
        else if (f[i] < 0)
        {
            f[i] = log(-f[i])/a;
            f[i] = -pow(exp(f[i] + 0.49*(sin(0.55 * f[i]) + sin(0.31*f[i]))), a);
        }
    }
}

//Reshape
void bbob2015::reshape(std::vector<std::vector<double>> &B, const std::vector<double> vector, int m, int n) const
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            B[i][j] = vector[j * m + i];
        }
    }
    return;
}

//Get an orthogonal basis.
void bbob2015::computeRotation(std::vector<std::vector<double>> &B, unsigned int seed, unsigned int _DIM) const
{
    double prod;
    unsigned int i, j, k; /*Loop over pairs of column vectors*/

    gauss(m_gvect, _DIM * _DIM, seed);
    reshape(B, m_gvect, _DIM, _DIM);
    /*1st coordinate is row, 2nd is column.*/

    for (i = 0; i < _DIM; i++)
    {
        for (j = 0; j < i; j++)
        {
            prod = 0;
            for (k = 0; k < _DIM; k++)
            {
                prod += B[k][i] * B[k][j];
            }
            for (k = 0; k < _DIM; k++)
            {
                B[k][i] -= prod * B[k][j];
            }
        }
        prod = 0;
        for (k = 0; k < _DIM; k++)
        {
            prod += B[k][i] * B[k][i];
        }
        for (k = 0; k < _DIM; k++)
        {
            B[k][i] /= sqrt(prod);
        }
    }
    return;
}

double bbob2015::myrand(void) const
{
    m_seed++;
    if (m_seed > 1e9)
        m_seed = 1;
    unif(m_uniftmp, 1, m_seed);
    return m_uniftmp[0];
}

double bbob2015::randn(void) const
{
    m_seedn++;
    if (m_seedn > 1e9)
        m_seedn = 1;
    gauss(m_uniftmp, 1, m_seedn);
    return m_uniftmp[0];
}

//Returns Gaussian model noisy value.
double bbob2015::FGauss(double Ftrue, double beta) const
{
    double Fval = Ftrue * exp(beta * randn());
    Fval += 1.01 * TOL;
    if (Ftrue < TOL)
    {
        Fval = Ftrue;
    }
    return Fval;
}

//Returns uniform model noisy value.
double bbob2015::FUniform(double Ftrue, double alpha, double beta) const
{
    double Fval = pow(myrand(), beta) * Ftrue * fmax(1., pow(1e9/(Ftrue+1e-99), alpha * myrand()));
    Fval += 1.01 * TOL;
    if (Ftrue < TOL)
    {
        Fval = Ftrue;
    }
    return Fval;
}

//Returns Cauchy model noisy value.
double bbob2015::FCauchy(double Ftrue, double alpha, double p) const
{
    double Fval;
    double tmp = randn()/fabs(randn()+1e-199);

    if (myrand() < p)
        Fval = Ftrue + alpha * fmax(0., 1e3 + tmp);
    else
        Fval = Ftrue + alpha * 1e3;

    Fval += 1.01 * TOL;
    if (Ftrue < TOL)
    {
        Fval = Ftrue;
    }
    return Fval;
}

//Compute the best function value.
double bbob2015::computeFopt(int _funcId, int _m_trialid) const
{
    unsigned int rseed, rrseed;
    if (_funcId == 4)
        rseed = 3;
    else if (_funcId == 18)
        rseed = 17;
    else if (_funcId == 101 || _funcId == 102 || _funcId == 103 || _funcId == 107 || _funcId == 108 || _funcId == 109)
        rseed = 1;
    else if (_funcId == 104 || _funcId == 105 || _funcId == 106 || _funcId == 110 || _funcId == 111 || _funcId == 112)
        rseed = 8;
    else if (_funcId == 113 || _funcId == 114 || _funcId == 115)
        rseed = 7;
    else if (_funcId == 116 || _funcId == 117 || _funcId == 118)
        rseed = 10;
    else if (_funcId == 119 || _funcId == 120 || _funcId == 121)
        rseed = 14;
    else if (_funcId == 122 || _funcId == 123 || _funcId == 124)
        rseed = 17;
    else if (_funcId == 125 || _funcId == 126 || _funcId == 127)
        rseed = 19;
    else if (_funcId == 128 || _funcId == 129 || _funcId == 130)
        rseed = 21;
    else
        rseed = _funcId;

    rrseed = rseed + 10000 * _m_trialid;
    gauss(m_gval, 1, rrseed);
    gauss(m_gval2, 1, rrseed + 1);
    return fmin(1000., fmax(-1000., (round(100.*100.*m_gval[0]/m_gval2[0])/100.)));
}

//Initialize benchmark helper variables.
void bbob2015::initbenchmarks(void) const
{
    unsigned int i;

    m_gval.resize(1);
    m_gval2.resize(1);
    m_gvect.resize(m_dim * m_dim);
    m_uniftmp.resize(2 * m_dim * m_dim);
    m_tmpvect.resize(m_dim);
    m_Xopt.resize(m_dim);

    m_tmx.resize(m_dim);
    m_rotation.resize(m_dim);
    m_rot2.resize(m_dim);
    m_linearTF.resize(m_dim);

    m_peaks21.resize(m_dim * NHIGHPEAKS21);
    m_rperm21.resize(fmax(m_dim, NHIGHPEAKS21 - 1));
    m_Xlocal21.resize(m_dim);
    m_arrScales21.resize(NHIGHPEAKS21);

    m_peaks22.resize(m_dim * NHIGHPEAKS22);
    m_rperm22.resize(fmax(m_dim, NHIGHPEAKS22 - 1));
    m_arrScales22.resize(NHIGHPEAKS22);
    m_Xlocal22.resize(m_dim);

    for (i = 0; i < m_dim; i++)
    {
        m_rotation[i].resize(m_dim);
        m_rot2[i].resize(m_dim);
        m_linearTF[i].resize(m_dim);
        m_Xlocal21[i].resize(NHIGHPEAKS21);
        m_Xlocal22[i].resize(NHIGHPEAKS22);
    }

    for (i = 0; i < NHIGHPEAKS21; i++)
        m_arrScales21[i].resize(m_dim);

    for (i = 0; i < NHIGHPEAKS22; i++)
        m_arrScales22[i].resize(m_dim);

    m_isInitDone = 0;

    return;
}

/*Noiseless functions testbed. All functions are ranged in [-5, 5]^DIM.*/

//Sphere function
bbob2015::TwoDoubles bbob2015::f1(const decision_vector &x) const
{
    unsigned int i, rseed; /*Loop over dim*/

    static unsigned int funcId = 1;
    double Fadd, r, Fval, Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        m_isInitDone = 1;
    }

    Fadd = m_Fopt;
    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        r = x[i] - m_Xopt[i];
        Ftrue += r * r;
    }
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/
    res.Ftrue = Ftrue;
    res.Fval = Fval;
    return res;
}

//Separable ellipsoid with monotone transformation, condition 1e6
bbob2015::TwoDoubles bbob2015::f2(const decision_vector &x) const
{
    unsigned int i, rseed; /*Loop over dim*/

    static double condition = 1e6;
    static unsigned int funcId = 2;
    double Fadd, Fval, Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        m_isInitDone = 1;
    }

    Fadd = m_Fopt;

    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = x[i] - m_Xopt[i];
    }

    monotoneTFosc(m_tmx);

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        Ftrue += pow(condition, ((double)i)/((double)(m_dim-1))) * m_tmx[i] * m_tmx[i];
    }
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Rastrigin with monotone transformation separable "condition" 10
bbob2015::TwoDoubles bbob2015::f3(const decision_vector &x) const
{
    unsigned int i, rseed; /*Loop over dim*/

    static unsigned int funcId = 3;
    static double condition = 10.;
    static double beta = 0.2;
    double tmp, tmp2, Fadd, Fval, Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        m_isInitDone = 1;
    }

    Fadd = m_Fopt;
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = x[i] - m_Xopt[i];
    }

    monotoneTFosc(m_tmx);
    for (i = 0; i < m_dim; i++)
    {
        tmp = ((double)i)/((double)(m_dim-1));
        if (m_tmx[i] > 0)
            m_tmx[i] = pow(m_tmx[i], 1 + beta * tmp * sqrt(m_tmx[i]));
        m_tmx[i] = pow(sqrt(condition), tmp) * m_tmx[i];
    }
    /* COMPUTATION core*/
    tmp = 0.;
    tmp2 = 0.;
    for (i = 0; i < m_dim; i++)
    {
        tmp += cos(2*M_PI*m_tmx[i]);
        tmp2 += m_tmx[i]*m_tmx[i];
    }
    Ftrue = 10 * (m_dim - tmp) + tmp2;
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Skew Rastrigin-Bueche, condition 10, skew-"condition" 100
bbob2015::TwoDoubles bbob2015::f4(const decision_vector &x) const
{
    unsigned int i, rseed; /*Loop over dim*/
    static unsigned int funcId = 4;
    static double condition = 10.;
    static double alpha = 100.;
    double tmp, tmp2, Fadd, Fval, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = 3 + 10000 * m_trialid; /* Not the same as before.*/
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        for (i = 0; i < m_dim; i += 2)
            m_Xopt[i] = fabs(m_Xopt[i]); /*Skew*/
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
            Fpen += tmp * tmp;
    }
    Fpen *= 1e2;
    Fadd += Fpen;

    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = x[i] - m_Xopt[i];
    }

    monotoneTFosc(m_tmx);
    for (i = 0; i < m_dim; i++)
    {
    if (i % 2 == 0 && m_tmx[i] > 0)
        m_tmx[i] = sqrt(alpha) * m_tmx[i];
    m_tmx[i] = pow(sqrt(condition), ((double)i)/((double)(m_dim-1))) * m_tmx[i];
    }
    /* COMPUTATION core*/
    tmp = 0.;
    tmp2 = 0.;
    for (i = 0; i < m_dim; i++)
    {
        tmp += cos(2*M_PI*m_tmx[i]);
        tmp2 += m_tmx[i]*m_tmx[i];
    }
    Ftrue = 10 * (m_dim - tmp) + tmp2;
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Linear slope
bbob2015::TwoDoubles bbob2015::f5(const decision_vector &x) const
{
    unsigned int i, rseed; /*Loop over dim*/
    static unsigned int funcId = 5;
    static double alpha = 100.;
    static double Fadd; /*Treatment is different from other functions.*/
    double tmp, Fval, Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        Fadd = m_Fopt;
        computeXopt(rseed, m_dim);
        for (i = 0; i < m_dim; i ++)
        {
            tmp = pow(sqrt(alpha), ((double)i)/((double)(m_dim-1)));
            if (m_Xopt[i] > 0)
            {
                m_Xopt[i] = 5.;
            }
            else if (m_Xopt[i] < 0)
            {
                m_Xopt[i] = -5.;
            }
            Fadd += 5. * tmp;
        }
        m_isInitDone = 1;
    }

    /* BOUNDARY HANDLING*/
    /* move "too" good coordinates back into domain*/
    for (i = 0; i < m_dim; i++)
    {
        if ((m_Xopt[i] == 5.) && (x[i] > 5))
            m_tmx[i] = 5.;
        else if ((m_Xopt[i] == -5.) && (x[i] < -5))
            m_tmx[i] = -5.;
        else
            m_tmx[i] = x[i];
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        if (m_Xopt[i] > 0)
            Ftrue -= pow(sqrt(alpha), ((double)i)/((double)(m_dim-1))) * m_tmx[i];
        else
            Ftrue += pow(sqrt(alpha), ((double)i)/((double)(m_dim-1))) * m_tmx[i];
    }
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Attractive sector function
bbob2015::TwoDoubles bbob2015::f6(const decision_vector &x) const
{
    unsigned int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 6;
    static double alpha = 100.;
    double Fadd, Fval, Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        static double condition = 10.;
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        computeRotation(m_rot2, rseed, m_dim);
        /* decouple scaling from function definition*/
        for (i = 0; i < m_dim; i ++)
        {
            for (j = 0; j < m_dim; j++)
            {
                m_linearTF[i][j] = 0.;
                for (k = 0; k < m_dim; k++)
                {
                    m_linearTF[i][j] += m_rotation[i][k] * pow(sqrt(condition), ((double)k)/((double)(m_dim-1))) * m_rot2[k][j];
                }
            }
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_linearTF[i][j] * (x[j] - m_Xopt[j]);
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        if (m_tmx[i] * m_Xopt[i] > 0)
            m_tmx[i] *= alpha;
        Ftrue += m_tmx[i] * m_tmx[i];
    }

    /*MonotoneTFosc...*/
    if (Ftrue > 0)
    {
        Ftrue = pow(exp(log(Ftrue)/0.1 + 0.49*(sin(log(Ftrue)/0.1) + sin(0.79*log(Ftrue)/0.1))), 0.1);
    }
    else if (Ftrue < 0)
    {
        Ftrue = -pow(exp(log(-Ftrue)/0.1 + 0.49*(sin(0.55 * log(-Ftrue)/0.1) + sin(0.31*log(-Ftrue)/0.1))), 0.1);
    }
    Ftrue = pow(Ftrue, 0.9);
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Step-ellipsoid condition 100
bbob2015::TwoDoubles bbob2015::f7(const decision_vector &x) const
{
    unsigned int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 7;
    static double condition = 100.;
    static double alpha = 10.;
    double x1, tmp, Fadd, Fval, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        computeRotation(m_rot2, rseed, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmpvect[i] = 0.;
        tmp = sqrt(pow(condition/10., ((double)i)/((double)(m_dim-1))));
        for (j = 0; j < m_dim; j++)
        {
            m_tmpvect[i] += tmp * m_rot2[i][j] * (x[j] - m_Xopt[j]);
        }
    }
    x1 = m_tmpvect[0];

    for (i = 0; i < m_dim; i++)
    {
        if (fabs(m_tmpvect[i]) > 0.5)
            m_tmpvect[i] = round(m_tmpvect[i]);
        else
            m_tmpvect[i] = round(alpha * m_tmpvect[i])/alpha;
    }

    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * m_tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        Ftrue += pow(condition, ((double)i)/((double)(m_dim-1))) * m_tmx[i] * m_tmx[i];
    }
    Ftrue = 0.1 * fmax(1e-4 * fabs(x1), Ftrue);

    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Rosenbrock, non-rotated
bbob2015::TwoDoubles bbob2015::f8(const decision_vector &x) const
{
    static unsigned int funcId = 8;

    unsigned int i, rseed; /*Loop over dim*/
    static double scales;
    double tmp, Fadd, Fval, Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;

        scales = fmax(1., sqrt((double)m_dim) / 8.);
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        for (i = 0; i < m_dim; i ++)
            m_Xopt[i] *= 0.75;
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = scales * (x[i] - m_Xopt[i]) + 1;
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim - 1; i++)
    {
        tmp = (m_tmx[i] * m_tmx[i] - m_tmx[i+1]);
        Ftrue += tmp * tmp;
    }
    Ftrue *= 1e2;
    for (i = 0; i < m_dim - 1; i ++)
    {
        tmp = (m_tmx[i] - 1.);
        Ftrue += tmp * tmp;
    }
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Rosenbrock, rotated
bbob2015::TwoDoubles bbob2015::f9(const decision_vector &x) const
{
    unsigned int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 9;
    double scales, tmp, Fadd, Fval, Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);

        computeRotation(m_rotation, rseed, m_dim);
        scales = fmax(1., sqrt((double)m_dim) / 8.);
        for (i = 0; i < m_dim; i ++)
        {
            for (j = 0; j < m_dim; j++)
                m_linearTF[i][j] = scales * m_rotation[i][j];
        }

        /*compute Xopt*/
        for (i = 0; i < m_dim; i++)
        {
            m_Xopt[i] = 0.;
            for (j = 0; j < m_dim; j++)
            {
                m_Xopt[i] += m_linearTF[j][i] * 0.5/scales/scales;
            }
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.5;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_linearTF[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim - 1; i++)
    {
        tmp = (m_tmx[i] * m_tmx[i] - m_tmx[i+1]);
        Ftrue += tmp * tmp;
    }
    Ftrue *= 1e2;
    for (i = 0; i < m_dim - 1; i ++)
    {
        tmp = (m_tmx[i] - 1.);
        Ftrue += tmp * tmp;
    }

    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Ellipsoid with monotone transformation, condition 1e6
bbob2015::TwoDoubles bbob2015::f10(const decision_vector &x) const
{
    unsigned int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 10;
    static double condition = 1e6;
    double Fadd, Fval, Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;
    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
    }

    monotoneTFosc(m_tmx);
    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        Ftrue += pow(condition, ((double)i)/((double)(m_dim-1))) * m_tmx[i] * m_tmx[i];
    }
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Discus (tablet) with monotone transformation, condition 1e6
bbob2015::TwoDoubles bbob2015::f11(const decision_vector &x) const
{
    unsigned int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 11;
    static double condition = 1e6;
    double Fadd, Fval, Ftrue;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;
    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
    }

    monotoneTFosc(m_tmx);

    /* COMPUTATION core*/
    Ftrue = condition * m_tmx[0] * m_tmx[0];
    for (i = 1; i < m_dim; i++)
    {
        Ftrue += m_tmx[i] * m_tmx[i];
    }
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/
    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Bent cigar with asymmetric space distortion, condition 1e6
bbob2015::TwoDoubles bbob2015::f12(const decision_vector &x) const
{
    unsigned int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 12;
    static double condition = 1e6;
    static double beta = 0.5;
    double Fadd, Fval, Ftrue;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed + 1000000, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;
    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmpvect[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmpvect[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
        if (m_tmpvect[i] > 0)
        {
            m_tmpvect[i] = pow(m_tmpvect[i], 1 + beta * ((double)i)/((double)(m_dim-1)) * sqrt(m_tmpvect[i]));
        }
    }

    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * m_tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    Ftrue = m_tmx[0] * m_tmx[0];
    for (i = 1; i < m_dim; i++)
    {
        Ftrue += condition * m_tmx[i] * m_tmx[i];
    }
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Sharp ridge
bbob2015::TwoDoubles bbob2015::f13(const decision_vector &x) const
{
    unsigned int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 13;
    static double condition = 10.;
    static double alpha = 100.;
    double Fadd, Fval, Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        computeRotation(m_rot2, rseed, m_dim);

        for (i = 0; i < m_dim; i++)
        {
            for (j = 0; j < m_dim; j++)
            {
                m_linearTF[i][j] = 0.;
                for (k = 0; k < m_dim; k++)
                {
                    m_linearTF[i][j] += m_rotation[i][k] * pow(sqrt(condition), ((double)k)/((double)(m_dim-1))) * m_rot2[k][j];
                }
            }
        }
        m_isInitDone = 1;
        }
    Fadd = m_Fopt;
    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_linearTF[i][j] * (x[j] - m_Xopt[j]);
        }
    }

    /* COMPUTATION core*/
    for (i = 1; i < m_dim; i++)
    {
        Ftrue += m_tmx[i] * m_tmx[i];
    }
    Ftrue = alpha * sqrt(Ftrue);
    Ftrue += m_tmx[0] * m_tmx[0];

    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Sum of different powers, between x^2 and x^6
bbob2015::TwoDoubles bbob2015::f14(const decision_vector &x) const
{
    unsigned int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 14;
    static double alpha = 4.;
    double Fadd, Fval, Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;
    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        Ftrue += pow(fabs(m_tmx[i]), 2. + alpha * ((double)i)/((double)(m_dim-1)));
    }
    Ftrue = sqrt(Ftrue);

    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}


//Rastrigin with asymmetric non-linear distortion, "condition" 10
bbob2015::TwoDoubles bbob2015::f15(const decision_vector &x) const
{
    unsigned int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 15;
    static double condition = 10.;
    static double beta = 0.2;
    double tmp = 0., tmp2 = 0., Fadd, Fval, Ftrue;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        computeRotation(m_rot2, rseed, m_dim);
        for (i = 0; i < m_dim; i++)
        {
            for (j = 0; j < m_dim; j++)
            {
                m_linearTF[i][j] = 0.;
                for (k = 0; k < m_dim; k++)
                {
                    m_linearTF[i][j] += m_rotation[i][k] * pow(sqrt(condition), ((double)k)/((double)(m_dim-1))) * m_rot2[k][j];
                }
            }
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;
    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmpvect[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmpvect[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
    }

    monotoneTFosc(m_tmpvect);
    for (i = 0; i < m_dim; i++)
    {
        if (m_tmpvect[i] > 0)
            m_tmpvect[i] = pow(m_tmpvect[i], 1 + beta * ((double)i)/((double)(m_dim-1)) * sqrt(m_tmpvect[i]));
    }
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_linearTF[i][j] * m_tmpvect[j];
        }
    }
    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        tmp += cos(2. * M_PI * m_tmx[i]);
        tmp2 += m_tmx[i] * m_tmx[i];
    }
    Ftrue = 10. * ((double)m_dim - tmp) + tmp2;
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Weierstrass, condition 100
bbob2015::TwoDoubles bbob2015::f16(const decision_vector &x) const
{
    unsigned int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 16;
    static double condition = 100.;
    static double aK[12];
    static double bK[12];
    static double F0;
    double tmp, Fadd, Fval, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        computeRotation(m_rot2, rseed, m_dim);
        for (i = 0; i < m_dim; i++)
        {
            for (j = 0; j < m_dim; j++)
            {
                m_linearTF[i][j] = 0.;
                for (k = 0; k < m_dim; k++)
                {
                    m_linearTF[i][j] += m_rotation[i][k] * pow(1./sqrt(condition), ((double)k)/((double)(m_dim-1))) * m_rot2[k][j];
                }
            }
        }

        F0 = 0.;
        for (i = 0; i < 12; i ++)
        {
            aK[i] = pow(0.5, (double)i);
            bK[i] = pow(3., (double)i);
            F0 += aK[i] * cos(2 * M_PI * bK[i] * 0.5);
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 10./(double)m_dim * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmpvect[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmpvect[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
    }

    monotoneTFosc(m_tmpvect);
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_linearTF[i][j] * m_tmpvect[j];
        }
    }
    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = 0.;
        for (j = 0; j < 12; j++)
        {
            tmp += cos(2 * M_PI * (m_tmx[i] + 0.5) * bK[j]) * aK[j];
        }
        Ftrue += tmp;
    }
    Ftrue = 10. * pow(Ftrue/(double)m_dim - F0, 3.);
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Schaffers F7 with asymmetric non-linear transformation, condition 10
bbob2015::TwoDoubles bbob2015::f17(const decision_vector &x) const
{
    unsigned int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 17;
    static double condition = 10.;
    static double beta = 0.5;
    double tmp, Fadd, Fval, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        computeRotation(m_rot2, rseed, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 10. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmpvect[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmpvect[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
        if (m_tmpvect[i] > 0)
            m_tmpvect[i] = pow(m_tmpvect[i], 1 + beta * ((double)i)/((double)(m_dim-1)) * sqrt(m_tmpvect[i]));
    }

    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        tmp = pow(sqrt(condition), ((double)i)/((double)(m_dim-1)));
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += tmp * m_rot2[i][j] * m_tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim - 1; i++)
    {
        tmp = m_tmx[i] * m_tmx[i] + m_tmx[i+1] * m_tmx[i+1];
        Ftrue += pow(tmp, 0.25) * (pow(sin(50 * pow(tmp, 0.1)), 2.) + 1.);
    }
    Ftrue = pow(Ftrue/(double)(m_dim - 1), 2.);
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Schaffers F7 with asymmetric non-linear transformation, condition 1000
bbob2015::TwoDoubles bbob2015::f18(const decision_vector &x) const
{
    unsigned int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 18;
    static double condition = 1e3;
    static double beta = 0.5;
    double tmp, Fadd, Fval, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = 17 + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        computeRotation(m_rot2, rseed, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;
    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 10. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmpvect[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmpvect[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
        if (m_tmpvect[i] > 0)
            m_tmpvect[i] = pow(m_tmpvect[i], 1. + beta * ((double)i)/((double)(m_dim-1)) * sqrt(m_tmpvect[i]));
    }

    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        tmp = pow(sqrt(condition), ((double)i)/((double)(m_dim-1)));
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += tmp * m_rot2[i][j] * m_tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim - 1; i++)
    {
        tmp = m_tmx[i] * m_tmx[i] + m_tmx[i+1] * m_tmx[i+1];
        Ftrue += pow(tmp, 0.25) * (pow(sin(50. * pow(tmp, 0.1)), 2.) + 1.);
    }
    Ftrue = pow(Ftrue/(double)(m_dim - 1), 2.);
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//F8F2 sum of Griewank-Rosenbrock 2-D blocks
bbob2015::TwoDoubles bbob2015::f19(const decision_vector &x) const
{
    unsigned int i, j, rseed; /*Loop over dim*/
    static unsigned int funcId = 19;
    double scales, F2, tmp = 0., tmp2, Fadd, Fval, Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);

        scales = fmax(1., sqrt((double)m_dim) / 8.);
        computeRotation(m_rotation, rseed, m_dim);
        for (i = 0; i < m_dim; i ++)
        {
            for (j = 0; j < m_dim; j++)
            {
                m_linearTF[i][j] = scales * m_rotation[i][j];
            }
        }
        /*compute Xopt*/
        for (i = 0; i < m_dim; i++)
        {
            m_Xopt[i] = 0.;
            for (j = 0; j < m_dim; j++)
            {
                m_Xopt[i] += m_linearTF[j][i] * 0.5/scales/scales;
            }
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;
    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.5;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_linearTF[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim - 1; i++)
    {
        tmp2 = m_tmx[i] * m_tmx[i] -m_tmx[i+1];
        F2 = 100. * tmp2 * tmp2;
        tmp2 = 1 - m_tmx[i];
        F2 += tmp2 * tmp2;
        tmp += F2 / 4000. - cos(F2);
    }
    Ftrue = 10. + 10. * tmp / (double)(m_dim - 1);

    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Schwefel with tridiagonal variable transformation
bbob2015::TwoDoubles bbob2015::f20(const decision_vector &x) const
{
    unsigned int i, rseed; /*Loop over dim*/
    static unsigned int funcId = 20;
    static double condition = 10.;
    double tmp, Fadd, Fval, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        unif(m_tmpvect, m_dim, rseed);
        //compute Xopt
        for (i = 0; i < m_dim; i++)
        {
            m_Xopt[i] = 0.5 * 4.2096874633;
            if (m_tmpvect[i] - 0.5 < 0)
                m_Xopt[i] *= -1.;
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmpvect[i] = 2. * x[i];
        if (m_Xopt[i] < 0.)
        m_tmpvect[i] *= -1.;
    }

    m_tmx[0] = m_tmpvect[0];
    for (i = 1; i < m_dim; i++)
    {
        m_tmx[i] = m_tmpvect[i] + 0.25 * (m_tmpvect[i-1] - 2. * fabs(m_Xopt[i-1]));
    }

    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] -= 2 * fabs(m_Xopt[i]);
        m_tmx[i] *= pow(sqrt(condition), ((double)i)/((double)(m_dim-1)));
        m_tmx[i] = 100. * (m_tmx[i] + 2 * fabs(m_Xopt[i]));
    }

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(m_tmx[i]) - 500.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 0.01 * Fpen;

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        Ftrue += m_tmx[i] * sin(sqrt(fabs(m_tmx[i])));
    }
    Ftrue = 0.01 * ((418.9828872724339) - Ftrue / (double)m_dim);

    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Gallagher with 101 Gaussian m_peaks, condition up to 1000, one global m_rotation
bbob2015::TwoDoubles bbob2015::f21(const decision_vector &x) const
{
    unsigned int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 21;
    static double fitvalues[2] = {1.1, 9.1};
    static double maxcondition = 1000.;
    static double arrCondition[NHIGHPEAKS21];
    static double peakvalues[NHIGHPEAKS21];
    static double a = 0.1;
    double tmp2, f = 0., Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    double fac = -0.5 / (double)m_dim;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeRotation(m_rotation, rseed, m_dim);
        m_peaks = m_peaks21;
        unif(m_peaks, NHIGHPEAKS21 - 1, rseed);
        m_rperm = m_rperm21;
        for (i = 0; i < NHIGHPEAKS21 - 1; i++)
            m_rperm[i] = i;

        std::sort(m_rperm.begin(), m_rperm.begin() + NHIGHPEAKS21 - 1, compare_doubles(*this));

        /* Random permutation*/
        arrCondition[0] = sqrt(maxcondition);
        peakvalues[0] = 10;
        for (i = 1; i < NHIGHPEAKS21; i++)
        {
            arrCondition[i] = pow(maxcondition, (double)(m_rperm[i-1])/((double)(NHIGHPEAKS21-2)));
            peakvalues[i] = (double)(i-1)/(double)(NHIGHPEAKS21-2) * (fitvalues[1] - fitvalues[0]) + fitvalues[0];
        }
        m_arrScales = m_arrScales21;
        for (i = 0; i < NHIGHPEAKS21; i++)
        {
            unif(m_peaks, m_dim, rseed + 1000 * i);
            for (j = 0; j < m_dim; j++)
                m_rperm[j] = j;

            std::sort(m_rperm.begin(), m_rperm.begin() + m_dim, compare_doubles(*this));

            for (j = 0; j < m_dim; j++)
            {
                m_arrScales[i][j] = pow(arrCondition[i], ((double)m_rperm[j])/((double)(m_dim-1)) - 0.5);
            }
        }

        unif(m_peaks, m_dim * NHIGHPEAKS21, rseed);
        m_Xlocal = m_Xlocal21;
        for (i = 0; i < m_dim; i++)
        {
            /*compute Xopt*/
            m_Xopt[i] = 0.8 * (10. * m_peaks[i] -5.);
            for (j = 0; j < NHIGHPEAKS21; j++)
            {
                m_Xlocal[i][j] = 0.;
                for (k = 0; k < m_dim; k++)
                {
                    m_Xlocal[i][j] += m_rotation[i][k] * (10. * m_peaks[j * m_dim + k] -5.);
                }
                if (j == 0)
                    m_Xlocal[i][j] *= 0.8;
            }
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < NHIGHPEAKS21; i++)
    {
        tmp2 = 0.;
        for (j = 0; j < m_dim; j++)
        {
            tmp = (m_tmx[j] - m_Xlocal[j][i]);
            tmp2 += m_arrScales[i][j] * tmp * tmp;
        }
        tmp2 = peakvalues[i] * exp(fac * tmp2);
        f = fmax(f, tmp2);
    }

    f = 10. - f;
    if (f > 0)
    {
        Ftrue = log(f)/a;
        Ftrue = pow(exp(Ftrue + 0.49*(sin(Ftrue) + sin(0.79*Ftrue))), a);
    }
    else if (f < 0)
    {
        Ftrue = log(-f)/a;
        Ftrue = -pow(exp(Ftrue + 0.49*(sin(0.55 * Ftrue) + sin(0.31*Ftrue))), a);
    }
    else
        Ftrue = f;

    Ftrue *= Ftrue;
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Gallagher with 21 Gaussian m_peaks, condition up to 1000, one global m_rotation
bbob2015::TwoDoubles bbob2015::f22(const decision_vector &x) const
{
    unsigned int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 22;
    static double fitvalues[2] = {1.1, 9.1};
    static double maxcondition = 1000.;
    static double arrCondition[NHIGHPEAKS22];
    static double peakvalues[NHIGHPEAKS22];
    static double a = 0.1;
    double tmp2, f = 0., Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    double fac = -0.5 / (double)m_dim;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeRotation(m_rotation, rseed, m_dim);
        m_peaks = m_peaks22;
        unif(m_peaks, NHIGHPEAKS22 - 1, rseed);
        m_rperm = m_rperm22;
        for (i = 0; i < NHIGHPEAKS22 - 1; i++)
        m_rperm[i] = i;

        std::sort(m_rperm.begin(), m_rperm.begin() + NHIGHPEAKS22 - 1, compare_doubles(*this));

        /* Random permutation*/
        arrCondition[0] = maxcondition;
        peakvalues[0] = 10;
        for (i = 1; i < NHIGHPEAKS22; i++)
        {
            arrCondition[i] = pow(maxcondition, (double)(m_rperm[i-1])/((double)(NHIGHPEAKS22-2)));
            peakvalues[i] = (double)(i-1)/(double)(NHIGHPEAKS22-2) * (fitvalues[1] - fitvalues[0]) + fitvalues[0];
        }
        m_arrScales = m_arrScales22;
        for (i = 0; i < NHIGHPEAKS22; i++)
        {
            unif(m_peaks, m_dim, rseed + 1000 * i);
            for (j = 0; j < m_dim; j++)
                m_rperm[j] = j;

            std::sort(m_rperm.begin(), m_rperm.begin() + m_dim, compare_doubles(*this));

            for (j = 0; j < m_dim; j++)
            {
                m_arrScales[i][j] = pow(arrCondition[i], ((double)m_rperm[j])/((double)(m_dim-1)) - 0.5);
            }
        }

        unif(m_peaks, m_dim * NHIGHPEAKS22, rseed);
        m_Xlocal = m_Xlocal22;
        for (i = 0; i < m_dim; i++)
        {
            //compute Xopt
            m_Xopt[i] = 0.8 * (9.8 * m_peaks[i] -4.9);
            for (j = 0; j < NHIGHPEAKS22; j++)
            {
                m_Xlocal[i][j] = 0.;
                for (k = 0; k < m_dim; k++)
                {
                    m_Xlocal[i][j] += m_rotation[i][k] * (9.8 * m_peaks[j * m_dim + k] -4.9);
                }
                if (j == 0)
                    m_Xlocal[i][j] *= 0.8;
            }
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < NHIGHPEAKS22; i++)
    {
        tmp2 = 0.;
        for (j = 0; j < m_dim; j++)
        {
            tmp = (m_tmx[j] - m_Xlocal[j][i]);
            tmp2 += m_arrScales[i][j] * tmp * tmp;
        }
        tmp2 = peakvalues[i] * exp(fac * tmp2);
        f = fmax(f, tmp2);
    }

    f = 10. - f;
    if (f > 0)
    {
        Ftrue = log(f)/a;
        Ftrue = pow(exp(Ftrue + 0.49*(sin(Ftrue) + sin(0.79*Ftrue))), a);
    }
    else if (f < 0)
    {
        Ftrue = log(-f)/a;
        Ftrue = -pow(exp(Ftrue + 0.49*(sin(0.55 * Ftrue) + sin(0.31*Ftrue))), a);
    }
    else
        Ftrue = f;

    Ftrue *= Ftrue;
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/
    /* free(m_Xopt); //Not used!*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Katsuura function
bbob2015::TwoDoubles bbob2015::f23(const decision_vector &x) const
{
    unsigned int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 23;
    static double condition = 100.;
    double Fadd = 0., Fpen = 0., tmp, Ftrue = 0., arr, prod = 1., tmp2, Fval;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        computeRotation(m_rot2, rseed, m_dim);

        for (i = 0; i < m_dim; i++)
        {
            for (j = 0; j < m_dim; j++)
            {
                m_linearTF[i][j] = 0.;
                for (k = 0; k < m_dim; k++)
                {
                    m_linearTF[i][j] += m_rotation[i][k] * pow(sqrt(condition), ((double)k)/(double)(m_dim - 1)) * m_rot2[k][j];
                }
            }
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_linearTF[i][j] * (x[j] - m_Xopt[j]);
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = 0.;
        for (j = 1; j < 33; j++)
        {
            tmp2 = pow(2., (double)j);
            arr = m_tmx[i] * tmp2;
            tmp += fabs(arr - round(arr)) / tmp2;
        }
        tmp = 1. + tmp * (double)(i + 1);
        prod *= tmp;
    }
    Ftrue = 10./(double)m_dim/(double)m_dim * (-1. + pow(prod, 10./pow((double)m_dim, 1.2)));
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}


//Lunacek bi-Rastrigin, condition 100 in PPSN 2008, Rastrigin part rotated and scaled
bbob2015::TwoDoubles bbob2015::f24(const decision_vector &x) const
{
    unsigned int i, j, k, rseed; /*Loop over dim*/
    static unsigned int funcId = 24;
    static double condition = 100.;
    static double mu1 = 2.5;
    double Fadd, Fpen = 0., tmp, Ftrue = 0., tmp2 = 0., tmp3 = 0., tmp4 = 0., Fval;
    double s = 1. - 0.5 / (sqrt((double)(m_dim + 20)) - 4.1);
    static double d = 1.;
    double mu2 = -sqrt((mu1 * mu1 - d) / s);
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = funcId + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        computeRotation(m_rot2, rseed, m_dim);
        gauss(m_tmpvect, m_dim, rseed);

        //compute Xopt
        for (i = 0; i < m_dim; i++)
        {
            m_Xopt[i] = 0.5 * mu1;
            if (m_tmpvect[i] < 0.)
                m_Xopt[i] *= -1.;
        }

        for (i = 0; i < m_dim; i++)
        {
            for (j = 0; j < m_dim; j++)
            {
                m_linearTF[i][j] = 0.;
                for (k = 0; k < m_dim; k++)
                {
                    m_linearTF[i][j] += m_rotation[i][k] * pow(sqrt(condition), ((double)k)/((double)(m_dim-1))) * m_rot2[k][j];
                }
            }
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 1e4 * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 2. * x[i];
        if (m_Xopt[i] < 0.)
            m_tmx[i] *= -1.;
    }

    /* COMPUTATION core*/
    tmp = 0.;
    for (i = 0; i < m_dim; i++)
    {
        tmp2 += (m_tmx[i] - mu1) * (m_tmx[i] - mu1);
        tmp3 += (m_tmx[i] - mu2) * (m_tmx[i] - mu2);
        tmp4 = 0.;
        for (j = 0; j < m_dim; j++)
        {
            tmp4 += m_linearTF[i][j] * (m_tmx[j] - mu1);
        }
        tmp += cos(2 * M_PI * tmp4);
    }
    Ftrue = fmin(tmp2, d * (double)m_dim + s * tmp3) + 10. * ((double)m_dim - tmp);
    Ftrue += Fadd;
    Fval = Ftrue; /* without noise*/

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}


/* Noisy functions testbed. All functions are ranged in [-5, 5]^m_dim.*/

//Sphere with moderate Gauss noise
bbob2015::TwoDoubles bbob2015::f101(const decision_vector &x) const
{
    unsigned int i; /*Loop over dim*/
    int rseed;

    static unsigned int funcId = 101;
    static unsigned int rrseed = 1;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = (x[i] - m_Xopt[i]);
        Ftrue += tmp * tmp;
    }

    Fval = FGauss(Ftrue, 0.01);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Sphere with moderate uniform noise
bbob2015::TwoDoubles bbob2015::f102(const decision_vector &x) const
{
    unsigned int i; /*Loop over dim*/
    int rseed;

    static unsigned int funcId = 102;
    static unsigned int rrseed = 1;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = (x[i] - m_Xopt[i]);
        Ftrue += tmp * tmp;
    }
    Fval = FUniform(Ftrue, 0.01 * (0.49 + 1./m_dim), 0.01);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Sphere with moderate Cauchy noise
bbob2015::TwoDoubles bbob2015::f103(const decision_vector &x) const
{
    unsigned int i; /*Loop over dim*/
    int rseed;

    static unsigned int funcId = 103;
    static unsigned int rrseed = 1;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = (x[i] - m_Xopt[i]);
        Ftrue += tmp * tmp;
    }
    Fval = FCauchy(Ftrue, 0.01, 0.05);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Rosenbrock non-rotated with moderate Gauss noise
bbob2015::TwoDoubles bbob2015::f104(const decision_vector &x) const
{
    unsigned int i; /*Loop over dim*/
    int rseed;

    static int funcId = 104;
    static int rrseed = 8;
    static double scales;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        scales = fmax(1., sqrt((double)m_dim) / 8.);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = scales * (x[i] - 0.75 * m_Xopt[i]) + 1;
    }

    /* COMPUTATION core*/
    Ftrue = 0.;
    for (i = 0; i < m_dim - 1; i++)
    {
        tmp = (m_tmx[i] * m_tmx[i] - m_tmx[i+1]);
        Ftrue += tmp * tmp;
    }
    Ftrue *= 1e2;
    for (i = 0; i < m_dim - 1; i ++)
    {
        tmp = (m_tmx[i] - 1);
        Ftrue += tmp * tmp;
    }

    Fval = FGauss(Ftrue, 0.01);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Rosenbrock non-rotated with moderate uniform noise
bbob2015::TwoDoubles bbob2015::f105(const decision_vector &x) const
{
    unsigned int i; /*Loop over dim*/
    int rseed;

    static int funcId = 105;
    static int rrseed = 8;
    static double scales;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        scales = fmax(1., sqrt((double)m_dim) / 8.);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = scales * (x[i] - 0.75 * m_Xopt[i]) + 1;
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim - 1; i++)
    {
        tmp = (m_tmx[i] * m_tmx[i] - m_tmx[i+1]);
        Ftrue += tmp * tmp;
    }
    Ftrue *= 1e2;
    for (i = 0; i < m_dim - 1; i ++)
    {
        tmp = (m_tmx[i] - 1);
        Ftrue += tmp * tmp;
    }

    Fval = FUniform(Ftrue, 0.01 * (0.49 + 1./m_dim), 0.01);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Rosenbrock non-rotated with moderate Cauchy noise
bbob2015::TwoDoubles bbob2015::f106(const decision_vector &x) const
{
    unsigned int i; /*Loop over dim*/
    int rseed;

    static int funcId = 106;
    static int rrseed = 8;
    static double scales;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        scales = fmax(1., sqrt((double)m_dim) / 8.);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = scales * (x[i] - 0.75 * m_Xopt[i]) + 1;
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim - 1; i++)
    {
        tmp = (m_tmx[i] * m_tmx[i] - m_tmx[i+1]);
        Ftrue += tmp * tmp;
    }
    Ftrue *= 1e2;
    for (i = 0; i < m_dim - 1; i ++)
    {
        tmp = (m_tmx[i] - 1);
        Ftrue += tmp * tmp;
    }

    Fval = FCauchy(Ftrue, 0.01, 0.05);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Sphere with Gauss noise
bbob2015::TwoDoubles bbob2015::f107(const decision_vector &x) const
{
    unsigned int i; /*Loop over dim*/
    int rseed;

    static int funcId = 107;
    static int rrseed = 1;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = (x[i] - m_Xopt[i]);
        Ftrue += tmp * tmp;
    }
    Fval = FGauss(Ftrue, 1.);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Sphere with uniform noise
bbob2015::TwoDoubles bbob2015::f108(const decision_vector &x) const
{
    unsigned int i; /*Loop over dim*/
    int rseed;

    static int funcId = 108;
    static int rrseed = 1;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = (x[i] - m_Xopt[i]);
        Ftrue += tmp * tmp;
    }
    Fval = FUniform(Ftrue, 0.49 + 1./m_dim, 1.);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Sphere with Cauchy noise
bbob2015::TwoDoubles bbob2015::f109(const decision_vector &x) const
{
    unsigned int i; /*Loop over dim*/
    int rseed;

    static int funcId = 109;
    static int rrseed = 1;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = (x[i] - m_Xopt[i]);
        Ftrue += tmp * tmp;
    }
    Fval = FCauchy(Ftrue, 1., 0.2);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Rosenbrock non-rotated with Gauss noise
bbob2015::TwoDoubles bbob2015::f110(const decision_vector &x) const
{
    unsigned int i; /*Loop over dim*/
    int rseed;

    static int funcId = 110;
    static int rrseed = 8;
    static double scales;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        scales = fmax(1., sqrt((double)m_dim) / 8.);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = scales * (x[i] - 0.75 * m_Xopt[i]) + 1;
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim - 1; i++)
    {
        tmp = (m_tmx[i] * m_tmx[i] - m_tmx[i+1]);
        Ftrue += tmp * tmp;
    }
    Ftrue *= 1e2;
    for (i = 0; i < m_dim - 1; i ++)
    {
        tmp = (m_tmx[i] - 1);
        Ftrue += tmp * tmp;
    }
    Fval = FGauss(Ftrue, 1.);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Rosenbrock non-rotated with moderate uniform noise
bbob2015::TwoDoubles bbob2015::f111(const decision_vector &x) const
{
    unsigned int i; /*Loop over dim*/
    int rseed;

    static int funcId = 111;
    static int rrseed = 8;
    static double scales;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        scales = fmax(1., sqrt((double)m_dim) / 8.);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = scales * (x[i] - 0.75 * m_Xopt[i]) + 1;
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim - 1; i++)
    {
        tmp = (m_tmx[i] * m_tmx[i] - m_tmx[i+1]);
        Ftrue += tmp * tmp;
    }
    Ftrue *= 1e2;
    for (i = 0; i < m_dim - 1; i ++)
    {
        tmp = (m_tmx[i] - 1);
        Ftrue += tmp * tmp;
    }
    Fval = FUniform(Ftrue, 0.49 + 1./m_dim, 1.);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Rosenbrock non-rotated with moderate Cauchy noise
bbob2015::TwoDoubles bbob2015::f112(const decision_vector &x) const
{
    unsigned int i; /*Loop over dim*/
    int rseed;

    static int funcId = 112;
    static int rrseed = 8;
    static double scales;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        scales = fmax(1., sqrt((double)m_dim) / 8.);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = scales * (x[i] - 0.75 * m_Xopt[i]) + 1;
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim - 1; i++)
    {
        tmp = (m_tmx[i] * m_tmx[i] - m_tmx[i+1]);
        Ftrue += tmp * tmp;
    }
    Ftrue *= 1e2;
    for (i = 0; i < m_dim - 1; i ++)
    {
        tmp = (m_tmx[i] - 1);
        Ftrue += tmp * tmp;
    }
    Fval = FCauchy(Ftrue, 1., 0.2);
    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Step-ellipsoid with gauss noise, condition 100
bbob2015::TwoDoubles bbob2015::f113(const decision_vector &x) const
{
    unsigned int i, j; /*Loop over dim*/
    int rseed;

    static int funcId = 113;
    static int rrseed = 7;
    static double condition = 100.;
    static double alpha = 10.;
    double x1, Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        computeRotation(m_rot2, rseed, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmpvect[i] = 0.;
        tmp = sqrt(pow(condition/10., ((double)i)/((double)(m_dim-1))));
        for (j = 0; j < m_dim; j++)
        {
            m_tmpvect[i] += tmp * m_rot2[i][j] * (x[j] - m_Xopt[j]);
        }
    }
    x1 = m_tmpvect[0];
    for (i = 0; i < m_dim; i++)
    {
        if (fabs(m_tmpvect[i]) > 0.5)
        {
            m_tmpvect[i] = round(m_tmpvect[i]);
        }
        else
        {
            m_tmpvect[i] = round(alpha * m_tmpvect[i])/alpha;
        }
    }

    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * m_tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        Ftrue += pow(condition, ((double)i)/((double)(m_dim-1))) * m_tmx[i] * m_tmx[i];
    }
    Ftrue = 0.1 * fmax(1e-4 * fabs(x1), Ftrue);
    Fval = FGauss(Ftrue, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Step-ellipsoid with uniform noise, condition 100
bbob2015::TwoDoubles bbob2015::f114(const decision_vector &x) const
{
    unsigned int i, j; /*Loop over dim*/
    int rseed;

    static int funcId = 114;
    static int rrseed = 7;
    static double condition = 100.;
    static double alpha = 10.;
    double x1, Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        computeRotation(m_rot2, rseed, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmpvect[i] = 0.;
        tmp = sqrt(pow(condition/10., ((double)i)/((double)(m_dim-1))));
        for (j = 0; j < m_dim; j++)
        {
            m_tmpvect[i] += tmp * m_rot2[i][j] * (x[j] - m_Xopt[j]);
        }
    }
    x1 = m_tmpvect[0];
    for (i = 0; i < m_dim; i++)
    {
        if (fabs(m_tmpvect[i]) > 0.5)
        {
            m_tmpvect[i] = round(m_tmpvect[i]);
        }
        else
        {
            m_tmpvect[i] = round(alpha * m_tmpvect[i])/alpha;
        }
    }

    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * m_tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        Ftrue += pow(condition, ((double)i)/((double)(m_dim-1))) * m_tmx[i] * m_tmx[i];
    }
    Ftrue = 0.1 * fmax(1e-4 * fabs(x1), Ftrue);
    Fval = FUniform(Ftrue, 0.49 + 1./m_dim, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Step-ellipsoid with Cauchy noise, condition 100
bbob2015::TwoDoubles bbob2015::f115(const decision_vector &x) const
{
    unsigned int i, j; /*Loop over dim*/
    int rseed;

    static int funcId = 115;
    static int rrseed = 7;
    static double condition = 100.;
    static double alpha = 10.;
    double x1, Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        computeRotation(m_rot2, rseed, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmpvect[i] = 0.;
        tmp = sqrt(pow(condition/10., ((double)i)/((double)(m_dim-1))));
        for (j = 0; j < m_dim; j++)
        {
            m_tmpvect[i] += tmp * m_rot2[i][j] * (x[j] - m_Xopt[j]);
        }
    }
    x1 = m_tmpvect[0];
    for (i = 0; i < m_dim; i++)
    {
        if (fabs(m_tmpvect[i]) > 0.5)
        {
            m_tmpvect[i] = round(m_tmpvect[i]);
        }
        else
        {
            m_tmpvect[i] = round(alpha * m_tmpvect[i])/alpha;
        }
    }

    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * m_tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        Ftrue += pow(condition, ((double)i)/((double)(m_dim-1))) * m_tmx[i] * m_tmx[i];
    }
    Ftrue = 0.1 * fmax(1e-4 * fabs(x1), Ftrue);
    Fval = FCauchy(Ftrue, 1., 0.2);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Ellipsoid with Gauss noise, monotone x-transformation, condition 1e4
bbob2015::TwoDoubles bbob2015::f116(const decision_vector &x) const
{
    unsigned int i, j; /*Loop over dim*/
    int rseed;

    static int funcId = 116;
    static int rrseed = 10;
    static double condition = 1e4;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
    }

    monotoneTFosc(m_tmx);
    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        Ftrue += pow(condition, ((double)i)/((double)(m_dim-1))) * m_tmx[i] * m_tmx[i];
    }

    Fval = FGauss(Ftrue, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Ellipsoid with uniform noise, monotone x-transformation, condition 1e4
bbob2015::TwoDoubles bbob2015::f117(const decision_vector &x) const
{
    unsigned int i, j; /*Loop over dim*/
    int rseed;

    static int funcId = 117;
    static int rrseed = 10;
    static double condition = 1e4;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
    }

    monotoneTFosc(m_tmx);
    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        Ftrue += pow(condition, ((double)i)/((double)(m_dim-1))) * m_tmx[i] * m_tmx[i];
    }

    Fval = FUniform(Ftrue, 0.49 + 1./m_dim, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Ellipsoid with Cauchy noise, monotone x-transformation, condition 1e4
bbob2015::TwoDoubles bbob2015::f118(const decision_vector &x) const
{
    unsigned int i, j; /*Loop over dim*/
    int rseed;

    static int funcId = 118;
    static int rrseed = 10;
    static double condition = 1e4;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
    }

    monotoneTFosc(m_tmx);
    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        Ftrue += pow(condition, ((double)i)/((double)(m_dim-1))) * m_tmx[i] * m_tmx[i];
    }

    Fval = FCauchy(Ftrue, 1., 0.2);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Sum of different powers with Gauss Noise, between x^2 and x^6
bbob2015::TwoDoubles bbob2015::f119(const decision_vector &x) const
{
    unsigned int i, j; /*Loop over dim*/
    int rseed;

    static int funcId = 119;
    static int rrseed = 14;
    static double alpha = 4.;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        Ftrue += pow(fabs(m_tmx[i]), 2 + alpha * ((double)i)/((double)(m_dim-1)));
    }
    Ftrue = sqrt(Ftrue);

    Fval = FGauss(Ftrue, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Sum of different powers with uniform noise, between x^2 and x^6
bbob2015::TwoDoubles bbob2015::f120(const decision_vector &x) const
{
    unsigned int i, j; /*Loop over dim*/
    int rseed;

    static int funcId = 120;
    static int rrseed = 14;
    static double alpha = 4.;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        Ftrue += pow(fabs(m_tmx[i]), 2 + alpha * ((double)i)/((double)(m_dim-1)));
    }
    Ftrue = sqrt(Ftrue);

    Fval = FUniform(Ftrue, 0.49 + 1./m_dim, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Sum of different powers with seldom Cauchy Noise, between x^2 and x^6
bbob2015::TwoDoubles bbob2015::f121(const decision_vector &x) const
{
    unsigned int i, j; /*Loop over dim*/
    int rseed;

    static int funcId = 121;
    static int rrseed = 14;
    static double alpha = 4.;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim; i++)
    {
        Ftrue += pow(fabs(m_tmx[i]), 2 + alpha * ((double)i)/((double)(m_dim-1)));
    }
    Ftrue = sqrt(Ftrue);

    Fval = FCauchy(Ftrue, 1., 0.2);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Schaffers F7 with Gauss noise, with asymmetric non-linear transformation, condition 10
bbob2015::TwoDoubles bbob2015::f122(const decision_vector &x) const
{
    unsigned int i, j; /*Loop over dim*/
    int rseed;

    static int funcId = 122;
    static int rrseed = 17;
    static double condition = 10.;
    static double beta = 0.5;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        computeRotation(m_rot2, rseed, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmpvect[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmpvect[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
        if (m_tmpvect[i] > 0)
            m_tmpvect[i] = pow(m_tmpvect[i], 1 + beta * ((double)i)/((double)(m_dim-1)) * sqrt(m_tmpvect[i]));
    }

    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        tmp = pow(sqrt(condition), ((double)i)/((double)(m_dim-1)));
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += tmp * m_rot2[i][j] * m_tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim - 1; i++)
    {
        tmp = m_tmx[i] * m_tmx[i] + m_tmx[i+1] * m_tmx[i+1];
        Ftrue += pow(tmp, 0.25) * (pow(sin(50. * pow(tmp, 0.1)), 2.) + 1.);
    }
    Ftrue = pow(Ftrue/(double)(m_dim - 1), 2.);

    Fval = FGauss(Ftrue, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Schaffers F7 with uniform noise, with asymmetric non-linear transformation, condition 10
bbob2015::TwoDoubles bbob2015::f123(const decision_vector &x) const
{
    unsigned int i, j; /*Loop over dim*/
    int rseed;

    static int funcId = 123;
    static int rrseed = 17;
    static double condition = 10.;
    static double beta = 0.5;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        computeRotation(m_rot2, rseed, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmpvect[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmpvect[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
        if (m_tmpvect[i] > 0)
            m_tmpvect[i] = pow(m_tmpvect[i], 1 + beta * ((double)i)/((double)(m_dim-1)) * sqrt(m_tmpvect[i]));
    }

    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        tmp = pow(sqrt(condition), ((double)i)/((double)(m_dim-1)));
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += tmp * m_rot2[i][j] * m_tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim - 1; i++)
    {
        tmp = m_tmx[i] * m_tmx[i] + m_tmx[i+1] * m_tmx[i+1];
        Ftrue += pow(tmp, 0.25) * (pow(sin(50. * pow(tmp, 0.1)), 2.) + 1.);
    }
    Ftrue = pow(Ftrue/(double)(m_dim - 1), 2.);

    Fval = FUniform(Ftrue, 0.49 + 1./m_dim, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Schaffers F7 with seldom Cauchy noise, with asymmetric non-linear transformation, condition 10
bbob2015::TwoDoubles bbob2015::f124(const decision_vector &x) const
{
    unsigned int i, j; /*Loop over dim*/
    int rseed;

    static int funcId = 124;
    static int rrseed = 17;
    static double condition = 10.;
    static double beta = 0.5;
    double Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeXopt(rseed, m_dim);
        computeRotation(m_rotation, rseed + 1000000, m_dim);
        computeRotation(m_rot2, rseed, m_dim);
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmpvect[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmpvect[i] += m_rotation[i][j] * (x[j] - m_Xopt[j]);
        }
        if (m_tmpvect[i] > 0)
            m_tmpvect[i] = pow(m_tmpvect[i], 1 + beta * ((double)i)/((double)(m_dim-1)) * sqrt(m_tmpvect[i]));
    }

    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        tmp = pow(sqrt(condition), ((double)i)/((double)(m_dim-1)));
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += tmp * m_rot2[i][j] * m_tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < m_dim - 1; i++)
    {
        tmp = m_tmx[i] * m_tmx[i] + m_tmx[i+1] * m_tmx[i+1];
        Ftrue += pow(tmp, 0.25) * (pow(sin(50. * pow(tmp, 0.1)), 2.) + 1.);
    }
    Ftrue = pow(Ftrue/(double)(m_dim - 1), 2.);

    Fval = FCauchy(Ftrue, 1., 0.2);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//F8F2 sum of Griewank-Rosenbrock 2-D blocks with Gauss noise
bbob2015::TwoDoubles bbob2015::f125(const decision_vector &x) const
{
    unsigned int i, j; /*Loop over dim*/
    int rseed;

    static int funcId = 125;
    static int rrseed = 19;
    static double scales;
    double F2, Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        scales = fmax(1., sqrt((double)m_dim) / 8.);
        computeRotation(m_rotation, rseed, m_dim);

        //compute Xopt
        for (i = 0; i < m_dim; i++)
        {
            m_Xopt[i] = 0.;
            for (j = 0; j < m_dim; j++)
            {
                m_Xopt[i] += m_rotation[j][i] * 0.5/scales;
            }
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.5;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += scales * m_rotation[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    tmp = 0.;
    for (i = 0; i < m_dim - 1; i++)
    {
        F2 = 100. * (m_tmx[i] * m_tmx[i] - m_tmx[i+1]) * (m_tmx[i] * m_tmx[i] - m_tmx[i+1]) + (1 - m_tmx[i]) * (1 - m_tmx[i]);
        tmp += F2 / 4000. - cos(F2);
    }
    Ftrue = 1. + 1. * tmp / (double)(m_dim - 1);

    Fval = FGauss(Ftrue, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//F8F2 sum of Griewank-Rosenbrock 2-D blocks with uniform noise
bbob2015::TwoDoubles bbob2015::f126(const decision_vector &x) const
{
    unsigned int i, j; /*Loop over dim*/
    int rseed;

    static int funcId = 126;
    static int rrseed = 19;
    static double scales;
    double F2, Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);

        scales = fmax(1., sqrt((double)m_dim) / 8.);
        computeRotation(m_rotation, rseed, m_dim);

        /*compute Xopt*/
        for (i = 0; i < m_dim; i++)
        {
            m_Xopt[i] = 0.;
            for (j = 0; j < m_dim; j++)
            {
                m_Xopt[i] += m_rotation[j][i] * 0.5/scales;
            }
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.5;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += scales * m_rotation[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    tmp = 0.;
    for (i = 0; i < m_dim - 1; i++)
    {
        F2 = 100. * (m_tmx[i] * m_tmx[i] - m_tmx[i+1]) * (m_tmx[i] * m_tmx[i] - m_tmx[i+1]) + (1 - m_tmx[i]) * (1 - m_tmx[i]);
        tmp += F2 / 4000. - cos(F2);
    }
    Ftrue = 1. + 1. * tmp / (double)(m_dim - 1);

    Fval = FUniform(Ftrue, 0.49 + 1./m_dim, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//F8F2 sum of Griewank-Rosenbrock 2-D blocks with seldom Cauchy noise
bbob2015::TwoDoubles bbob2015::f127(const decision_vector &x) const
{
    unsigned int i, j; /*Loop over dim*/
    int rseed;

    static int funcId = 127;
    static int rrseed = 19;
    static double scales;
    double F2, Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);

        scales = fmax(1., sqrt((double)m_dim) / 8.);
        computeRotation(m_rotation, rseed, m_dim);

        /*compute Xopt*/
        for (i = 0; i < m_dim; i++)
        {
            m_Xopt[i] = 0.;
            for (j = 0; j < m_dim; j++)
            {
                m_Xopt[i] += m_rotation[j][i] * 0.5/scales;
            }
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.5;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += scales * m_rotation[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    tmp = 0.;
    for (i = 0; i < m_dim - 1; i++)
    {
        F2 = 100. * (m_tmx[i] * m_tmx[i] - m_tmx[i+1]) * (m_tmx[i] * m_tmx[i] - m_tmx[i+1]) + (1 - m_tmx[i]) * (1 - m_tmx[i]);
        tmp += F2 / 4000. - cos(F2);
    }
    Ftrue = 1. + 1. * tmp / (double)(m_dim - 1);

    Fval = FCauchy(Ftrue, 1., 0.2);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Gallagher with 101 Gaussian peaks with Gauss noise, condition up to 1000, one global m_rotation
bbob2015::TwoDoubles bbob2015::f128(const decision_vector &x) const
{
    unsigned int i, j, k; /*Loop over dim*/
    int rseed;

    static int funcId = 128;
    static int rrseed = 21;
    static double fitvalues[2] = {1.1, 9.1};
    static double maxcondition = 1000.;
    static double arrCondition[NHIGHPEAKS21];
    static double peakvalues[NHIGHPEAKS21];
    static double a = 0.1;
    double tmp2, f = 0., Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    double fac = -0.5 / (double)m_dim;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeRotation(m_rotation, rseed, m_dim);
        m_peaks = m_peaks21;
        unif(m_peaks, NHIGHPEAKS21 - 1, rseed);
        m_rperm = m_rperm21;
        for (i = 0; i < NHIGHPEAKS21 - 1; i++)
        m_rperm[i] = i;

        std::sort(m_rperm.begin(), m_rperm.begin() + NHIGHPEAKS21 - 1, compare_doubles(*this));

        /* Random permutation*/

        arrCondition[0] = sqrt(maxcondition);
        peakvalues[0] = 10;
        for (i = 1; i < NHIGHPEAKS21; i++)
        {
            arrCondition[i] = pow(maxcondition, (double)(m_rperm[i-1])/((double)(NHIGHPEAKS21-2)));
            peakvalues[i] = (double)(i-1)/(double)(NHIGHPEAKS21-2) * (fitvalues[1] - fitvalues[0]) + fitvalues[0];
        }
        m_arrScales = m_arrScales21;
        for (i = 0; i < NHIGHPEAKS21; i++)
        {
            unif(m_peaks, m_dim, rseed + 1000 * i);
            for (j = 0; j < m_dim; j++)
            m_rperm[j] = j;

            std::sort(m_rperm.begin(), m_rperm.begin() + m_dim, compare_doubles(*this));

            for (j = 0; j < m_dim; j++)
            {
                m_arrScales[i][j] = pow(arrCondition[i], ((double)m_rperm[j])/((double)(m_dim-1)) - 0.5);
            }
        }

        unif(m_peaks, m_dim * NHIGHPEAKS21, rseed);
        m_Xlocal = m_Xlocal21;
        for (i = 0; i < m_dim; i++)
        {
            /*compute Xopt*/
            m_Xopt[i] = 0.8 * (10. * m_peaks[i] -5.);
            for (j = 0; j < NHIGHPEAKS21; j++)
            {
                m_Xlocal[i][j] = 0.;
                for (k = 0; k < m_dim; k++)
                {
                    m_Xlocal[i][j] += m_rotation[i][k] * (10. * m_peaks[j * m_dim + k] -5.);
                }
                if (j == 0)
                    m_Xlocal[i][j] *= 0.8;
            }
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < NHIGHPEAKS21; i++)
    {
        tmp2 = 0.;
        for (j = 0; j < m_dim; j++)
        {
            tmp2 += m_arrScales[i][j] * (m_tmx[j] - m_Xlocal[j][i]) * (m_tmx[j] - m_Xlocal[j][i]);
        }
        tmp2 = peakvalues[i] * exp(fac * tmp2);
        f = fmax(f, tmp2);
    }

    f = 10 - f;
    /*monotoneTFosc*/
    if (f > 0)
    {
        Ftrue = log(f)/a;
        Ftrue = pow(exp(Ftrue + 0.49*(sin(Ftrue) + sin(0.79*Ftrue))), a);
    }
    else if (f < 0)
    {
        Ftrue = log(-f)/a;
        Ftrue = -pow(exp(Ftrue + 0.49*(sin(0.55 * Ftrue) + sin(0.31*Ftrue))), a);
    }
    else
        Ftrue = f;

    Ftrue *= Ftrue;

    Fval = FGauss(Ftrue, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Gallagher with 101 Gaussian peaks with uniform noise, condition up to 1000, one global m_rotation
bbob2015::TwoDoubles bbob2015::f129(const decision_vector &x) const
{
    unsigned int i, j, k; /*Loop over dim*/
    int rseed;

    static int funcId = 129;
    static int rrseed = 21;
    static double fitvalues[2] = {1.1, 9.1};
    static double maxcondition = 1000.;
    static double arrCondition[NHIGHPEAKS21];
    static double peakvalues[NHIGHPEAKS21];
    static double a = 0.1;
    double tmp2, f = 0., Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    double fac = -0.5 / (double)m_dim;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeRotation(m_rotation, rseed, m_dim);
        m_peaks = m_peaks21;
        unif(m_peaks, NHIGHPEAKS21 - 1, rseed);
        m_rperm = m_rperm21;
        for (i = 0; i < NHIGHPEAKS21 - 1; i++)
        m_rperm[i] = i;

        std::sort(m_rperm.begin(), m_rperm.begin() + NHIGHPEAKS21 - 1, compare_doubles(*this));

        /* Random permutation*/

        arrCondition[0] = sqrt(maxcondition);
        peakvalues[0] = 10;
        for (i = 1; i < NHIGHPEAKS21; i++)
        {
            arrCondition[i] = pow(maxcondition, (double)(m_rperm[i-1])/((double)(NHIGHPEAKS21-2)));
            peakvalues[i] = (double)(i-1)/(double)(NHIGHPEAKS21-2) * (fitvalues[1] - fitvalues[0]) + fitvalues[0];
        }
        m_arrScales = m_arrScales21;
        for (i = 0; i < NHIGHPEAKS21; i++)
        {
            unif(m_peaks, m_dim, rseed + 1000 * i);
            for (j = 0; j < m_dim; j++)
            m_rperm[j] = j;

            std::sort(m_rperm.begin(), m_rperm.begin() + m_dim, compare_doubles(*this));

            for (j = 0; j < m_dim; j++)
            {
                m_arrScales[i][j] = pow(arrCondition[i], ((double)m_rperm[j])/((double)(m_dim-1)) - 0.5);
            }
        }

        unif(m_peaks, m_dim * NHIGHPEAKS21, rseed);
        m_Xlocal = m_Xlocal21;
        for (i = 0; i < m_dim; i++)
        {
            /*compute Xopt*/
            m_Xopt[i] = 0.8 * (10. * m_peaks[i] -5.);
            for (j = 0; j < NHIGHPEAKS21; j++)
            {
                m_Xlocal[i][j] = 0.;
                for (k = 0; k < m_dim; k++)
                {
                    m_Xlocal[i][j] += m_rotation[i][k] * (10. * m_peaks[j * m_dim + k] -5.);
                }
                if (j == 0)
                    m_Xlocal[i][j] *= 0.8;
            }
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < NHIGHPEAKS21; i++)
    {
        tmp2 = 0.;
        for (j = 0; j < m_dim; j++)
        {
            tmp2 += m_arrScales[i][j] * (m_tmx[j] - m_Xlocal[j][i]) * (m_tmx[j] - m_Xlocal[j][i]);
        }
        tmp2 = peakvalues[i] * exp(fac * tmp2);
        f = fmax(f, tmp2);
    }

    f = 10 - f;
    /*monotoneTFosc*/
    if (f > 0)
    {
        Ftrue = log(f)/a;
        Ftrue = pow(exp(Ftrue + 0.49*(sin(Ftrue) + sin(0.79*Ftrue))), a);
    }
    else if (f < 0)
    {
        Ftrue = log(-f)/a;
        Ftrue = -pow(exp(Ftrue + 0.49*(sin(0.55 * Ftrue) + sin(0.31*Ftrue))), a);
    }
    else
        Ftrue = f;

    Ftrue *= Ftrue;

    Fval = FUniform(Ftrue, 0.49 + 1./m_dim, 1.);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

//Gallagher with 101 Gaussian peaks with Cauchy noise, condition up to 1000, one global m_rotation
bbob2015::TwoDoubles bbob2015::f130(const decision_vector &x) const
{
    unsigned int i, j, k; /*Loop over dim*/
    int rseed;

    static int funcId = 130;
    static int rrseed = 21;
    static double fitvalues[2] = {1.1, 9.1};
    static double maxcondition = 1000.;
    static double arrCondition[NHIGHPEAKS21];
    static double peakvalues[NHIGHPEAKS21];
    static double a = 0.1;
    double tmp2, f = 0., Fadd, Fval, tmp, Fpen = 0., Ftrue = 0.;
    double fac = -0.5 / (double)m_dim;
    bbob2015::TwoDoubles res;

    if (!m_isInitDone)
    {
        rseed = rrseed + 10000 * m_trialid;
        /*INITIALIZATION*/
        m_Fopt = computeFopt(funcId, m_trialid);
        computeRotation(m_rotation, rseed, m_dim);
        m_peaks = m_peaks21;
        unif(m_peaks, NHIGHPEAKS21 - 1, rseed);
        m_rperm = m_rperm21;
        for (i = 0; i < NHIGHPEAKS21 - 1; i++)
        m_rperm[i] = i;

        std::sort(m_rperm.begin(), m_rperm.begin() + NHIGHPEAKS21 - 1, compare_doubles(*this));

        /* Random permutation*/

        arrCondition[0] = sqrt(maxcondition);
        peakvalues[0] = 10;
        for (i = 1; i < NHIGHPEAKS21; i++)
        {
            arrCondition[i] = pow(maxcondition, (double)(m_rperm[i-1])/((double)(NHIGHPEAKS21-2)));
            peakvalues[i] = (double)(i-1)/(double)(NHIGHPEAKS21-2) * (fitvalues[1] - fitvalues[0]) + fitvalues[0];
        }
        m_arrScales = m_arrScales21;
        for (i = 0; i < NHIGHPEAKS21; i++)
        {
            unif(m_peaks, m_dim, rseed + 1000 * i);
            for (j = 0; j < m_dim; j++)
            m_rperm[j] = j;

            std::sort(m_rperm.begin(), m_rperm.begin() + m_dim, compare_doubles(*this));

            for (j = 0; j < m_dim; j++)
            {
                m_arrScales[i][j] = pow(arrCondition[i], ((double)m_rperm[j])/((double)(m_dim-1)) - 0.5);
            }
        }

        unif(m_peaks, m_dim * NHIGHPEAKS21, rseed);
        m_Xlocal = m_Xlocal21;
        for (i = 0; i < m_dim; i++)
        {
            /*compute Xopt*/
            m_Xopt[i] = 0.8 * (10. * m_peaks[i] -5.);
            for (j = 0; j < NHIGHPEAKS21; j++)
            {
                m_Xlocal[i][j] = 0.;
                for (k = 0; k < m_dim; k++)
                {
                    m_Xlocal[i][j] += m_rotation[i][k] * (10. * m_peaks[j * m_dim + k] -5.);
                }
                if (j == 0)
                    m_Xlocal[i][j] *= 0.8;
            }
        }
        m_isInitDone = 1;
    }
    Fadd = m_Fopt;

    /* BOUNDARY HANDLING*/
    for (i = 0; i < m_dim; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
    Fadd += 100. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < m_dim; i++)
    {
        m_tmx[i] = 0.;
        for (j = 0; j < m_dim; j++)
        {
            m_tmx[i] += m_rotation[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < NHIGHPEAKS21; i++)
    {
        tmp2 = 0.;
        for (j = 0; j < m_dim; j++)
        {
            tmp2 += m_arrScales[i][j] * (m_tmx[j] - m_Xlocal[j][i]) * (m_tmx[j] - m_Xlocal[j][i]);
        }
        tmp2 = peakvalues[i] * exp(fac * tmp2);
        f = fmax(f, tmp2);
    }

    f = 10 - f;
    /*monotoneTFosc*/
    if (f > 0)
    {
        Ftrue = log(f)/a;
        Ftrue = pow(exp(Ftrue + 0.49*(sin(Ftrue) + sin(0.79*Ftrue))), a);
    }
    else if (f < 0)
    {
        Ftrue = log(-f)/a;
        Ftrue = -pow(exp(Ftrue + 0.49*(sin(0.55 * Ftrue) + sin(0.31*Ftrue))), a);
    }
    else
        Ftrue = f;

    Ftrue *= Ftrue;

    Fval = FCauchy(Ftrue, 1., 0.2);

    Ftrue += Fadd;
    Fval += Fadd;

    res.Fval = Fval;
    res.Ftrue = Ftrue;
    return res;
}

void bbob2015::set_seed(unsigned int seed) const {
    m_seed = seed;
    m_seedn = seed;
    // As the problem is now muted we must reset the caches that contain the evaluations w.r.t. the old seed
    reset_caches();
}

base_ptr bbob2015::clone() const
{
    return base_ptr(new bbob2015(*this));
}
}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::bbob2015)
