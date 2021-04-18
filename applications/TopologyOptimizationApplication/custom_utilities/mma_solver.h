// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

#if !defined(KRATOS_MMA_SOLVER_H_INCLUDED)
#define  KRATOS_MMA_SOLVER_H_INCLUDED


#include <iostream>
#include <string>
#include <algorithm>
#include <iomanip>   


// External includes
#include <pybind11/pybind11.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
#include "containers/array_1d.h"
#include <vector>

// Application includes
#include "topology_optimization_application.h"
#include "spatial_containers/spatial_containers.h" // For kd-tree
#include "custom_utilities/filter_function.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{


///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Solution utility to filter results.
/** Detail class definition.

 */

class MMACalculateNewDensities
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(MMACalculateNewDensities);

    /// Default constructor.
    MMACalculateNewDensities(std::string FilterFunctionType)
    {
        		// Set precision for output
		std::cout.precision(12);

		// Set type of weighting function

		// Type 1: Gaussian function
		std::string linear("linear");
		if(FilterFunctionType.compare(linear)==0)
			mFilterFunctionType = 1;
    }

    /// Destructor.
    virtual ~MMACalculateNewDensities()
    {
    }
    
    ///@name Operators
	///@{

        int n = 10;
        int m= 2;
        int iter = 0;
        double xmamieps = 0.000001;
        double epsimin = std::sqrt(n+m)*0.0000000001;
        double raa0 = 0.00001;
        double albefa = 0.1;
        double asyminit = 0.5;
	    double asymdec=0.7;
	    double asyminc =1.2;
        double move = 0.5;
        double a = 1;
	    double c = 1;
	    double d = 1;
	    double y = 1;
        double z = 1;
	    double lam = 1;
	    double mu = 1; 
        double s = 1;
	    double low = 1;
	    double upp = 1;
	    double alpha = 1;
	    double beta = 1;
	    double p0 = 1;
	    double q0 = 1;
	    double pij = 1;
	    double qij = 1;
        double b = 1;
	    double grad = 1;
	    double hess = 1;
/* 	    double xold1 = 0;
	    double xold2 = 1; */


	///@}
	///@name Operations

    double UpdateMMA(double xval, const double dfdx, const double gx, const double dgdx, const double xmin, const double xmax, const double xold1, const double xold2)
    {
        KRATOS_TRY;
        std::cout << "  Hallo aus dem Update"<< std::endl;

        // Generate the subproblem
        GenSub(xval, dfdx, gx, dgdx, xmin, xmax, xold1, xold2);

        // Update xolds
/*         double old = xold2;
        xold2 = xold1; */
/*         std::copy_n(xval, n, xold1.data()); */


        // Solve the dual with an interior point method
        SolveDIP(xval);

        // Solve the dual with a steepest ascent method
        // SolveDSA(xval);
        return xval;


        KRATOS_CATCH("");
    }
    	///@name Access
	///@{


	///@}
	///@name Inquiry
	///@{


	///@}
	///@name Input and output
	///@{


	///@}
	///@name Friends
	///@{


	///@}

	protected:
	///@name Protected static Member Variables
	///@{


	///@}
	///@name Protected member Variables
	///@{


	///@}
	///@name Protected Operators
	///@{


	///@}
	///@name Protected Operations
	///@{


	///@}
	///@name Protected  Access
	///@{


	///@}
	///@name Protected Inquiry
	///@{


	///@}
	///@name Protected LifeCycle
	///@{


	///@}

	private:
	///@name Static Member Variables
	///@{

    void SolveDIP(double x) 
    {
        KRATOS_TRY
        std::cout << "  Hallo aus dem Solve"<< std::endl;

        
        lam = c / 2.0;
        mu = 1.0;

        const double tol = epsimin; // 1.0e-9*sqrt(m+n);
        double epsi = 1.0;
        double err = 1.0;
        int loop;

        while (epsi > tol) 
        {

            loop = 0;
            while (err > 0.9 * epsi && loop < 100) 
            {
                loop++;

                // Set up Newton system
                XYZofLAMBDA(x);
                DualGrad(x);
                grad = -1.0 * grad - epsi / lam;
                DualHess(x);

                // Solve Newton system
                if (m > 1) {
                    /* Factorize(hess.data(), m);
                    Solve(hess.data(), grad.data(), m); */

                    s = grad;
                    
                } else if (m > 0) {
                    s = grad / hess;
                }

                // Get the full search direction
                s = -mu + epsi / lam - s * mu / lam;

                // Perform linesearch and update lam and mu
                DualLineSearch();

                XYZofLAMBDA(x);

                // Compute KKT res
                err = DualResidual(x, epsi);
            }
            epsi = epsi * 0.1;
        }
        
        KRATOS_CATCH("");
    }


    void SolveDSA(double x) 
    {
        KRATOS_TRY

        lam = 1.0;

        const double tol = epsimin; // 1.0e-9*sqrt(m+n);
        double err = 1.0;
        int loop = 0;

        while (err > tol && loop < 500) 
        {
            loop++;
            XYZofLAMBDA(x);
            DualGrad(x);
            double theta = 1.0;
            err = 0.0;
            lam = std::max(0.0, lam + theta * grad);
            err += grad * grad;
            
            err = std::sqrt(err);
        }
        
        KRATOS_CATCH("");
        
    }

    double DualResidual(double x, double epsi) 
    {
        KRATOS_TRY

        double res = 0;


        res = -b - a * z - y + mu;


/*             res[j + m] = mu[j] * lam[j] - epsi;
            for (int i = 0; i < n; i++) {
                res[j] += pij[i * m + j] / (upp[i] - x[i]) + qij[i * m + j] / (x[i] - low[i]);
            } */
        

        double nrI = 0.0;
        if (nrI < std::abs(res)) 
        {
            nrI = std::abs(res);
        }

/*         delete[] res; */

        return nrI;

        KRATOS_CATCH("");
    }

    void DualLineSearch() 
    {

        KRATOS_TRY

        double theta = 1.005;

        if (theta < -1.01 * s / lam) 
        {
            theta = -1.01 * s / lam;
        }
            
        if (theta < -1.01 * s / mu) 
        {
            theta = -1.01 * s / mu;
        }
        theta = 1.0 / theta;


        lam = lam + theta * s;
        mu = mu + theta * s;

        KRATOS_CATCH("");
        
    }


    void DualHess(double x) 
    {
        KRATOS_TRY

        double df2 = 0;
        double PQ = 0;
        #ifdef MMA_WITH_OPENMP
        #pragma omp parallel for
        #endif
        double pjlam = p0;
        double qjlam = q0;
        pjlam += pij * lam;
        qjlam += qij * lam;
        PQ = pij / pow(upp - x, 2.0) - qij / pow(x - low, 2.0);
        
        df2 = -1.0 / (2.0 * pjlam / pow(upp - x, 3.0) + 2.0 * qjlam / pow(x - low, 3.0));
        double xp = (sqrt(pjlam) * low + sqrt(qjlam) * upp) / (sqrt(pjlam) + sqrt(qjlam));
        if (xp < alpha) 
        {
            df2 = 0.0;
        }
        if (xp > beta) 
        {
            df2 = 0.0;
        }
        

        // Create the matrix/matrix/matrix product: PQ^T * diag(df2) * PQ
        double tmp = 0;
        #ifdef MMA_WITH_OPENMP
        #pragma omp parallel for
        #endif
        tmp = 0.0;
        tmp += PQ  * df2;
    
        hess = 0.0;
        hess += tmp * PQ;
            
            
        

        double lamai = 0.0;

        if (lam < 0.0) 
        {
            lam = 0.0;
        }
        lamai += lam * a;

        if (lam > c) 
        {
            hess += -1.0;
        }
        hess += -mu / lam;


        if (lamai > 0.0)
        {
            hess += -10.0 * a * a;
                
        }

        // pos def check
        double HessTrace = 0.0;
        HessTrace += hess;
        
        double HessCorr = 1e-4 * HessTrace / m;

        if (-1.0 * HessCorr < 1.0e-7) 
        {
            HessCorr = -1.0e-7;
        }


        hess += HessCorr;


/*         delete df2;
        delete PQ;
        delete tmp; */

        KRATOS_CATCH("");
    }

    void DualGrad(double x) 
    {
        KRATOS_TRY

        grad = -b - a * z - y;
        grad += pij / (upp - x) + qij / (x - low);

        KRATOS_CATCH("");
            
    }

    void XYZofLAMBDA(double x) 
    
    {
        KRATOS_TRY

        double lamai = 0.0;
        if (lam < 0.0) 
        {
            lam = 0;
        }
        y = std::max(0.0, lam - c); // Note y=(lam-c)/d - however d is fixed at one !!
            lamai += lam * a;
        
        z = std::max(0.0, 10.0 * (lamai - 1.0)); // SINCE a0 = 1.0

        #ifdef MMA_WITH_OPENMP
        #pragma omp parallel for
        #endif
            double pjlam = p0;
            double qjlam = q0;

            pjlam += pij * lam;
            qjlam += qij * lam;
            
            x = (sqrt(pjlam) * low + sqrt(qjlam) * upp) / (sqrt(pjlam) + sqrt(qjlam));
            if (x < alpha) 
            {
                x = alpha;
            }
            if (x > beta) 
            {
                x = beta;
            }
        

        KRATOS_CATCH("");
    }
    
    
    void GenSub(double xval, const double dfdx, const double gx, const double dgdx, const double xmin, const double xmax, const double xold1, const double xold2)
    {
        KRATOS_TRY
            // Forward the iterator
            iter++;
            std::cout << "  Hallo aus dem Subproblem"<< std::endl;


            // Set asymptotes
            if (iter < 3) {
                #ifdef MMA_WITH_OPENMP
                #pragma omp parallel for
                #endif
                low = xval - asyminit * (xmax - xmin);
                upp = xval + asyminit * (xmax - xmin);
    
            } else {
                #ifdef MMA_WITH_OPENMP
                #pragma omp parallel for
                #endif
                double zzz = (xval - xold1) * (xold1 - xold2);
                double gamma;
                    if (zzz < 0.0) {
                        gamma = asymdec;
                    } else if (zzz > 0.0) {
                        gamma = asyminc;
                    } else {
                        gamma = 1.0;
                    }
                    low = xval - gamma * (xold1 - low);
                    upp= xval + gamma * (upp - xold1);

                    double xmami = std::max(xmamieps, xmax - xmin);
                    // double xmami = xmax[i] - xmin[i];
                    low = std::max(low, xval - 100.0 * xmami);
                    low = std::min(low, xval - 1.0e-5 * xmami);
                    upp = std::max(upp, xval + 1.0e-5 * xmami);
                    upp = std::min(upp, xval + 100.0 * xmami);

                    double xmi = xmin - 1.0e-6;
                    double xma = xmax + 1.0e-6;
                    if (xval < xmi) {
                        low = xval - (xma - xval) / 0.9;
                        upp = xval + (xma - xval) / 0.9;
                    }
                    if (xval > xma) {
                        low = xval - (xval - xmi) / 0.9;
                        upp = xval + (xval - xmi) / 0.9;
                    }
                }

            // Set bounds and the coefficients for the approximation
            // double raa0 = 0.5*1e-6;
            #ifdef MMA_WITH_OPENMP
            #pragma omp parallel for
            #endif
            // Compute bounds alpha and beta
            alpha = std::max(xmin, low + albefa * (xval - low));
            alpha = std::max(alpha, xval - move * (xmax - xmin));
            alpha = std::min(alpha, xmax);
            beta = std::min(xmax, upp - albefa * (upp - xval));
            beta = std::min(beta, xval + move * (xmax - xmin));
            beta = std::max(beta, xmin);

            // Objective function
            {
                double dfdxp = std::max(0.0, dfdx);
                double dfdxm = std::max(0.0, -1.0 * dfdx);
                double xmamiinv = 1.0 / std::max(xmamieps, xmax - xmin);
                double pq = 0.001 * std::abs(dfdx) + raa0 * xmamiinv;
                p0 = std::pow(upp - xval, 2.0) * (dfdxp + pq);
                q0 = std::pow(xval - low, 2.0) * (dfdxm + pq);
            }

            // Constraints
                double dgdxp = std::max(0.0, dgdx);
                double dgdxm = std::max(0.0, -1.0 * dgdx);
                double xmamiinv = 1.0 / std::max(xmamieps, xmax - xmin);
                double pq = 0.001 * std::abs(dgdx) + raa0 * xmamiinv;
                pij = std::pow(upp - xval, 2.0) * (dgdxp + pq);
                qij= std::pow(xval - low, 2.0) * (dgdxm + pq);
        
            

            // The constant for the constraints
                b = -gx;
                b += pij / (upp - xval) + qij / (xval - low);

        KRATOS_CATCH("");
    }
/* 
    void Factorize(double *K, int n) 
    {

            for (int s = 0; s < n - 1; s++) {
                for (int i = s + 1; i < n; i++) {
                    K[i * n + s] = K[i * n + s] / K[s * n + s];
                    for (int j = s + 1; j < n; j++) {
                        K[i * n + j] = K[i * n + j] - K[i * n + s] * K[s * n + j];
                    }
                }
            }
    } */

/*     void Solve(double K, double x, int n) 
    {


        double a = 0.0;

        a = a - K * x;
            
        x = x + a;
        

      x[n - 1] = x[n - 1] / K[(n - 1) * n + (n - 1)]; 

        for (int i = n - 2; i >= 0; i--) {
            double a = x[i];
            for (int j = i + 1; j < n; j++) {
                a = a - K[i * n + j] * x[j];
            }
            x[i] = a / K[i * n + i];
        }
    }
 */

	///@}
	///@name Member Variables
	///@{
        	
	unsigned int mFilterFunctionType;


	///@}
	///@name Private Operators
	///@{


	///@}
	///@name Private Operations
	///@{


	///@}
	///@name Private  Access
	///@{


	///@}
	///@name Private Inquiry
	///@{


	///@}
	///@name Un accessible methods
	///@{

	/// Assignment operator.
	//FilterFunction& operator=(FilterFunction const& rOther);

	/// Copy constructor.
	//FilterFunction(FilterFunction const& rOther);


	///@}

}; // Class TopologyFilteringUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif	/* KRATOS_MMA_SOLVER_H_INCLUDED */