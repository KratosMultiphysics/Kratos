//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela, Ruben Zorrilla
//
//

#ifndef KRATOS_TAIT_EQUATION_UTILITIES_H
#define	KRATOS_TAIT_EQUATION_UTILITIES_H

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "custom_utilities/tait_equation_utility.h"

namespace Kratos
{

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

///This class provides the utilities needed for the implementation of the Tait model
class TaitEquationUtility
{
public:

    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    
    /** 
     * This function calculates the specific volume given the coefficients vectors of the Tait Equation
     * @param bm  is a vector of size 10 containing the coefficients of the tait equation for the molten state: 
     *            bm = [0.0,b1m, b2m, b3m, b4m, b5, b6, b7, b8, b9 ]
     * @param bs  is a vector of size 10 containing the coefficients of the tait equation for the molten state: 
     *            [0.0,b1s, b2s, b3s, b4s, b5, b6, b7, b8, b9 ]
     * @param T is the temperature in Kelvin
     * @param P is the ABSOLUTE pressure in Pascals
    */
    static double CalculateV(const Vector& bm, const Vector& bs, const double T, const double P)
    {
        KRATOS_DEBUG_ERROR_IF(bm.size() != 10) << "expected size for TAIT_PARAMETERS_MOLTEN_STATE is 10. Got instead a size of " << bm.size() << " provided vector is " << bm << std::endl;
        KRATOS_DEBUG_ERROR_IF(bs.size() != 10) << "expected size for TAIT_PARAMETERS_SOLID_STATE  is 10. Got instead a size of " << bs.size() << " provided vector is " << bs << std::endl;
        KRATOS_DEBUG_ERROR_IF(T < 0.0) << "absolute temperature cannot be below the absolute zero. Current value of T = " << T << std::endl;
//         KRATOS_DEBUG_ERROR_IF(P < 0.0) << "absolute pressure cannot be below zero. Current value of P =  " << P << std::endl;

        const double Tt = bm[5] + bm[6]*P; //glass transition gradient 
        
        if(T > Tt )
        {
            const double V0 = bm[1] + bm[2]*(T-bm[5]);
            const double Vt = 0.0;
            const double B  = bm[3]+bm[4]*(T-bm[5]);
            return V0*(1.0 - 0.0894*std::log(1.0+P/B)) - Vt;
            
        }
        else
        {
            const double V0 = bs[1] + bs[2]*(T-bs[5]);
            const double Vt = bs[7]*std::exp(bs[8]*(T-bs[5]) - bs[9]*P);
            const double B  = bs[3]+bs[4]*(T-bs[5]);
            return V0*(1.0 - 0.0894*std::log(1.0+P/B)) - Vt;
        }
    }
 
    /** 
     * This function calculates the density using the Tait Equation
     * @param bm  is a vector of size 10 containing the coefficients of the tait equation for the molten state: 
     *            bm = [0.0,b1m, b2m, b3m, b4m, b5, b6, b7, b8, b9 ]
     * @param bs  is a vector of size 10 containing the coefficients of the tait equation for the molten state: 
     *            [0.0,b1s, b2s, b3s, b4s, b5, b6, b7, b8, b9 ]
     * @param T is the temperature in Kelvin
     * @param P is the ABSOLUTE pressure in Pascals
    */
    static double CalculateRho(const Vector& bm, const Vector& bs, const double T, const double P)
    {
        const double V = CalculateV(bm,bs,T,P);
        return 1.0/V;
    }
    
    /** 
     * This function calculates the specific volume given the coefficients vectors of the Tait Equation
     * @param bm  is a vector of size 10 containing the coefficients of the tait equation for the molten state: 
     *            bm = [0.0,b1m, b2m, b3m, b4m, b5, b6, b7, b8, b9 ]
     * @param bs  is a vector of size 10 containing the coefficients of the tait equation for the molten state: 
     *            [0.0,b1s, b2s, b3s, b4s, b5, b6, b7, b8, b9 ]
     * @param T is the temperature in Kelvin
     * @param P is the ABSOLUTE pressure in Pascals
     * @param dP is the pressure increment used in the calculation of the bulk_modulus by perturbation
    */
    static double CalculateBulkModulus(const Vector& bm, const Vector& bs, const double T, const double P, const double dP)
    {
        const double rho1 = CalculateRho(bm,bs,T,P);
        const double rho2 = CalculateRho(bm,bs,T,P+dP);
        
        const double bulk_modulus = dP/(rho2-rho1); //-(V*dP)/(Vperturbed-V); //note that we define the bulk modulus as k=Dp/Drho (that is, the inverse of what is on the paper)
        return bulk_modulus;
    }
    ///@} // Operators

private:

    ///@name Auxiliary Data types
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ///@} // Member variables
    ///@name Private Operations
    ///@{

    ///@} // Private Operations

};

///@} // Kratos classes

///@}

} // namespace Kratos.


#endif	/* KRATOS_TAIT_EQUATION_UTILITIES_H */
