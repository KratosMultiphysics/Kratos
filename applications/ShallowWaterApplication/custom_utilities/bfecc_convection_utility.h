//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


#ifndef KRATOS_BFECC_CONVECTION_UTILITY_H_INCLUDED
#define KRATOS_BFECC_CONVECTION_UTILITY_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/variables.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/binbased_fast_point_locator.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

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

/// Explicit convection utility
/** Convection of scalars and vectors for shallow water equations using BFECC correction
*/
template<std::size_t TDim> class BFECCConvectionUtility
{

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(BFECCConvectionUtility<TDim>);

    ///@}
    ///@name Life Cycle
    ///@{

    BFECCConvectionUtility(ModelPart& rThisModelPart, Parameters ThisParameters = Parameters());


    ~BFECCConvectionUtility() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method sets at Tn+1 the material value from Tn, using the velocity variable to convect
     * @param rVar The variable to convect
     * @param rConvVar The velocity variable
     */
    template<class TVarType, class TType>
    void Convect(const TVarType& rVar, const Variable<array_1d<double,3>>& rConvVar);

    /**
     * @brief This method updates the search structure if the mesh has been modified
     */
    void UpdateSearchDatabase();

    /**
     * @brief This method copy the value from Tn to Tn+1 if the variable is fixed
     * @param The variable to reset if is fixed
     */
    template<class TVarType>
    void ResetBoundaryConditions(const TVarType& rVar);

    /**
     * This method copies the variable from Tn+1 to Tn
     * @param rVar The variable to copy to the previous time step
     */
    template<class TVarType>
    void CopyVariableToPreviousTimeStep(const TVarType& rVar);

    ///@}
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

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    BinBasedFastPointLocator<TDim> mSearchStructure;
    int mMaxResults;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    bool RK2Convect(
        const double Dt,
        array_1d<double,3>& rPosition,
        const array_1d<double,3>& rInitialVelocity,
        Vector& rN,
        Element::Pointer& pElement,
        typename BinBasedFastPointLocator<TDim>::ResultIteratorType& rResultBegin,
        const int VelocitySign,
        const Variable<array_1d<double,3> >& rConvVar);

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

};

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_BFECC_CONVECTION_UTILITY_H_INCLUDED  defined