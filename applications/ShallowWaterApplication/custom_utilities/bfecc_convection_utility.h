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

template<std::size_t TDim> class BFECCConvectionUtility
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(BFECCConvectionUtility<TDim>);


    BFECCConvectionUtility(ModelPart& rThisModelPart, Parameters ThisParameters = Parameters());


    ~BFECCConvectionUtility() = default;


    template<class TVarType, class TType>
    void Convect(const TVarType& rVar, const Variable<array_1d<double,3>>& conv_var);


    void UpdateSearchDatabase();


    template<class TVarType>
    void ResetBoundaryConditions(const TVarType& rVar);


    template<class TVarType>
    void CopyVariableToPreviousTimeStep(const TVarType& rVar);


private:

    ModelPart& mrModelPart;
    BinBasedFastPointLocator<TDim> mSearchStructure;
    int mMaxResults;


    bool RK2Convect(
        const double Dt,
        array_1d<double,3>& rPosition,
        const array_1d<double,3>& rInitialVelocity,
        Vector& rN,
        Element::Pointer& pElement,
        typename BinBasedFastPointLocator<TDim>::ResultIteratorType& rResultBegin,
        const int VelocitySign,
        const Variable<array_1d<double,3> >& rConvVar);

};

} // namespace Kratos.

#endif // KRATOS_BFECC_CONVECTION_UTILITY_H_INCLUDED  defined