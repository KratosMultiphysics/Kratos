// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_FLUX_CONDITION_H_INCLUDED )
#define  KRATOS_FLUX_CONDITION_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "includes/condition.h"
#include "geometries/geometry.h"
#include "includes/variables.h"

namespace Kratos
{

///@addtogroup ConvectionDiffusionApplication
///@{

namespace FluxConditionInternals
{
///@name Auxiliary data structure to hold FEM data
///@{

template< unsigned int TNodeNumber >
class IntegrationData
{
public:

    IntegrationData(
        Geometry< Node<3> >& rGeometry,
        const Variable<double>& rFluxVar
        ):
        mGaussPoint(0),
        mNodalFluxes(TNodeNumber,0.0)
    {
        NumGauss = rGeometry.IntegrationPointsNumber(GeometryData::GI_GAUSS_2);
        Vector DetJ = ZeroVector(NumGauss);
        rGeometry.DeterminantOfJacobian(DetJ,GeometryData::GI_GAUSS_2);

        mShapeFunctionValues.resize(NumGauss,TNodeNumber);
        mIntegrationWeights.resize(NumGauss);
        
        noalias(mShapeFunctionValues) = rGeometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
            
        const auto& IntegrationPoints = rGeometry.IntegrationPoints(GeometryData::GI_GAUSS_2);
            
        for (unsigned int g = 0; g < NumGauss; g++)
        {
            mIntegrationWeights[g] = DetJ[g]*IntegrationPoints[g].Weight();
        }

        for (unsigned int i = 0; i < TNodeNumber; i++)
        {
            mNodalFluxes[i] = rGeometry[i].FastGetSolutionStepValue(rFluxVar);
        }
    }

    void SetCurrentGaussPoint(unsigned int g)
    {
        mGaussPoint = g;
    }

    double N(unsigned int i) const
    {
        return mShapeFunctionValues(mGaussPoint,i);
    }

    double NodalFlux(unsigned int i) const
    {
        return mNodalFluxes[i];
    }

    double GaussPointFlux() const
    {
        double flux = mNodalFluxes[0]*mShapeFunctionValues(mGaussPoint,0);
        for (unsigned int i = 1; i < TNodeNumber; i++)
        {
            flux += mNodalFluxes[i]*mShapeFunctionValues(mGaussPoint,i);
        }
        return flux;
    }

    double IntegrationWeight() const
    {
        return mIntegrationWeights[mGaussPoint];
    }

    unsigned int NumGauss;

private:

    unsigned int mGaussPoint;

    array_1d<double,TNodeNumber> mNodalFluxes;

    Matrix mShapeFunctionValues;

    Vector mIntegrationWeights;
};

///@}

}

///@}
///@name Kratos Classes
///@{

/// A basic Neumann condition for convection-diffusion problems.
/** It applies a flux condition based on the nodal values of the
 *  variable defined as SurfaceSourceVariable by the
 *  CONVECTION_DIFFUSION_SETTINGS variable in the given ProcessInfo.
 */
template< unsigned int TNodeNumber >
class FluxCondition: public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FluxCondition
    KRATOS_CLASS_POINTER_DEFINITION(FluxCondition);

    typedef Condition::MatrixType MatrixType;
    typedef Condition::VectorType VectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    FluxCondition(
        IndexType NewId,
        Geometry< Node<3> >::Pointer pGeometry);

    FluxCondition(
        IndexType NewId,
        Geometry< Node<3> >::Pointer pGeometry,
        Properties::Pointer pProperties);

    /// Destructor.
    ~FluxCondition() override;
    
    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties) const override;
        
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(
        DofsVectorType& ConditionalDofList,
        ProcessInfo& CurrentProcessInfo) override;

    GeometryData::IntegrationMethod GetIntegrationMethod() override;

    void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                                             std::vector<array_1d<double, 3 > >& rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                             std::vector<double>& rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                                             std::vector<array_1d<double, 6 > >& rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                             std::vector<Vector>& rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                             std::vector<Matrix>& rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;
    
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

protected:
    
    ///@name Protected Operations
    ///@{

    void AddIntegrationPointRHSContribution(
        VectorType& rRightHandSideVector,
        const FluxConditionInternals::IntegrationData<TNodeNumber>& rData);

    void CalculateNormal(array_1d<double,3>& An );

    ///@}

private:
      
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    FluxCondition();

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;
    
    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    FluxCondition& operator=(FluxCondition const& rOther);

    /// Copy constructor.
    FluxCondition(FluxCondition const& rOther);

    ///@}

}; // Class FluxCondition

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FLUX_CONDITION_H_INCLUDED  defined
