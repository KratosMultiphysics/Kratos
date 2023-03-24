// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Mohamed Nabi
//
//
//


#pragma once

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/T_condition.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/condition_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template<unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) TMicroClimateFluxCondition : public TCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TMicroClimateFluxCondition);
    
    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using TCondition<TDim,TNumNodes>::mThisIntegrationMethod;
    
    bool mIsInitialised = false;
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    TMicroClimateFluxCondition() : TCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    TMicroClimateFluxCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : TCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    TMicroClimateFluxCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : TCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

    // Destructor
    ~TMicroClimateFluxCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;
 
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct ElementVariables
    {
        double albedoCoefficient;
        double firstCoverStorageCoefficient;
        double secondCoverStorageCoefficient;
        double thirdCoverStorageCoefficient;
        double buildEnvironmentRadiation;
        double minimalStorage;
        double maximalStorage;

        double NormalFlux;
        double IntegrationCoefficient;
        array_1d<double,TNumNodes> Np;
        array_1d<double,TNumNodes> TVector;
        double roughnessTemperature = 0.0;
        double netRadiation = 0.0;
        double waterStorage = 0.0;

        double previousRoughnessTemperature;
        double previousStorage;
        double previousRadiation;

        array_1d<double, TNumNodes> leftHandSideFlux;
        array_1d<double, TNumNodes> rightHandSideFlux;
        //Matrix GradNpT;
        Vector detJContainer;
        Matrix NContainer;
        //GeometryType::ShapeFunctionsGradientsType DN_DXContainer;
        // needed for updated Lagrangian:
        double detJ;
        BoundedMatrix<double, TNumNodes, TNumNodes> TMatrix;
        array_1d<double, TNumNodes> righthandside;
        array_1d<double, TNumNodes> lefthandside;
    };
    
    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                
    void Initialize(const ProcessInfo& rCurrentProcessInfo);

    void InitializeElementVariables(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo);

    //void CalculateRHS(VectorType& rRightHandSideVector,
    //                  const ProcessInfo& CurrentProcessInfo) override;
    
    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);


    //virtual void CalculateIntegrationCoefficient(double& rIntegrationCoefficient,
    //    const Matrix& Jacobian,
    //    const double& Weight);

    virtual void CalculateIntegrationCoefficient(double& rIntegrationCoefficient,
        const Matrix& Jacobian,
        const double& Weight);

    void CalculateRoughness(const ProcessInfo& CurrentProcessInfo,
        ElementVariables& rVariables);

    void CalculateNodalFluxes(const ProcessInfo& CurrentProcessInfo,
        ElementVariables& rVariables);

    //void CalculateAll(MatrixType& rLeftHandSideMatrix,
    //    VectorType& rRightHandSideVector,
    //    const ProcessInfo& CurrentProcessInfo);

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo);

    void InitializeProperties(ElementVariables& rVariables);

    void CalculateKinematics(ElementVariables& rVariables, unsigned int PointNumber);

    double CalculateIntegrationCoefficient(const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        unsigned int PointNumber, double detJ);
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Serialization
    
    friend class Serializer;
    
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
    }
    
}; // class TMicroClimateFluxCondition.

} // namespace Kratos.