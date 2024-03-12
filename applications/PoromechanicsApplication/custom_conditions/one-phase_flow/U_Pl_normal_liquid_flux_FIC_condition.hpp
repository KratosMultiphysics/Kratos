//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


#if !defined(KRATOS_U_PL_NORMAL_LIQUID_FLUX_FIC_CONDITION_H_INCLUDED )
#define  KRATOS_U_PL_NORMAL_LIQUID_FLUX_FIC_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/one-phase_flow/U_Pl_condition.hpp"
#include "custom_conditions/one-phase_flow/U_Pl_normal_liquid_flux_condition.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPlNormalLiquidFluxFICCondition : public UPlNormalLiquidFluxCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPlNormalLiquidFluxFICCondition );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPlCondition<TDim,TNumNodes>::mThisIntegrationMethod;
    typedef typename UPlNormalLiquidFluxCondition<TDim,TNumNodes>::NormalLiquidFluxVariables NormalLiquidFluxVariables;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPlNormalLiquidFluxFICCondition() : UPlNormalLiquidFluxCondition<TDim,TNumNodes>() {}

    // Constructor 1
    UPlNormalLiquidFluxFICCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPlNormalLiquidFluxCondition<TDim,TNumNodes>(NewId, pGeometry) {}

    // Constructor 2
    UPlNormalLiquidFluxFICCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPlNormalLiquidFluxCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties)
    {
        mThisIntegrationMethod = this->GetIntegrationMethod();
    }

    // Destructor
    ~UPlNormalLiquidFluxFICCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct NormalLiquidFluxFICVariables
    {
        double DtPressureCoefficient;
        double ElementLength;
        double BiotModulusInverse;

        array_1d<double,TNumNodes> DtPressureVector;
        BoundedMatrix<double,TNumNodes,TNumNodes> PMatrix;
    };

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateElementLength(double& rElementLength, const GeometryType& Geom);


    void CalculateAndAddLHSStabilization(MatrixType& rLeftHandSideMatrix, NormalLiquidFluxVariables& rVariables, NormalLiquidFluxFICVariables& rFICVariables);

    void CalculateAndAddBoundaryMassMatrix(MatrixType& rLeftHandSideMatrix, NormalLiquidFluxVariables& rVariables, NormalLiquidFluxFICVariables& rFICVariables);


    void CalculateAndAddRHSStabilization(VectorType& rRightHandSideVector, NormalLiquidFluxVariables& rVariables, NormalLiquidFluxFICVariables& rFICVariables);

    void CalculateAndAddBoundaryMassFlow(VectorType& rRightHandSideVector, NormalLiquidFluxVariables& rVariables, NormalLiquidFluxFICVariables& rFICVariables);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        typedef UPlNormalLiquidFluxCondition<TDim,TNumNodes> BaseCondition;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseCondition )
    }

    void load(Serializer& rSerializer) override
    {
        typedef UPlNormalLiquidFluxCondition<TDim,TNumNodes> BaseCondition;
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseCondition )
    }

}; // class UPlNormalLiquidFluxFICCondition.

} // namespace Kratos.

#endif // KRATOS_U_PL_NORMAL_LIQUID_FLUX_FIC_CONDITION_H_INCLUDED defined
