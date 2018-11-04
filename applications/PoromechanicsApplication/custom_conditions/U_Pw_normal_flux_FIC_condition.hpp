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


#if !defined(KRATOS_U_PW_NORMAL_FLUX_FIC_CONDITION_H_INCLUDED )
#define  KRATOS_U_PW_NORMAL_FLUX_FIC_CONDITION_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_conditions/U_Pw_condition.hpp"
#include "custom_conditions/U_Pw_normal_flux_condition.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPwNormalFluxFICCondition : public UPwNormalFluxCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( UPwNormalFluxFICCondition );
    
    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPwCondition<TDim,TNumNodes>::mThisIntegrationMethod;
    typedef typename UPwNormalFluxCondition<TDim,TNumNodes>::NormalFluxVariables NormalFluxVariables;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPwNormalFluxFICCondition() : UPwNormalFluxCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    UPwNormalFluxFICCondition( IndexType NewId, GeometryType::Pointer pGeometry ) : UPwNormalFluxCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    UPwNormalFluxFICCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPwNormalFluxCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

    // Destructor
    ~UPwNormalFluxFICCondition() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const override;
 
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:   
    
    struct NormalFluxFICVariables
    {
        double DtPressureCoefficient;
        double ElementLength;
        double BiotModulusInverse;
        
        array_1d<double,TNumNodes> DtPressureVector;
        BoundedMatrix<double,TNumNodes,TNumNodes> PMatrix;
    };
    
    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo ) override;
    
    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo ) override;
    
    void CalculateElementLength(double& rElementLength, const GeometryType& Geom);
    
    
    void CalculateAndAddLHSStabilization(MatrixType& rLeftHandSideMatrix, NormalFluxVariables& rVariables, NormalFluxFICVariables& rFICVariables);

    void CalculateAndAddBoundaryMassMatrix(MatrixType& rLeftHandSideMatrix, NormalFluxVariables& rVariables, NormalFluxFICVariables& rFICVariables);
    

    void CalculateAndAddRHSStabilization(VectorType& rRightHandSideVector, NormalFluxVariables& rVariables, NormalFluxFICVariables& rFICVariables);

    void CalculateAndAddBoundaryMassFlow(VectorType& rRightHandSideVector, NormalFluxVariables& rVariables, NormalFluxFICVariables& rFICVariables);
    
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
    
}; // class UPwNormalFluxFICCondition.

} // namespace Kratos.

#endif // KRATOS_U_PW_NORMAL_FLUX_FIC_CONDITION_H_INCLUDED defined 
