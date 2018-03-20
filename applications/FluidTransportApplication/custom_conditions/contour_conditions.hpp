//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Albert Puigferrat Perez
//                   Ignasi de Pouplana
//

#if !defined(KRATOS_CONTOUR_CONDITIONS_INCLUDED )
#define  KRATOS_CONTOUR_CONDITIONS_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
// #include "custom_conditions/U_Pw_condition.hpp"
// #include "custom_conditions/U_Pw_face_load_condition.hpp"
// #include "custom_utilities/element_utilities.hpp"
// #include "custom_utilities/condition_utilities.hpp"
// #include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(FLUID_TRANSPORT_APPLICATION) ContourConditions : public UPwFaceLoadCondition<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( ContourConditions );
    
    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPwCondition<TDim,TNumNodes>::mThisIntegrationMethod;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    ContourConditions() : UPwFaceLoadCondition<TDim,TNumNodes>() {}
    
    // Constructor 1
    ContourConditions( IndexType NewId, GeometryType::Pointer pGeometry ) : UPwFaceLoadCondition<TDim,TNumNodes>(NewId, pGeometry) {}
    
    // Constructor 2
    ContourConditions( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : UPwFaceLoadCondition<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

    // Destructor
    virtual ~ContourConditions() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties ) const;
 
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:   
    
    struct NormalFluxVariables
    {
        double NormalFlux;
        double IntegrationCoefficient;
        array_1d<double,TNumNodes> N;
        array_1d<double,TNumNodes> PVector;
    };
    
    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                                    
    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo );
    
    void CalculateAndAddRHS(VectorType& rRightHandSideVector, NormalFluxVariables& rVariables);
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Serialization
    
    friend class Serializer;
    
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
    }
    
}; // class ContourConditions.

} // namespace Kratos.

#endif // KRATOS_CONTOUR_CONDITIONS_INCLUDED defined 
