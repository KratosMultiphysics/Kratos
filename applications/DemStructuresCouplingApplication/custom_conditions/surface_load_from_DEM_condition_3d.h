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

#if !defined(KRATOS_DEM_SURFACE_LOAD_FROM_DEM_CONDITION_3D_H_INCLUDED )
#define  KRATOS_SURFACE_LOAD_FROM_DEM_CONDITION_3D_H_INCLUDED

// Project includes
#include "includes/define.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "custom_conditions/surface_load_condition_3d.h"

// Application includes
#include "dem_structures_coupling_application_variables.h"


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

class KRATOS_API(DEM_STRUCTURES_COUPLING_APPLICATION)  SurfaceLoadFromDEMCondition3D
    : public SurfaceLoadCondition3D
{
public:

    ///@name Type Definitions
    ///@{

    typedef SurfaceLoadCondition3D BaseType;

    // Counted pointer of SurfaceLoadFromDEMCondition3D
    KRATOS_CLASS_POINTER_DEFINITION( SurfaceLoadFromDEMCondition3D );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    SurfaceLoadFromDEMCondition3D();

    // Constructor using an array of nodes
    SurfaceLoadFromDEMCondition3D(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        );

    // Constructor using an array of nodes with properties
    SurfaceLoadFromDEMCondition3D(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        );

    // Destructor
    ~SurfaceLoadFromDEMCondition3D() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // Name Operations
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * returns the used integration method.
     */
    GeometryData::IntegrationMethod GetIntegrationMethod() override;

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

    /**
     * This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix: The LHS
     * @param rRightHandSideVector: The RHS
     * @param rCurrentProcessInfo: The current process info instance
     * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        ) override;


    virtual void InterpolateSurfaceLoad(array_1d<double,3>& r_surface_load,
                                        const Matrix& n_container,
                                        const unsigned int& number_of_nodes,
                                        const unsigned int& g_point);

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

private:
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

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
    ///@name Private LifeCycle
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SurfaceLoadCondition3D );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SurfaceLoadCondition3D );
    }


}; // class SurfaceLoadFromDEMCondition3D.

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_SURFACE_LOAD_FROM_DEM_CONDITION_3D_H_INCLUDED  defined
