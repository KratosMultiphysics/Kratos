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

// System includes
#if !defined(KRATOS_LINE_LOAD_FROM_DEM_CONDITION_2D_H_INCLUDED )
#define  KRATOS_LINE_LOAD_FROM_DEM_CONDITION_2D_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "includes/define.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "custom_conditions/line_load_condition.h"

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

/**
 * @class LineLoadFromDEMCondition2D
 */
template<std::size_t TDim>
class KRATOS_API(DEM_STRUCTURES_COUPLING_APPLICATION) LineLoadFromDEMCondition2D
    : public LineLoadCondition<TDim>
{
public:
    ///@name Type Definitions
    ///@{

    /// We define the base class LineLoadCondition
    typedef LineLoadCondition<TDim> BaseType;

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;

    /// Counted pointer of LineLoadFromDEMCondition2D
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( LineLoadFromDEMCondition2D );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    LineLoadFromDEMCondition2D()
        : LineLoadCondition<TDim>() {}

    // Constructor using an array of nodes
    LineLoadFromDEMCondition2D( IndexType NewId, GeometryType::Pointer pGeometry )
        : LineLoadCondition<TDim>( NewId, pGeometry ) {}

    // Constructor using an array of nodes with properties
    LineLoadFromDEMCondition2D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : LineLoadCondition<TDim>( NewId, pGeometry, pProperties ) {}

    // Destructor
    ~LineLoadFromDEMCondition2D() override {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new condition pointer
     * @param NewId the ID of the new condition
     * @param ThisNodes the nodes of the new condition
     * @param pProperties the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new condition pointer
     * @param NewId the ID of the new condition
     * @param pGeom the geometry to be employed
     * @param pProperties the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new condition pointer and clones the previous condition data
     * @param NewId the ID of the new condition
     * @param ThisNodes the nodes of the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone (
        IndexType NewId,
        NodesArrayType const& ThisNodes
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

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix: The LHS
     * @param rRightHandSideVector: The RHS
     * @param rCurrentProcessInfo: The current process info instance
     * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        ) override;

    virtual void InterpolateLineLoad(array_1d<double,3>& r_surface_load,
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


    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
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
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LineLoadCondition<TDim> );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LineLoadCondition<TDim> );
    }


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //LineLoadFromDEMCondition2D& operator=(const LineLoadFromDEMCondition2D& rOther);

    /// Copy constructor.
    //LineLoadFromDEMCondition2D(const LineLoadFromDEMCondition2D& rOther);


    ///@}

}; // Class LineLoadFromDEMCondition2D

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_LINE_LOAD_FROM_DEM_CONDITION_2D_H_INCLUDED  defined


