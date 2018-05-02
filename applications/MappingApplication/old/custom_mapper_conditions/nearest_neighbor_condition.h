//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_NEAREST_NEIGHBOR_CONDITION_H)
#define KRATOS_NEAREST_NEIGHBOR_CONDITION_H

// System includes

// External includes

// Project includes
#include "includes/condition.h"


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

class NearestNeighborCondition : public Condition {
public:

    ///@name Type Definitions
    ///@{

    typedef Condition BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of NearestNeighborCondition
    KRATOS_CLASS_POINTER_DEFINITION(NearestNeighborCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Constructor using an array of nodes
    */
    NearestNeighborCondition(IndexType NewId, const NodesArrayType& ThisNodes);

    /**
    * Constructor using Geometry
    */
    NearestNeighborCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    /**
    * Constructor using Properties
    */
    NearestNeighborCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /**
    * Copy Constructor
    */
    NearestNeighborCondition(NearestNeighborCondition const& rOther);

    /**
    * Destructor
    */
    virtual ~NearestNeighborCondition();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    NearestNeighborCondition & operator=(NearestNeighborCondition const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /**
    * CONDITIONS inherited from this class have to implement next
    * Create and Clone methods: MANDATORY
    */

    /**
    * creates a new condition pointer
    * @param NewId: the ID of the new condition
    * @param ThisNodes: the nodes of the new condition
    * @param pProperties: the properties assigned to the new condition
    * @return a Pointer to the new condition
    */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    /**
    * creates a new condition pointer
    * @param NewId: the ID of the new condition
    * @param pGeom: the geometry to be employed
    * @param pProperties: the properties assigned to the new condition
    * @return a Pointer to the new condition
    */
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    /**
    * creates a new condition pointer and clones the previous condition data
    * @param NewId: the ID of the new condition
    * @param ThisNodes: the nodes of the new condition
    * @param pProperties: the properties assigned to the new condition
    * @return a Pointer to the new condition
    */
    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    /**
    * this determines the condition equation ID vector for all condition
    * DOFs
    * @param rResult: the condition equation ID vector
    * @param rCurrentProcessInfo: the current process info instance
    */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override;

    /**
    * CONDITIONS inherited from this class have to implement next
    * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
    * they can be managed internally with a private method to do the same calculations
    * only once: MANDATORY
    */

    /**
    * this is called during the assembling process in order
    * to calculate all condition contributions to the global system
    * matrix and the right hand side
    * @param rLeftHandSideMatrix: the condition left hand side matrix
    * @param rRightHandSideVector: the condition right hand side
    * @param rCurrentProcessInfo: the current process info instance
    */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    /**
    * this is called during the assembling process in order
    * to calculate the condition left hand side matrix only
    * @param rLeftHandSideMatrix: the condition left hand side matrix
    * @param rCurrentProcessInfo: the current process info instance
    */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;

    /**
    * this is called during the assembling process in order
    * to calculate the condition right hand side vector only
    * @param rRightHandSideVector: the condition right hand side vector
    * @param rCurrentProcessInfo: the current process info instance
    */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    /**
    * This method provides the place to perform checks on the completeness of the input
    * and the compatibility with the problem options as well as the contitutive laws selected
    * It is designed to be called only once (or anyway, not often) typically at the beginning
    * of the calculations, so to verify that nothing is missing from the input
    * or that no common error is found.
    * @param rCurrentProcessInfo
    * this method is: MANDATORY
    */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    //TODO bring "CleanMemory" back => VertexMorphing

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


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

    NearestNeighborCondition() {} // Default Constructor for Serializer

    ///@}
    ///@name Protected Operations
    ///@{

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
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

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

}; // Class NearestNeighborCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_NEAREST_NEIGHBOR_CONDITION_H  defined
