//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//

// System includes
#if !defined(KRATOS_ACOUSTIC_STRUCTURE_COUPLING_CONDITION_H_INCLUDED )
#define  KRATOS_ACOUSTIC_STRUCTURE_COUPLING_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/condition.h"
#include "includes/variables.h"
#include "mor_application_variables.h"

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
 * @class AcousticStructureCouplingCondition
 * @ingroup StructuralMechanicsApplication
 * @brief This class is the responsible to add the contributions of the RHS and LHS of the line loads of the structure
 * @details It allows to consider different types of pressure and line loads
 * @tparam TDim The dimension of the condition
 * @author Vicente Mataix Ferrandiz
 */
template<std::size_t TDim, bool TIsMapping>
class AcousticStructureCouplingCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// We define the base class BaseLoadCondition
    typedef Condition BaseType;

    /// Dfinition of the index type
    typedef BaseType::IndexType IndexType;

    /// Definition of the size type
    typedef BaseType::SizeType SizeType;

    /// Definition of the node type
    typedef BaseType::NodeType NodeType;

    /// Definition of the properties type
    typedef BaseType::PropertiesType PropertiesType;

    /// Definition of the geometry type with given NodeType
    typedef BaseType::GeometryType GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef BaseType::NodesArrayType NodesArrayType;

    /// Counted pointer of AcousticStructureCouplingCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( AcousticStructureCouplingCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor using an array of nodes
    AcousticStructureCouplingCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        );

    // Constructor using an array of nodes with properties
    AcousticStructureCouplingCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        );

    /// Destructor.
    ~AcousticStructureCouplingCondition() override;

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
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The vector containing the equation id
     * @param rCurrentProcessInfo The current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
     * @param rElementalDofList The vector containing the dof of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetDofList(
        DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets on rValues the nodal displacements
     * @param rValues The values of displacements
     * @param Step The step to be computed
     */
    void GetValuesVector(
        Vector& rValues,
        int Step = 0
        ) const override;

    /**
     * @brief Sets on rValues the nodal velocities
     * @param rValues The values of velocities
     * @param Step The step to be computed
     */
    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0
        ) const override;

    /**
     * @brief Sets on rValues the nodal accelerations
     * @param rValues The values of accelerations
     * @param Step The step to be computed
     */
    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0
        ) const override;

    /**
     * @brief This function provides a more general interface to the element.
     * @details It is designed so that rLHSvariables and rRHSvariables are passed to the element thus telling what is the desired output
     * @param rLeftHandSideMatrices container with the output left hand side matrices
     * @param rLHSVariables paramter describing the expected LHSs
     * @param rRightHandSideVectors container for the desired RHS output
     * @param rRHSVariables parameter describing the expected RHSs
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental right hand side vector only
      * @param rRightHandSideVector the elemental right hand side vector
      * @param rCurrentProcessInfo the current process info instance
      */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental mass matrix
      * @param rMassMatrix the elemental mass matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental damping matrix
      * @param rDampingMatrix the elemental damping matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfoO
        ) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rCurrentProcessInfo The current process info instance
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) override;

    /**
     * @brief Check if Rotational Dof existant
     * @return Trues if exists, false otherwise
     */
    virtual bool HasRotDof() const
    {
        if( !TIsMapping )
            return (GetGeometry()[0].HasDofFor(ROTATION_Z) && GetGeometry().size() == 2); //Z????
        else
            return (this->GetValue(MAPPING_NODES)[0]->HasDofFor(ROTATION_Z) && GetGeometry().size() == 2);
    }

    /**
     * @brief This method computes the DoF block size
     * @return The size of the DoF block
     */
    unsigned int GetBlockSize() const
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        if( HasRotDof() ) { // if it has rotations
            if(dim == 2)
                return 4;
            else if(dim == 3)
                return 7;
            else
                KRATOS_ERROR << "The condition only works for 2D and 3D elements";
        } else {
            return dim+1;
        }
    }

    /**
     * @brief This method computed the equation system size
     * @return The size
     */
    unsigned int GetSystemSize()
    {
        // unsigned int size = 0;
        if( !TIsMapping ) {
            return 0;
        } else {
            const SizeType number_of_nodes = GetGeometry().size();
            auto& mapped_nodes = this->GetValue(MAPPING_NODES);
            const SizeType n_mapped_nodes = mapped_nodes.size();
            const unsigned int dim = GetGeometry().WorkingSpaceDimension();

            if( dim == 2 ) {
                return n_mapped_nodes * (this->GetBlockSize() - 1) + number_of_nodes;
            } else if( dim == 3 ) {
                return n_mapped_nodes * 3 + number_of_nodes;
            } else {
                KRATOS_ERROR << "The condition only works for 2D and 3D elements";
            }
        }
        // return size;
    }

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
    std::string Info() const override
    {
        std::string s;
        TIsMapping ? s = " with " : s = " without ";
        std::stringstream buffer;
        buffer << "AcousticStructureCouplingCondition #" << Id() << s << "mapping";
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        std::string s;
        TIsMapping ? s = " with " : s = " without ";
        rOStream << "AcousticStructureCouplingCondition #" << Id() << s << "mapping";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
        // rOStream << "\n\n\t Is mapping condition: " << TIsMapping << std::endl;
    }

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
     * @brief This functions calculates LHS and mass contribution
     * @param rLeftHandSideMatrix: The LHS / mass matrix
     * @param rRightHandSideVector: The RHS
     * @param rCurrentProcessInfo: The current process info instance
     * @param CalculateMassMatrixFlag: The flag to set if to compute mass. otherwise stiffness
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateMassMatrixFlag,
        const bool CalculateVectorFlag
        ) ;

    /**
     * This functions computes the integration weight to consider
     * @param IntegrationPoints: The array containing the integration points
     * @param PointNumber: The id of the integration point considered
     * @param detJ: The determinant of the jacobian of the element
     */
    virtual double GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const SizeType PointNumber,
        const double detJ
        ) const;

    /**
     * @brief This method provides the local axis
     * @param rLocalAxis The local axis
     * @param rJacobian The jacobian matrix
     */
    void GetLocalAxis1(
        array_1d<double, 3>& rLocalAxis,
        const Matrix& rJacobian
        ) const;

    /**
     * @brief This method provides the local axis
     * @param rLocalAxis The local axis
     */
    void GetLocalAxis2(
        array_1d<double, 3>& rLocalAxis,
        const Matrix& rJacobian
        ) const;

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

    // A protected default constructor necessary for serialization
    AcousticStructureCouplingCondition() {};

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, AcousticStructureCouplingCondition );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, AcousticStructureCouplingCondition );
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //AcousticStructureCouplingCondition& operator=(const AcousticStructureCouplingCondition& rOther);

    /// Copy constructor.
    //AcousticStructureCouplingCondition(const AcousticStructureCouplingCondition& rOther);


    ///@}

}; // Class AcousticStructureCouplingCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<std::size_t TDim, bool TIsMapping>
inline std::istream& operator >> (std::istream& rIStream,
        AcousticStructureCouplingCondition<TDim,TIsMapping>& rThis);
/// output stream function
template<std::size_t TDim, bool TIsMapping>
inline std::ostream& operator << (std::ostream& rOStream,
        const AcousticStructureCouplingCondition<TDim,TIsMapping>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_ACOUSTIC_STRUCTURE_COUPLING_CONDITION_H_INCLUDED  defined


