//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined (KRATOS_RVE_PERIODICITY_UTILITY_H_INCLUDED)
#define KRATOS_RVE_PERIODICITY_UTILITY_H_INCLUDED

// System includes

// External includes
#include <tuple>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{
///@name Kratos Classes
///@{

/**
 * @class RVEPeriodicityUtility
 * @ingroup StructuralMechanicsApplication
 * @brief This defines a class to define periodic BC to a RVE
 * @details It uses MPC in order to set the periodic BC
 * @author Riccardo Rossi
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) RVEPeriodicityUtility
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RVEPeriodicityUtility
    KRATOS_CLASS_POINTER_DEFINITION(RVEPeriodicityUtility);

    /// Definition of the index type
    typedef std::size_t IndexType;

    /// Definition of the data tuple type
    typedef std::tuple<std::vector<IndexType>, std::vector<double>, Vector> DataTupletype;

    /// The DoF type definition
    typedef Dof<double> DofType;

    /// The DoF pointer vector type definition
    typedef std::vector< DofType::Pointer > DofPointerVectorType;

    /// Definition of the node
    typedef Node<3> NodeType;

    /// Definition of the component of variable type
    typedef Variable<double> DoubleVariableType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RVEPeriodicityUtility(
        ModelPart& rDestinationModelPart,
        std::size_t EchoLevel = 0
        ) : mrModelPart(rDestinationModelPart)
          , mEchoLevel(EchoLevel)  { }

    /// Copy constructor.
    RVEPeriodicityUtility(RVEPeriodicityUtility const& rOther) = delete;

    /// Destructor.
    virtual ~RVEPeriodicityUtility() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    RVEPeriodicityUtility &operator=(RVEPeriodicityUtility const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function assign a pairing condition between two modelparts which contain two flat faces, parallel to  each other and separated by a distance rDistance.
     *  @details Note that this function should be called multiple times to pair the different faces in a box.
     *  @param rMasterModelPart master part to be paired
     *  @param rSlaveModelPart slave in the pairing
     *  @param rStrainTensor strain tensor which will be used in computing the pairing conditions the condition to be guaranteed will be that :    uslave = umaster + rStrainTensor * rDirection
     *  @param rDirection  a node with coordinates Xs on the slave, will be paired to the corresponding point with coordinates Xm on the master
     *         Xm will be computed as      Xm = Xs - rDirection
     */
    void AssignPeriodicity(ModelPart& rMasterModelPart,
                           ModelPart& rSlaveModelPart,
                           const Matrix& rStrainTensor,
                           const Vector& rDirection);

    /** this function finalizes the computation of the pairings. It can be called ONLY ONCE
     * @param rVariable is the value to which the pairing condition will be applied (needs to be a Variable with components)
     */
    void Finalize(const Variable<array_1d<double, 3>>& rVariable);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "RVEPeriodicityUtility";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "RVEPeriodicityUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

  protected:
    ///@name Protected Operations
    ///@{

    /**
    * @brief This method appends the weights and ids to construct the MPC
    * @param rAux The auxiliar map containing the ids and tuples
    * @param MasterId The id of the master dof
    * @param MasterWeight The constribution of the master dof in the MPC
    * @param rFinalMastersIds The resulting vector of ids
    * @param rFinalMastersWeights The resulting vector of weights
    * @param rFinalT The resulting vector of constants (rigid displacements)
    */
    void AppendIdsAndWeights(
        std::map<IndexType, DataTupletype>& rAux,
        const IndexType MasterId,
        const double MasterWeight,
        std::vector<IndexType>& rFinalMastersIds,
        std::vector<double>& rFinalMastersWeights,
        Vector& rFinalT
        );

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart; /// The model part where to apply the constraints

    std::map<IndexType, DataTupletype> mAuxPairings; /// This map contains the pairings

    std::size_t mEchoLevel = 0; /// The echo level of the utility

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method helps to create the constraints
     * @param rConstraintId The new constraint ID
     * @param rVar The variable to be set with the constraint
     * @param pSlaveNode The pointer to the slave node
     * @param rMasterIds The id of the master dof
     * @param rRelationMatrix The relation matrix between master slave dofs
     * @param rTranslationVector The rigid motion vector of the dofs
     */
    MasterSlaveConstraint::Pointer  GenerateConstraint(
        IndexType& rConstraintId,
        const DoubleVariableType& rVar,
        NodeType::Pointer pSlaveNode,
        const std::vector<IndexType>& rMasterIds,
        const Matrix& rRelationMatrix,
        const Vector& rTranslationVector
        );

    ///@}

}; // Class RVEPeriodicityUtility

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream& rIStream,
                                RVEPeriodicityUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream &operator<<(std::ostream& rOStream,
                                const RVEPeriodicityUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RVE_PERIODICITY_UTILITY_H_INCLUDED defined
