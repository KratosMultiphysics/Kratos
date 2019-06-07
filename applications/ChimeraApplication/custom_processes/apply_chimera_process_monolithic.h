// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//					 Navaneeth K Narayanan
//					 Rishith Ellath Meethal
// ==============================================================================
//

#if !defined(KRATOS_CUSTOM_APPLY_CHIMERA_MONOLITHIC_H_INCLUDED)
#define KRATOS_CUSTOM_APPLY_CHIMERA_MONOLITHIC_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include "omp.h"

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "containers/model.h"
#include "geometries/geometry_data.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "includes/linear_master_slave_constraint.h"

// Application includes
#include "chimera_application_variables.h"
#include "custom_processes/custom_calculate_signed_distance_process.h"
#include "custom_hole_cutting_process.h"
#include "utilities/binbased_fast_point_locator.h"

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

template <int TDim>
class ApplyChimeraProcessMonolithic : public Process
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    typedef ProcessInfo::Pointer ProcessInfoPointerType;
    typedef typename BinBasedFastPointLocator<TDim>::Pointer BinBasedPointLocatorPointerType;
    typedef Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>> VariableComponentType;
    typedef std::size_t IndexType;
    typedef Kratos::Variable<double> VariableType;
    typedef std::vector<IndexType> ConstraintIdsVectorType;
    typedef typename ModelPart::MasterSlaveConstraintType MasterSlaveConstraintType;
    typedef typename ModelPart::MasterSlaveConstraintContainerType MasterSlaveConstraintContainerType;
    typedef std::vector<MasterSlaveConstraintContainerType> MasterSlaveContainerVectorType;
    KRATOS_CLASS_POINTER_DEFINITION(ApplyChimeraProcessMonolithic);

    ///@}
    ///@name Life Cycle
    ///@{

    ApplyChimeraProcessMonolithic(ModelPart &rMainModelPart, Parameters iParameters);

    /// Destructor.
    virtual ~ApplyChimeraProcessMonolithic() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitializeSolutionStep() override;

    void ExecuteFinalizeSolutionStep() override;

    virtual std::string Info() const override
    {
        return "ApplyChimeraProcessMonolithic";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "ApplyChimeraProcessMonolithic" << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const override
    {
        KRATOS_INFO("\nNumber of slave nodes :: ") << std::endl;
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
    BinBasedPointLocatorPointerType mpBinLocatorForBackground; // Template argument 3 stands for 3D case
    BinBasedPointLocatorPointerType mpBinLocatorForPatch;

    //for monolithic
    CustomHoleCuttingProcess::Pointer mpHoleCuttingProcess;
    typename CustomCalculateSignedDistanceProcess<TDim>::Pointer mpCalculateDistanceProcess;
    ModelPart &mrMainModelPart;
    double mOverlapDistance;
    int mNumberOfLevels;
    std::vector<int> mLevelTable;
    Parameters mParameters;
    std::string m_background_model_part_name;
    std::string m_patch_boundary_model_part_name;
    std::string m_domain_boundary_model_part_name;
    std::string m_patch_inside_boundary_model_part_name;
    std::string m_patch_model_part_name;

    IndexType mNumberOfConstraintsAdded;

    std::unordered_map<IndexType, ConstraintIdsVectorType> mNodeIdToConstraintIdsMap;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void DoChimeraLoop();

    void FormulateChimera(int MainDomainOrNot);

    void CreateConstraintIds(std::vector<int> &rIdVector, const IndexType NumberOfConstraintsRequired);

    void SetOverlapDistance(double distance);

    void ApplyMpcConstraint(ModelPart &rBoundaryModelPart, BinBasedPointLocatorPointerType &pBinLocator);

    /**
     * @brief Computes the bounding box of the modelpart given. The low and high points (brute force way)
     * @param rModelPart Modelpart for which the bounding box is to be computed.
     * @param rLowPoint The lowest point in the modelpart (returned)
     * @param rHighPoint The highest point in the modelpart (returned)
     */
    void GetBoundingBox(ModelPart &rModelPart, std::vector<double> &rLowPoint, std::vector<double> &rHighPoint);

    /**
     * @brief Checks if two given modelparts (A and B) have bounding box overlaps
     *                  Here in Chimera A is usually the background and B is the patch
     * @param rModelPartA ModelPartA
     * @param rModelPartB ModelPartB
     * @return bool if the bounding boxes intersect or not.
     */
    bool BoundingBoxTest(ModelPart &rModelPartA, ModelPart &rModelPartB);

    /**
     * @brief Extracts the boundary of the modelpart when there is an internal boundary (surface)
     * @param rModelPart The modelpart for which boundary is to be extracted.
     * @param rInsideBoundary the internal boundary of the rModelPart
     * @param rExtractedBoundaryModelPart The result. That is the extracted boundary of the modelpart
     */
    void FindOutsideBoundaryOfModelPartGivenInside(ModelPart &rModelPart, ModelPart &rInsideBoundary, ModelPart &rExtractedBoundaryModelPart);

    /**
     * @brief Applies the master-slave constraint between the given master and slave nodes with corresponding variable.
     * @param rMasterSlaveContainer The container to which the constraint to be added (useful to so OpenMP loop)
     * @param rCloneConstraint The prototype of constraint which is to be added.
     * @param ConstraintId The ID of the constraint to be added.
     * @param rMasterNode The Master node of the constraint.
     * @param rMasterVariable The variable for the master node.
     * @param rSlaveNode The Slave node of the constraint.
     * @param rSlaveVariable The variable for the slave node.
     * @param Weight The weight of the Master node.
     * @param Constant The constant of the master slave relation.
     */
    template <typename TVariableType>
    inline void AddMasterSlaveRelation(MasterSlaveConstraintContainerType &rMasterSlaveContainer,
                                const LinearMasterSlaveConstraint& rCloneConstraint,
                                unsigned int ConstraintId,
                                Node<3> &rMasterNode,
                                TVariableType &rMasterVariable,
                                Node<3> &rSlaveNode,
                                TVariableType &rSlaveVariable,
                                const double Weight,
                                const double Constant = 0.0)
    {
        rSlaveNode.Set(SLAVE);
        ModelPart::MasterSlaveConstraintType::Pointer p_new_constraint = rCloneConstraint.Create(ConstraintId, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant);
        p_new_constraint->Set(TO_ERASE);
        mNodeIdToConstraintIdsMap[rSlaveNode.Id()].push_back(ConstraintId);
        rMasterSlaveContainer.insert(rMasterSlaveContainer.begin(), p_new_constraint);
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ApplyChimeraProcessMonolithic &operator=(ApplyChimeraProcessMonolithic const &rOther);

    ///@}
}; // Class ApplyChimeraProcessMonolithic

} // namespace Kratos.

#endif //  KRATOS_CUSTOM_APPLY_CHIMERA_MONOLITHIC_H_INCLUDED defined
