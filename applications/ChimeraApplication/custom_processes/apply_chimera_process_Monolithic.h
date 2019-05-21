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

// External includes
#include "includes/kratos_flags.h"
#include "utilities/binbased_fast_point_locator.h"

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_flags.h"
#include "includes/model_part.h"
#include "containers/model.h"
#include "geometries/geometry_data.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "elements/distance_calculation_element_simplex.h"

// Application includes
#include "chimera_application_variables.h"
#include "custom_processes/custom_calculate_signed_distance_process.h"
#include "custom_hole_cutting_process.h"
#include "custom_utilities/vtk_output.hpp"

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

/// Short class definition.

template <std::size_t TDim>
class ApplyChimeraProcessMonolithic : public Process
{
  public:
    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ApplyChimeraProcessMonolithic
    KRATOS_CLASS_POINTER_DEFINITION(ApplyChimeraProcessMonolithic);
    typedef ProcessInfo::Pointer ProcessInfoPointerType;
    typedef typename BinBasedFastPointLocator<TDim>::Pointer BinBasedPointLocatorPointerType;
    typedef Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>> VariableComponentType;
    typedef std::size_t IndexType;
    typedef Kratos::Variable<double> VariableType;
    typedef std::vector<IndexType> ConstraintIdsVectorType;


    ///@}
    ///@name Life Cycle
    ///@{

    ApplyChimeraProcessMonolithic(ModelPart &MainModelPart, Parameters rParameters) : Process(Flags()), mrMainModelPart(MainModelPart), m_parameters(rParameters)
    {
        Parameters default_parameters(R"(
            {
               	"chimera"   :   [
									[{
										"model_part_name":"GENERIC_background",
										"model_part_inside_boundary_name" :"GENERIC_domainboundary",
										"overlap_distance":0.045
									}],
									[{
										"model_part_name":"GENERIC_patch_1_1",
										"model_part_inside_boundary_name":"GENERIC_structure_1_1",
										"overlap_distance":0.045
									}],
									[{
										"model_part_name":"GENERIC_patch_2_1",
										"model_part_inside_boundary_name":"GENERIC_strcuture2_1",
										"overlap_distance":0.045
									}]
								]
            })");

        NumberOfLevels = m_parameters.size();
        for (int i = 0; i < NumberOfLevels; i++)
            LevelTable.push_back(m_parameters[i].size());

        ProcessInfoPointerType info = mrMainModelPart.pGetProcessInfo();
        this->pHoleCuttingProcess = CustomHoleCuttingProcess::Pointer(new CustomHoleCuttingProcess());
        this->pCalculateDistanceProcess = typename CustomCalculateSignedDistanceProcess<TDim>::Pointer(new CustomCalculateSignedDistanceProcess<TDim>());

        mNumberOfConstraintsAdded = 0;
    }

    /// Destructor.
    virtual ~ApplyChimeraProcessMonolithic()
    {
        Clear();
    }

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    virtual void Execute() override
    {
    }

    virtual void Clear()
    {
        KRATOS_INFO("Monolithic Chimera process is cleared") << std::endl;
    }

    void ExecuteBeforeSolutionLoop() override
    {
    }

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        // Actual execution of the functionality of this class
        for (ModelPart::ElementsContainerType::iterator it = mrMainModelPart.ElementsBegin(); it != mrMainModelPart.ElementsEnd(); ++it)
            it->SetValue(SPLIT_ELEMENT, false);
        DoChimeraLoop();
        KRATOS_CATCH("");
    }

    void ExecuteFinalizeSolutionStep() override
    {
        Clear();
        //for multipatch
        for (ModelPart::ElementsContainerType::iterator it = mrMainModelPart.ElementsBegin(); it != mrMainModelPart.ElementsEnd(); ++it)
        {
            it->Set(VISITED, false);
            it->SetValue(SPLIT_ELEMENT, false);
        }

        ModelPart::MasterSlaveConstraintContainerType::iterator it = mrMainModelPart.MasterSlaveConstraintsEnd() - 1;
        int constraintId = (*it).Id();

        for (int i = 1; i <= constraintId; ++i)
            mrMainModelPart.RemoveMasterSlaveConstraintFromAllLevels(i);
    }

    void ExecuteBeforeOutputStep() override
    {
    }

    void ExecuteAfterOutputStep() override
    {
    }

    void ExecuteFinalize() override
    {
    }

    void ApplyMpcConstraint(ModelPart &rBoundaryModelPart, BinBasedPointLocatorPointerType &pBinLocator)
    {
        //loop over nodes and find the triangle in which it falls, than do interpolation
        //array_1d<double, TDim + 1> N;
        Vector N;
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        const unsigned int n_boundary_nodes = rBoundaryModelPart.Nodes().size();
        ConstraintIdsVectorType ConstrainIdsForTheNode;
        std::size_t counter = 0;
        std::size_t removed_counter = 0;

        #pragma omp parallel for firstprivate(removed_counter) reduction(+:counter)
        for (unsigned int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn)
        {
            ModelPart::NodesContainerType::iterator iparticle = rBoundaryModelPart.NodesBegin() + i_bn;
            Node<3>::Pointer p_boundary_node = *(iparticle.base());

            bool NodeCoupled = false;
            if ((p_boundary_node)->IsDefined(VISITED))
                NodeCoupled = (p_boundary_node)->Is(VISITED);

            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
            Element::Pointer pElement;
            bool is_found = false;
            is_found = pBinLocator->FindPointOnMesh(p_boundary_node->Coordinates(), N, pElement, result_begin, max_results);

            if (NodeCoupled && is_found)
            {
                ConstrainIdsForTheNode = mNodeIdToConstraintIdsMap[p_boundary_node->Id()];
                for (auto const& constraint_id : ConstrainIdsForTheNode)
                {
                    mrMainModelPart.RemoveMasterSlaveConstraintFromAllLevels(constraint_id);
                    removed_counter++;
                }
                p_boundary_node->Set(VISITED, false);
            }

            // Initialise the boundary nodes dofs to 0 at ever time steps
            p_boundary_node->FastGetSolutionStepValue(VELOCITY_X, 0) = 0.0;
            p_boundary_node->FastGetSolutionStepValue(VELOCITY_Y, 0) = 0.0;
            if (TDim == 3)
                p_boundary_node->FastGetSolutionStepValue(VELOCITY_Z, 0) = 0.0;
            p_boundary_node->FastGetSolutionStepValue(PRESSURE, 0) = 0.0;

            if (is_found == true)
            {
                Geometry<Node<3>> &geom = pElement->GetGeometry();
                for (std::size_t i = 0; i < geom.size(); i++)
                {
                    //Interpolation of velocity
                    p_boundary_node->FastGetSolutionStepValue(VELOCITY_X, 0) += geom[i].GetDof(VELOCITY_X).GetSolutionStepValue(0) * N[i];
                    p_boundary_node->FastGetSolutionStepValue(VELOCITY_Y, 0) += geom[i].GetDof(VELOCITY_Y).GetSolutionStepValue(0) * N[i];

                    //Define master slave relation for velocity
                    AddMasterSlaveRelationWithNodesAndVariableComponents(geom[i], VELOCITY_X, *p_boundary_node, VELOCITY_X, N[i]);
                    AddMasterSlaveRelationWithNodesAndVariableComponents(geom[i], VELOCITY_Y, *p_boundary_node, VELOCITY_Y, N[i]);
                    if (TDim == 3)
                    {
                        //Interpolation of velocity
                        p_boundary_node->FastGetSolutionStepValue(VELOCITY_Z, 0) += geom[i].GetDof(VELOCITY_Z).GetSolutionStepValue(0) * N[i];
                        AddMasterSlaveRelationWithNodesAndVariableComponents(geom[i], VELOCITY_Z, *p_boundary_node, VELOCITY_Z, N[i]);
                    }

                    //Interpolation of pressure
                    p_boundary_node->FastGetSolutionStepValue(PRESSURE, 0) += geom[i].GetDof(PRESSURE).GetSolutionStepValue(0) * N[i];
                    //Defining master slave relation for pressure
                    AddMasterSlaveRelationWithNodesAndVariable(geom[i], PRESSURE, *p_boundary_node, PRESSURE, N[i]);
                    counter++;

                } // end of loop over host element nodes

                // Setting the buffer 1 same buffer 0
                p_boundary_node->FastGetSolutionStepValue(VELOCITY_X, 1) = p_boundary_node->FastGetSolutionStepValue(VELOCITY_X, 0);
                p_boundary_node->FastGetSolutionStepValue(VELOCITY_Y, 1) = p_boundary_node->FastGetSolutionStepValue(VELOCITY_Y, 0);
                if (TDim == 3)
                    p_boundary_node->FastGetSolutionStepValue(VELOCITY_Z, 1) = p_boundary_node->FastGetSolutionStepValue(VELOCITY_Z, 0);
                p_boundary_node->FastGetSolutionStepValue(PRESSURE, 1) = p_boundary_node->FastGetSolutionStepValue(PRESSURE, 0);
            }
            p_boundary_node->Set(VISITED, true);

        } // end of loop over boundary nodes

        counter /= TDim + 1;
        KRATOS_INFO("Pressure nodes from") << rBoundaryModelPart.Name() << " is coupled" << counter << std::endl;
        KRATOS_INFO("number of constriants created for this combination") << counter * 9 << std::endl;
        KRATOS_INFO("number of constraints removed for this combination") << removed_counter << std::endl;
    }

    void GetBoundingBox(ModelPart& rModelPart, std::vector<double>& rLowPoint, std::vector<double>& rHighPoint)
    {
        double rLowPoint0 = 1e10;
        double rLowPoint1 = 1e10;
        double rLowPoint2 = 1e10;

        double rHighPoint0 = -1e10;
        double rHighPoint1 = -1e10;
        double rHighPoint2 = -1e10;


        const unsigned int num_nodes = rModelPart.Nodes().size();
        #pragma omp parallel for reduction(min:rLowPoint0) reduction(min:rLowPoint1) reduction(min:rLowPoint2) reduction(max:rHighPoint0) reduction(max:rHighPoint1) reduction(max:rHighPoint2)
        for (unsigned int i_node=0; i_node<num_nodes; ++i_node)
        {
            ModelPart::NodesContainerType::iterator it_node = rModelPart.NodesBegin() + i_node;
            rLowPoint0 = std::min(it_node->X(), rLowPoint0);
            rLowPoint1 = std::min(it_node->Y(), rLowPoint1);
            rLowPoint2 = std::min(it_node->Z(), rLowPoint2);

            rHighPoint0 = std::max(it_node->X(), rHighPoint0);
            rHighPoint1 = std::max(it_node->Y(), rHighPoint1);
            rHighPoint2 = std::max(it_node->Z(), rHighPoint2);
        }

        rHighPoint[0] = rHighPoint0;
        rHighPoint[1] = rHighPoint1;
        rHighPoint[2] = rHighPoint2;

        rLowPoint[0] = rLowPoint0;
        rLowPoint[1] = rLowPoint1;
        rLowPoint[2] = rLowPoint2;
    }

    bool BoundingBoxTest(ModelPart &A, ModelPart &B) //background A and Patch B
    {
        std::vector<double> min_cornerA(3), max_cornerA(3), min_cornerB(3), max_cornerB(3);

        GetBoundingBox(A, min_cornerA, max_cornerA);
        GetBoundingBox(B, min_cornerB, max_cornerB);

        KRATOS_INFO("Bounding box of Background") << min_cornerA[0] << "::" << max_cornerA[0] << std::endl;
        KRATOS_INFO("Bounding box of patch") << min_cornerB[0] << "::" << max_cornerB[0] << std::endl;

        ProcessInfo &CurrentProcessInfo = A.GetProcessInfo();
        int idim = CurrentProcessInfo.GetValue(DOMAIN_SIZE);

        for (int i = 0; i < idim; i++)
        {
            if (min_cornerA[i] > max_cornerB[i])
                return false;
            if (max_cornerA[i] < min_cornerB[i])
                return false;
        }
        return true;
    }

    void FindOutsideBoundaryOfModelPartGivenInside(ModelPart &rModelPart, ModelPart &rInsideBoundary, ModelPart &rExtractedBoundaryModelPart)
    {
        std::size_t n_nodes = rModelPart.ElementsBegin()->GetGeometry().size();

        if (n_nodes == 3)
            this->pHoleCuttingProcess->ExtractOutsideBoundaryMesh(rInsideBoundary, rModelPart, rExtractedBoundaryModelPart);
        else if (n_nodes == 4)
            this->pHoleCuttingProcess->ExtractOutsideSurfaceMesh(rInsideBoundary, rModelPart, rExtractedBoundaryModelPart);
    }

    void DoChimeraLoop() //selecting patch and background combination for chimera method
    {
        for (ModelPart::ElementsContainerType::iterator it = mrMainModelPart.ElementsBegin(); it != mrMainModelPart.ElementsEnd(); ++it)
        {
            if (!it->Is(VISITED)) //for multipatch
                it->Set(ACTIVE, true);
        }

        for (ModelPart::NodesContainerType::iterator it = mrMainModelPart.NodesBegin(); it != mrMainModelPart.NodesEnd(); ++it)
            it->Set(VISITED, false);

        int MainDomainOrNot = 1;

        for (int BG_i = 0; BG_i < NumberOfLevels; BG_i++) // Iteration for selecting background
        {
            for (int BG_j = 0; BG_j < LevelTable[BG_i]; BG_j++) //TODO change the names
            {
                for (int patch_i = BG_i + 1; patch_i < NumberOfLevels; patch_i++) // Iteration for selecting patch
                {
                    for (int patch_j = 0; patch_j < LevelTable[patch_i]; patch_j++)
                    {
                        m_background_model_part_name = m_parameters[BG_i][BG_j]["model_part_name"].GetString();
                        m_domain_boundary_model_part_name = m_parameters[BG_i][BG_j]["model_part_inside_boundary_name"].GetString();
                        m_patch_model_part_name = m_parameters[patch_i][patch_j]["model_part_name"].GetString();
                        m_patch_inside_boundary_model_part_name = m_parameters[patch_i][patch_j]["model_part_inside_boundary_name"].GetString();

                        double mesh_size_1 = m_parameters[BG_i][BG_j]["overlap_distance"].GetDouble();
                        double mesh_size_2 = m_parameters[patch_i][patch_j]["overlap_distance"].GetDouble();

                        if (mesh_size_1 > mesh_size_2)
                            m_overlap_distance = mesh_size_1;
                        else
                            m_overlap_distance = mesh_size_2;

                        KRATOS_INFO("Formulating Chimera for the combination \n background::") << m_background_model_part_name << "  \t Patch::" << m_patch_model_part_name << std::endl;

                        MainDomainOrNot = 1;
                        if (BG_i == 0) // a check to identify computational Domain boundary
                            MainDomainOrNot = -1;

                        FormulateChimera(MainDomainOrNot);
                    }
                }
            }
        }
        KRATOS_INFO("End of chimera loop") << std::endl;
        int number_of_constraints = mrMainModelPart.MasterSlaveConstraints().size();
        KRATOS_INFO("Total number of constraints created so far") << number_of_constraints << std::endl;
    }

    //Apply Chimera with or without overlap

    void FormulateChimera(int MainDomainOrNot)
    {
        ModelPart &rBackgroundModelPart = mrMainModelPart.GetSubModelPart(m_background_model_part_name);
        ModelPart &rPatchModelPart = mrMainModelPart.GetSubModelPart(m_patch_model_part_name);
        ModelPart &rDomainBoundaryModelPart = mrMainModelPart.GetSubModelPart(m_domain_boundary_model_part_name);
        ModelPart &rPatchInsideBoundaryModelPart = mrMainModelPart.GetSubModelPart(m_patch_inside_boundary_model_part_name);

        this->pBinLocatorForBackground = BinBasedPointLocatorPointerType(new BinBasedFastPointLocator<TDim>(rBackgroundModelPart));
        this->pBinLocatorForPatch = BinBasedPointLocatorPointerType(new BinBasedFastPointLocator<TDim>(rPatchModelPart));
        this->pBinLocatorForBackground->UpdateSearchDatabase();
        this->pBinLocatorForPatch->UpdateSearchDatabase();

        const double epsilon = 1e-12;
        if (m_overlap_distance < epsilon)
            KRATOS_THROW_ERROR("", "Overlap distance should be a positive number \n", "");

        if (m_overlap_distance > epsilon)
        {
            Model &current_model = mrMainModelPart.GetModel();
            ModelPart &pHoleModelPart = current_model.CreateModelPart("HoleModelpart");
            ModelPart &pHoleBoundaryModelPart = current_model.CreateModelPart("HoleBoundaryModelPart");
            ModelPart &pModifiedPatchBoundaryModelPart = current_model.CreateModelPart("ModifiedPatchBoundary");
            ModelPart &pModifiedPatchModelPart = current_model.CreateModelPart("ModifiedPatch");
            bool BBoxOverlapTest = BoundingBoxTest(rBackgroundModelPart, rPatchModelPart); // true if they dont overlap

            if (BBoxOverlapTest)
            {
                KRATOS_INFO("Bounding boxes overlap , So finding the modified patch boundary") << std::endl;
                this->pCalculateDistanceProcess->CalculateSignedDistance(rPatchModelPart, rDomainBoundaryModelPart);
                this->pHoleCuttingProcess->RemoveOutOfDomainPatchAndReturnModifiedPatch(rPatchModelPart, rPatchInsideBoundaryModelPart, pModifiedPatchModelPart, pModifiedPatchBoundaryModelPart, MainDomainOrNot);
            }
            else
            {
                KRATOS_INFO("Bounding boxes does NOT overlap , So finding outside boundary of patch using the inside boundary") << std::endl;
                FindOutsideBoundaryOfModelPartGivenInside(rPatchModelPart, rPatchInsideBoundaryModelPart, pModifiedPatchBoundaryModelPart);
            }

            this->pCalculateDistanceProcess->CalculateSignedDistance(rBackgroundModelPart, pModifiedPatchBoundaryModelPart);
            this->pHoleCuttingProcess->CreateHoleAfterDistance(rBackgroundModelPart, pHoleModelPart, pHoleBoundaryModelPart, m_overlap_distance);

            //for multipatch
            for (ModelPart::ElementsContainerType::iterator it = pHoleModelPart.ElementsBegin(); it != pHoleModelPart.ElementsEnd(); ++it)
                it->Set(VISITED, true);

            KRATOS_INFO("Formulate Chimera: Calculated Nodal area and nodal mass") << std::endl;

            IndexType num_constraints_required = (TDim+1) * (pModifiedPatchBoundaryModelPart.Nodes().size() +  pHoleBoundaryModelPart.Nodes().size());

            CreateConstraintIds(num_constraints_required);
            mNumberOfConstraintsAdded = 0;

            ApplyMpcConstraint(pModifiedPatchBoundaryModelPart, pBinLocatorForBackground);
            ApplyMpcConstraint(pHoleBoundaryModelPart, pBinLocatorForPatch);

            KRATOS_INFO("Patch boundary coupled with background & HoleBoundary  coupled with patch using nearest element approach") << std::endl;
            KRATOS_INFO("Formulate Chimera: Appplied MPCs") << std::endl;

            current_model.DeleteModelPart("HoleModelpart");
            current_model.DeleteModelPart("HoleBoundaryModelPart");
            current_model.DeleteModelPart("ModifiedPatchBoundary");
            current_model.DeleteModelPart("ModifiedPatch");
        }
        KRATOS_INFO("End of Formulate Chimera") << std::endl;
    }

    void CreateConstraintIds(const IndexType NumberOfConstraintsRequired)
    {
        IndexType max_constraint_id = 0;
        // Get current maximum constraint ID
        if (mrMainModelPart.MasterSlaveConstraints().size() != 0)
        {
            mrMainModelPart.MasterSlaveConstraints().Sort();
            ModelPart::MasterSlaveConstraintContainerType::iterator it = mrMainModelPart.MasterSlaveConstraintsEnd() - 1;
            max_constraint_id = (*it).Id();
        }
        else
            max_constraint_id = 1;

        // Now create a vector size NumberOfConstraintsRequired
        mConstraintsIdVector.resize(NumberOfConstraintsRequired);
        std::iota (std::begin(mConstraintsIdVector), std::end(mConstraintsIdVector), max_constraint_id+1); // Fill with consecutive integers
    }

    void SetOverlapDistance(double distance)
    {
        this->m_overlap_distance = distance;
    }

    void AddMasterSlaveRelationWithNodesAndVariableComponents(Node<3> &MasterNode, VariableComponentType &MasterVariable, Node<3> &SlaveNode, VariableComponentType &SlaveVariable, double weight, double constant = 0.0)
    {
        SlaveNode.Set(SLAVE);
        #pragma omp critical
        {
            mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", mConstraintsIdVector[mNumberOfConstraintsAdded], MasterNode, MasterVariable, SlaveNode, SlaveVariable, weight, constant);
            mNodeIdToConstraintIdsMap[SlaveNode.Id()].push_back(mConstraintsIdVector[mNumberOfConstraintsAdded]);
            ++mNumberOfConstraintsAdded;
        }
    }

    // Functions with use two variables
    void AddMasterSlaveRelationWithNodesAndVariable(Node<3> &MasterNode, VariableType &MasterVariable, Node<3> &SlaveNode, VariableType &SlaveVariable, double weight, double constant = 0.0)
    {
        SlaveNode.Set(SLAVE);
        #pragma omp critical
        {
            mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", mConstraintsIdVector[mNumberOfConstraintsAdded], MasterNode, MasterVariable, SlaveNode, SlaveVariable, weight, constant);
            mNodeIdToConstraintIdsMap[SlaveNode.Id()].push_back(mConstraintsIdVector[mNumberOfConstraintsAdded]);
            ++mNumberOfConstraintsAdded;
        }
    }
    /**
		Activates the constraint set or deactivates
		@arg isActive true/false
		*/
    virtual std::string Info() const override
    {
        return "ApplyChimeraProcessMonolithic";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "ApplyChimeraProcessMonolithic";
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

    //ModelPart &mrBackGroundModelPart;
    //ModelPart &mrPatchSurfaceModelPart;
    BinBasedPointLocatorPointerType pBinLocatorForBackground; // Template argument 3 stands for 3D case
    BinBasedPointLocatorPointerType pBinLocatorForPatch;

    //for monolithic
    CustomHoleCuttingProcess::Pointer pHoleCuttingProcess;
    typename CustomCalculateSignedDistanceProcess<TDim>::Pointer pCalculateDistanceProcess;
    ModelPart &mrMainModelPart;
    double m_overlap_distance;
    int NumberOfLevels;
    std::vector<int> LevelTable;
    Parameters m_parameters;
    std::string m_background_model_part_name;
    std::string m_patch_boundary_model_part_name;
    std::string m_domain_boundary_model_part_name;
    std::string m_patch_inside_boundary_model_part_name;
    std::string m_patch_model_part_name;

    std::vector<int> mConstraintsIdVector;
    IndexType mNumberOfConstraintsAdded;

    ConstraintIdsVectorType VectorOfConstraintIds;
    std::unordered_map<IndexType, ConstraintIdsVectorType> mNodeIdToConstraintIdsMap;
    // epsilon
    //static const double epsilon;

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
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ApplyChimeraProcessMonolithic &operator=(ApplyChimeraProcessMonolithic const &rOther);

    ///@}
}; // Class ApplyChimeraProcessMonolithic

} // namespace Kratos.

#endif //  KRATOS_CUSTOM_APPLY_CHIMERA_MONOLITHIC_H_INCLUDED defined
