//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_SELECT_MESH_ELEMENTS_FOR_FLUIDS_CUT_PFEM_PROCESS_H_INCLUDED)
#define KRATOS_SELECT_MESH_ELEMENTS_FOR_FLUIDS_CUT_PFEM_PROCESS_H_INCLUDED

// External includes

// System includes

// Project includes
#include "custom_processes/select_mesh_elements_for_fluids_process.hpp"

/// VARIABLES used:
// Data:
// Flags:    (checked) TO_ERASE, BOUNDARY, NEW_ENTITY
//           (set)
//           (modified)
//           (reset)
//(set):=(set in this process)

namespace Kratos
{

    ///@name Kratos Classes
    ///@{

    /// Refine Mesh Elements Process 2D and 3D
    /** The process labels the elements to be refined in the mesher
    it applies a size constraint to elements that must be refined.

*/
    class SelectMeshElementsForFluidsCutPfemProcess
        : public SelectMeshElementsForFluidsProcess
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of Process
        KRATOS_CLASS_POINTER_DEFINITION(SelectMeshElementsForFluidsCutPfemProcess);

        typedef ModelPart::ConditionType ConditionType;
        typedef ModelPart::PropertiesType PropertiesType;
        typedef ConditionType::GeometryType GeometryType;
        typedef GlobalPointersVector<Node> NodeWeakPtrVectorType;
        typedef std::size_t SizeType;

        typedef SelectMeshElementsForFluidsProcess BaseType;
        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        SelectMeshElementsForFluidsCutPfemProcess(ModelPart &rModelPart,
                                                  MesherUtilities::MeshingParameters &rRemeshingParameters,
                                                  int EchoLevel)
            : BaseType(rModelPart, rRemeshingParameters, EchoLevel)
        {
            KRATOS_INFO("SelectMeshElementsForFluidsCutPfemProcess") << " activated " << std::endl;
        }

        /// Destructor.
        virtual ~SelectMeshElementsForFluidsCutPfemProcess() {}

        ///@}
        ///@name Operators
        ///@{

        /// This operator is provided to call the process as a function and simply calls the Execute method.
        void operator()()
        {
            Execute();
        }

        ///@}
        ///@name Operations
        ///@{

        /// Execute method is used to execute the Process algorithms.
        void Execute() override
        {
            KRATOS_TRY

            if (mEchoLevel > 1)
            {
                std::cout << " [ SELECT MESH ELEMENTS FOR CUT PFEM ANALYSIS in PfemFluid: (" << mrRemesh.OutMesh.GetNumberOfElements() << ") " << std::endl;
                std::cout << "MODEL PART InNumberOfElements " << mrRemesh.InMesh.GetNumberOfElements() << std::endl;
                std::cout << "MODEL PART InNumberOfPoints " << mrRemesh.InMesh.GetNumberOfPoints() << std::endl;
                std::cout << "MODEL PART OutNumberOfElements " << mrRemesh.OutMesh.GetNumberOfElements() << std::endl;
                std::cout << "MODEL PART OutNumberOfPoints " << mrRemesh.OutMesh.GetNumberOfPoints() << std::endl;
            }
            const int &OutNumberOfElements = mrRemesh.OutMesh.GetNumberOfElements();
            mrRemesh.PreservedElements.clear();
            mrRemesh.PreservedElements.resize(OutNumberOfElements, false);
            std::fill(mrRemesh.PreservedElements.begin(), mrRemesh.PreservedElements.end(), 0);
            mrRemesh.MeshElementsSelectedFlag = true;

            MesherUtilities MesherUtils;
            double ModelPartVolume = MesherUtils.ComputeModelPartVolume(mrModelPart);
            double CriticalVolume = 0.05 * ModelPartVolume / double(mrModelPart.Elements().size());

            mrRemesh.Info->NumberOfElements = 0;

            const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
            double currentTime = rCurrentProcessInfo[TIME];
            double deltaTime = rCurrentProcessInfo[DELTA_TIME];
            int number_of_slivers = 0;

            bool refiningBox = false;
            for (SizeType index = 0; index < mrRemesh.UseRefiningBox.size(); index++)
            {
                if (mrRemesh.UseRefiningBox[index] == true && currentTime > mrRemesh.RefiningBoxInitialTime[index] && currentTime < mrRemesh.RefiningBoxFinalTime[index])
                {
                    refiningBox = true;
                }
            }

            if (mrRemesh.ExecutionOptions.IsNot(MesherUtilities::SELECT_TESSELLATION_ELEMENTS))
            {
                for (int el = 0; el < OutNumberOfElements; el++)
                {
                    mrRemesh.PreservedElements[el] = 1;
                    mrRemesh.Info->NumberOfElements += 1;
                }
            }
            else
            {
                if (mEchoLevel > 1)
                    std::cout << " Start Element Selection " << OutNumberOfElements << std::endl;

                ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();
                const SizeType nds = element_begin->GetGeometry().size();
                const SizeType dimension = element_begin->GetGeometry().WorkingSpaceDimension();
                int *OutElementList = mrRemesh.OutMesh.GetElementList();
                ModelPart::NodesContainerType &rNodes = mrModelPart.Nodes();

                int el = 0;
                int number = 0;

                // #pragma omp parallel for reduction(+:number) private(el)
                for (el = 0; el < OutNumberOfElements; el++)
                {
                    Geometry<Node> vertices;
                    double meanMeshSize = mrRemesh.Refine->CriticalRadius; // this must be inside because if there is a refined zone, each element has a different critical radius
                    double distance_tolerance = 0.05 * meanMeshSize;
                    SizeType numfreesurf = 0;
                    SizeType numboundary = 0;
                    SizeType numrigid = 0;
                    SizeType numInletNodes = 0;
                    SizeType numisolated = 0;
                    bool noremesh = false;
                    std::vector<double> normVelocityP;
                    normVelocityP.resize(nds, false);
                    SizeType checkedNodes = 0;
                    bool increaseAlfa = false;
                    SizeType previouslyFreeSurfaceNodes = 0;
                    SizeType previouslyIsolatedNodes = 0;
                    SizeType sumPreviouslyIsolatedFreeSurf = 0;
                    SizeType sumIsolatedFreeSurf = 0;
                    std::vector<array_1d<double, 3>> nodesCoordinates;
                    nodesCoordinates.resize(nds);
                    std::vector<array_1d<double, 3>> nodesVelocities;
                    nodesVelocities.resize(nds);
                    SizeType isolatedNodesInTheElement = 0;
                    double rigidNodeLocalMeshSize = 0;
                    double rigidNodeMeshCounter = 0;

                    for (SizeType pn = 0; pn < nds; pn++)
                    {
                        if (OutElementList[el * nds + pn] <= 0)
                            std::cout << " ERROR: something is wrong: nodal id < 0 " << el << std::endl;

                        if ((SizeType)OutElementList[el * nds + pn] > mrRemesh.NodalPreIds.size())
                        {
                            std::cout << " ERROR: something is wrong: node out of bounds " << std::endl;
                            break;
                        }
                        vertices.push_back(rNodes(OutElementList[el * nds + pn]));

                        if (vertices.back().IsNot(RIGID) && vertices.back().IsNot(SOLID) && vertices.back().GetSolutionStepValue(DISTANCE) > distance_tolerance)
                        {
                            isolatedNodesInTheElement += vertices.back().FastGetSolutionStepValue(ISOLATED_NODE);
                        }
                        // check flags on nodes
                        if (vertices.back().Is(ISOLATED))
                        {
                            numisolated++;
                        }
                        if (vertices.back().Is(PFEMFlags::PREVIOUS_FREESURFACE))
                        {
                            previouslyFreeSurfaceNodes++;
                        }
                        if (vertices.back().Is(PFEMFlags::PREVIOUS_ISOLATED))
                        {
                            previouslyIsolatedNodes++;
                        }
                        if (vertices.back().Is(BOUNDARY))
                        {
                            numboundary++;
                        }
                        if (vertices.back().GetValue(NO_MESH))
                        {
                            noremesh = true;
                        }

                        if (vertices.back().Is(RIGID) || vertices.back().Is(SOLID) || (vertices.back().GetSolutionStepValue(DISTANCE) < 0.0 && dimension == 3))
                        // if (vertices.back().Is(RIGID) || vertices.back().Is(SOLID))
                        {
                            if (vertices.back().Is(RIGID) || (vertices.back().GetSolutionStepValue(DISTANCE) < 0.0 && dimension == 3))
                            {
                                rigidNodeLocalMeshSize += vertices.back().FastGetSolutionStepValue(NODAL_H_WALL);
                                rigidNodeMeshCounter += 1.0;
                            }

                            numrigid++;
                        }

                        if (vertices.back().IsNot(RIGID) && vertices.back().Is(BOUNDARY) && vertices.back().GetSolutionStepValue(DISTANCE) > distance_tolerance)
                        {
                            numfreesurf++;
                            const array_1d<double, 3> &velocityP0 = vertices.back().FastGetSolutionStepValue(VELOCITY, 0);
                            normVelocityP[pn] = norm_2(velocityP0);
                            nodesVelocities[pn] = velocityP0;
                            checkedNodes++;
                        }
                        else if (vertices.back().Is(ISOLATED))
                        {
                            const array_1d<double, 3> &velocityP0 = vertices.back().FastGetSolutionStepValue(VELOCITY, 0);
                            normVelocityP[pn] = norm_2(velocityP0);
                            nodesVelocities[pn] = velocityP0;
                            checkedNodes++;
                        }
                        if (vertices.back().Is(INLET))
                        {
                            numInletNodes++;
                        }

                        if (refiningBox == true && vertices.back().IsNot(RIGID) && vertices.back().GetSolutionStepValue(DISTANCE) > distance_tolerance)
                        {
                            if (dimension == 2)
                            {
                                MesherUtils.DefineMeshSizeInTransitionZones2D(mrRemesh, currentTime, vertices.back().Coordinates(), meanMeshSize, increaseAlfa);
                            }
                            else if (dimension == 3)
                            {
                                MesherUtils.DefineMeshSizeInTransitionZones3D(mrRemesh, currentTime, vertices.back().Coordinates(), meanMeshSize, increaseAlfa);
                            }
                            CriticalVolume = 0.05 * (std::pow(meanMeshSize, 3) / (6.0 * std::sqrt(2)));
                        }

                        if (dimension == 3)
                        {
                            nodesCoordinates[pn] = vertices.back().Coordinates();
                        }
                    }

                    double Alpha = mrRemesh.AlphaParameter; //*nds;

                    if (rigidNodeMeshCounter > 0)
                    {
                        this->SetMeshSizeInBoundaryZones(meanMeshSize, previouslyFreeSurfaceNodes, numfreesurf, currentTime, deltaTime, rigidNodeLocalMeshSize, rigidNodeMeshCounter);
                    }

                    if (refiningBox == true)
                    {
                        this->IncreaseAlphaForRefininedZones(Alpha, increaseAlfa, nds, numfreesurf, numrigid, numisolated);
                    }

                    sumIsolatedFreeSurf = numisolated + numfreesurf;
                    sumPreviouslyIsolatedFreeSurf = previouslyFreeSurfaceNodes + previouslyIsolatedNodes;

                    this->ModifyAlpha(Alpha, dimension, nds, numfreesurf, numrigid, numisolated, numInletNodes, previouslyIsolatedNodes, previouslyFreeSurfaceNodes);

                    bool accepted = false;

                    accepted = MesherUtils.AlphaShape(Alpha, vertices, dimension, meanMeshSize);

                    if (numrigid == nds || noremesh == true)
                    {
                        accepted = false;
                    }

                    if (accepted == true && (numfreesurf == nds || sumIsolatedFreeSurf == nds || sumPreviouslyIsolatedFreeSurf == nds))
                    {
                        if (dimension == 2)
                        {
                            // this is to avoid the formation of isolated elements with different velocity fields. They can give convergence problems
                            if ((numfreesurf == nds || sumIsolatedFreeSurf == nds) && numrigid == 0)
                            {
                                if (checkedNodes == nds)
                                {
                                    this->ControlSkewedElements2D(accepted, normVelocityP, nodesVelocities);
                                }
                            }
                            Geometry<Node> *triangle = new Triangle2D3<Node>(vertices);
                            double elementArea = triangle->Area();
                            if (elementArea < CriticalVolume)
                            {
                                accepted = false;
                                // std::cout << "ATTENTION!!! ELEMENT WITH AREA = " << elementArea << " versus critical area " << CriticalVolume << std::endl;
                            }
                            delete triangle;
                        }
                        else if (dimension == 3)
                        {
                            if ((numfreesurf == nds || sumIsolatedFreeSurf == nds || previouslyIsolatedNodes == nds || previouslyFreeSurfaceNodes == nds) && numrigid == 0 && isolatedNodesInTheElement > 0)
                            {
                                if (checkedNodes == nds)
                                {
                                    this->ControlSkewedElements3D(accepted, normVelocityP, nodesVelocities);
                                }
                            }
                        }
                    }

                    // // to control that the element has a good shape
                    if (dimension == 3 && accepted && numrigid < 3 &&
                        (previouslyIsolatedNodes == 4 || previouslyFreeSurfaceNodes == 4 || sumIsolatedFreeSurf == 4 || numfreesurf == 4 || numisolated == 4 || (numrigid == 2 && isolatedNodesInTheElement > 1)))
                    {
                        this->ControlSliverElements(accepted, number_of_slivers, vertices, nodesCoordinates, CriticalVolume);
                    }

                    if (accepted)
                    {
                        number += 1;
                        mrRemesh.PreservedElements[el] = number;
                    }
                }
                mrRemesh.Info->NumberOfElements = number;
            }

            if (mEchoLevel > 1)
            {
                std::cout << "Number of Preserved Fluid Elements " << mrRemesh.Info->NumberOfElements << " (slivers detected: " << number_of_slivers << ") " << std::endl;
                std::cout << "TOTAL removed nodes " << mrRemesh.Info->RemovedNodes << std::endl;
            }
            if (mrRemesh.ExecutionOptions.IsNot(MesherUtilities::KEEP_ISOLATED_NODES))
            {
                ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();
                const SizeType nds = (*element_begin).GetGeometry().size();

                int *OutElementList = mrRemesh.OutMesh.GetElementList();

                ModelPart::NodesContainerType &rNodes = mrModelPart.Nodes();

                // check engaged nodes
                for (int el = 0; el < OutNumberOfElements; el++)
                {
                    if (mrRemesh.PreservedElements[el])
                    {
                        for (SizeType pn = 0; pn < nds; pn++)
                        {
                            // set vertices
                            rNodes[OutElementList[el * nds + pn]].Set(BLOCKED);
                        }
                    }
                }

                int count_release = 0;
                for (ModelPart::NodesContainerType::iterator i_node = rNodes.begin(); i_node != rNodes.end(); i_node++)
                {
                    if (i_node->IsNot(BLOCKED))
                    {
                        if (!(i_node->Is(FREE_SURFACE) || i_node->Is(RIGID) || i_node->GetSolutionStepValue(DISTANCE) < 0.0))
                        {
                            i_node->Set(TO_ERASE);
                            if (mEchoLevel > 0)
                                std::cout << " NODE " << i_node->Id() << " RELEASE " << std::endl;
                            if (i_node->Is(BOUNDARY))
                                std::cout << " ERROR: node " << i_node->Id() << " IS BOUNDARY RELEASE " << std::endl;
                        }
                    }
                    if (i_node->Is(TO_ERASE))
                    {
                        count_release++;
                    }

                    i_node->Reset(BLOCKED);
                }

                if (mEchoLevel > 0)
                    std::cout << "   fluid NUMBER OF RELEASED NODES " << count_release << std::endl;
            }
            else
            {
                ModelPart::NodesContainerType &rNodes = mrModelPart.Nodes();
                for (ModelPart::NodesContainerType::iterator i_node = rNodes.begin(); i_node != rNodes.end(); i_node++)
                {
                    i_node->Reset(BLOCKED);
                }
            }

            mrRemesh.InputInitializedFlag = false;
            mMesherUtilities.SetNodes(mrModelPart, mrRemesh);
            mrRemesh.InputInitializedFlag = true;

            if (mEchoLevel > 1)
            {
                std::cout << "   Generated_Elements :" << OutNumberOfElements << std::endl;
                std::cout << "   Passed_AlphaShape  :" << mrRemesh.Info->NumberOfElements << std::endl;
                std::cout << "   SELECT MESH ELEMENTS ]; " << std::endl;
            }

            KRATOS_CATCH("")
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
            return "SelectMeshElementsForFluidsCutPfemProcess";
        }

        /// Print information about this object.
        void PrintInfo(std::ostream &rOStream) const override
        {
            rOStream << "SelectMeshElementsForFluidsCutPfemProcess";
        }

        /// Print object's data.
        void PrintData(std::ostream &rOStream) const override
        {
        }

        ///@}
        ///@name Friends
        ///@{

        ///@}

    private:
        ///@name Static Member Variables
        ///@{

        ///@}
        ///@name Static Member Variables
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
        ///@name Un accessible methods
        ///@{

        /// Assignment operator.
        SelectMeshElementsForFluidsCutPfemProcess &operator=(SelectMeshElementsForFluidsCutPfemProcess const &rOther);

        /// this function is a private function

        /// Copy constructor.
        // Process(Process const& rOther);

        ///@}

    }; // Class Process

    ///@}

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// input stream function
    inline std::istream &operator>>(std::istream &rIStream,
                                    SelectMeshElementsForFluidsCutPfemProcess &rThis);

    /// output stream function
    inline std::ostream &operator<<(std::ostream &rOStream,
                                    const SelectMeshElementsForFluidsCutPfemProcess &rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

} // namespace Kratos.

#endif // KRATOS_SELECT_MESH_ELEMENTS_FOR_FLUIDS_CUT_PFEM_PROCESS_H_INCLUDED defined
