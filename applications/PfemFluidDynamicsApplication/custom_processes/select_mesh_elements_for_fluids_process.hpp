//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_SELECT_MESH_ELEMENTS_FOR_FLUIDS_PROCESS_H_INCLUDED)
#define KRATOS_SELECT_MESH_ELEMENTS_FOR_FLUIDS_PROCESS_H_INCLUDED

// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesher_utilities.hpp"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "custom_processes/mesher_process.hpp"

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
    class SelectMeshElementsForFluidsProcess
        : public MesherProcess
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of Process
        KRATOS_CLASS_POINTER_DEFINITION(SelectMeshElementsForFluidsProcess);

        typedef ModelPart::ConditionType ConditionType;
        typedef ModelPart::PropertiesType PropertiesType;
        typedef ConditionType::GeometryType GeometryType;
        typedef GlobalPointersVector<Node> NodeWeakPtrVectorType;
        typedef std::size_t SizeType;
        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        SelectMeshElementsForFluidsProcess(ModelPart &rModelPart,
                                           MesherUtilities::MeshingParameters &rRemeshingParameters,
                                           int EchoLevel)
            : mrModelPart(rModelPart),
              mrRemesh(rRemeshingParameters)
        {
            mEchoLevel = EchoLevel;
        }

        /// Destructor.
        virtual ~SelectMeshElementsForFluidsProcess() {}

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
                std::cout << " [ SELECT MESH ELEMENTS in PfemFluid: (" << mrRemesh.OutMesh.GetNumberOfElements() << ") " << std::endl;
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
                    SizeType numfreesurf = 0;
                    SizeType numboundary = 0;
                    SizeType numrigid = 0;
                    SizeType numInletNodes = 0;
                    SizeType numisolated = 0;
                    bool noremesh = false;
                    std::vector<double> normVelocityP;
                    normVelocityP.resize(nds, false);
                    SizeType checkedNodes = 0;
                    SizeType countIsolatedWallNodes = 0;
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

                        if (vertices.back().IsNot(RIGID) && vertices.back().IsNot(SOLID))
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

                        if (vertices.back().Is(RIGID) || vertices.back().Is(SOLID))
                        {
                            if (vertices.back().Is(RIGID))
                            {
                                rigidNodeLocalMeshSize += vertices.back().FastGetSolutionStepValue(NODAL_H_WALL);
                                rigidNodeMeshCounter += 1.0;
                            }

                            numrigid++;

                            NodeWeakPtrVectorType &rN = vertices.back().GetValue(NEIGHBOUR_NODES);
                            bool localIsolatedWallNode = true;
                            for (SizeType i = 0; i < rN.size(); i++)
                            {
                                if (rN[i].IsNot(RIGID))
                                {
                                    localIsolatedWallNode = false;
                                }
                            }
                            if (localIsolatedWallNode == true)
                            {
                                countIsolatedWallNodes++;
                            }
                        }

                        if (vertices.back().IsNot(RIGID) && vertices.back().Is(BOUNDARY))
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

                        if (refiningBox == true && vertices.back().IsNot(RIGID))
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
                        SetMeshSizeInBoundaryZones(meanMeshSize, previouslyFreeSurfaceNodes, numfreesurf, currentTime, deltaTime, rigidNodeLocalMeshSize, rigidNodeMeshCounter);
                    }

                    if (refiningBox == true)
                    {
                        IncreaseAlphaForRefininedZones(Alpha, increaseAlfa, nds, numfreesurf, numrigid, numisolated);
                    }

                    sumIsolatedFreeSurf = numisolated + numfreesurf;
                    sumPreviouslyIsolatedFreeSurf = previouslyFreeSurfaceNodes + previouslyIsolatedNodes;

                    ModifyAlpha(Alpha, dimension, nds, numfreesurf, numrigid, numisolated, numInletNodes, previouslyIsolatedNodes, previouslyFreeSurfaceNodes);

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
                                    ControlSkewedElements2D(accepted, normVelocityP, nodesVelocities);
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
                                    ControlSkewedElements3D(accepted, normVelocityP, nodesVelocities);
                                }
                            }
                        }
                    }

                    // // to control that the element has a good shape
                    if (dimension == 3 && accepted && numrigid < 3 &&
                        (previouslyIsolatedNodes == 4 || previouslyFreeSurfaceNodes == 4 || sumIsolatedFreeSurf == 4 || numfreesurf == 4 || numisolated == 4 || (numrigid == 2 && isolatedNodesInTheElement > 1)))
                    {
                        ControlSliverElements(accepted, number_of_slivers, vertices, nodesCoordinates, CriticalVolume);
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
                        if (!(i_node->Is(FREE_SURFACE) || i_node->Is(RIGID)))
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
            return "SelectMeshElementsForFluidsProcess";
        }

        /// Print information about this object.
        void PrintInfo(std::ostream &rOStream) const override
        {
            rOStream << "SelectMeshElementsForFluidsProcess";
        }

        /// Print object's data.
        void PrintData(std::ostream &rOStream) const override
        {
        }

        ///@}
        ///@name Friends
        ///@{

        ///@}

    protected:
        ///@name Static Member Variables
        ///@{

        ///@}
        ///@name Static Member Variables
        ///@{
        ModelPart &mrModelPart;
        MesherUtilities::MeshingParameters &mrRemesh;
        MesherUtilities mMesherUtilities;
        int mEchoLevel;

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

        void SetMeshSizeInBoundaryZones(double &meanMeshSize,
                                        bool previouslyFreeSurfaceNodes,
                                        bool numfreesurf,
                                        double currentTime,
                                        double deltaTime,
                                        SizeType rigidNodeLocalMeshSize,
                                        SizeType rigidNodeMeshCounter)
        {
            KRATOS_TRY
            const double rigidWallMeshSize = rigidNodeLocalMeshSize / rigidNodeMeshCounter;
            const double ratio = rigidWallMeshSize / meanMeshSize;
            double tolerance = 1.8;
            if (currentTime < 10 * deltaTime)
            {
                tolerance = 1.5;
            }
            if (ratio > tolerance && numfreesurf == 0 && previouslyFreeSurfaceNodes == 0)
            {
                meanMeshSize *= 0.5;
                meanMeshSize += 0.5 * rigidWallMeshSize;
            }

            KRATOS_CATCH("")
        }

        void IncreaseAlphaForRefininedZones(double &Alpha,
                                            bool increaseAlfa,
                                            SizeType nds,
                                            SizeType numfreesurf,
                                            SizeType numrigid,
                                            SizeType numisolated)
        {
            KRATOS_TRY
            if (increaseAlfa == true)
            {
                if (numfreesurf < nds && numisolated == 0)
                {
                    Alpha *= 1.2;
                }
                else if (numfreesurf == 0 && numrigid == 0 && numisolated == 0)
                {
                    Alpha *= 1.4;
                }
                else if (numfreesurf == 0 && numrigid > (0.5 * nds) && numisolated == 0)
                {
                    Alpha *= 5.0;
                }
                else if (numfreesurf == 0 && numrigid > 0 && numisolated == 0)
                {
                    Alpha *= 1.8;
                }
            }

            KRATOS_CATCH("")
        }

        void ModifyAlpha(double &Alpha,
                         const SizeType dimension,
                         SizeType nds,
                         SizeType numfreesurf,
                         SizeType numrigid,
                         SizeType numisolated,
                         SizeType numInletNodes,
                         SizeType previouslyIsolatedNodes,
                         SizeType previouslyFreeSurfaceNodes)
        {
            KRATOS_TRY
            if (dimension == 2)
            {
                if (numrigid == 0 && numfreesurf == 0 && numisolated == 0 && previouslyIsolatedNodes == 0 && previouslyFreeSurfaceNodes == 0)
                {
                    Alpha *= 1.5;
                }
                else if (numfreesurf == 0 && numisolated == 0 && previouslyIsolatedNodes == 0 && previouslyFreeSurfaceNodes == 0)
                {
                    Alpha *= 1.25;
                }
                else if (numisolated == 0 && previouslyIsolatedNodes == 0 && numfreesurf < nds && previouslyFreeSurfaceNodes < nds)
                {
                    Alpha *= 1.125;
                }
                else
                {
                    Alpha *= 0.975;
                }
            }
            else if (dimension == 3)
            {
                if (numrigid == 0 && numfreesurf == 0 && numisolated == 0 && previouslyIsolatedNodes == 0 && previouslyFreeSurfaceNodes == 0)
                {
                    Alpha *= 1.5;
                }
                else if (numfreesurf == 0 && numisolated == 0 && previouslyIsolatedNodes == 0 && previouslyFreeSurfaceNodes == 0)
                {
                    Alpha *= 1.25;
                }
                else if (numisolated == 0 && previouslyIsolatedNodes == 0 && numfreesurf < nds && previouslyFreeSurfaceNodes < nds)
                {
                    Alpha *= 1.05;
                }
                else
                {
                    Alpha *= 0.95;
                }

                if (numrigid == nds)
                {
                    Alpha *= 0.95;
                }
                if (mrRemesh.ExecutionOptions.Is(MesherUtilities::REFINE_WALL_CORNER))
                {
                    if (numrigid == 3 && numfreesurf == 0 && numisolated == 0)
                    {
                        Alpha *= 1.1;
                    }
                    if (numrigid == 2 && numfreesurf == 0 && numisolated == 0)
                    {
                        Alpha *= 1.05;
                    }
                }
            }

            if (numInletNodes > 0)
            {
                Alpha *= 1.5;
            }

            KRATOS_CATCH("")
        }

        void ControlSkewedElements2D(bool &accepted,
                                     std::vector<double> &normVelocityP,
                                     std::vector<array_1d<double, 3>> &nodesVelocities)
        {
            KRATOS_TRY

            const double maxValue = 1.5;
            const double minValue = 1.0 / maxValue;
            if (normVelocityP[0] / normVelocityP[1] > maxValue || normVelocityP[0] / normVelocityP[1] < minValue ||
                normVelocityP[0] / normVelocityP[2] > maxValue || normVelocityP[0] / normVelocityP[2] < minValue ||
                normVelocityP[1] / normVelocityP[2] > maxValue || normVelocityP[1] / normVelocityP[2] < minValue)
            {
                accepted = false;
            }
            else
            {
                const double cosAngle01 = (nodesVelocities[0][0] * nodesVelocities[1][0] + nodesVelocities[0][1] * nodesVelocities[1][1]) /
                                          (std::sqrt(std::pow(nodesVelocities[0][0], 2) + std::pow(nodesVelocities[0][1], 2)) *
                                           std::sqrt(std::pow(nodesVelocities[1][0], 2) + std::pow(nodesVelocities[1][1], 2)));
                const double cosAngle02 = (nodesVelocities[0][0] * nodesVelocities[2][0] + nodesVelocities[0][1] * nodesVelocities[2][1]) /
                                          (std::sqrt(std::pow(nodesVelocities[0][0], 2) + std::pow(nodesVelocities[0][1], 2)) *
                                           std::sqrt(std::pow(nodesVelocities[2][0], 2) + std::pow(nodesVelocities[2][1], 2)));
                const double cosAngle12 = (nodesVelocities[1][0] * nodesVelocities[2][0] + nodesVelocities[1][1] * nodesVelocities[2][1]) /
                                          (std::sqrt(std::pow(nodesVelocities[1][0], 2) + std::pow(nodesVelocities[1][1], 2)) *
                                           std::sqrt(std::pow(nodesVelocities[2][0], 2) + std::pow(nodesVelocities[2][1], 2)));

                if (std::abs(cosAngle01) < 0.95 || std::abs(cosAngle02) < 0.95 || std::abs(cosAngle12) < 0.95)
                {
                    accepted = false;
                }
            }
            KRATOS_CATCH("")
        }

        void ControlSkewedElements3D(bool &accepted,
                                     std::vector<double> &normVelocityP,
                                     std::vector<array_1d<double, 3>> &nodesVelocities)
        {
            KRATOS_TRY
            const double maxValue = 2.5;
            const double minValue = 1.0 / maxValue;
            if (normVelocityP[0] / normVelocityP[1] < minValue || normVelocityP[0] / normVelocityP[2] < minValue || normVelocityP[0] / normVelocityP[3] < minValue ||
                normVelocityP[0] / normVelocityP[1] > maxValue || normVelocityP[0] / normVelocityP[2] > maxValue || normVelocityP[0] / normVelocityP[3] > maxValue ||
                normVelocityP[1] / normVelocityP[2] < minValue || normVelocityP[1] / normVelocityP[3] < minValue ||
                normVelocityP[1] / normVelocityP[2] > maxValue || normVelocityP[1] / normVelocityP[3] > maxValue ||
                normVelocityP[2] / normVelocityP[3] < minValue ||
                normVelocityP[2] / normVelocityP[3] > maxValue)
            {
                accepted = false;
            }
            else
            {
                const double cosAngle01 = (nodesVelocities[0][0] * nodesVelocities[1][0] + nodesVelocities[0][1] * nodesVelocities[1][1] + nodesVelocities[0][1] * nodesVelocities[1][2]) /
                                          (std::sqrt(std::pow(nodesVelocities[0][0], 2) + std::pow(nodesVelocities[0][1], 2) + std::pow(nodesVelocities[0][2], 2)) *
                                           std::sqrt(std::pow(nodesVelocities[1][0], 2) + std::pow(nodesVelocities[1][1], 2) + std::pow(nodesVelocities[1][2], 2)));
                const double cosAngle02 = (nodesVelocities[0][0] * nodesVelocities[2][0] + nodesVelocities[0][1] * nodesVelocities[2][1] + nodesVelocities[0][1] * nodesVelocities[2][2]) /
                                          (std::sqrt(std::pow(nodesVelocities[0][0], 2) + std::pow(nodesVelocities[0][1], 2) + std::pow(nodesVelocities[0][2], 2)) *
                                           std::sqrt(std::pow(nodesVelocities[2][0], 2) + std::pow(nodesVelocities[2][1], 2) + std::pow(nodesVelocities[2][2], 2)));
                const double cosAngle03 = (nodesVelocities[0][0] * nodesVelocities[3][0] + nodesVelocities[0][1] * nodesVelocities[3][1] + nodesVelocities[0][1] * nodesVelocities[3][2]) /
                                          (std::sqrt(std::pow(nodesVelocities[0][0], 2) + std::pow(nodesVelocities[0][1], 2) + std::pow(nodesVelocities[0][2], 2)) *
                                           std::sqrt(std::pow(nodesVelocities[3][0], 2) + std::pow(nodesVelocities[3][1], 2) + std::pow(nodesVelocities[3][2], 2)));
                const double cosAngle12 = (nodesVelocities[1][0] * nodesVelocities[2][0] + nodesVelocities[1][1] * nodesVelocities[2][1] + nodesVelocities[1][1] * nodesVelocities[2][2]) /
                                          (std::sqrt(std::pow(nodesVelocities[1][0], 2) + std::pow(nodesVelocities[1][1], 2) + std::pow(nodesVelocities[1][2], 2)) *
                                           std::sqrt(std::pow(nodesVelocities[2][0], 2) + std::pow(nodesVelocities[2][1], 2) + std::pow(nodesVelocities[2][2], 2)));
                const double cosAngle13 = (nodesVelocities[1][0] * nodesVelocities[3][0] + nodesVelocities[1][1] * nodesVelocities[3][1] + nodesVelocities[1][1] * nodesVelocities[3][2]) /
                                          (std::sqrt(std::pow(nodesVelocities[1][0], 2) + std::pow(nodesVelocities[1][1], 2) + std::pow(nodesVelocities[1][2], 2)) *
                                           std::sqrt(std::pow(nodesVelocities[3][0], 2) + std::pow(nodesVelocities[3][1], 2) + std::pow(nodesVelocities[3][2], 2)));
                const double cosAngle23 = (nodesVelocities[2][0] * nodesVelocities[3][0] + nodesVelocities[2][1] * nodesVelocities[3][1] + nodesVelocities[2][1] * nodesVelocities[3][2]) /
                                          (std::sqrt(std::pow(nodesVelocities[2][0], 2) + std::pow(nodesVelocities[2][1], 2) + std::pow(nodesVelocities[2][2], 2)) *
                                           std::sqrt(std::pow(nodesVelocities[3][0], 2) + std::pow(nodesVelocities[3][1], 2) + std::pow(nodesVelocities[3][2], 2)));

                if (std::abs(cosAngle01) < 0.85 || std::abs(cosAngle02) < 0.85 || std::abs(cosAngle03) < 0.85 || std::abs(cosAngle12) < 0.85 || std::abs(cosAngle13) < 0.85 || std::abs(cosAngle23) < 0.85)
                {
                    accepted = false;
                    // std::cout << "The angle between the velocity vectors is too big" << std::endl;
                }
            }
            KRATOS_CATCH("")
        }

        void ControlSliverElements(bool &accepted,
                                   int &number_of_slivers,
                                   Geometry<Node> vertices,
                                   std::vector<array_1d<double, 3>> nodesCoordinates,
                                   double CriticalVolume)
        {
            KRATOS_TRY
            Geometry<Node> *tetrahedron = new Tetrahedra3D4<Node>(vertices);
            double Volume = tetrahedron->Volume();

            // a1 slope x for plane on the first triangular face of the tetrahedra (nodes A,B,C)
            // b1 slope y for plane on the first triangular face of the tetrahedra (nodes A,B,C)
            // c1 slope z for plane on the first triangular face of the tetrahedra (nodes A,B,C)
            const double a1 = (nodesCoordinates[1][1] - nodesCoordinates[0][1]) * (nodesCoordinates[2][2] - nodesCoordinates[0][2]) - (nodesCoordinates[2][1] - nodesCoordinates[0][1]) * (nodesCoordinates[1][2] - nodesCoordinates[0][2]);
            const double b1 = (nodesCoordinates[1][2] - nodesCoordinates[0][2]) * (nodesCoordinates[2][0] - nodesCoordinates[0][0]) - (nodesCoordinates[2][2] - nodesCoordinates[0][2]) * (nodesCoordinates[1][0] - nodesCoordinates[0][0]);
            const double c1 = (nodesCoordinates[1][0] - nodesCoordinates[0][0]) * (nodesCoordinates[2][1] - nodesCoordinates[0][1]) - (nodesCoordinates[2][0] - nodesCoordinates[0][0]) * (nodesCoordinates[1][1] - nodesCoordinates[0][1]);
            // a2 slope x for plane on the second triangular face of the tetrahedra (nodes A,B,D)
            // b2 slope y for plane on the second triangular face of the tetrahedra (nodes A,B,D)
            // c2 slope z for plane on the second triangular face of the tetrahedra (nodes A,B,D)
            const double a2 = (nodesCoordinates[1][1] - nodesCoordinates[0][1]) * (nodesCoordinates[3][2] - nodesCoordinates[0][2]) - (nodesCoordinates[3][1] - nodesCoordinates[0][1]) * (nodesCoordinates[1][2] - nodesCoordinates[0][2]);
            const double b2 = (nodesCoordinates[1][2] - nodesCoordinates[0][2]) * (nodesCoordinates[3][0] - nodesCoordinates[0][0]) - (nodesCoordinates[3][2] - nodesCoordinates[0][2]) * (nodesCoordinates[1][0] - nodesCoordinates[0][0]);
            const double c2 = (nodesCoordinates[1][0] - nodesCoordinates[0][0]) * (nodesCoordinates[3][1] - nodesCoordinates[0][1]) - (nodesCoordinates[3][0] - nodesCoordinates[0][0]) * (nodesCoordinates[1][1] - nodesCoordinates[0][1]);
            // a3 slope x for plane on the third triangular face of the tetrahedra (nodes B,C,D)
            // b3 slope y for plane on the third triangular face of the tetrahedra (nodes B,C,D)
            // c3 slope z for plane on the third triangular face of the tetrahedra (nodes B,C,D)
            const double a3 = (nodesCoordinates[1][1] - nodesCoordinates[2][1]) * (nodesCoordinates[3][2] - nodesCoordinates[2][2]) - (nodesCoordinates[3][1] - nodesCoordinates[2][1]) * (nodesCoordinates[1][2] - nodesCoordinates[2][2]);
            const double b3 = (nodesCoordinates[1][2] - nodesCoordinates[2][2]) * (nodesCoordinates[3][0] - nodesCoordinates[2][0]) - (nodesCoordinates[3][2] - nodesCoordinates[2][2]) * (nodesCoordinates[1][0] - nodesCoordinates[2][0]);
            const double c3 = (nodesCoordinates[1][0] - nodesCoordinates[2][0]) * (nodesCoordinates[3][1] - nodesCoordinates[2][1]) - (nodesCoordinates[3][0] - nodesCoordinates[2][0]) * (nodesCoordinates[1][1] - nodesCoordinates[2][1]);
            // a4 slope x for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
            // b4 slope y for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
            // c4 slope z for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
            const double a4 = (nodesCoordinates[0][1] - nodesCoordinates[2][1]) * (nodesCoordinates[3][2] - nodesCoordinates[2][2]) - (nodesCoordinates[3][1] - nodesCoordinates[2][1]) * (nodesCoordinates[0][2] - nodesCoordinates[2][2]);
            const double b4 = (nodesCoordinates[0][2] - nodesCoordinates[2][2]) * (nodesCoordinates[3][0] - nodesCoordinates[2][0]) - (nodesCoordinates[3][2] - nodesCoordinates[2][2]) * (nodesCoordinates[0][0] - nodesCoordinates[2][0]);
            const double c4 = (nodesCoordinates[0][0] - nodesCoordinates[2][0]) * (nodesCoordinates[3][1] - nodesCoordinates[2][1]) - (nodesCoordinates[3][0] - nodesCoordinates[2][0]) * (nodesCoordinates[0][1] - nodesCoordinates[2][1]);

            const double cosAngle12 = (a1 * a2 + b1 * b2 + c1 * c2) / (std::sqrt(std::pow(a1, 2) + std::pow(b1, 2) + std::pow(c1, 2)) * std::sqrt(std::pow(a2, 2) + std::pow(b2, 2) + std::pow(c2, 2)));
            const double cosAngle13 = (a1 * a3 + b1 * b3 + c1 * c3) / (std::sqrt(std::pow(a1, 2) + std::pow(b1, 2) + std::pow(c1, 2)) * std::sqrt(std::pow(a3, 2) + std::pow(b3, 2) + std::pow(c3, 2)));
            const double cosAngle14 = (a1 * a4 + b1 * b4 + c1 * c4) / (std::sqrt(std::pow(a1, 2) + std::pow(b1, 2) + std::pow(c1, 2)) * std::sqrt(std::pow(a4, 2) + std::pow(b4, 2) + std::pow(c4, 2)));
            const double cosAngle23 = (a3 * a2 + b3 * b2 + c3 * c2) / (std::sqrt(std::pow(a3, 2) + std::pow(b3, 2) + std::pow(c3, 2)) * std::sqrt(std::pow(a2, 2) + std::pow(b2, 2) + std::pow(c2, 2)));
            const double cosAngle24 = (a4 * a2 + b4 * b2 + c4 * c2) / (std::sqrt(std::pow(a4, 2) + std::pow(b4, 2) + std::pow(c4, 2)) * std::sqrt(std::pow(a2, 2) + std::pow(b2, 2) + std::pow(c2, 2)));
            const double cosAngle34 = (a4 * a3 + b4 * b3 + c4 * c3) / (std::sqrt(std::pow(a4, 2) + std::pow(b4, 2) + std::pow(c4, 2)) * std::sqrt(std::pow(a3, 2) + std::pow(b3, 2) + std::pow(c3, 2)));

            if (std::abs(cosAngle12) > 0.999 || std::abs(cosAngle13) > 0.999 || std::abs(cosAngle14) > 0.999 || std::abs(cosAngle23) > 0.999 || std::abs(cosAngle24) > 0.999 || std::abs(cosAngle34) > 0.999) // if two faces are coplanar, I will erase the element (which is probably a sliver)
            {
                accepted = false;
                number_of_slivers++;
            }
            else if (Volume <= CriticalVolume)
            {
                accepted = false;
                number_of_slivers++;
            }
            delete tetrahedron;
            KRATOS_CATCH("")
        }

    private:
        /// Assignment operator.
        SelectMeshElementsForFluidsProcess &operator=(SelectMeshElementsForFluidsProcess const &rOther);

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
                                    SelectMeshElementsForFluidsProcess &rThis);

    /// output stream function
    inline std::ostream &operator<<(std::ostream &rOStream,
                                    const SelectMeshElementsForFluidsProcess &rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

} // namespace Kratos.

#endif // KRATOS_SELECT_MESH_ELEMENTS_FOR_FLUIDS_PROCESS_H_INCLUDED defined
