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

///VARIABLES used:
//Data:
//Flags:    (checked) TO_ERASE, BOUNDARY, NEW_ENTITY
//          (set)
//          (modified)
//          (reset)
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

        typedef GlobalPointersVector<Node<3>> NodeWeakPtrVectorType;
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
            mrRemesh.PreservedElements.resize(OutNumberOfElements);
            std::fill(mrRemesh.PreservedElements.begin(), mrRemesh.PreservedElements.end(), 0);
            mrRemesh.MeshElementsSelectedFlag = true;

            MesherUtilities MesherUtils;
            double ModelPartVolume = MesherUtils.ComputeModelPartVolume(mrModelPart);
            double CriticalVolume = 0.05 * ModelPartVolume / double(mrModelPart.Elements().size());

            mrRemesh.Info->NumberOfElements = 0;

            const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
            double currentTime = rCurrentProcessInfo[TIME];
            double timeInterval = rCurrentProcessInfo[DELTA_TIME];
            bool firstMesh = false;
            if (currentTime < 2 * timeInterval)
            {
                firstMesh = true;
            }

            bool box_side_element = false;
            bool wrong_added_node = false;

            int number_of_slivers = 0;

            bool refiningBox = mrRemesh.UseRefiningBox;
            if (!(mrRemesh.UseRefiningBox == true && currentTime > mrRemesh.RefiningBoxInitialTime && currentTime < mrRemesh.RefiningBoxFinalTime))
            {
                refiningBox = false;
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
                const unsigned int nds = element_begin->GetGeometry().size();
                const unsigned int dimension = element_begin->GetGeometry().WorkingSpaceDimension();
                int *OutElementList = mrRemesh.OutMesh.GetElementList();
                ModelPart::NodesContainerType &rNodes = mrModelPart.Nodes();

                int el = 0;
                int number = 0;

                //#pragma omp parallel for reduction(+:number) private(el)
                for (el = 0; el < OutNumberOfElements; el++)
                {
                    Geometry<Node<3>> vertices;
                    double MeanMeshSize = mrRemesh.Refine->CriticalRadius; // this must be inside because if there is a refined zone, each element has a different critical radius
                    unsigned int numfreesurf = 0;
                    unsigned int numboundary = 0;
                    unsigned int numrigid = 0;
                    unsigned int numinlet = 0;
                    unsigned int numisolated = 0;
                    bool noremesh = false;
                    std::vector<double> normVelocityP;
                    normVelocityP.resize(nds);
                    unsigned int checkedNodes = 0;
                    box_side_element = false;
                    unsigned int countIsolatedWallNodes = 0;
                    bool increaseAlfa = false;
                    unsigned int previouslyFreeSurfaceNodes = 0;
                    unsigned int previouslyIsolatedNodes = 0;
                    unsigned int sumPreviouslyIsolatedFreeSurf = 0;
                    unsigned int sumIsolatedFreeSurf = 0;
                    std::vector<array_1d<double, 3>> nodesCoordinates;
                    nodesCoordinates.resize(nds);
                    std::vector<array_1d<double, 3>> nodesVelocities;
                    nodesVelocities.resize(nds);
                    unsigned int isolatedNodesInTheElement = 0;
                    for (unsigned int pn = 0; pn < nds; pn++)
                    {
                        if (OutElementList[el * nds + pn] <= 0)
                            std::cout << " ERROR: something is wrong: nodal id < 0 " << el << std::endl;

                        if ((unsigned int)OutElementList[el * nds + pn] > mrRemesh.NodalPreIds.size())
                        {
                            wrong_added_node = true;
                            std::cout << " ERROR: something is wrong: node out of bounds " << std::endl;
                            break;
                        }
                        vertices.push_back(rNodes(OutElementList[el * nds + pn]));

                        if (vertices.back().IsNot(RIGID) && vertices.back().IsNot(SOLID))
                        {
                            isolatedNodesInTheElement += vertices.back().FastGetSolutionStepValue(ISOLATED_NODE);
                            // if (isolatedNode==true)
                            // std::cout << el << "node " << vertices.back().Id() << " " << vertices.back().X() << " " << vertices.back().Y() << std::endl;
                        }
                        // check flags on nodes
                        if (vertices.back().Is(ISOLATED))
                        {
                            numisolated++;
                        }
                        if (vertices.back().Is(BOUNDARY))
                        {
                            numboundary++;
                        }
                        if (vertices.back().GetValue(NO_MESH))
                        {
                            noremesh = true;
                        }

                        previouslyFreeSurfaceNodes += vertices.back().FastGetSolutionStepValue(FREESURFACE);       //it is 1 if it was free-surface (set in build_mesh_boundary_for_fluids)
                        previouslyIsolatedNodes += vertices.back().FastGetSolutionStepValue(PREVIOUS_FREESURFACE); //it is 1 if it was isolated (set in build_mesh_boundary_for_fluids)

                        if (vertices.back().Is(RIGID) || vertices.back().Is(SOLID))
                        {
                            numrigid++;

                            NodeWeakPtrVectorType &rN = vertices.back().GetValue(NEIGHBOUR_NODES);
                            bool localIsolatedWallNode = true;
                            for (unsigned int i = 0; i < rN.size(); i++)
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
                            numinlet++;
                        }

                        if (refiningBox == true && vertices.back().IsNot(RIGID))
                        {
                            if (dimension == 2)
                            {
                                SetAlphaForRefinedZones2D(MeanMeshSize, increaseAlfa, vertices.back().X(), vertices.back().Y());
                            }
                            else if (dimension == 3)
                            {
                                SetAlphaForRefinedZones3D(MeanMeshSize, increaseAlfa, vertices.back().X(), vertices.back().Y(), vertices.back().Z());
                            }
                        }
                        if (dimension == 3)
                        {
                            nodesCoordinates[pn] = vertices.back().Coordinates();
                        }
                    }

                    if (box_side_element || wrong_added_node)
                    {
                        std::cout << " ,,,,,,,,,,,,,,,,,,,,,,,,,,,,, Box_Side_Element " << std::endl;
                        continue;
                    }

                    double Alpha = mrRemesh.AlphaParameter; //*nds;

                    if (refiningBox == true)
                    {

                        IncreaseAlphaForRefininedZones(Alpha, increaseAlfa, nds, numfreesurf, numrigid, numisolated);

                        if (dimension == 3)
                        {
                            Alpha *= 1.1;
                        }
                    }

                    sumIsolatedFreeSurf = numisolated + numfreesurf;
                    sumPreviouslyIsolatedFreeSurf = previouslyFreeSurfaceNodes + previouslyIsolatedNodes;

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
                        // else if (numfreesurf < nds && numisolated < nds && previouslyIsolatedNodes < 3 && previouslyFreeSurfaceNodes < nds && sumPreviouslyIsolatedFreeSurf < nds && sumIsolatedFreeSurf < nds)
                        // {
                        //     Alpha *= 1.05;
                        // }

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
                    if (firstMesh == true)
                    {
                        Alpha *= 1.15;
                    }

                    if (numinlet > 0)
                    {
                        Alpha *= 1.5;
                    }

                    bool accepted = false;

                    accepted = MesherUtils.AlphaShape(Alpha, vertices, dimension, MeanMeshSize);

                    if (numrigid == nds || noremesh == true)
                    {
                        accepted = false;
                    }

                    if (accepted == true && (numfreesurf == nds || sumIsolatedFreeSurf == nds || sumPreviouslyIsolatedFreeSurf == nds) && firstMesh == false)
                    {
                        if (dimension == 2)
                        {
                            // this is to avoid the formation of isolated elements with different velocity fields. They can give convergence problems
                            if ((numfreesurf == nds || sumIsolatedFreeSurf == nds) && numrigid == 0)
                            {
                                if (checkedNodes == nds)
                                {
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
                                        double cosAngle01 = (nodesVelocities[0][0] * nodesVelocities[1][0] + nodesVelocities[0][1] * nodesVelocities[1][1]) /
                                                            (sqrt(pow(nodesVelocities[0][0], 2) + pow(nodesVelocities[0][1], 2)) *
                                                             sqrt(pow(nodesVelocities[1][0], 2) + pow(nodesVelocities[1][1], 2)));
                                        double cosAngle02 = (nodesVelocities[0][0] * nodesVelocities[2][0] + nodesVelocities[0][1] * nodesVelocities[2][1]) /
                                                            (sqrt(pow(nodesVelocities[0][0], 2) + pow(nodesVelocities[0][1], 2)) *
                                                             sqrt(pow(nodesVelocities[2][0], 2) + pow(nodesVelocities[2][1], 2)));
                                        double cosAngle12 = (nodesVelocities[1][0] * nodesVelocities[2][0] + nodesVelocities[1][1] * nodesVelocities[2][1]) /
                                                            (sqrt(pow(nodesVelocities[1][0], 2) + pow(nodesVelocities[1][1], 2)) *
                                                             sqrt(pow(nodesVelocities[2][0], 2) + pow(nodesVelocities[2][1], 2)));

                                        if (fabs(cosAngle01) < 0.95 || fabs(cosAngle02) < 0.95 || fabs(cosAngle12) < 0.95)
                                        {
                                            accepted = false;
                                            // std::cout << isolatedNodesInTheElement << " isolatedNodesInTheElement The angle between the velocity vectors is too big" << std::endl;
                                        }
                                    }
                                }
                            }
                            Geometry<Node<3>> *triangle = new Triangle2D3<Node<3>>(vertices);
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
                                        double cosAngle01 = (nodesVelocities[0][0] * nodesVelocities[1][0] + nodesVelocities[0][1] * nodesVelocities[1][1] + nodesVelocities[0][1] * nodesVelocities[1][2]) /
                                                            (sqrt(pow(nodesVelocities[0][0], 2) + pow(nodesVelocities[0][1], 2) + pow(nodesVelocities[0][2], 2)) *
                                                             sqrt(pow(nodesVelocities[1][0], 2) + pow(nodesVelocities[1][1], 2) + pow(nodesVelocities[1][2], 2)));
                                        double cosAngle02 = (nodesVelocities[0][0] * nodesVelocities[2][0] + nodesVelocities[0][1] * nodesVelocities[2][1] + nodesVelocities[0][1] * nodesVelocities[2][2]) /
                                                            (sqrt(pow(nodesVelocities[0][0], 2) + pow(nodesVelocities[0][1], 2) + pow(nodesVelocities[0][2], 2)) *
                                                             sqrt(pow(nodesVelocities[2][0], 2) + pow(nodesVelocities[2][1], 2) + pow(nodesVelocities[2][2], 2)));
                                        double cosAngle03 = (nodesVelocities[0][0] * nodesVelocities[3][0] + nodesVelocities[0][1] * nodesVelocities[3][1] + nodesVelocities[0][1] * nodesVelocities[3][2]) /
                                                            (sqrt(pow(nodesVelocities[0][0], 2) + pow(nodesVelocities[0][1], 2) + pow(nodesVelocities[0][2], 2)) *
                                                             sqrt(pow(nodesVelocities[3][0], 2) + pow(nodesVelocities[3][1], 2) + pow(nodesVelocities[3][2], 2)));
                                        double cosAngle12 = (nodesVelocities[1][0] * nodesVelocities[2][0] + nodesVelocities[1][1] * nodesVelocities[2][1] + nodesVelocities[1][1] * nodesVelocities[2][2]) /
                                                            (sqrt(pow(nodesVelocities[1][0], 2) + pow(nodesVelocities[1][1], 2) + pow(nodesVelocities[1][2], 2)) *
                                                             sqrt(pow(nodesVelocities[2][0], 2) + pow(nodesVelocities[2][1], 2) + pow(nodesVelocities[2][2], 2)));
                                        double cosAngle13 = (nodesVelocities[1][0] * nodesVelocities[3][0] + nodesVelocities[1][1] * nodesVelocities[3][1] + nodesVelocities[1][1] * nodesVelocities[3][2]) /
                                                            (sqrt(pow(nodesVelocities[1][0], 2) + pow(nodesVelocities[1][1], 2) + pow(nodesVelocities[1][2], 2)) *
                                                             sqrt(pow(nodesVelocities[3][0], 2) + pow(nodesVelocities[3][1], 2) + pow(nodesVelocities[3][2], 2)));
                                        double cosAngle23 = (nodesVelocities[2][0] * nodesVelocities[3][0] + nodesVelocities[2][1] * nodesVelocities[3][1] + nodesVelocities[2][1] * nodesVelocities[3][2]) /
                                                            (sqrt(pow(nodesVelocities[2][0], 2) + pow(nodesVelocities[2][1], 2) + pow(nodesVelocities[2][2], 2)) *
                                                             sqrt(pow(nodesVelocities[3][0], 2) + pow(nodesVelocities[3][1], 2) + pow(nodesVelocities[3][2], 2)));

                                        if (fabs(cosAngle01) < 0.85 || fabs(cosAngle02) < 0.85 || fabs(cosAngle03) < 0.85 || fabs(cosAngle12) < 0.85 || fabs(cosAngle13) < 0.85 || fabs(cosAngle23) < 0.85)
                                        {
                                            accepted = false;
                                            // std::cout << "The angle between the velocity vectors is too big" << std::endl;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // // to control that the element has a good shape
                    if (dimension == 3 && accepted && numrigid < 3 &&
                        (previouslyIsolatedNodes == 4 || previouslyFreeSurfaceNodes == 4 || sumIsolatedFreeSurf == 4 || numfreesurf == 4 || numisolated == 4 || (numrigid == 2 && isolatedNodesInTheElement > 1)))
                    {
                        Geometry<Node<3>> *tetrahedron = new Tetrahedra3D4<Node<3>>(vertices);
                        double Volume = tetrahedron->Volume();

                        double a1 = 0; //slope x for plane on the first triangular face of the tetrahedra (nodes A,B,C)
                        double b1 = 0; //slope y for plane on the first triangular face of the tetrahedra (nodes A,B,C)
                        double c1 = 0; //slope z for plane on the first triangular face of the tetrahedra (nodes A,B,C)
                        a1 = (nodesCoordinates[1][1] - nodesCoordinates[0][1]) * (nodesCoordinates[2][2] - nodesCoordinates[0][2]) - (nodesCoordinates[2][1] - nodesCoordinates[0][1]) * (nodesCoordinates[1][2] - nodesCoordinates[0][2]);
                        b1 = (nodesCoordinates[1][2] - nodesCoordinates[0][2]) * (nodesCoordinates[2][0] - nodesCoordinates[0][0]) - (nodesCoordinates[2][2] - nodesCoordinates[0][2]) * (nodesCoordinates[1][0] - nodesCoordinates[0][0]);
                        c1 = (nodesCoordinates[1][0] - nodesCoordinates[0][0]) * (nodesCoordinates[2][1] - nodesCoordinates[0][1]) - (nodesCoordinates[2][0] - nodesCoordinates[0][0]) * (nodesCoordinates[1][1] - nodesCoordinates[0][1]);
                        double a2 = 0; //slope x for plane on the second triangular face of the tetrahedra (nodes A,B,D)
                        double b2 = 0; //slope y for plane on the second triangular face of the tetrahedra (nodes A,B,D)
                        double c2 = 0; //slope z for plane on the second triangular face of the tetrahedra (nodes A,B,D)
                        a2 = (nodesCoordinates[1][1] - nodesCoordinates[0][1]) * (nodesCoordinates[3][2] - nodesCoordinates[0][2]) - (nodesCoordinates[3][1] - nodesCoordinates[0][1]) * (nodesCoordinates[1][2] - nodesCoordinates[0][2]);
                        b2 = (nodesCoordinates[1][2] - nodesCoordinates[0][2]) * (nodesCoordinates[3][0] - nodesCoordinates[0][0]) - (nodesCoordinates[3][2] - nodesCoordinates[0][2]) * (nodesCoordinates[1][0] - nodesCoordinates[0][0]);
                        c2 = (nodesCoordinates[1][0] - nodesCoordinates[0][0]) * (nodesCoordinates[3][1] - nodesCoordinates[0][1]) - (nodesCoordinates[3][0] - nodesCoordinates[0][0]) * (nodesCoordinates[1][1] - nodesCoordinates[0][1]);
                        double a3 = 0; //slope x for plane on the third triangular face of the tetrahedra (nodes B,C,D)
                        double b3 = 0; //slope y for plane on the third triangular face of the tetrahedra (nodes B,C,D)
                        double c3 = 0; //slope z for plane on the third triangular face of the tetrahedra (nodes B,C,D)
                        a3 = (nodesCoordinates[1][1] - nodesCoordinates[2][1]) * (nodesCoordinates[3][2] - nodesCoordinates[2][2]) - (nodesCoordinates[3][1] - nodesCoordinates[2][1]) * (nodesCoordinates[1][2] - nodesCoordinates[2][2]);
                        b3 = (nodesCoordinates[1][2] - nodesCoordinates[2][2]) * (nodesCoordinates[3][0] - nodesCoordinates[2][0]) - (nodesCoordinates[3][2] - nodesCoordinates[2][2]) * (nodesCoordinates[1][0] - nodesCoordinates[2][0]);
                        c3 = (nodesCoordinates[1][0] - nodesCoordinates[2][0]) * (nodesCoordinates[3][1] - nodesCoordinates[2][1]) - (nodesCoordinates[3][0] - nodesCoordinates[2][0]) * (nodesCoordinates[1][1] - nodesCoordinates[2][1]);
                        double a4 = 0; //slope x for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
                        double b4 = 0; //slope y for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
                        double c4 = 0; //slope z for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
                        a4 = (nodesCoordinates[0][1] - nodesCoordinates[2][1]) * (nodesCoordinates[3][2] - nodesCoordinates[2][2]) - (nodesCoordinates[3][1] - nodesCoordinates[2][1]) * (nodesCoordinates[0][2] - nodesCoordinates[2][2]);
                        b4 = (nodesCoordinates[0][2] - nodesCoordinates[2][2]) * (nodesCoordinates[3][0] - nodesCoordinates[2][0]) - (nodesCoordinates[3][2] - nodesCoordinates[2][2]) * (nodesCoordinates[0][0] - nodesCoordinates[2][0]);
                        c4 = (nodesCoordinates[0][0] - nodesCoordinates[2][0]) * (nodesCoordinates[3][1] - nodesCoordinates[2][1]) - (nodesCoordinates[3][0] - nodesCoordinates[2][0]) * (nodesCoordinates[0][1] - nodesCoordinates[2][1]);

                        double cosAngle12 = (a1 * a2 + b1 * b2 + c1 * c2) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
                        double cosAngle13 = (a1 * a3 + b1 * b3 + c1 * c3) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)));
                        double cosAngle14 = (a1 * a4 + b1 * b4 + c1 * c4) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)));
                        double cosAngle23 = (a3 * a2 + b3 * b2 + c3 * c2) / (sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
                        double cosAngle24 = (a4 * a2 + b4 * b2 + c4 * c2) / (sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
                        double cosAngle34 = (a4 * a3 + b4 * b3 + c4 * c3) / (sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)) * sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)));

                        if (fabs(cosAngle12) > 0.999 || fabs(cosAngle13) > 0.999 || fabs(cosAngle14) > 0.999 || fabs(cosAngle23) > 0.999 || fabs(cosAngle24) > 0.999 || fabs(cosAngle34) > 0.999) // if two faces are coplanar, I will erase the element (which is probably a sliver)
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
                    }

                    // // to control that the element has a good shape
                    // if (accepted && (numfreesurf > 0 || numrigid == nds))
                    //     {
                    //         Geometry<Node<3>> *tetrahedron = new Tetrahedra3D4<Node<3>>(vertices);

                    //         double Volume = tetrahedron->Volume();
                    //         double CriticalVolume = 0.01 * mrRemesh.Refine->MeanVolume;
                    //         if(Volume==0){
                    //             std::cout<<" !!!!! Volume==0 ";
                    //             array_1d<double, 3> CoorDifference = vertices[0].Coordinates() - vertices[1].Coordinates();
                    //             double SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
                    //             double meanLength = sqrt(SquaredLength) / 6.0;
                    //             CoorDifference = vertices[0].Coordinates() - vertices[2].Coordinates();
                    //             SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
                    //             meanLength += sqrt(SquaredLength) / 6.0;
                    //             CoorDifference = vertices[0].Coordinates() - vertices[3].Coordinates();
                    //             SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
                    //             meanLength += sqrt(SquaredLength) / 6.0;
                    //             CoorDifference = vertices[1].Coordinates() - vertices[2].Coordinates();
                    //             SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
                    //             meanLength += sqrt(SquaredLength) / 6.0;
                    //             CoorDifference = vertices[1].Coordinates() - vertices[3].Coordinates();
                    //             SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
                    //             meanLength += sqrt(SquaredLength) / 6.0;
                    //             CoorDifference = vertices[2].Coordinates() - vertices[3].Coordinates();
                    //             SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
                    //             meanLength += sqrt(SquaredLength) / 6.0;
                    //             Volume = pow(meanLength, 3) * sqrt(2) / 12.0;
                    //             std::cout<<" now volume is  "<<Volume<<std::endl;
                    //         }
                    //         if (CriticalVolume == 0)
                    //         {
                    //             array_1d<double, 3> CoorDifference = vertices[0].Coordinates() - vertices[1].Coordinates();
                    //             double SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
                    //             double meanLength = sqrt(SquaredLength) / 6.0;
                    //             CoorDifference = vertices[0].Coordinates() - vertices[2].Coordinates();
                    //             SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
                    //             meanLength += sqrt(SquaredLength) / 6.0;
                    //             CoorDifference = vertices[0].Coordinates() - vertices[3].Coordinates();
                    //             SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
                    //             meanLength += sqrt(SquaredLength) / 6.0;
                    //             CoorDifference = vertices[1].Coordinates() - vertices[2].Coordinates();
                    //             SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
                    //             meanLength += sqrt(SquaredLength) / 6.0;
                    //             CoorDifference = vertices[1].Coordinates() - vertices[3].Coordinates();
                    //             SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
                    //             meanLength += sqrt(SquaredLength) / 6.0;
                    //             CoorDifference = vertices[2].Coordinates() - vertices[3].Coordinates();
                    //             SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
                    //             meanLength += sqrt(SquaredLength) / 6.0;
                    //             double regularTetrahedronVolume = pow(meanLength, 3) * sqrt(2) / 12.0;
                    //             CriticalVolume = 0.00001 * regularTetrahedronVolume;
                    //         }

                    //         if (fabs(Volume) < CriticalVolume)
                    //         {
                    //             accepted = false;
                    //             number_of_slivers++;
                    //         }
                    //         delete tetrahedron;
                    //     }
                    //}

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
                const unsigned int nds = (*element_begin).GetGeometry().size();

                int *OutElementList = mrRemesh.OutMesh.GetElementList();

                ModelPart::NodesContainerType &rNodes = mrModelPart.Nodes();

                //check engaged nodes
                for (int el = 0; el < OutNumberOfElements; el++)
                {
                    if (mrRemesh.PreservedElements[el])
                    {
                        for (unsigned int pn = 0; pn < nds; pn++)
                        {
                            //set vertices
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

    private:
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
        void SetAlphaForRefinedZones2D(double &MeanMeshSize, bool &increaseAlfa, double coorX, double coorY)
        {

            KRATOS_TRY
            array_1d<double, 3> RefiningBoxMinimumPoint = mrRemesh.RefiningBoxMinimumPoint;
            array_1d<double, 3> RefiningBoxMaximumPoint = mrRemesh.RefiningBoxMaximumPoint;
            array_1d<double, 3> minExternalPoint = mrRemesh.RefiningBoxMinExternalPoint;
            array_1d<double, 3> minInternalPoint = mrRemesh.RefiningBoxMinInternalPoint;
            array_1d<double, 3> maxExternalPoint = mrRemesh.RefiningBoxMaxExternalPoint;
            array_1d<double, 3> maxInternalPoint = mrRemesh.RefiningBoxMaxInternalPoint;
            double distance = 2.0 * mrRemesh.Refine->CriticalRadius;
            double seperation = 0;
            double coefficient = 0;
            if (coorX > RefiningBoxMinimumPoint[0] && coorY > RefiningBoxMinimumPoint[1] &&
                coorX < RefiningBoxMaximumPoint[0] && coorY < RefiningBoxMaximumPoint[1])
            {
                MeanMeshSize = mrRemesh.RefiningBoxMeshSize;
            }
            else if (coorX < RefiningBoxMinimumPoint[0] && coorX > (minExternalPoint[0] - distance) && coorY > minExternalPoint[1] && coorY < maxExternalPoint[1])
            {
                seperation = coorX - RefiningBoxMinimumPoint[0];
                coefficient = fabs(seperation) / (distance + MeanMeshSize);
                MeanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
                increaseAlfa = true;
            }
            else if (coorY < RefiningBoxMinimumPoint[1] && coorY > (minExternalPoint[1] - distance) && coorX > minExternalPoint[0] && coorX < maxExternalPoint[0])
            {
                seperation = coorY - RefiningBoxMinimumPoint[1];
                coefficient = fabs(seperation) / (distance + MeanMeshSize);
                MeanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
                increaseAlfa = true;
            }
            else if (coorX > RefiningBoxMaximumPoint[0] && coorX < (maxExternalPoint[0] + distance) && coorY > minExternalPoint[1] && coorY < maxExternalPoint[1])
            {
                seperation = coorX - RefiningBoxMaximumPoint[0];
                coefficient = fabs(seperation) / (distance + MeanMeshSize);
                MeanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
                increaseAlfa = true;
            }
            else if (coorY > RefiningBoxMaximumPoint[1] && coorY < (maxExternalPoint[1] + distance) && coorX > minExternalPoint[0] && coorX < maxExternalPoint[0])
            {
                seperation = coorY - RefiningBoxMaximumPoint[1];
                coefficient = fabs(seperation) / (distance + MeanMeshSize);
                MeanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
                increaseAlfa = true;
            }

            KRATOS_CATCH("")
        }

        void SetAlphaForRefinedZones3D(double &MeanMeshSize, bool &increaseAlfa, double coorX, double coorY, double coorZ)
        {

            KRATOS_TRY
            array_1d<double, 3> RefiningBoxMinimumPoint = mrRemesh.RefiningBoxMinimumPoint;
            array_1d<double, 3> RefiningBoxMaximumPoint = mrRemesh.RefiningBoxMaximumPoint;
            array_1d<double, 3> minExternalPoint = mrRemesh.RefiningBoxMinExternalPoint;
            array_1d<double, 3> minInternalPoint = mrRemesh.RefiningBoxMinInternalPoint;
            array_1d<double, 3> maxExternalPoint = mrRemesh.RefiningBoxMaxExternalPoint;
            array_1d<double, 3> maxInternalPoint = mrRemesh.RefiningBoxMaxInternalPoint;
            double distance = 2.0 * mrRemesh.Refine->CriticalRadius;
            double seperation = 0;
            double coefficient = 0;

            if (coorX > RefiningBoxMinimumPoint[0] && coorX < RefiningBoxMaximumPoint[0] &&
                coorY > RefiningBoxMinimumPoint[1] && coorY < RefiningBoxMaximumPoint[1] &&
                coorZ > RefiningBoxMinimumPoint[2] && coorZ < RefiningBoxMaximumPoint[2])
            {
                MeanMeshSize = mrRemesh.RefiningBoxMeshSize;
            }
            else if (coorX < RefiningBoxMinimumPoint[0] && coorX > (minExternalPoint[0] - distance) && coorY > minExternalPoint[1] && coorY < maxExternalPoint[1] && coorZ > minExternalPoint[2] && coorZ < maxExternalPoint[2])
            {
                seperation = coorX - RefiningBoxMinimumPoint[0];
                coefficient = fabs(seperation) / (distance + MeanMeshSize);
                MeanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
                increaseAlfa = true;
            }
            else if (coorY < RefiningBoxMinimumPoint[1] && coorY > (minExternalPoint[1] - distance) && coorX > minExternalPoint[0] && coorX < maxExternalPoint[0] && coorZ > minExternalPoint[2] && coorZ < maxExternalPoint[2])
            {
                seperation = coorY - RefiningBoxMinimumPoint[1];
                coefficient = fabs(seperation) / (distance + MeanMeshSize);
                MeanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
                increaseAlfa = true;
            }
            else if (coorZ < RefiningBoxMinimumPoint[2] && coorZ > (minExternalPoint[2] - distance) && coorX > minExternalPoint[0] && coorX < maxExternalPoint[0] && coorY > minExternalPoint[1] && coorY < maxExternalPoint[1])
            {
                seperation = coorZ - RefiningBoxMinimumPoint[2];
                coefficient = fabs(seperation) / (distance + MeanMeshSize);
                MeanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
                increaseAlfa = true;
            }
            else if (coorX > RefiningBoxMaximumPoint[0] && coorX < (maxExternalPoint[0] + distance) && coorY > minExternalPoint[1] && coorY < maxExternalPoint[1] && coorZ > minExternalPoint[2] && coorZ < maxExternalPoint[2])
            {
                seperation = coorX - RefiningBoxMaximumPoint[0];
                coefficient = fabs(seperation) / (distance + MeanMeshSize);
                MeanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
                increaseAlfa = true;
            }
            else if (coorY > RefiningBoxMaximumPoint[1] && coorY < (maxExternalPoint[1] + distance) && coorX > minExternalPoint[0] && coorX < maxExternalPoint[0] && coorZ > minExternalPoint[2] && coorZ < maxExternalPoint[2])
            {
                seperation = coorY - RefiningBoxMaximumPoint[1];
                coefficient = fabs(seperation) / (distance + MeanMeshSize);
                MeanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
                increaseAlfa = true;
            }
            else if (coorZ > RefiningBoxMaximumPoint[2] && coorZ < (maxExternalPoint[2] + distance) && coorX > minExternalPoint[0] && coorX < maxExternalPoint[0] && coorY > minExternalPoint[1] && coorY < maxExternalPoint[1])
            {
                seperation = coorZ - RefiningBoxMaximumPoint[2];
                coefficient = fabs(seperation) / (distance + MeanMeshSize);
                MeanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
                increaseAlfa = true;
            }

            KRATOS_CATCH("")
        }

        void IncreaseAlphaForRefininedZones(double &Alpha,
                                            bool increaseAlfa,
                                            unsigned int nds,
                                            unsigned int numfreesurf,
                                            unsigned int numrigid,
                                            unsigned int numisolated)
        {
            KRATOS_TRY

            if (increaseAlfa == true)
            {
                if (numfreesurf < nds && numisolated == 0)
                {
                    Alpha *= 1.275;
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
            if (numfreesurf < (0.5 * nds) && (numrigid < (0.5 * nds) && numfreesurf > 0))
            {
                if (numisolated > 0)
                {
                    Alpha *= 1.0;
                }
                else if (numfreesurf == 0)
                {
                    Alpha *= 1.1;
                }
                else
                {
                    Alpha *= 1.05;
                }
            }
            KRATOS_CATCH("")
        }

        /// Assignment operator.
        SelectMeshElementsForFluidsProcess &operator=(SelectMeshElementsForFluidsProcess const &rOther);

        /// this function is a private function

        /// Copy constructor.
        //Process(Process const& rOther);

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
