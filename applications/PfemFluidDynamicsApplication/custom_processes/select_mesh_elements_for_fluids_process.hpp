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
//StepData: NODAL_H, CONTACT_FORCE
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
                        checkedNodes++;
                    }
                    else if (vertices.back().Is(ISOLATED))
                    {
                        checkedNodes++;
                        const array_1d<double, 3> &velocityP0 = vertices.back().FastGetSolutionStepValue(VELOCITY, 0);
                        normVelocityP[pn] = norm_2(velocityP0);
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

                if (dimension == 2)
                {
                    if ((numfreesurf == nds || (numisolated + numfreesurf) == nds) && firstMesh == false)
                    {
                        if (checkedNodes == nds)
                        {
                            const double maxValue = 1.5;
                            const double minValue = 1.0 / maxValue;
                            if (normVelocityP[0] / normVelocityP[1] > maxValue || normVelocityP[0] / normVelocityP[1] < minValue ||
                                normVelocityP[0] / normVelocityP[2] > maxValue || normVelocityP[0] / normVelocityP[2] < minValue ||
                                normVelocityP[1] / normVelocityP[2] > maxValue || normVelocityP[1] / normVelocityP[2] < minValue)
                            {
                                Alpha *= 0;
                            }
                        }
                        else
                        {
                            KRATOS_INFO("ATTENTION!!! CHECKED NODES= ") << checkedNodes << " and the nodes are " << nds << std::endl;
                            Alpha *= 0;
                        }
                    }

                    if (numrigid == 0 && numfreesurf == 0 && numisolated == 0)
                    {
                        Alpha *= 1.75;
                    }
                    else if (numrigid > 0 && numfreesurf == 0 && numisolated == 0){
                        Alpha *= 1.1;
                    }else
                    {
                        Alpha *= 1.04;
                    }
                }
                else if (dimension == 3)
                {
                    // if (numfreesurf == nds || (numisolated + numfreesurf) == nds)
                    // {
                    //     if (checkedNodes == nds)
                    //     {
                    //         const double maxValue = 1.5;
                    //         const double minValue = 1.0 / maxValue;
                    //         if (normVelocityP[0] / normVelocityP[1] < minValue || normVelocityP[0] / normVelocityP[2] < minValue || normVelocityP[0] / normVelocityP[3] < minValue ||
                    //             normVelocityP[0] / normVelocityP[1] > maxValue || normVelocityP[0] / normVelocityP[2] < maxValue || normVelocityP[0] / normVelocityP[3] > maxValue ||
                    //             normVelocityP[1] / normVelocityP[2] < minValue || normVelocityP[1] / normVelocityP[3] < minValue ||
                    //             normVelocityP[1] / normVelocityP[2] > maxValue || normVelocityP[1] / normVelocityP[3] < maxValue ||
                    //             normVelocityP[2] / normVelocityP[3] < minValue ||
                    //             normVelocityP[2] / normVelocityP[3] > maxValue)
                    //         {
                    //             Alpha *= 0;
                    //         }
                    //     }
                    //     else
                    //     {
                    //         std::cout << "ATTENTION!!! CHECKED NODES= " << checkedNodes << " and the nodes are " << nds << std::endl;
                    //         Alpha *= 0;
                    //     }
                    // }

                    if (numrigid == 0 && numfreesurf == 0 && numisolated == 0)
                    {
                        Alpha *= 1.75;
                    }
                    else
                    {
                        Alpha *= 1.125;
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
                if (firstMesh == true)
                {
                    Alpha *= 1.15;
                }

                if (numinlet > 0)
                {
                    Alpha *= 1.5;
                }

                bool accepted = false;
                MesherUtilities MesherUtils;

                    accepted = MesherUtils.AlphaShape(Alpha, vertices, dimension, MeanMeshSize);


                if (numrigid == nds || noremesh == true)
                {
                    accepted = false;
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
                // }

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
