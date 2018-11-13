//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velázquez
//

#if !defined(KRATOS_STRESS_TO_NODES_PROCESS)
#define KRATOS_STRESS_TO_NODES_PROCESS

#include <fstream>
#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"
#include "fem_to_dem_application_variables.h"
#include "includes/kratos_flags.h"
#include "includes/define.h"

namespace Kratos
{

class StressToNodesProcess : public Process
{

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  protected:
    struct NodeStresses
    {
        Vector EffectiveStressVector;
        int NElems;
        double Damage;

        NodeStresses()
        {
            EffectiveStressVector = ZeroVector(6); // 6 components to allow 3D case
            // 2D: sigma = [Sxx,Syy, 0, Sxy, 0, 0]
            // 3D: sigma = [Sxx,Syy,Szz,Sxy,Syz,Sxz]
            NElems = 0;
            Damage = 0.0;
        }
    };
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  public:
    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(StressToNodesProcess);

    typedef ModelPart::ElementsContainerType ElementsArrayType;

    // Constructor
    StressToNodesProcess(ModelPart &r_model_part, int Dimension) : mr_model_part(r_model_part)
    {
        mNNodes = mr_model_part.NumberOfNodes();
        mNElements = mr_model_part.NumberOfElements();
        mDimension = Dimension;
    }

    // Destructor
    virtual ~StressToNodesProcess() {}

    // --------------------------------------------------------------------
    void Execute()
    {
        int max_id;
        this->ObtainMaximumNodeId(max_id);
        NodeStresses *NodeStressesVector = new NodeStresses[max_id];
        this->StressExtrapolationAndSmoothing(NodeStressesVector);
        delete[] NodeStressesVector;
    }

    // --------------------------------------------------------------------
    void StressExtrapolationAndSmoothing(NodeStresses *pNodeStressesVector)
    {

        Vector GaussPointsStresses;
        // Loop over elements to extrapolate the stress to the nodes
        for (ElementsArrayType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it) {
            bool condition_is_active = true;
            if ((*it)->IsDefined(ACTIVE)) {
                condition_is_active = (*it)->Is(ACTIVE);
            }
            if (condition_is_active) {
                if ((*it)->GetGeometry().PointsNumber() == 3) {
                    GaussPointsStresses = (*it)->GetValue(STRESS_VECTOR);
                    for (unsigned int i = 0; i < 3; i++) {
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].EffectiveStressVector[0] += GaussPointsStresses[0];
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].EffectiveStressVector[1] += GaussPointsStresses[1];
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].EffectiveStressVector[3] += GaussPointsStresses[2];
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].NElems += 1;
                    }
                } else {
                    GaussPointsStresses = (*it)->GetValue(STRESS_VECTOR);

                    for (unsigned int i = 0; i < 4; i++) {
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].EffectiveStressVector[0] += GaussPointsStresses[0];
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].EffectiveStressVector[1] += GaussPointsStresses[1];
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].EffectiveStressVector[2] += GaussPointsStresses[2];
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].EffectiveStressVector[3] += GaussPointsStresses[3];
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].EffectiveStressVector[4] += GaussPointsStresses[4];
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].EffectiveStressVector[5] += GaussPointsStresses[5];
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].NElems += 1;
                    }
                }
            }
        }

        // Ponderate over the elements coincident on that node
        for (unsigned int i = 0; i < mNNodes; i++) {
            pNodeStressesVector[i].EffectiveStressVector = pNodeStressesVector[i].EffectiveStressVector / pNodeStressesVector[i].NElems;
        }

        // Loop to compute the max eq. stress in order to normalize
        for (ModelPart::NodeIterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); ++it) {
            int Id = (*it).Id();
            const Vector nodal_stress = pNodeStressesVector[Id - 1].EffectiveStressVector;

            // Compute the norm of the vector -> Put VonMises stress
            double &norm = it->GetSolutionStepValue(EQUIVALENT_NODAL_STRESS);
            norm = this->CalculateStressInvariant(nodal_stress);
        }

        // Loop to compute the max eq. stress and normalize
        double MaxEqStress = 0.0;
        for (ModelPart::NodeIterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); ++it) {
            double &norm = it->GetSolutionStepValue(EQUIVALENT_NODAL_STRESS);
            if (norm > MaxEqStress)
                MaxEqStress = norm;
        }

        for (ModelPart::NodeIterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); ++it) {
            double &norm = it->GetSolutionStepValue(EQUIVALENT_NODAL_STRESS);
            norm /= (MaxEqStress * 1.0e4);
            // norm /= (MaxEqStress);
        }
    }

    // --------------------------------------------------------------------
    double CalculateStressInvariant(const Vector &StressVector) // Returns Von Mises stress
    {
        const double I1 = StressVector[0] + StressVector[1] + StressVector[2];

        Vector Deviator = StressVector;
        const double Pmean = I1 / 3.0;

        Deviator[0] -= Pmean;
        Deviator[1] -= Pmean;
        Deviator[2] -= Pmean;

        const double J2 = 0.5 * (Deviator[0] * Deviator[0] + Deviator[1] * Deviator[1] + Deviator[2] * Deviator[2]) +
                    (Deviator[3] * Deviator[3] + Deviator[4] * Deviator[4] + Deviator[5] * Deviator[5]);

        return std::sqrt(3.0 * J2);
    }

    void ObtainMaximumNodeId(int &rmax_id)
    {
        int aux = 0;
        int id;

        for (ModelPart::NodeIterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); ++it) {
            id = (*it).Id();
            if (id > aux)
                aux = id;
        }
        rmax_id = aux;
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  protected:
    // Member Variables
    ModelPart &mr_model_part;
    unsigned int mNNodes, mNElements, mDimension;

}; // Class

} // namespace Kratos

#endif /* KRATOS_STRESS_TO_NODES_PROCESS defined */