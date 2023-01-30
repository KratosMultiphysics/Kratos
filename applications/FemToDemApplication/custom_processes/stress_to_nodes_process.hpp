//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#if !defined(KRATOS_STRESS_TO_NODES_PROCESS)
#define KRATOS_STRESS_TO_NODES_PROCESS

#include "processes/process.h"

namespace Kratos
{

class StressToNodesProcess : public Process
{
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

  public:
    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(StressToNodesProcess);

    typedef ModelPart::ElementsContainerType ElementsArrayType;

    // Constructor
    StressToNodesProcess(ModelPart &r_model_part, int Dimension) : mrModelPart(r_model_part)
    {
        mNNodes = mrModelPart.NumberOfNodes();
        mNElements = mrModelPart.NumberOfElements();
        mDimension = Dimension;
    }

    // Destructor
    virtual ~StressToNodesProcess() {}

    void Execute()
    {
        int max_id;
        this->ObtainMaximumNodeId(max_id);
        NodeStresses *NodeStressesVector = new NodeStresses[max_id];
        this->StressExtrapolationAndSmoothing(NodeStressesVector);
        delete[] NodeStressesVector;
    }

    void StressExtrapolationAndSmoothing(NodeStresses *pNodeStressesVector)
    {
        Vector gauss_point_stress;
        double gauss_point_damage;
        // Loop over elements to extrapolate the stress to the nodes
        for (ElementsArrayType::ptr_iterator it = mrModelPart.Elements().ptr_begin(); it != mrModelPart.Elements().ptr_end(); ++it) {
            auto& r_geometry = (*it)->GetGeometry();
            gauss_point_stress = (*it)->GetValue(STRESS_VECTOR);
            gauss_point_damage = (*it)->GetValue(DAMAGE_ELEMENT);

            if (r_geometry.PointsNumber() == 3) { // 2D version
                for (unsigned int i = 0; i < 3; i++) {
                    const int node_id = r_geometry.GetPoint(i).Id();
                    pNodeStressesVector[node_id - 1].EffectiveStressVector[0] += gauss_point_stress[0];
                    pNodeStressesVector[node_id - 1].EffectiveStressVector[1] += gauss_point_stress[1];
                    pNodeStressesVector[node_id - 1].EffectiveStressVector[3] += gauss_point_stress[2];
                    pNodeStressesVector[node_id - 1].Damage += gauss_point_damage;
                    pNodeStressesVector[node_id - 1].NElems += 1;
                }
            } else { // 3D version
                for (unsigned int i = 0; i < 4; i++) {
                    const int node_id = r_geometry.GetPoint(i).Id();
                    pNodeStressesVector[node_id - 1].EffectiveStressVector[0] += gauss_point_stress[0];
                    pNodeStressesVector[node_id - 1].EffectiveStressVector[1] += gauss_point_stress[1];
                    pNodeStressesVector[node_id - 1].EffectiveStressVector[2] += gauss_point_stress[2];
                    pNodeStressesVector[node_id - 1].EffectiveStressVector[3] += gauss_point_stress[3];
                    pNodeStressesVector[node_id - 1].EffectiveStressVector[4] += gauss_point_stress[4];
                    pNodeStressesVector[node_id - 1].EffectiveStressVector[5] += gauss_point_stress[5];
                    pNodeStressesVector[node_id - 1].Damage += gauss_point_damage;
                    pNodeStressesVector[node_id - 1].NElems += 1;
                }
            }
        }

        // Ponderate over the elements coincident on that node
        for (unsigned int i = 0; i < mNNodes; i++) {
            pNodeStressesVector[i].EffectiveStressVector /= pNodeStressesVector[i].NElems;
            pNodeStressesVector[i].Damage /= pNodeStressesVector[i].NElems;
        }

        // Loop to compute the equivalent streses at each node
        for (ModelPart::NodeIterator it = mrModelPart.NodesBegin(); it != mrModelPart.NodesEnd(); ++it) {
            int Id = (*it).Id();
            const Vector& r_nodal_stress = pNodeStressesVector[Id - 1].EffectiveStressVector;

            // Compute the norm of the vector -> Put VonMises stress
            double &norm = it->GetSolutionStepValue(EQUIVALENT_NODAL_STRESS);
            norm = this->CalculateVonMisesStress(r_nodal_stress);
        }

        // Loop to compute the max eq. stress to normalize the inidicator
        double max_equivalent_stress = 0.0;
        for (ModelPart::NodeIterator it = mrModelPart.NodesBegin(); it != mrModelPart.NodesEnd(); ++it) {
            const double &norm = it->GetSolutionStepValue(EQUIVALENT_NODAL_STRESS);
            if (norm > max_equivalent_stress)
                max_equivalent_stress = norm;
        }

        for (ModelPart::NodeIterator it = mrModelPart.NodesBegin(); it != mrModelPart.NodesEnd(); ++it) {
            double &norm = it->GetSolutionStepValue(EQUIVALENT_NODAL_STRESS);
            int Id = (*it).Id();
            double nodal_damage = pNodeStressesVector[Id - 1].Damage;
            // if (nodal_damage >= 0.9 && norm < max_equivalent_stress / 1.0e2) norm = max_equivalent_stress / 1.0e2;
            norm /= (max_equivalent_stress * 1.0e4);
        }
    }

    // --------------------------------------------------------------------
    double CalculateVonMisesStress(const Vector& rStressVector) // Returns Von Mises stress
    {
        const double I1 = rStressVector[0] + rStressVector[1] + rStressVector[2];
        Vector Deviator = rStressVector;
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

        for (ModelPart::NodeIterator it = mrModelPart.NodesBegin(); it != mrModelPart.NodesEnd(); ++it) {
            id = (*it).Id();
            if (id > aux)
                aux = id;
        }
        rmax_id = aux;
    }

  protected:
    // Member Variables
    ModelPart& mrModelPart;
    unsigned int mNNodes, mNElements, mDimension;

}; // Class StressToNodesProcess

} // namespace Kratos

#endif /* KRATOS_STRESS_TO_NODES_PROCESS defined */