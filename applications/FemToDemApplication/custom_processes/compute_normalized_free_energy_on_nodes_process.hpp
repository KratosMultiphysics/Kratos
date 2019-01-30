//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#if !defined(KRATOS_COMPUTE_NORMALIZED_FREE_ENERGY_ON_NODES_PROCESS)
#define KRATOS_COMPUTE_NORMALIZED_FREE_ENERGY_ON_NODES_PROCESS

#include "processes/process.h"

namespace Kratos
{
class ComputeNormalizedFreeEnergyOnNodesProcess : public Process
{
  protected:
    struct NodeNormalizedFreeEnergy
    {
        int NElems;
        double NormalizedFreeEnergy;

		NodeNormalizedFreeEnergy()
        {
            NElems = 0;
            NormalizedFreeEnergy = 0.0;
        }
    };

  public:
    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeNormalizedFreeEnergyOnNodesProcess);

    typedef ModelPart::ElementsContainerType ElementsArrayType;

    // Constructor
    ComputeNormalizedFreeEnergyOnNodesProcess(ModelPart &r_model_part) : mrModelPart(r_model_part)
    {
        mNNodes = mrModelPart.NumberOfNodes();
    }

    // Destructor
    virtual ~ComputeNormalizedFreeEnergyOnNodesProcess() {}

    void Execute()
    {
        int max_id;
        this->ObtainMaximumNodeId(max_id);
        NodeNormalizedFreeEnergy *NodeNormalizedFreeEnergyVector = new NodeNormalizedFreeEnergy[max_id];
        this->NormalizedFreeEnergyExtrapolation(NodeNormalizedFreeEnergyVector);
        delete[] NodeNormalizedFreeEnergyVector;
    }

    void NormalizedFreeEnergyExtrapolation(NodeNormalizedFreeEnergy *pNodeNormalizedFreeEnergyVector)
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
                    // const int node_id = r_geometry.GetPoint(i).Id();
                    // pNodeStressesVector[node_id - 1].EffectiveStressVector[0] += gauss_point_stress[0];
                    // pNodeStressesVector[node_id - 1].EffectiveStressVector[1] += gauss_point_stress[1];
                    // pNodeStressesVector[node_id - 1].EffectiveStressVector[3] += gauss_point_stress[2];
                    // pNodeStressesVector[node_id - 1].Damage += gauss_point_damage;
                    // pNodeStressesVector[node_id - 1].NElems += 1;
                }
            } else { // 3D version
                for (unsigned int i = 0; i < 4; i++) {
                    // const int node_id = r_geometry.GetPoint(i).Id();
                    // pNodeStressesVector[node_id - 1].EffectiveStressVector[0] += gauss_point_stress[0];
                    // pNodeStressesVector[node_id - 1].EffectiveStressVector[1] += gauss_point_stress[1];
                    // pNodeStressesVector[node_id - 1].EffectiveStressVector[2] += gauss_point_stress[2];
                    // pNodeStressesVector[node_id - 1].EffectiveStressVector[3] += gauss_point_stress[3];
                    // pNodeStressesVector[node_id - 1].EffectiveStressVector[4] += gauss_point_stress[4];
                    // pNodeStressesVector[node_id - 1].EffectiveStressVector[5] += gauss_point_stress[5];
                    // pNodeStressesVector[node_id - 1].Damage += gauss_point_damage;
                    // pNodeStressesVector[node_id - 1].NElems += 1;
                }
            }
        }

        // Ponderate over the elements coincident on that node
        for (unsigned int i = 0; i < mNNodes; i++) {
            pNodeNormalizedFreeEnergyVector[i].NormalizedFreeEnergy /= pNodeNormalizedFreeEnergyVector[i].NElems;
        }
    }

	double CalculateCharacteristicLength(Geometry<Node<3>>& rGeometry)
	{

	}

    double CalculateNormalizedFreeEnergy(
        const Vector& rStrainVector, 
        const Vector& rStressVector, 
        const double damage, 
        const Properties& rMatProps,
        const Geometry<Node<3>>& rGeometry)
    {
        const double fracture_energy_tension = rMatProps[FRAC_ENERGY_T];
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
    unsigned int mNNodes;

}; // Class ComputeNormalizedFreeEnergyOnNodesProcess

} // namespace Kratos

#endif /* KRATOS_COMPUTE_NORMALIZED_FREE_ENERGY_ON_NODES_PROCESS defined */