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
            auto& r_mat_properties = (*it)->GetProperties();

            if (r_geometry.PointsNumber() == 3) { // 2D version
                for (unsigned int i = 0; i < 3; i++) {


                }
            } else { // 3D version
                for (unsigned int i = 0; i < 4; i++) {


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
        const double yield_tension = rMatProps[YIELD_STRESS_T];
        const double yield_compression = rMatProps[YIELD_STRESS_C];
        const double ratio = yield_compression / yield_tension;
        const double fracture_energy_compression = fracture_energy_tension * std::pow(ratio, 2.0);
        const double density = rMatProps[DENSITY];


    }

    void ComputeTensionFactor2D(const Vector& rStressVector, double& rTensionFactor)
    {
        Vector principal_stress_vector;
        this->CalculatePrincipalStresses2D(rStressVector, principal_stress_vector);
        double SumA = 0.0, SumB = 0.0, SumC = 0.0;
        for (unsigned int cont = 0; cont < 2; cont++) {
            SumA += std::abs(principal_stress_vector[cont]);
            SumB += 0.5 * (principal_stress_vector[cont]  + std::abs(principal_stress_vector[cont]));
            SumC += 0.5 * (-principal_stress_vector[cont] + std::abs(principal_stress_vector[cont]));
        }
        rTensionFactor = SumB / SumA;
    }

    void CalculatePrincipalStresses2D(const Vector& rStressVector, Vector& rPrincipalStressVector)
    {
        rPrincipalStressVector.resize(2);
        rPrincipalStressVector[0] = 0.5 * (rStressVector[0] + rStressVector[1]) + 
            std::sqrt(std::pow(0.5 * (rStressVector[0] - rStressVector[1]), 2) + std::pow(rStressVector[2], 2));
    	rPrincipalStressVector[1] = 0.5 * (rStressVector[0] + rStressVector[1]) - 
		std::sqrt(std::pow(0.5 * (rStressVector[0] - rStressVector[1]), 2) + std::pow(rStressVector[2], 2));
    }

    void CalculatePrincipalStresses3D(const Vector& rStressVector, Vector& rPrincipalStressVector)
    {
        rPrincipalStressVector.resize(3);
        const double I1 = this->CalculateI1Invariant(rStressVector);
        const double I2 = this->CalculateI2Invariant(rStressVector);
        const double I3 = this->CalculateI3Invariant(rStressVector);
        const double II1 = I1 * I1;

        const double Num = (2.0 * II1 - 9.0 * I2) * I1 + 27.0 * I3;
        const double Denom = (II1 - 3.0 * I2);

        if (Denom != 0.0) {
            double phi = Num / (2.0 * Denom * std::sqrt(Denom));

            if (std::abs(phi) > 1.0) {
                if (phi > 0.0)
                    phi = 1.0;
                else
                    phi = -1.0;
            }

            const double acosphi = std::acos(phi);
            phi = acosphi / 3.0;

            const double aux1 = 2.0 / 3.0 * std::sqrt(II1 - 3.0 * I2);
            const double aux2 = I1 / 3.0;

            rPrincipalStressVector[0] = aux2 + aux1 * std::cos(phi);
            rPrincipalStressVector[1] = aux2 + aux1 * std::cos(phi - 2.09439510239);
            rPrincipalStressVector[2] = aux2 + aux1 * std::cos(phi - 4.18879020478);
        } else {
            rPrincipalStressVector = ZeroVector(3);
        }
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

    double CalculateI1Invariant(const Vector& rStressVector)
    {
        return rStressVector[0] + rStressVector[1] + rStressVector[2];
    }

    double CalculateI2Invariant(const Vector& rStressVector)
    {
        return (rStressVector[0] + rStressVector[2]) * rStressVector[1] + rStressVector[0] * rStressVector[2] +
            -rStressVector[3] * rStressVector[3] - rStressVector[4] * rStressVector[4] - rStressVector[5] * rStressVector[5];
    }

    double CalculateI3Invariant(const Vector& rStressVector)
    {
        return (rStressVector[1] * rStressVector[2] - rStressVector[4] * rStressVector[4]) * rStressVector[0] -
            rStressVector[1] * rStressVector[5] * rStressVector[5] - rStressVector[2] * rStressVector[3] * rStressVector[3] +
            2.0 * rStressVector[3] * rStressVector[4] * rStressVector[5];
    }

  protected:
    // Member Variables
    ModelPart& mrModelPart;
    unsigned int mNNodes;

}; // Class ComputeNormalizedFreeEnergyOnNodesProcess

} // namespace Kratos

#endif /* KRATOS_COMPUTE_NORMALIZED_FREE_ENERGY_ON_NODES_PROCESS defined */