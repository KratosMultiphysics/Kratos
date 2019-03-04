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
    ComputeNormalizedFreeEnergyOnNodesProcess(ModelPart &r_model_part, int Dimension) : mrModelPart(r_model_part)
    {
        mNNodes = mrModelPart.NumberOfNodes();
        mDimension = Dimension;
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

            const Vector& r_strain_vector = (*it)->GetValue(STRAIN_VECTOR);
            const Vector& r_stress_vector = (*it)->GetValue(STRESS_VECTOR);
            const double damage = (*it)->GetValue(DAMAGE_ELEMENT);

            // Compute Normalized Free Energy on that element
            const double normalized_free_energy = this->CalculateNormalizedFreeEnergy(r_strain_vector,
                                                                                      r_stress_vector,
                                                                                      damage,
                                                                                      r_mat_properties,
                                                                                      r_geometry);

            for (unsigned int i = 0; i < r_geometry.PointsNumber(); i++) {
                const int node_id = r_geometry.GetPoint(i).Id();
                pNodeNormalizedFreeEnergyVector[node_id - 1].NormalizedFreeEnergy += normalized_free_energy;
                pNodeNormalizedFreeEnergyVector[node_id - 1].NElems += 1;
            }
        }

        // Ponderate over the elements coincident on that node
        for (unsigned int i = 0; i < mNNodes; i++) {
            pNodeNormalizedFreeEnergyVector[i].NormalizedFreeEnergy /= pNodeNormalizedFreeEnergyVector[i].NElems;
        }

        // Loop over nodes to assign the variable
        for (ModelPart::NodeIterator it = mrModelPart.NodesBegin(); it != mrModelPart.NodesEnd(); ++it) {
            int Id = (*it).Id();
            const double nodal_free_energy = pNodeNormalizedFreeEnergyVector[Id - 1].NormalizedFreeEnergy;

            double &norm = it->GetSolutionStepValue(EQUIVALENT_NODAL_STRESS);
            norm = nodal_free_energy;
        }
    }

    double CalculateNormalizedFreeEnergy(
        const Vector& rStrainVector, 
        const Vector& rStressVector, 
        const double Damage, 
        const Properties& rMatProps,
        Geometry<Node<3>>& rGeometry)
    {
        const double fracture_energy_tension = rMatProps[FRAC_ENERGY_T];
        const double yield_tension = rMatProps[YIELD_STRESS_T];
        const double yield_compression = rMatProps[YIELD_STRESS_C];
        const double ratio = yield_compression / yield_tension;
        const double fracture_energy_compression = fracture_energy_tension * std::pow(ratio, 2.0);
        const double density = rMatProps[DENSITY];
        double normalized_free_energy;

        if (mDimension == 2) { // 2D version
            const double characteristic_length = this->CalculateCharacteristicLength2D(rGeometry);
            const double g_t = fracture_energy_tension / characteristic_length;
            const double g_c = fracture_energy_compression / characteristic_length;
            const double r   = this->ComputeTensionFactor2D(rStressVector);
            normalized_free_energy = inner_prod(rStrainVector, rStressVector);
            normalized_free_energy *= (1.0 - Damage);
            normalized_free_energy /= (2.0 * density);
            normalized_free_energy *= (r / g_t + (1.0 - r) / g_c);
            return normalized_free_energy;
        } else if (mDimension == 3) { // 3D version
            const double characteristic_length = this->CalculateCharacteristicLength3D(rGeometry);
            const double g_t = fracture_energy_tension / characteristic_length;
            const double g_c = fracture_energy_compression / characteristic_length;
            const double r   = this->ComputeTensionFactor3D(rStressVector);
            normalized_free_energy = inner_prod(rStrainVector, rStressVector);
            normalized_free_energy *= (1.0 - Damage);
            normalized_free_energy /= (2.0 * density);
            normalized_free_energy *= (r / g_t + (1.0 - r) / g_c);
            return normalized_free_energy;
        } else {
            KRATOS_ERROR << "Dimension is wrong..." << std::endl;
        }
    }

	double CalculateCharacteristicLength2D(const Geometry<Node<3>>& rGeometry)
	{
        Vector node_1_coordinates = ZeroVector(3);
        Vector node_2_coordinates = ZeroVector(3);
        Vector node_3_coordinates = ZeroVector(3);
        node_1_coordinates[0] = rGeometry[0].X0();
        node_1_coordinates[1] = rGeometry[0].Y0();
        node_1_coordinates[2] = rGeometry[0].Z0();
        node_2_coordinates[0] = rGeometry[1].X0();
        node_2_coordinates[1] = rGeometry[1].Y0();
        node_2_coordinates[2] = rGeometry[1].Z0();
        node_3_coordinates[0] = rGeometry[2].X0();
        node_3_coordinates[1] = rGeometry[2].Y0();
        node_3_coordinates[2] = rGeometry[2].Z0();

        const double length1 = MathUtils<double>::Norm3(node_1_coordinates-node_2_coordinates);
        const double length2 = MathUtils<double>::Norm3(node_2_coordinates-node_3_coordinates);
        const double length3 = MathUtils<double>::Norm3(node_3_coordinates-node_1_coordinates);
        return (length1 + length2 + length3) / 3.0;
	}

	double CalculateCharacteristicLength3D(Geometry<Node<3>>& rGeometry)
	{
    	Vector lengths = ZeroVector(6);
        auto& r_edges = rGeometry.Edges();
        for (unsigned int edge = 0; edge < 6; edge++) { // Loop over edges
            const double X1 = r_edges[edge][0].X0();
            const double X2 = r_edges[edge][1].X0();
            const double Y1 = r_edges[edge][0].Y0();
            const double Y2 = r_edges[edge][1].Y0();
            const double Z1 = r_edges[edge][0].Z0();
            const double Z2 = r_edges[edge][1].Z0();
            lengths[edge] = std::sqrt(std::pow((X1 - X2), 2.0) + std::pow((Y1 - Y2), 2.0) + std::pow((Z1 - Z2), 2.0));
        }
        return (lengths[0] + lengths[1] + lengths[2] + lengths[3] + lengths[4] + lengths[5]) / 6.0;
	}

    double ComputeTensionFactor2D(const Vector& rStressVector)
    {
        Vector principal_stress_vector;
        this->CalculatePrincipalStresses2D(rStressVector, principal_stress_vector);
        double SumA = 0.0, SumB = 0.0, SumC = 0.0;
        for (unsigned int cont = 0; cont < 2; cont++) {
            SumA += std::abs(principal_stress_vector[cont]);
            SumB += 0.5 * (principal_stress_vector[cont]  + std::abs(principal_stress_vector[cont]));
            SumC += 0.5 * (-principal_stress_vector[cont] + std::abs(principal_stress_vector[cont]));
        }
        return SumB / SumA;
    }

    double ComputeTensionFactor3D(const Vector& rStressVector)
    {
        Vector principal_stress_vector;
        this->CalculatePrincipalStresses3D(rStressVector, principal_stress_vector);
        double SumA = 0.0, SumB = 0.0, SumC = 0.0;
        for (unsigned int cont = 0; cont < 3; cont++) {
            SumA += std::abs(principal_stress_vector[cont]);
            SumB += 0.5 * (principal_stress_vector[cont]  + std::abs(principal_stress_vector[cont]));
            SumC += 0.5 * (-principal_stress_vector[cont] + std::abs(principal_stress_vector[cont]));
        }
        return SumB / SumA;
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

        const double numerator = (2.0 * II1 - 9.0 * I2) * I1 + 27.0 * I3;
        const double denominator = (II1 - 3.0 * I2);

        if (denominator != 0.0) {
            double phi = numerator / (2.0 * denominator * std::sqrt(denominator));

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
    unsigned int mNNodes, mDimension;

}; // Class ComputeNormalizedFreeEnergyOnNodesProcess

} // namespace Kratos

#endif /* KRATOS_COMPUTE_NORMALIZED_FREE_ENERGY_ON_NODES_PROCESS defined */