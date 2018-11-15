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

#if !defined(KRATOS_DAMAGE_TO_NODES_PROCESS)
#define KRATOS_DAMAGE_TO_NODES_PROCESS

#include "processes/process.h"

namespace Kratos
{

class DamageToNodesProcess : public Process
{

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  protected:
    struct NodeStresses
    {
        int NElems;
        double Damage;

        NodeStresses()
        {
            NElems = 0;
            Damage = 0.0;
        }
    };
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  public:
    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(DamageToNodesProcess);

    typedef ModelPart::ElementsContainerType ElementsArrayType;

    // Constructor
    DamageToNodesProcess(ModelPart &r_model_part, int Dimension) : mr_model_part(r_model_part)
    {
        mNNodes = mr_model_part.NumberOfNodes();
        mNElements = mr_model_part.NumberOfElements();
        mDimension = Dimension;
    }

    // Destructor
    virtual ~DamageToNodesProcess() {}

    // --------------------------------------------------------------------
    void Execute()
    {
        int max_id;
        this->ObtainMaximumNodeId(max_id);

        NodeStresses *NodeStressesVector = new NodeStresses[max_id];

        this->DamageExtrapolationAndSmoothing(NodeStressesVector);

        delete[] NodeStressesVector;
    }

    // --------------------------------------------------------------------
    void DamageExtrapolationAndSmoothing(NodeStresses *pNodeStressesVector)
    {
        double damage;
        // Loop over elements to extrapolate the stress to the nodes
        for (ElementsArrayType::ptr_iterator it = mr_model_part.Elements().ptr_begin(); it != mr_model_part.Elements().ptr_end(); ++it)
        {
            bool condition_is_active = true;
            if ((*it)->IsDefined(ACTIVE))
            {
                condition_is_active = (*it)->Is(ACTIVE);
            }

            if (condition_is_active)
            {
                if ((*it)->GetGeometry().PointsNumber() == 3)
                {
                    damage = (*it)->GetValue(DAMAGE_ELEMENT);

                    for (int i = 0; i < 3; i++)
                    {
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].NElems += 1;
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].Damage += damage;
                    }
                }
                else
                {
                    damage = (*it)->GetValue(DAMAGE_ELEMENT);

                    for (int i = 0; i < 4; i++)
                    {
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].NElems += 1;
                        pNodeStressesVector[(*it)->GetGeometry().GetPoint(i).Id() - 1].Damage += damage;
                    }
                }
            }
        }

        // Ponderate over the elements coincident on that node
        for (unsigned int i = 0; i < mNNodes; i++)
        {
            pNodeStressesVector[i].Damage = pNodeStressesVector[i].Damage / pNodeStressesVector[i].NElems;
        }

        // Loop to compute the max eq. stress in order to normalize
        for (ModelPart::NodeIterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); ++it)
        {
            const int Id = (*it).Id();
            double &damage = it->GetSolutionStepValue(NODAL_DAMAGE);
            damage = pNodeStressesVector[Id - 1].Damage;
        }
    }

    void ObtainMaximumNodeId(int &rmax_id)
    {
        int aux = 0;
        int id;

        for (ModelPart::NodeIterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); ++it)
        {
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