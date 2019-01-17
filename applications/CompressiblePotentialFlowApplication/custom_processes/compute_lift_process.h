//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#ifndef KRATOS_COMPUTE_LIFT_PROCESS_H
#define KRATOS_COMPUTE_LIFT_PROCESS_H


#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "processes/process.h"
#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"

#include <string>
#include <iostream>
#include <sstream>

#include <boost/functional/hash.hpp> //TODO: remove this dependence when Kratos has en internal one
#include <unordered_map> //TODO: remove this dependence when Kratos has en internal one
#include <utility>

namespace Kratos
{

class ComputeLiftProcess: public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(ComputeLiftProcess);

    typedef ModelPart::ElementType ElementType;
    typedef ModelPart::ConditionType ConditionType;
    typedef Node < 3 > NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    typedef Element::WeakPointer ElementWeakPointerType;
    
    typedef Element::Pointer ElementPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for ComputeLiftProcess Process
    ComputeLiftProcess(ModelPart& rModelPart,
                Vector& rResultForce                             
                    ):
        Process(),
        mrModelPart(rModelPart),
        mrResultForce(rResultForce)
    {
    }

    /// Destructor.
    ~ComputeLiftProcess() override {}


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


    /// Check elements to make sure that their jacobian is positive and conditions to ensure that their face normals point outwards
    void Execute() override
    {
        KRATOS_TRY;
        double Cl=0;
        double Cd=0;
        double Rz=0;

        for(auto it=mrModelPart.ConditionsBegin(); it!=mrModelPart.ConditionsEnd(); ++it)
        {
            if (it->IsNot(BOUNDARY)){
                GeometryType &rGeom = it->GetGeometry();
                WeakPointerVector<Element> ElementCandidates;
                GetElementCandidates(ElementCandidates, rGeom);

                std::vector<IndexType> NodeIds, ElementNodeIds;
                GetSortedIds(NodeIds, rGeom);
                Element::WeakPointer pElem=FindParentElement(NodeIds, ElementNodeIds, ElementCandidates);

                std::vector<double> cp;
                Vector normal;
                IndexType elem_id = pElem.lock()->Id();
                mrModelPart.pGetElement(elem_id,0)->GetValueOnIntegrationPoints(PRESSURE,cp,mrModelPart.GetProcessInfo());
                double cpressure=cp[0];
                normal=it->GetValue(NORMAL);
                                
                Cl += cpressure*normal(1);
                Cd += cpressure*normal(0);
                Rz += cpressure*normal(2);
            }
        }
        mrResultForce[0]=Cd;
        mrResultForce[1]=Cl;
        mrResultForce[2]=Rz;
                 
        KRATOS_CATCH("");
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
        return "ComputeLiftProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeLiftProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}


private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    Vector& mrResultForce;
    Flags mrOptions;



    void GetElementCandidates(WeakPointerVector<Element> &ElementCandidates, GeometryType &rGeom)
    {
        for (SizeType i = 0; i < rGeom.WorkingSpaceDimension(); i++)
        {
            if (rGeom[i].Has(NEIGHBOUR_ELEMENTS)) {
                WeakPointerVector<Element> &rNodeElementCandidates = rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);
                for (SizeType j = 0; j < rNodeElementCandidates.size(); j++)
                    ElementCandidates.push_back(rNodeElementCandidates(j));
            } else
                std::cout << "fuck" << std::endl; 
        }
    }

    void GetSortedIds(std::vector<IndexType> &Ids,
                        const GeometryType &rGeom)
    {
        Ids.resize(rGeom.PointsNumber());
        for (SizeType i = 0; i < Ids.size(); i++)
            Ids[i] = rGeom[i].Id();
        std::sort(Ids.begin(), Ids.end());
    }

    Element::WeakPointer FindParentElement(std::vector<IndexType> &NodeIds,
                            std::vector<IndexType> &ElementNodeIds,
                            WeakPointerVector<Element> ElementCandidates)
    {
        for (SizeType i = 0; i < ElementCandidates.size(); i++)
        {
            
            GeometryType &rElemGeom = ElementCandidates[i].GetGeometry();
            
            GetSortedIds(ElementNodeIds, rElemGeom);

            if (std::includes(ElementNodeIds.begin(), ElementNodeIds.end(), NodeIds.begin(), NodeIds.end()))
            {
                return ElementCandidates(i);
            }
        }
        Element::WeakPointer void_ptr;
        return void_ptr;
    }


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ComputeLiftProcess& operator=(ComputeLiftProcess const& rOther);

    /// Copy constructor.
    ComputeLiftProcess(ComputeLiftProcess const& rOther);


    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ComputeLiftProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ComputeLiftProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



} // namespace Kratos


#endif // KRATOS_ComputeLift_LEVEL_SET_PROCESS_H
