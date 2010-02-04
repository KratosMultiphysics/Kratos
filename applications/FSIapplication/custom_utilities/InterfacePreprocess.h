/*
 * File:   InterfacePreprocess.h
 * Author: jcotela
 *
 * Created on 19 January 2010, 11:44
 */

#if !defined(KRATOS_INTERFACE_PREPROCESS_MAPPER_H_INCLUDED )
#define  KRATOS_INTERFACE_PREPROCESS_MAPPER_H_INCLUDED

#include "includes/model_part.h"
#include "geometries/triangle_3d_3.h"
#include "includes/define.h"

namespace Kratos {

    class InterfacePreprocess {
        // Creates Model Parts containing the interface, to be used by AdvancedNMPointsMapper
    public:

        void GenerateInterfacePart(const ModelPart& rOriginPart, ModelPart& InterfacePart)
        // Generate a new ModelPart containing only the interface
        // It will contain only triangular conditions, regardless of what was used while meshing
        {

            // Store pointers to all interface nodes
            for (
                    ModelPart::NodesContainerType::const_iterator node_it = rOriginPart.NodesBegin();
                    node_it != rOriginPart.NodesEnd();
                    node_it++)
            {
                if (node_it->FastGetSolutionStepValue(IS_INTERFACE) == 1.0) {
                    InterfacePart.Nodes().push_back( *(node_it.base()) );
                }
            }

            // Generate triangular Conditions from original interface conditions
            ModelPart::ConditionsContainerType aux;

            for (
                    ModelPart::ConditionsContainerType::const_iterator cond_it = rOriginPart.ConditionsBegin();
                    cond_it != rOriginPart.ConditionsEnd();
                    cond_it++)
            {
                if (
                        ((*cond_it).GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE) == 1.0) &&
                        ((*cond_it).GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE) == 1.0) &&
                        ((*cond_it).GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE) == 1.0))
                {
                    aux.push_back( *(cond_it.base()) );
                }
            }

            GenerateTriangularConditions(aux,InterfacePart.Conditions());

        }

        void GenerateTriangularConditions(const ModelPart::ConditionsContainerType& rOriginConds,
                ModelPart::ConditionsContainerType& rTriConds)
        // Create a set of linear trinagular conditions from a generic condition set
        {
            // Define a condition to use as reference for all new triangle conditions
            const Condition& rCondition = KratosComponents<Condition>::Get("Face3D3N"); /*"Condition3D"*/

            // Required information for new conditions: Id, geometry and properties
            Condition::IndexType TriId = 1; // Id

            // Loop over origin conditions and create a set of triangular ones
            for (
                    ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin();
                    it != rOriginConds.end();
                    it++)
            {

                if (it->GetGeometry().PointsNumber() == 3)
                {
                    rTriConds.push_back( rCondition.Create(TriId++, it->GetGeometry(), it->pGetProperties()));
                } else if (it->GetGeometry().PointsNumber() == 4) {
//                    Triangle3D3< Node < 3 > >::Pointer Tri1(new Triangle3D3< Node<3> >(it->GetGeometry()[0], it->GetGeometry()[1], it->GetGeometry()[2]));
//                    Triangle3D3< Node < 3 > >::Pointer Tri2(new Triangle3D3< Node<3> >(it->GetGeometry()[2], it->GetGeometry()[3], it->GetGeometry()[0]));
                    Triangle3D3< Node<3> > Tri1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(2));
                    Triangle3D3< Node<3> > Tri2(it->GetGeometry()(2), it->GetGeometry()(3), it->GetGeometry()(0));

                    rTriConds.push_back(rCondition.Create(TriId++, Tri1, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri2, it->pGetProperties()));
                } else if (it->GetGeometry().PointsNumber() == 6) {
                    Triangle3D3< Node<3> > Tri1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(5));
                    Triangle3D3< Node<3> > Tri2(it->GetGeometry()(1), it->GetGeometry()(2), it->GetGeometry()(3));
                    Triangle3D3< Node<3> > Tri3(it->GetGeometry()(1), it->GetGeometry()(3), it->GetGeometry()(5));
                    Triangle3D3< Node<3> > Tri4(it->GetGeometry()(3), it->GetGeometry()(4), it->GetGeometry()(5));

                    rTriConds.push_back(rCondition.Create(TriId++, Tri1, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri2, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri3, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri4, it->pGetProperties()));
                } else if (it->GetGeometry().PointsNumber() == 8) {
                    Triangle3D3< Node<3> > Tri1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(7));
                    Triangle3D3< Node<3> > Tri2(it->GetGeometry()(1), it->GetGeometry()(5), it->GetGeometry()(7));
                    Triangle3D3< Node<3> > Tri3(it->GetGeometry()(1), it->GetGeometry()(3), it->GetGeometry()(5));
                    Triangle3D3< Node<3> > Tri4(it->GetGeometry()(1), it->GetGeometry()(2), it->GetGeometry()(3));
                    Triangle3D3< Node<3> > Tri5(it->GetGeometry()(3), it->GetGeometry()(4), it->GetGeometry()(5));
                    Triangle3D3< Node<3> > Tri6(it->GetGeometry()(5), it->GetGeometry()(6), it->GetGeometry()(7));

                    rTriConds.push_back(rCondition.Create(TriId++, Tri1, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri2, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri3, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri4, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri5, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri6, it->pGetProperties()));
                } else if (it->GetGeometry().PointsNumber() == 9) {
                    Triangle3D3< Node<3> > Tri1(it->GetGeometry()(0), it->GetGeometry()(1), it->GetGeometry()(8));
                    Triangle3D3< Node<3> > Tri2(it->GetGeometry()(1), it->GetGeometry()(2), it->GetGeometry()(3));
                    Triangle3D3< Node<3> > Tri3(it->GetGeometry()(1), it->GetGeometry()(3), it->GetGeometry()(8));
                    Triangle3D3< Node<3> > Tri4(it->GetGeometry()(8), it->GetGeometry()(3), it->GetGeometry()(4));
                    Triangle3D3< Node<3> > Tri5(it->GetGeometry()(8), it->GetGeometry()(4), it->GetGeometry()(5));
                    Triangle3D3< Node<3> > Tri6(it->GetGeometry()(5), it->GetGeometry()(6), it->GetGeometry()(7));
                    Triangle3D3< Node<3> > Tri7(it->GetGeometry()(5), it->GetGeometry()(7), it->GetGeometry()(8));
                    Triangle3D3< Node<3> > Tri8(it->GetGeometry()(0), it->GetGeometry()(8), it->GetGeometry()(7));

                    rTriConds.push_back(rCondition.Create(TriId++, Tri1, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri2, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri3, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri4, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri5, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri6, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri7, it->pGetProperties()));
                    rTriConds.push_back(rCondition.Create(TriId++, Tri8, it->pGetProperties()));
                }
            }
        }
    };
}

#endif
