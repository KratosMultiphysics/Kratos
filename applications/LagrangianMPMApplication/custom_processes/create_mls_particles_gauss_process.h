//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Zhiming Guo
//                   Riccardo Rossi
//






#if !defined(KRATOS_CREATE_MLS_PARTICLE_GAUSS_H_INCLUDED )
#define  KRATOS_CREATE_MLS_PARTICLE_GAUSS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "lagrangian_mpm_application_variables.h"
#include "lagrangian_mpm_application.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
//#include "custom_elements/MLSparticle.h"
#include "geometries/point_2d.h"
#include "geometries/geometry.h"

#include "includes/variables.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"



#include "spatial_containers/spatial_search.h"
#include "spatial_containers/bins_dynamic.h"

#include "utilities/spatial_containers_configure.h"



#include "includes/constitutive_law.h"

#include "includes/properties.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
typedef  ModelPart::NodesContainerType NodesContainerType;
typedef  ModelPart::ElementsContainerType ElementsContainerType;
//typedef  ElementsContainerType temp_container;

typedef SpatialSearch::DistanceType                       DistanceType;
typedef std::vector<double>                               RadiusArrayType;
typedef std::vector<DistanceType>                         VectorDistanceType;
typedef SpatialSearch::ResultNodesContainerType           ResultNodesContainerType;
typedef std::vector<ResultNodesContainerType>             VectorResultNodesContainerType;
typedef GeometryData::IntegrationMethod IntegrationMethod;
typedef Geometry<Node<3> >::GeometryType GeometryType;
///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class CreateMLSParticleGauss
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CalculateNodalAreaProcess
    KRATOS_CLASS_POINTER_DEFINITION(CreateMLSParticleGauss);

    ///@}
    ///@name Life Cycle
    ///@{

    CreateMLSParticleGauss(ModelPart& rModelPart,
                           const std::string ElementName,
                           const std::string ConditionName,
                           unsigned int PropertyId)
        : mrModelPart(rModelPart), mElementName(ElementName), mConditionName(ConditionName), mPropertyId(PropertyId)
    {
        mDomainSize = rModelPart.GetProcessInfo()[DOMAIN_SIZE];
    }

    /// Destructor.
    virtual ~CreateMLSParticleGauss()
    {
    }


    ///@}
    ///@name Operators
    ///@{bins_dynamic.h

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    virtual void Execute()
    {
        KRATOS_TRY;

        //create a copy of the pointers to the elements and conditions in the model, and then remove all elements and conditions
        ModelPart::ElementsContainerType old_elements = mrModelPart.Elements();
        ModelPart::ConditionsContainerType old_conditions = mrModelPart.Conditions();

	
        for(auto iel = mrModelPart.ElementsBegin(); iel!=mrModelPart.ElementsEnd(); iel++)
            iel->Set(TO_ERASE, true);
        mrModelPart.RemoveElements( TO_ERASE );

        for(auto icond = mrModelPart.ConditionsBegin(); icond!=mrModelPart.ConditionsEnd(); icond++)
            icond->Set(TO_ERASE, true);
        mrModelPart.RemoveConditions( TO_ERASE );

        const Element& reference_element = KratosComponents<Element>::Get(mElementName.c_str());
//         const Condition& reference_condition = KratosComponents<Condition>::Get(mConditionName.c_str());


        if(mDomainSize == 2)
        {

            double tot_area = 0.0;

            //create a search structure for the nodes
            typedef Node<3>                                         NodeType;
            typedef Node<3>::Pointer                                NodePointerType;
            typedef ModelPart::NodesContainerType::ContainerType    NodesContainerType;
            ModelPart::ElementsContainerType                        temp_container;
            ModelPart::ConditionsContainerType                      temp_condition_container;

            BinsDynamic<3, NodeType, NodesContainerType> bins(mrModelPart.NodesArray().begin(),mrModelPart.NodesArray().end());

            unsigned int idcount = mrModelPart.Nodes().size();//0;
            double area;

            double search_radius = 0.0;
            double effective_radius = 0.0;
            double mlength = 0.0;

            Properties::Pointer properties = mrModelPart.GetMesh().pGetProperties(mPropertyId);

            ProcessInfo& rCurrentProcessInfo  =mrModelPart.GetProcessInfo();

            for(auto iel =old_elements.begin(); iel!=old_elements.end(); iel++)
            {
                unsigned int MaximumNumberOfResults = 1000;
                std::vector<NodePointerType>                            results(MaximumNumberOfResults);
                std::vector<double>                                     distances(MaximumNumberOfResults);

                //calculating shape functions values
                Geometry< Node<3> >& geom = iel->GetGeometry();

                double lmax = 0.0;
                array_1d<double,3> side;
                for(unsigned int i=0; i<geom.size(); i++)
                {
                    for(unsigned int j=i+1; j<geom.size(); j++)
                    {
                        side = geom[j] - geom[i];
                        double l = norm_2(side);
                        if(l > lmax)
                            lmax = l;
                    }
                }

                area = geom.Area();
                tot_area += area;

                search_radius = 1.0 * lmax; //2.0*lmax;
                effective_radius = 0.5 * lmax;
		        mlength = geom.Length();
		        
		        
		        

				if(iel->Id()==1)
				{
					KRATOS_WATCH(search_radius);
					KRATOS_WATCH(effective_radius);
				}

                //create an auxiliary point
                Node<3> aux_node(0,0.0,0.0,0.0);


                //get gauss points
//                 mThisIntegrationMethod = iel->GetGeometry().GetDefaultIntegrationMethod();

                const Element::GeometryType::IntegrationPointsArrayType& integration_points = iel->GetGeometry().IntegrationPoints( GeometryData::GI_GAUSS_2);
                Vector detJ0( integration_points.size() );

                const Matrix& Ncontainer = iel->GetGeometry().ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );

                iel->GetGeometry().DeterminantOfJacobian( detJ0 );

                for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
                {
                    //getting informations for integration
                    double IntegrationWeight = integration_points[PointNumber].Weight();

                    double Agauss = IntegrationWeight*detJ0[PointNumber];
                    //KRATOS_WATCH(Agauss);

                    idcount += 1;

                    //compute the coordinates of the gauss point as  the sum of N[i]*xnode
                    array_1d<double,3> coordinates_of_gauss = ZeroVector(3);
                    for(unsigned int k=0; k<geom.size(); k++)
                    {
                        coordinates_of_gauss += Ncontainer(PointNumber,k)*geom[k].Coordinates();
                    }

                    //TODO: do not create a node here!!!
                    //Node<3> search_coordinates(0,coordinates_of_gauss[0],coordinates_of_gauss[1],coordinates_of_gauss[2]);
                    aux_node.Coordinates() = coordinates_of_gauss;
                    unsigned int nresults = bins.SearchInRadius(aux_node,search_radius,results.begin(),distances.begin(),MaximumNumberOfResults);
                    //KRATOS_WATCH(nresults);
                    //once you have the list of points use them to build a geometry
                    Geometry< Node<3> >::Pointer pgeom = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >());
                    unsigned int counter = 0;
                    for(auto it = results.begin(); it != results.end(); ++it)
                    {
                        if(counter++ < nresults)  pgeom->push_back( *(it.base()));
                    }
                    //generating a new element
                    Element::Pointer p_element = reference_element.Create(idcount, pgeom, properties);
                    mrModelPart.AddElement(p_element);

                    //assignging variables needed to the element
                    p_element->SetValue(GAUSS_POINT_COORDINATES, coordinates_of_gauss);
                    p_element->SetValue(GAUSS_AREA, Agauss);
                    p_element->SetValue(EFFECTIVE_RADIUS, effective_radius);
                    p_element->SetValue(SEARCH_RADIUS, search_radius);


                    //set the constitutive law
                    std::vector< ConstitutiveLaw::Pointer > constitutive_law_vec(1);
                    constitutive_law_vec[0] = properties->GetValue(CONSTITUTIVE_LAW);
                    p_element->SetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vec, rCurrentProcessInfo);
                }
            }
            
            
            for(ModelPart::NodesContainerType::iterator inode = mrModelPart.NodesBegin();
                inode!=mrModelPart.NodesEnd(); inode++)
        {
			//std::cout<<" ID NODE "<<inode->Id()<<std::endl;
			//std::cout<<" mlength "<<mlength<<std::endl;
			inode->FastGetSolutionStepValue(NODAL_RADIUS) = mlength;
			
		}
            
            
            KRATOS_WATCH(mrModelPart);

//***************************************************************************************************************************



            //for(ModelPart::NodesContainerType::iterator inode = mrModelPart.NodesBegin();
                //inode!=mrModelPart.NodesEnd(); inode++)
            //{
                //unsigned int MaximumNumberOfResults = 1000;
                //std::vector<NodePointerType>                            results(MaximumNumberOfResults);
                //std::vector<double>                                     distances(MaximumNumberOfResults);
                ////double search_radius = 4.8;
                ////double effective_radius = 0.7;
                ////double search_radius = 0.0003;
                ////double effective_radius_lm = effective_radius/3;
                 ////inode->IsFixed(DISPLACEMENT_X) &&

                //if(inode->IsFixed(DISPLACEMENT_X) && inode->IsFixed(DISPLACEMENT_Y))
                //{
                    //inode->Free(DISPLACEMENT_X);
                    //inode->Free(DISPLACEMENT_Y);

                    //idcount += 1;


                    //array_1d<double,3> coordinates_of_node = ZeroVector(3);



                    //unsigned int center_id = inode->Id();


                    ////coordinates_of_node = inode->Coordinates();
                    ////Node<3> search_coordinates(0,coordinates_of_node[0],coordinates_of_node[1],coordinates_of_node[2]);
                    ////KRATOS_WATCH("--1--");

                    ////KRATOS_WATCH(coordinates_of_node);

                    //unsigned int nresults = bins.SearchInRadius(*inode,search_radius,results.begin(),distances.begin(),MaximumNumberOfResults);


                    ////once you have the list of points use them to build a geometry
                    //Geometry< Node<3> >::Pointer pcgeom = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >());

                    ////first of all add the "central" node
                    //pcgeom->push_back( *(inode.base()) );
                     ////KRATOS_WATCH((inode)->Id())

                    //unsigned int counter=0;
                    //for(NodesContainerType::iterator ir = results.begin(); ir != results.end(); ++ir)
                    //{

                        //if(counter++ < nresults)
                        //{
                            //if(center_id != (*ir)->Id() )
                            //{
                                ////KRATOS_WATCH((*ir)->Id())
                                //pcgeom->push_back( *(ir.base()));
                                ////KRATOS_WATCH(*it);
                            //}
                        //}

                    //}

                    ////KRATOS_WATCH("--2--");
                    //Condition::Pointer p_condition = Condition::Pointer(new LagrangeMultiplierCondition2D0(idcount, pcgeom, properties));
                    //mrModelPart.AddCondition(p_condition);

                    ////KRATOS_WATCH("--3--");

                    ////temp_condition_container.push_back(p_condition);
                    ////p_condition->SetValue(TEMP_POS, coordinates_of_node);
                    ////p_condition->SetValue(NODAL_AREA, Agauss);
                    //p_condition->SetValue(SEARCH_RADIUS, search_radius);
                    //p_condition->SetValue(EFFECTIVE_RADIUS, effective_radius);
                    //p_condition->SetValue(CENTER_ID, center_id);
                   //// KRATOS_WATCH(*p_condition);

                //}
			//}
//***************************************************************************************************************************
                /*if(inode->Y() >= 1.01)
                {
                    //if(inode->IsFixed(DISPLACEMENT_X) && inode->IsFixed(DISPLACEMENT_Y))
                    //{
                        //inode->Free(DISPLACEMENT_X);
                        //inode->Free(DISPLACEMENT_Y);

                        idcount += 1;


                        array_1d<double,3> coordinates_of_node = ZeroVector(3);



                        unsigned int center_id = inode->Id();


                        //coordinates_of_node = inode->Coordinates();
                        //Node<3> search_coordinates(0,coordinates_of_node[0],coordinates_of_node[1],coordinates_of_node[2]);
                        //KRATOS_WATCH("--1--");

                        //KRATOS_WATCH(coordinates_of_node);

                        unsigned int nresults = bins.SearchInRadius(*inode,search_radius,results.begin(),distances.begin(),MaximumNumberOfResults);


                        //once you have the list of points use them to build a geometry
                        Geometry< Node<3> >::Pointer pcgeom = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >());

                        //first of all add the "central" node
                        pcgeom->push_back( *(inode.base()) );
                         //KRATOS_WATCH((inode)->Id())

                        unsigned int counter=0;
                        for(NodesContainerType::iterator ir = results.begin(); ir != results.end(); ++ir)
                        {

                            if(counter++ < nresults)
                            {
                                if(center_id != (*ir)->Id() )
                                {
                                    //KRATOS_WATCH((*ir)->Id())
                                    pcgeom->push_back( *(ir.base()));
                                    //KRATOS_WATCH(*it);
                                }
                            }

                        }

                        //KRATOS_WATCH("--2--");
                        Condition::Pointer p_condition = Condition::Pointer(new LagrangeMultiplierCondition2D(idcount, pcgeom, properties));
                        mrModelPart.AddCondition(p_condition);

                        //KRATOS_WATCH("--3--");

                        p_condition->SetValue(SEARCH_RADIUS, search_radius);
                        p_condition->SetValue(EFFECTIVE_RADIUS, effective_radius);
                        p_condition->SetValue(CENTER_ID, center_id);
                       // KRATOS_WATCH(*p_condition);
                    //}


                }*/
                //if((inode->Y()) == 0)
                /*if(inode->IsFixed(DISPLACEMENT_X) && inode->IsFixed(DISPLACEMENT_Y))
                {



                    inode->Free(DISPLACEMENT_X);
                    inode->Free(DISPLACEMENT_Y);

                    idcount += 1;


                    array_1d<double,3> coordinates_of_node = ZeroVector(3);



                    unsigned int center_id = inode->Id();


                    //coordinates_of_node = inode->Coordinates();
                    //Node<3> search_coordinates(0,coordinates_of_node[0],coordinates_of_node[1],coordinates_of_node[2]);
                    KRATOS_WATCH("--1--");

                    //KRATOS_WATCH(coordinates_of_node);

                    unsigned int nresults = bins.SearchInRadius(*inode,search_radius,results.begin(),distances.begin(),MaximumNumberOfResults);


                    //once you have the list of points use them to build a geometry
                    Geometry< Node<3> >::Pointer pcgeom = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >());

                    //first of all add the "central" node
                    pcgeom->push_back( *(inode.base()) );
                     //KRATOS_WATCH((inode)->Id())

                    unsigned int counter=0;
                    for(NodesContainerType::iterator ir = results.begin(); ir != results.end(); ++ir)
                    {

                        if(counter++ < nresults)
                        {
                            if(center_id != (*ir)->Id() )
                            {
                                //KRATOS_WATCH((*ir)->Id())
                                pcgeom->push_back( *(ir.base()));
                                //KRATOS_WATCH(*it);
                            }
                        }

                    }

                    KRATOS_WATCH("--2--");
                    Condition::Pointer p_condition = Condition::Pointer(new LagrangeMultiplierCondition2D0(idcount, pcgeom, properties));
                    mrModelPart.AddCondition(p_condition);

                    KRATOS_WATCH("--3--");

                    //temp_condition_container.push_back(p_condition);
                    //p_condition->SetValue(TEMP_POS, coordinates_of_node);
                    //p_condition->SetValue(NODAL_AREA, Agauss);
                    p_condition->SetValue(EFFECTIVE_RADIUS, effective_radius);
                    p_condition->SetValue(SEARCH_RADIUS, search_radius);
                   // KRATOS_WATCH(*p_condition);

                }*/


/*

                if((inode->Y()) == 0.0)
                {
                    inode->Free(DISPLACEMENT_X);
                    inode->Free(DISPLACEMENT_Y);

                    idcount += 1;


                    array_1d<double,3> coordinates_of_node = ZeroVector(3);



                    unsigned int center_id = inode->Id();


                    //coordinates_of_node = inode->Coordinates();
                    //Node<3> search_coordinates(0,coordinates_of_node[0],coordinates_of_node[1],coordinates_of_node[2]);
                    KRATOS_WATCH("--1--");

                    //KRATOS_WATCH(coordinates_of_node);

                    unsigned int nresults = bins.SearchInRadius(*inode,search_radius,results.begin(),distances.begin(),MaximumNumberOfResults);


                    //once you have the list of points use them to build a geometry
                    Geometry< Node<3> >::Pointer pcgeom = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >());

                    //first of all add the "central" node
                    pcgeom->push_back( *(inode.base()) );
                     //KRATOS_WATCH((inode)->Id())

                    unsigned int counter=0;
                    for(NodesContainerType::iterator ir = results.begin(); ir != results.end(); ++ir)
                    {

                        if(counter++ < nresults)
                        {
                            if(center_id != (*ir)->Id() )
                            {
                                //KRATOS_WATCH((*ir)->Id())
                                pcgeom->push_back( *(ir.base()));
                                //KRATOS_WATCH(*it);
                            }
                        }

                    }

                    KRATOS_WATCH("--2--");
                    Condition::Pointer p_condition = Condition::Pointer(new LagrangeMultiplierCondition2DY(idcount, pcgeom, properties));
                    mrModelPart.AddCondition(p_condition);

                    KRATOS_WATCH("--3--");

                    //temp_condition_container.push_back(p_condition);
                    //p_condition->SetValue(TEMP_POS, coordinates_of_node);
                    //p_condition->SetValue(NODAL_AREA, Agauss);
                    p_condition->SetValue(EFFECTIVE_RADIUS, effective_radius);
                   // KRATOS_WATCH(*p_condition);

                }

                if((inode->X()) == 0.0)
                {
                    inode->Free(DISPLACEMENT_X);
                    inode->Free(DISPLACEMENT_Y);

                    idcount += 1;


                    array_1d<double,3> coordinates_of_node = ZeroVector(3);



                    unsigned int center_id = inode->Id();


                    //coordinates_of_node = inode->Coordinates();
                    //Node<3> search_coordinates(0,coordinates_of_node[0],coordinates_of_node[1],coordinates_of_node[2]);
                    KRATOS_WATCH("--1--");

                    //KRATOS_WATCH(coordinates_of_node);

                    unsigned int nresults = bins.SearchInRadius(*inode,search_radius,results.begin(),distances.begin(),MaximumNumberOfResults);


                    //once you have the list of points use them to build a geometry
                    Geometry< Node<3> >::Pointer pcgeom = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >());

                    //first of all add the "central" node
                    pcgeom->push_back( *(inode.base()) );
                     //KRATOS_WATCH((inode)->Id())

                    unsigned int counter=0;
                    for(NodesContainerType::iterator ir = results.begin(); ir != results.end(); ++ir)
                    {

                        if(counter++ < nresults)
                        {
                            if(center_id != (*ir)->Id() )
                            {
                                //KRATOS_WATCH((*ir)->Id())
                                pcgeom->push_back( *(ir.base()));
                                //KRATOS_WATCH(*it);
                            }
                        }

                    }

                    KRATOS_WATCH("--2--");
                    Condition::Pointer p_condition = Condition::Pointer(new LagrangeMultiplierCondition2DX(idcount, pcgeom, properties));
                    mrModelPart.AddCondition(p_condition);

                    KRATOS_WATCH("--3--");

                    //temp_condition_container.push_back(p_condition);
                    //p_condition->SetValue(TEMP_POS, coordinates_of_node);
                    //p_condition->SetValue(NODAL_AREA, Agauss);
                    p_condition->SetValue(EFFECTIVE_RADIUS, effective_radius);
                   // KRATOS_WATCH(*p_condition);

                }*/
                //KRATOS_WATCH(inode->FastGetSolutionStepValue(FORCE_Y));
                //if(inode->GetSolutionStepValue(FORCE_Y) != 0)
                /*if(inode->X() == 7.0)
                //if(inode->Y() > 1.0 && inode->X() == 0.5)
                {
                    idcount += 1;
                    Geometry< Node<3> >::Pointer pfgeom = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >());
                    pfgeom->push_back( *(inode.base()));
                    Condition::Pointer pf_condition = Condition::Pointer(new PointForce2D(idcount, pfgeom, properties));
                    mrModelPart.AddCondition(pf_condition);

                    //KRATOS_WATCH(*p_condition);
                    //temp_condition_container.push_back(p_condition);
                }*/

            

               //temp_condition_container.swap(mrModelPart.Conditions());

//////////////////////////////////////////////////////////////////////////////////////////////////


     }
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
    virtual std::string Info() const
    {
        return "CreateMLSParticleGauss";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "CreateMLSParticleGauss";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    ModelPart& mrModelPart;
    const std::string mElementName;
    const std::string mConditionName;
    unsigned int mDomainSize;
    unsigned int mPropertyId;




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

    /// Assignment operator.
    CreateMLSParticleGauss& operator=(CreateMLSParticleGauss const& rOther);

    /// Copy constructor.
    //CalculateNodalAreaProcess(CalculateNodalAreaProcess const& rOther);


    ///@}

}; // Class CreateMLSParticleGauss

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CreateMLSParticleGauss& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CreateMLSParticleGauss& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif //KRATOS_CREATE_MLS_PARTICLE_GAUSS_H_INCLUDED   defined


