/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-10-31 17:51:34 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_CREATE_SPH_PARTICLE_H_INCLUDED )
#define  KRATOS_CREATE_SPH_PARTICLE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "meshless_application.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "custom_elements/SPHparticle.h"
#include "geometries/point_3d.h"

#include "includes/variables.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
typedef  ModelPart::NodesContainerType NodesContainerType;
typedef  ModelPart::ElementsContainerType ElementsContainerType;


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
class CreateSPHParticle
        : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CalculateNodalAreaProcess
    KRATOS_CLASS_POINTER_DEFINITION(CreateSPHParticle);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /// avg_elems ------ expected number of neighbour elements per node.,
    /// avg_nodes ------ expected number of neighbour Nodes
    /// the better the guess for the quantities above the less memory occupied and the fastest the algorithm
    CreateSPHParticle(ModelPart& model_part, unsigned int domain_size, float element_choice)
        : mr_model_part(model_part), mdomain_size(domain_size), element_choice(element_choice)
    {
    }

    /// Destructor.
    virtual ~CreateSPHParticle()
    {
    }


    ///@}
    ///@name Operators
    ///@{

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

        //set to zero the nodal area and nodal mass
        for (ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin();in!=mr_model_part.NodesEnd(); in++)
        {
            in->FastGetSolutionStepValue(NODAL_AREA) = 0.00;
            in->FastGetSolutionStepValue(NODAL_MASS) = 0.00;
        }





        if(mdomain_size == 2)
        {
            double area = 0.0;

            for(ModelPart::ElementsContainerType::iterator i = mr_model_part.ElementsBegin();
                i!=mr_model_part.ElementsEnd(); i++)
            {
                //calculating shape functions values
                Geometry< Node<3> >& geom = i->GetGeometry();

                area = GeometryUtils::CalculateVolume2D(geom);

                area *= 0.333333333333333333333333333;


                geom[0].FastGetSolutionStepValue(NODAL_AREA) += area;
                geom[1].FastGetSolutionStepValue(NODAL_AREA) += area;
                geom[2].FastGetSolutionStepValue(NODAL_AREA) += area;
            }
        }


        //*****************************************************//
        //STUFF TO ASSIGN NODAL AREA TO ALL BOUNDARY NODES

        double aux1=0;
        double aux2=0;

        for (ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin();in!=mr_model_part.NodesEnd(); in++)
        {
            //            if (in->FastGetSolutionStepValue(IS_STRUCTURE)== 1.0){
            aux1 += in->FastGetSolutionStepValue(NODAL_AREA);
            aux2 += 1;

            //            }
        }

        double average_nodal_area_for_Boundary = aux1/aux2;

        // CHAPUZA TO ASSIGN EQUAL NODAL AREA ( SO MASS ) TO ALL PARTICLES
        for (ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin();in!=mr_model_part.NodesEnd(); in++)
        {
            //            if (in->FastGetSolutionStepValue(IS_STRUCTURE)== 1.0){
            in->FastGetSolutionStepValue(NODAL_AREA)= average_nodal_area_for_Boundary;

            //            }
        }

        //****************************************************************//






        mr_model_part.GetCommunicator().AssembleCurrentData(NODAL_AREA);

        for (ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin();in!=mr_model_part.NodesEnd(); in++)
        {
            in->FastGetSolutionStepValue(NODAL_MASS) = (in->FastGetSolutionStepValue(NODAL_AREA))*(in->FastGetSolutionStepValue(DENSITY));
        }

        double total_mass=0;
        double inter_node_counter = 0;
        double min_x = 1000;
        double max_x = -1000;
        double min_y = 1000;
        double max_y = -1000;
        double density=0;

        for (ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin();in!=mr_model_part.NodesEnd(); in++)
        {
            //total_mass += in->FastGetSolutionStepValue(NODAL_MASS);
            density = in->FastGetSolutionStepValue(DENSITY);
            if ( (in->FastGetSolutionStepValue(IS_WET)== 1.0) && (in->FastGetSolutionStepValue(IS_STRUCTURE)== 0.0) ){
                inter_node_counter += 1;
            }
            if ( (in->FastGetSolutionStepValue(IS_WET)== 1.0)) {
                if (in->X0()>max_x){max_x=in->X0();}

                if (in->Y0()>max_y){max_y=in->Y0();}

                if (in->X0()<min_x){min_x=in->X0();}

                if (in->Y0()<min_y){min_y=in->Y0();}
            }

        }

        double total_area = (max_x - min_x) * (max_y - min_y);
        KRATOS_WATCH(total_area);
        KRATOS_WATCH(max_x);
        KRATOS_WATCH(max_y);
        KRATOS_WATCH(min_x);
        KRATOS_WATCH(min_y);
        total_mass =  total_area * density;
        double mass_for_int_node = total_mass / inter_node_counter;
        for (ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin();in!=mr_model_part.NodesEnd(); in++)
        {

            in->FastGetSolutionStepValue(NODAL_MASS) = mass_for_int_node ;

        }


        KRATOS_WATCH(total_mass);










        mr_model_part.GetCommunicator().AssembleCurrentData(NODAL_MASS);

        mr_model_part.Elements().clear();

        Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(0);

        std::string ElementName;
        if (mdomain_size == 2){

            if (element_choice==1){ElementName = std::string("SPHparticlePoly");}
            else if (element_choice==2){ElementName = std::string("SPHparticlePolyPresSpiky");}
            else if (element_choice==3){ElementName = std::string("SPHparticleQuintic");}
            else if (element_choice==4){ElementName = std::string("SPHparticleC2");}
            else if (element_choice==5){ElementName = std::string("SPHparticlePolyPresQuad");}
            else if (element_choice==6){ElementName = std::string("SPHparticleGaus");}
            else {KRATOS_WATCH("Element does not exist!!!");}

        }
        else{
            KRATOS_WATCH("DOMAIN taken as 3D")
            ElementName = std::string("SPHparticle");}



        const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);
        //        KRATOS_WATCH(&rReferenceElement);
        //        KRATOS_WATCH(&(*properties));

        //        SPHparticle * refelement;

        int idcount =0;
        for (ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin();in!=mr_model_part.NodesEnd(); in++)
        {
            idcount += 1;


            Point3D< Node<3> > geom(*(in.base()));
            //            geom[0]=in->X();
            //            geom[1]=in->Y();
            //            geom[2]=in->Z();

            Element::Pointer p_element = rReferenceElement.Create(idcount, geom, properties);
            (mr_model_part.Elements()).push_back(p_element);


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
        return "CreateSPHParticle";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "CreateSPHParticle";
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
    ModelPart& mr_model_part;
    unsigned int mdomain_size;
    float element_choice;


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
    CreateSPHParticle& operator=(CreateSPHParticle const& rOther);

    /// Copy constructor.
    //CalculateNodalAreaProcess(CalculateNodalAreaProcess const& rOther);


    ///@}

}; // Class CalculateNodalAreaProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CreateSPHParticle& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CreateSPHParticle& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif //KRATOS_CREATE_SPH_PARTICLE_H_INCLUDED   defined


