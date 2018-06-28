// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Riccardo Rossi
//

#if !defined(KRATOS_FACE_HEAT_UTILITIES_INCLUDED )
#define  KRATOS_FACE_HEAT_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "convection_diffusion_application.h"

namespace Kratos
{
class FaceHeatUtilities
{
public:

    //**********************************************************************************************
    //**********************************************************************************************
    void ApplyFaceHeat(
        ModelPart::ConditionsContainerType& conditions ,
        double face_heat_source
    )
    {
        KRATOS_TRY
        for(ModelPart::ConditionsContainerType::iterator iii = conditions.begin(); iii != conditions.end(); iii++)
        {
            Geometry< Node<3> >& geom = iii->GetGeometry();

            for(unsigned int k = 0; k<geom.size(); k++)
                geom[k].FastGetSolutionStepValue(FACE_HEAT_FLUX) = face_heat_source;

        }
        std::cout << "Conditions are generated" << std::endl;

        KRATOS_CATCH("");
    }

    void GenerateModelPart(	ModelPart& OriginModelPart, ModelPart& DestinationModelPart, Element const& rReferenceElement, Condition const& rReferenceBoundaryCondition)
    {
        KRATOS_TRY;

        //assigning the nodes to the new model part
        DestinationModelPart.rProperties() = OriginModelPart.rProperties();
        DestinationModelPart.Nodes().clear();
        DestinationModelPart.Nodes() = OriginModelPart.Nodes();

        //generating the elements
        int id = 1;
        Properties::Pointer properties = OriginModelPart.GetMesh().pGetProperties(1);
        for(ModelPart::ElementsContainerType::iterator iii = OriginModelPart.ElementsBegin(); iii != OriginModelPart.ElementsEnd(); iii++)
        {
            Geometry< Node<3> >& geom = iii->GetGeometry();
            Properties::Pointer properties = iii->pGetProperties();
            Element::Pointer p_element = rReferenceElement.Create(id, geom ,properties);
            DestinationModelPart.Elements().push_back(p_element);
            id = id + 1;
        }
        std::cout << "Elements are generated" << std::endl;

        //generating the conditions
        id = 1;
        for(ModelPart::ConditionsContainerType::iterator iii = OriginModelPart.ConditionsBegin(); iii != OriginModelPart.ConditionsEnd(); iii++)
        {
            Geometry< Node<3> >& geom = iii->GetGeometry();
            double nfree_surf = 0;
            for(unsigned int k = 0; k<geom.size(); k++)
                nfree_surf += geom[k].FastGetSolutionStepValue(IS_FREE_SURFACE);

            if(nfree_surf > 1)
            {
                Properties::Pointer properties = iii->pGetProperties();
                Condition::Pointer p_condition = rReferenceBoundaryCondition.Create(id, geom,properties);
                DestinationModelPart.Conditions().push_back(p_condition);
                id = id + 1;
            }
        }
        std::cout << "Conditions are generated" << std::endl;

        KRATOS_CATCH("");
    }





    void ConditionModelPart(ModelPart& temperature_model_part, ModelPart& full_model_part, const int TDim)
    {
        KRATOS_TRY;


        temperature_model_part.Conditions().clear();
        int nd=TDim;
        Properties::Pointer properties = full_model_part.GetMesh().pGetProperties(1);


        //int id=full_model_part.Conditions().size();
        //int n_int=0.0;


        for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; im != full_model_part.ElementsEnd() ; ++im)
        {
            if(nd==2)
            {

                int n_int=im->GetGeometry()[0].FastGetSolutionStepValue(MATERIAL_VARIABLE);
                for(int j=1; j<nd+1; j++) n_int+= im->GetGeometry()[j].FastGetSolutionStepValue(MATERIAL_VARIABLE);

                if (n_int==3)
                {
                }
                else
                {


                    n_int=im->GetGeometry()[1].FastGetSolutionStepValue(MATERIAL_VARIABLE) + im->GetGeometry()[2].FastGetSolutionStepValue(MATERIAL_VARIABLE);

                    if(n_int==2)
                    {

                        Condition::NodesArrayType temp1;
                        temp1.reserve(2);
                        temp1.push_back(im->GetGeometry()(1));
                        temp1.push_back(im->GetGeometry()(2));
                        Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp1) );
                        int id = (im->Id()-1)*3;
                        Condition::Pointer p_cond(new ThermalFace2D(id, cond, properties));
                        temperature_model_part.Conditions().push_back(p_cond);

                    }

                    n_int= im->GetGeometry()[2].FastGetSolutionStepValue(MATERIAL_VARIABLE) + im->GetGeometry()[0].FastGetSolutionStepValue(MATERIAL_VARIABLE);


                    if(n_int==2)
                    {

                        Condition::NodesArrayType temp1;
                        temp1.reserve(2);
                        temp1.push_back(im->GetGeometry()(2));
                        temp1.push_back(im->GetGeometry()(0));
                        Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp1) );
                        int id = (im->Id()-1)*3+1;
                        Condition::Pointer p_cond(new ThermalFace2D(id, cond, properties));
                        temperature_model_part.Conditions().push_back(p_cond);
                    }


                    n_int= im->GetGeometry()[0].FastGetSolutionStepValue(MATERIAL_VARIABLE) + im->GetGeometry()[1].FastGetSolutionStepValue(MATERIAL_VARIABLE) ;

                    if(n_int==2)
                    {

                        Condition::NodesArrayType temp1;
                        temp1.reserve(2);
                        temp1.push_back(im->GetGeometry()(0));
                        temp1.push_back(im->GetGeometry()(1));
                        Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp1) );
                        int id = (im->Id()-1)*3+2;

                        Condition::Pointer p_cond(new ThermalFace2D(id, cond, properties));
                        temperature_model_part.Conditions().push_back(p_cond);

                    }

                }

            }

            else
            {

                int n_int=im->GetGeometry()[0].FastGetSolutionStepValue(MATERIAL_VARIABLE);
                for(int j=1; j<nd+1; j++) n_int+= im->GetGeometry()[j].FastGetSolutionStepValue(MATERIAL_VARIABLE);

                if (n_int==4)
                {

                }
                else
                {
                    n_int=0.0;
                    n_int=im->GetGeometry()[1].FastGetSolutionStepValue(MATERIAL_VARIABLE);
                    n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(MATERIAL_VARIABLE);
                    n_int+=im->GetGeometry()[3].FastGetSolutionStepValue(MATERIAL_VARIABLE);
                    if(n_int==3)
                    {

                        Condition::NodesArrayType temp;
                        temp.reserve(3);
                        temp.push_back(im->GetGeometry()(1));
                        temp.push_back(im->GetGeometry()(2));
                        temp.push_back(im->GetGeometry()(3));
                        Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Triangle3D3< Node<3> >(temp) );
                        int id = (im->Id()-1)*4;
                        Condition::Pointer p_cond = Condition::Pointer(new ThermalFace3D(id, cond, properties) );
                        temperature_model_part.Conditions().push_back(p_cond);
                    }
                    n_int=0.0;
                    n_int=im->GetGeometry()[0].FastGetSolutionStepValue(MATERIAL_VARIABLE);
                    n_int+=im->GetGeometry()[3].FastGetSolutionStepValue(MATERIAL_VARIABLE);
                    n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(MATERIAL_VARIABLE);
                    if(n_int==3)
                    {

                        Condition::NodesArrayType temp;
                        temp.reserve(3);
                        temp.push_back(im->GetGeometry()(0));
                        temp.push_back(im->GetGeometry()(3));
                        temp.push_back(im->GetGeometry()(2));
                        Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Triangle3D3< Node<3> >(temp) );
                        int id = (im->Id()-1)*4;
                        Condition::Pointer p_cond = Condition::Pointer(new ThermalFace3D(id, cond, properties) );
                        temperature_model_part.Conditions().push_back(p_cond);
                    }
                    n_int=0.0;
                    n_int=im->GetGeometry()[0].FastGetSolutionStepValue(MATERIAL_VARIABLE);
                    n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(MATERIAL_VARIABLE);
                    n_int+=im->GetGeometry()[3].FastGetSolutionStepValue(MATERIAL_VARIABLE);
                    if(n_int==3)
                    {

                        Condition::NodesArrayType temp;
                        temp.reserve(3);
                        temp.push_back(im->GetGeometry()(0));
                        temp.push_back(im->GetGeometry()(1));
                        temp.push_back(im->GetGeometry()(3));
                        Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Triangle3D3< Node<3> >(temp) );
                        int id = (im->Id()-1)*4;
                        Condition::Pointer p_cond = Condition::Pointer(new ThermalFace3D(id, cond, properties) );
                        temperature_model_part.Conditions().push_back(p_cond);
                    }


                    n_int=0.0;
                    n_int=im->GetGeometry()[0].FastGetSolutionStepValue(MATERIAL_VARIABLE);
                    n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(MATERIAL_VARIABLE);
                    n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(MATERIAL_VARIABLE);
                    if(n_int==3)
                    {

                        Condition::NodesArrayType temp;
                        temp.reserve(3);
                        temp.push_back(im->GetGeometry()(0));
                        temp.push_back(im->GetGeometry()(2));
                        temp.push_back(im->GetGeometry()(1));
                        Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Triangle3D3< Node<3> >(temp) );
                        int id = (im->Id()-1)*4;
                        Condition::Pointer p_cond = Condition::Pointer(new ThermalFace3D(id, cond, properties) );
                        temperature_model_part.Conditions().push_back(p_cond);

                    }

                }

            }
        }


        KRATOS_CATCH("");
    }


private:

};

}  // namespace Kratos.

#endif // KRATOS_FACE_HEAT_UTILITIES_INCLUDED  defined 


