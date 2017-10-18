//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti
//

#if !defined(KRATOS_RAD_FACE_UTILITIES_INCLUDED )
#define  KRATOS_RAD_FACE_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "radiation_application.h"




namespace Kratos
{
class RadFaceUtilities
{
public:

 
    void ConditionModelPart(ModelPart& temperature_model_part, ModelPart& full_model_part, const int TDim)
    {
        KRATOS_TRY;


        temperature_model_part.Conditions().clear();
        int nd=TDim;
        Properties::Pointer properties = full_model_part.GetMesh().pGetProperties(1);



        for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; im != full_model_part.ElementsEnd() ; ++im)
        {
	   if(nd==2)
            {

                int n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY);
                for(int j=1; j<nd+1; j++) n_int+= im->GetGeometry()[j].FastGetSolutionStepValue(IS_BOUNDARY);
				/*
                if (n_int==3)
                {
                }
                else
                {


                    n_int=im->GetGeometry()[1].FastGetSolutionStepValue(IS_BOUNDARY) + im->GetGeometry()[2].FastGetSolutionStepValue(IS_BOUNDARY);

                    if(n_int==2)
                    {

                        Condition::NodesArrayType temp1;
                        temp1.reserve(2);
                        temp1.push_back(im->GetGeometry()(1));
                        temp1.push_back(im->GetGeometry()(2));
                        Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp1) );
                        int id = (im->Id()-1)*3;
                        Condition::Pointer p_cond(new RadFace2D(id, cond, properties));
                        temperature_model_part.Conditions().push_back(p_cond);

                    }

                    n_int= im->GetGeometry()[2].FastGetSolutionStepValue(IS_BOUNDARY) + im->GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY);


                    if(n_int==2)
                    {

                        Condition::NodesArrayType temp1;
                        temp1.reserve(2);
                        temp1.push_back(im->GetGeometry()(2));
                        temp1.push_back(im->GetGeometry()(0));
                        Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp1) );
                        int id = (im->Id()-1)*3+1;
                        Condition::Pointer p_cond(new RadFace2D(id, cond, properties));
                        temperature_model_part.Conditions().push_back(p_cond);
                    }


                    n_int= im->GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY) + im->GetGeometry()[1].FastGetSolutionStepValue(IS_BOUNDARY) ;

                    if(n_int==2)
                    {

                        Condition::NodesArrayType temp1;
                        temp1.reserve(2);
                        temp1.push_back(im->GetGeometry()(0));
                        temp1.push_back(im->GetGeometry()(1));
                        Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp1) );
                        int id = (im->Id()-1)*3+2;

                        Condition::Pointer p_cond(new RadFace2D(id, cond, properties));
                        temperature_model_part.Conditions().push_back(p_cond);

                    }

                }*/

            }

            else
            {


		KRATOS_WATCH(nd);
                int n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY);


                for(int j=1; j<nd+1; j++) n_int+= im->GetGeometry()[j].FastGetSolutionStepValue(IS_BOUNDARY);

                if (n_int==4)
                {

                }
                else
                {
                    n_int=0.0;
                    n_int=im->GetGeometry()[1].FastGetSolutionStepValue(IS_BOUNDARY);
                    n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(IS_BOUNDARY);
                    n_int+=im->GetGeometry()[3].FastGetSolutionStepValue(IS_BOUNDARY);
                    if(n_int==3)
                    {

                        Condition::NodesArrayType temp;
                        temp.reserve(3);
                        temp.push_back(im->GetGeometry()(1));
                        temp.push_back(im->GetGeometry()(2));
                        temp.push_back(im->GetGeometry()(3));
                        Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Triangle3D3< Node<3> >(temp) );
                        int id = (im->Id()-1)*4;
                        Condition::Pointer p_cond = Condition::Pointer(new RadFace3D(id, cond, properties) );
                        temperature_model_part.Conditions().push_back(p_cond);
                    }
                    n_int=0.0;
                    n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY);
                    n_int+=im->GetGeometry()[3].FastGetSolutionStepValue(IS_BOUNDARY);
                    n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(IS_BOUNDARY);
                    if(n_int==3)
                    {

                        Condition::NodesArrayType temp;
                        temp.reserve(3);
                        temp.push_back(im->GetGeometry()(0));
                        temp.push_back(im->GetGeometry()(3));
                        temp.push_back(im->GetGeometry()(2));
                        Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Triangle3D3< Node<3> >(temp) );
                        int id = (im->Id()-1)*4;
                        Condition::Pointer p_cond = Condition::Pointer(new RadFace3D(id, cond, properties) );
                        temperature_model_part.Conditions().push_back(p_cond);
                    }
                    n_int=0.0;
                    n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY);
                    n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(IS_BOUNDARY);
                    n_int+=im->GetGeometry()[3].FastGetSolutionStepValue(IS_BOUNDARY);
                    if(n_int==3)
                    {

                        Condition::NodesArrayType temp;
                        temp.reserve(3);
                        temp.push_back(im->GetGeometry()(0));
                        temp.push_back(im->GetGeometry()(1));
                        temp.push_back(im->GetGeometry()(3));
                        Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Triangle3D3< Node<3> >(temp) );
                        int id = (im->Id()-1)*4;
                        Condition::Pointer p_cond = Condition::Pointer(new RadFace3D(id, cond, properties) );
                        temperature_model_part.Conditions().push_back(p_cond);
                    }


                    n_int=0.0;
                    n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY);
                    n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(IS_BOUNDARY);
                    n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(IS_BOUNDARY);
                    if(n_int==3)
                    {

                        Condition::NodesArrayType temp;
                        temp.reserve(3);
                        temp.push_back(im->GetGeometry()(0));
                        temp.push_back(im->GetGeometry()(2));
                        temp.push_back(im->GetGeometry()(1));
                        Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Triangle3D3< Node<3> >(temp) );
                        int id = (im->Id()-1)*4;
                        Condition::Pointer p_cond = Condition::Pointer(new RadFace3D(id, cond, properties) );
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


