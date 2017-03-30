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





#if !defined(COMPUTE_SHAPE_FUNCTIONS_PROCESS)
#define COMPUTE_SHAPE_FUNCTIONS_PROCESS



/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/variables.h"


#include "lagrangian_mpm_application_variables.h"
#include "utilities/math_utils.h"
#include "custom_utilities/compute_mls_shape_functions_utility.h"
#include <cmath>
#include <algorithm>



namespace Kratos
{



class ComputeMLSShapeFunctionsProcess : public Process
{

public:

    ComputeMLSShapeFunctionsProcess( ModelPart& model_part ) : Process(), mr_model_part( model_part )
    {
        /////////////

    }

    ~ComputeMLSShapeFunctionsProcess() {}



    virtual void Execute()

    {

        KRATOS_TRY;

        //const unsigned int dim = mr_model_part.GetProcessInfo()[DOMAIN_SIZE];
        const unsigned int dim = 2;

        for(ModelPart::ElementsContainerType::iterator iel = mr_model_part.ElementsBegin();
                iel!=mr_model_part.ElementsEnd(); iel++)
        {

            //calculating shape functions values

            Geometry< Node<3> >& geom = iel->GetGeometry();

            Vector& Ng = iel->GetValue(SHAPE_FUNCTIONS);

            Matrix& DN_DX = iel->GetValue(SHAPE_FUNCTIONS_DERIVATIVES);

            const array_1d<double,3>& xg = iel->GetValue(GAUSS_POINT_COORDINATES);

            const double h = iel->GetValue(EFFECTIVE_RADIUS);

            Ng = zero_vector<double> (geom.size());

            DN_DX = ZeroMatrix(geom.size(), dim);

            //KRATOS_WATCH(this->GetValue(TEMP_POS));

            Matrix node_coordinates(geom.size(),dim);

            for ( unsigned int i = 0; i < geom.size(); i++ )

            {

                for(unsigned int k=0; k<node_coordinates.size2(); k++)

                {

                    node_coordinates(i,k) = geom[i].Coordinates()[k];

                }

            }

            //Vector Ng( number_of_nodes);

            //Matrix DN_DX( number_of_nodes, 2 );

            LinearMLSKernel::ComputeMLSKernel(Ng, DN_DX, node_coordinates, xg, h);


        }


        for(ModelPart::ConditionsContainerType::iterator icon = mr_model_part.ConditionsBegin();
                icon!=mr_model_part.ConditionsEnd(); icon++)
        {


            //calculating shape functions values
            Geometry< Node<3> >& geom = icon->GetGeometry();

            Vector& Nc = icon->GetValue(SHAPE_FUNCTIONS);
            Matrix& DN_DXc = icon->GetValue(SHAPE_FUNCTIONS_DERIVATIVES);

            const double h = icon->GetValue(EFFECTIVE_RADIUS);




            Nc = zero_vector<double> (geom.size());
            DN_DXc = ZeroMatrix(geom.size(), dim);

            //KRATOS_WATCH(this->GetValue(TEMP_POS));
            Matrix node_coordinates(geom.size(),dim);
            for ( unsigned int i = 0; i < geom.size(); i++ )
            {
                for(unsigned int k=0; k<node_coordinates.size2(); k++)
                {
                    node_coordinates(i,k) = geom[i].Coordinates()[k];
                    //node_coordinates(i,k) = geom[i].GetInitialPosition()[k];
                }
            }
        //KRATOS_WATCH("bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb");
            //const array_1d<double,dim> xc = row(node_coordinates,0);
            const array_1d<double,dim> xc = icon->GetGeometry()[0].Coordinates();
        //KRATOS_WATCH("rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr");
            //Vector Ng( number_of_nodes);
            //Matrix DN_DX( number_of_nodes, 2 );
            LinearMLSKernel::ComputeMLSKernel(Nc, DN_DXc, node_coordinates, xc, h);
            //LinearULFLME::ComputeLMEShapef(Nc, DN_DXc, node_coordinates, xg, h);
            //KRATOS_WATCH(xc);

            //KRATOS_WATCH(node_coordinates);

        }

        KRATOS_CATCH("");
    }





private:

    ModelPart& mr_model_part;


};


}

#endif //COMPUTE_SHAPE_FUNCTIONS_PROCESS

