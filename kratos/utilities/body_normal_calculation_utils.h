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
//


#include "includes/model_part.h"

#if !defined(KRATOS_BODY_NORMAL_CALCULATION_UTILS )
#define  KRATOS_BODY_NORMAL_CALCULATION_UTILS


/* System includes */


/* External includes */


/* Project includes */
//#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"


namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */


/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/** Short class definition.
Detail class definition.

  \URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

	\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

	  \URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

		\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


			\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

			  \URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

				\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

				  \URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


*/
class BodyNormalCalculationUtils
{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */


    /** Destructor.
    */

    /*@} */
    /**@name Operators
    */
    /*@{ */


    /*@} */
    /**@name Operations */
    /*@{ */

    //***********************************************************************
    //***********************************************************************
    ///this function calculates the "area normal" (vector oriented as the normal
    ///with a dimension proportional to the area. This is done basing on the volume discretization.
    /**the concept is to consider the volume as filled by a fluid at constant pressure. The effect of the
     * pressure cancels on the internal nodes, but results in a net force being applied
     * on the boundary nodes. Such net force is proportional to the area of boundary
     * and oriented outwards, and we take it as definition of normal.
     * NOTE: this works only for simplex meshes, that is triangles and tetras.
     * @param ModelPart the model part with the discretization of the volume
     * @param dimension the working space dimension
     */
    void CalculateBodyNormals(
        ModelPart& r_model_part,
        int dimension)
    {
        KRATOS_TRY

        ModelPart::ElementsContainerType& rElements = r_model_part.Elements();

        //resetting the normals - only for the nodes on which we will do the calculate
        array_1d<double,3> zero = ZeroVector(3);
	
	for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin();
                it !=r_model_part.NodesEnd(); it++)
        {
            noalias(it->FastGetSolutionStepValue(NORMAL)) = zero;
        }

//         for(ModelPart::ElementsContainerType::iterator it =  rElements.begin();
//                 it !=rElements.end(); it++)
//         {
//             Element::GeometryType& rNodes = it->GetGeometry();
//             for(unsigned int in = 0; in<rNodes.size(); in++)
//                 noalias((rNodes[in]).FastGetSolutionStepValue(NORMAL)) = zero;
//         }


        //calculating the normals and storing on the conditions
        array_1d<double,3> An;
        if(dimension == 2)
        {
            BoundedMatrix<double,3,2> DN_DX;
            array_1d<double,3> N;
            double Volume;
            for(ModelPart::ElementsContainerType::iterator it =  rElements.begin(); it !=rElements.end(); it++)
            {
                Element::GeometryType& geom = it->GetGeometry();
                GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

                for(unsigned int i = 0; i<geom.size(); i++)
                {
                    array_1d<double,3>& normal = geom[i].FastGetSolutionStepValue(NORMAL);
                    for(unsigned int j=0; j<2; j++)
                    {
                        normal[j] += Volume*DN_DX(i,j);
                    }
                }
            }
        }
        else if(dimension == 3)
        {
            BoundedMatrix<double,4,3> DN_DX;
            array_1d<double,4> N;
            double Volume;
            for(ModelPart::ElementsContainerType::iterator it =  rElements.begin(); it !=rElements.end(); it++)
            {
                Element::GeometryType& geom = it->GetGeometry();
                GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

                for(unsigned int i = 0; i<geom.size(); i++)
                {
                    array_1d<double,3>& normal = geom[i].FastGetSolutionStepValue(NORMAL);
                    for(unsigned int j=0; j<3; j++)
                    {
                        normal[j] += Volume*DN_DX(i,j);
                    }
                }
            }
        }

        r_model_part.GetCommunicator().AssembleCurrentData(NORMAL);
// 			r_model_part.GetCommunicator().AssembleCurrentData(PRESSURE);


        KRATOS_CATCH("")

    }



    /*@} */
    /**@name Acces */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */
    //this function adds the Contribution of one of the geometries
    //to the corresponding nodes

    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Acces */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    //BodyNormalCalculationUtils(void);

    //BodyNormalCalculationUtils(BodyNormalCalculationUtils& rSource);


    /*@} */

}; /* Class ClassName */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_BODY_NORMAL_CALCULATION_UTILS defined */

