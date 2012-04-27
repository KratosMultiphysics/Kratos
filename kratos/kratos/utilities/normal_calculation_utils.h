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


/* *********************************************************
*
*   Last Modified by:    $Author: pooyan $
*   Date:                $Date: 2008-11-13 12:12:17 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/
#include "includes/model_part.h"

#if !defined(KRATOS_NORMAL_CALCULATION_UTILS )
#define  KRATOS_NORMAL_CALCULATION_UTILS


/* System includes */


/* External includes */


/* Project includes */
#include "utilities/math_utils.h"


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

/// Tool to evaluate the normals on nodes based on the normals of a set of surface conditions
class NormalCalculationUtils
{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
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

    /// Calculates the "area normal" (vector oriented as the normal with a dimension proportional to the area).
    /** This is done on the base of the Conditions provided which should be
      * understood as the surface elements of the area of interest.
      * @param rConditions A set of conditions defining the "skin" of a model
      * @param dimension Spatial dimension (2 or 3)
      * @note This function is not recommended for distributed (MPI) runs, as
      * the user has to ensure that the calculated normals are assembled between
      * processes. The overload of this function that takes a ModelPart is
      * preferable in ths case, as it performs the required communication.
      */
    void CalculateOnSimplex(ConditionsArrayType& rConditions,
                            int dimension)
    {
        KRATOS_TRY

        //resetting the normals
        array_1d<double,3> zero = Vector(3);
        noalias(zero) = ZeroVector(3);

        for(ConditionsArrayType::iterator it =  rConditions.begin();
                it !=rConditions.end(); it++)
        {
            Element::GeometryType& rNodes = it->GetGeometry();
            for(unsigned int in = 0; in<rNodes.size(); in++)
                noalias((rNodes[in]).GetSolutionStepValue(NORMAL)) = zero;
        }


        //calculating the normals and storing on the conditions
        array_1d<double,3> An;
        if(dimension == 2)
        {
            for(ConditionsArrayType::iterator it =  rConditions.begin();
                    it !=rConditions.end(); it++)
            {
                CalculateNormal2D(it,An);
            }
        }
        else if(dimension == 3)
        {
            array_1d<double,3> v1;
            array_1d<double,3> v2;
            for(ConditionsArrayType::iterator it =  rConditions.begin();
                    it !=rConditions.end(); it++)
            {
                //calculate the normal on the given condition
                CalculateNormal3D(it,An,v1,v2);
            }
        }

        //adding the normals to the nodes
        for(ConditionsArrayType::iterator it =  rConditions.begin();
                it !=rConditions.end(); it++)
        {
            Geometry<Node<3> >& pGeometry = (it)->GetGeometry();
            double coeff = 1.00/pGeometry.size();
            const array_1d<double,3>& An = it->GetValue(NORMAL);

            for(unsigned int i = 0; i<pGeometry.size(); i++)
            {
                noalias(pGeometry[i].FastGetSolutionStepValue(NORMAL)) += coeff * An;
            }
        }


        KRATOS_CATCH("")

    }

    /// Calculates the area normal (vector oriented as the normal with a dimension proportional to the area).
    /** This is done on the base of the Conditions provided which should be
      * understood as the surface elements of the area of interest.
      * @param rModelPart ModelPart of the problem. Must have a set of conditions defining the "skin" of the domain
      * @param dimension Spatial dimension (2 or 3)
      * @note Use this fuction instead of its overload taking a Conditions array for MPI applications,
      * as it will take care of communication between partitions.
      */
    void CalculateOnSimplex(ModelPart& rModelPart,
                            int Dimension)
    {
        this->CalculateOnSimplex(rModelPart.Conditions(),Dimension);
        rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);
    }

    /// Calculates the area normal (vector oriented as the normal with a dimension proportional to the area) using only nodes marked with a flag variable.
    /** This function is equivalent to other implementations of CalculateOnSimplex, but instead of using all conditions in the array, it only uses
      * those that contain a value of rVariable != Zero. This is useful in problems where a part of the boundary is a slip condition, as it provides
      * more reasonable values for the normals on the border between this area and other parts of the boundary. This function is safe to use in MPI.
      * @param rModelPart ModelPart of the problem. Must have a set of conditions defining the "skin" of the domain.
      * @param Dimension Spatial dimension (2 or 3).
      * @param rVariable The Kratos::Variable used to indicate which parts of the boundary will be used to calculate the normals.
      * @param Zero The 'off' value for the flag. Conditions where rVariable == Zero will be skipped for normal calculation.
      */
    template< class TValueType >
    void CalculateOnSimplex(ModelPart& rModelPart,
                            int Dimension,
                            Variable<TValueType>& rVariable,
                            const TValueType Zero)
    {
        KRATOS_TRY;

        // Reset normals
        const array_1d<double,3> ZeroNormal(3,0.0);

        for ( ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond )
        {
            Condition::GeometryType& rGeom = itCond->GetGeometry();
            for ( Condition::GeometryType::iterator itNode = rGeom.begin(); itNode != rGeom.end(); ++itNode)
                itNode->GetValue(NORMAL) = ZeroNormal;
        }

        // Calculate new condition normals, using only conditions with rVariable == rValue
        array_1d<double,3> An(3,0.0);

        if ( Dimension == 2 )
        {
            for ( ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond )
            {
                if ( itCond->GetValue(rVariable) != Zero )
                    CalculateNormal2D(itCond,An);
            }
        }
        else if ( Dimension == 3 )
        {
            array_1d<double,3> v1(3,0.0);
            array_1d<double,3> v2(3,0.0);

            for ( ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond )
            {
                if ( itCond->GetValue(rVariable) != Zero )
                    CalculateNormal3D(itCond,An,v1,v2);
            }
        }

        // Transfer normals to nodes
        for ( ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); ++itCond )
        {
            Condition::GeometryType& rGeom = itCond->GetGeometry();
            const double Coef = 1.0 / rGeom.PointsNumber();
            const array_1d<double,3>& rNormal = itCond->GetValue(NORMAL);
            for ( Condition::GeometryType::iterator itNode = rGeom.begin(); itNode != rGeom.end(); ++itNode)
                noalias( itNode->FastGetSolutionStepValue(NORMAL) ) += rNormal * Coef;
        }

        // For MPI: correct values on partition boundaries
        rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);

        KRATOS_CATCH("");
    }

    /// Calculates the area normal (vector oriented as the normal with a dimension proportional to the area) using only nodes marked with a flag variable.
    /** This function is equivalent to other implementations of CalculateOnSimplex, but instead of using all conditions in the array, it only uses
      * those that contain a value of rVariable != Zero. This is useful in problems where a part of the boundary is a slip condition, as it provides
      * more reasonable values for the normals on the border between this area and other parts of the boundary. This function is safe to use in MPI.
      * @param rModelPart ModelPart of the problem. Must have a set of conditions defining the "skin" of the domain.
      * @param Dimension Spatial dimension (2 or 3).
      * @param rVariable The Kratos::Variable used to indicate which parts of the boundary will be used to calculate the normals. Conditions where rVariable == Zero will be skipped.
      */
    template< class TValueType >
    void CalculateOnSimplex(ModelPart& rModelPart,
                            int Dimension,
                            Variable<TValueType>& rVariable)
    {
        CalculateOnSimplex(rModelPart,Dimension,rVariable,TValueType());
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
    static void CalculateNormal2D(ConditionsArrayType::iterator it, array_1d<double,3>& An)
    {
        Geometry<Node<3> >& pGeometry = (it)->GetGeometry();

        An[0] =    pGeometry[1].Y() - pGeometry[0].Y();
        An[1] = - (pGeometry[1].X() - pGeometry[0].X());
        An[2] =    0.00;

        array_1d<double,3>& normal = (it)->GetValue(NORMAL);
        noalias(normal) = An;

        // 				(it)->SetValue(NORMAL,An);
    }

    static void CalculateNormal3D(ConditionsArrayType::iterator it, array_1d<double,3>& An,
                                  array_1d<double,3>& v1,array_1d<double,3>& v2 )
    {
        Geometry<Node<3> >& pGeometry = (it)->GetGeometry();

        v1[0] = pGeometry[1].X() - pGeometry[0].X();
        v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
        v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

        v2[0] = pGeometry[2].X() - pGeometry[0].X();
        v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
        v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

        MathUtils<double>::CrossProduct(An,v1,v2);
        An *= 0.5;

        array_1d<double,3>& normal = (it)->GetValue(NORMAL);
        noalias(normal) = An;
        // 				noalias((it)->GetValue(NORMAL)) = An;
    }
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

    //NormalCalculationUtils(void);

    //NormalCalculationUtils(NormalCalculationUtils& rSource);


    /*@} */

}; /* Class ClassName */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_NORMAL_CALCULATION_UTILS  defined */

