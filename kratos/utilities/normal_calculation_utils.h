//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//
#include "includes/model_part.h"

#if !defined(KRATOS_NORMAL_CALCULATION_UTILS )
#define  KRATOS_NORMAL_CALCULATION_UTILS


/* System includes */


/* External includes */


/* Project includes */
#include "utilities/math_utils.h"
#include "includes/deprecated_variables.h"


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
                if (it->GetGeometry().PointsNumber() == 2)
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
                if (it->GetGeometry().PointsNumber() == 3)
                    CalculateNormal3D(it,An,v1,v2);
            }
        }

        //adding the normals to the nodes
        for(ConditionsArrayType::iterator it =  rConditions.begin();
                it !=rConditions.end(); it++)
        {
            Geometry<Node<3> >& pGeometry = (it)->GetGeometry();
            double coeff = 1.00/pGeometry.size();
	    const array_1d<double,3>& normal = it->GetValue(NORMAL);
            for(unsigned int i = 0; i<pGeometry.size(); i++)
            {
                noalias(pGeometry[i].FastGetSolutionStepValue(NORMAL)) += coeff * normal;
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
    
    
    /**this function swaps the normal of all of the conditions in a model part
	 * This is done by swapping the two first nodes in the geometry and is thus appropriate for simplicial elements
	 * @param rModelPart ModelPart of the problem. Must have a set of conditions defining the "skin" of the domain
	 */    
    void SwapNormals(ModelPart& rModelPart)
	{
		KRATOS_TRY
        //KRATOS_WATCH(rModelPart.ConditionsBegin()->GetGeometry()[0].Id())
		for(ModelPart::ConditionsContainerType::iterator it=rModelPart.ConditionsBegin(); it!=rModelPart.ConditionsEnd(); it++)
		{
            Node<3>::Pointer paux = it->GetGeometry()(0);
            it->GetGeometry()(0) = it->GetGeometry()(1);
            it->GetGeometry()(1) = paux;
			//it->GetGeometry()(0).swap(it->GetGeometry()(1));
		}
        //KRATOS_WATCH(rModelPart.ConditionsBegin()->GetGeometry()[0].Id())
		KRATOS_CATCH("")
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

	for(ModelPart::NodesContainerType::iterator it =  rModelPart.NodesBegin();
                it !=rModelPart.NodesEnd(); it++)
        {
            noalias(it->FastGetSolutionStepValue(NORMAL)) = ZeroNormal;
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

    /// Calculates the area normal (vector oriented as the normal with a dimension proportional to the area) using only nodes marked with a flag variable and 
    /** detecting corners. Corners are defined as nodes that recieves more than 2 normals from their neighbor conditions with a difference in angle greater than Alpha . 
      *  This function is equivalent to other implementations of CalculateOnSimplex, but instead of using all conditions in the array, it only uses
      * those that contain a value of rVariable != Zero. This is useful in problems where a part of the boundary is a slip condition, as it provides
      * more reasonable values for the normals on the border between this area and other parts of the boundary. This function is safe to use in MPI.
      * @param rModelPart ModelPart of the problem. Must have a set of conditions defining the "skin" of the domain.
      * @param Dimension Spatial dimension (2 or 3).
      * @param rVariable The Kratos::Variable used to indicate which parts of the boundary will be used to calculate the normals. Conditions where rVariable == Zero will be skipped.
      * @param rAlpha the maximum angle to distinguish normals.
      */

    template< class TValueType >
    void CalculateOnSimplex(ModelPart& rModelPart,
                            int Dimension,
                            Variable<TValueType>& rVariable,
                            const TValueType Zero,const double rAlpha)
    {
        KRATOS_TRY;

        // Reset normals
        const array_1d<double,3> ZeroNormal(3,0.0);

	for(ModelPart::NodesContainerType::iterator it =  rModelPart.NodesBegin();
                it !=rModelPart.NodesEnd(); it++)
        {
            noalias(it->FastGetSolutionStepValue(NORMAL)) = ZeroNormal;
            it->FastGetSolutionStepValue(NODAL_PAUX) = 0.0;	    
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
        
        
      //loop over nodes to set normals
	for(ModelPart::NodesContainerType::iterator it =  rModelPart.NodesBegin();
                it !=rModelPart.NodesEnd(); it++)
        {
	  std::vector< array_1d<double,3> > N_Mat; 
	  N_Mat.reserve(10);
	  double nodal_area = 0.0;
	 
	  WeakPointerVector<Condition >& ng_cond = it->GetValue(NEIGHBOUR_CONDITIONS);
	  
	  if(ng_cond.size() != 0){	  
	    for(WeakPointerVector<Condition >::iterator ic = ng_cond.begin(); ic!=ng_cond.end(); ic++)
	    {
		Condition::GeometryType& pGeom = ic->GetGeometry();
		const array_1d<double,3>&  rNormal = ic->GetValue(NORMAL);
		const double Coef = 1.0 / pGeom.PointsNumber();
		double norm_normal = norm_2( rNormal );
		
		if(norm_normal != 0.0)
		{
		  nodal_area += Coef * norm_normal;
		  
		  if(N_Mat.size() == 0.0)
		      N_Mat.push_back( rNormal * Coef );
		  else{
			
			int added = 0;
			for(unsigned int ii=0; ii<N_Mat.size();++ii)
			  {
			  const array_1d<double,3>& temp_normal = N_Mat[ii];
			  double norm_temp = norm_2( temp_normal );
			  
			  double cos_alpha=temp_normal[0]*rNormal[0] + temp_normal[1]*rNormal[1] +temp_normal[2]*rNormal[2];
			  cos_alpha /= (norm_temp*norm_normal);
			  
			  if( cos_alpha > cos(0.017453293*rAlpha) ){
			    N_Mat[ii] += rNormal * Coef;
			    added = 1;}		      
			  }
			  
			if(!added)
			    N_Mat.push_back( rNormal*Coef );
			  
		    }
		  }
	      }
	  }
	  //compute NORMAL and mark 
	  array_1d<double,3> sum_Normal(3,0.0);	  
	  
	  for(unsigned int ii=0; ii<N_Mat.size(); ++ii){
	    sum_Normal += N_Mat[ii];
	  }
	  
	  noalias( it->FastGetSolutionStepValue(NORMAL) ) = sum_Normal;
	  it->FastGetSolutionStepValue(NODAL_PAUX) = nodal_area;
	  //assign IS_SLIP = 0 for vertices
	  if(N_Mat.size() == 2){
// 	    it->SetValue(IS_SLIP,0);	  
	    it->FastGetSolutionStepValue(IS_SLIP)=20.0;}	
	  else if(N_Mat.size() == 3)
	    it->FastGetSolutionStepValue(IS_SLIP)=30.0;	
	  else if(N_Mat.size() == 1)
	    it->FastGetSolutionStepValue(IS_SLIP)=10.0;	
	  
	   
	  
	}
               
        // For MPI: correct values on partition boundaries
        rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);
        rModelPart.GetCommunicator().AssembleCurrentData(NODAL_PAUX);
	
        KRATOS_CATCH("");
    }
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	template< class TValueType >
    void CalculateOnSimplexLowMemory(ModelPart& rModelPart,
                            int Dimension,
                            Variable<TValueType>& rVariable,
                            const TValueType Zero,const double rAlpha)
    {
        KRATOS_TRY;

        // Reset normals
        const array_1d<double,3> ZeroNormal(3,0.0);

	for(ModelPart::NodesContainerType::iterator it =  rModelPart.NodesBegin();
                it !=rModelPart.NodesEnd(); it++)
        {
            noalias(it->GetValue(NORMAL)) = ZeroNormal;
            it->GetValue(NODAL_PAUX) = 0.0;	    
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
        
        
      //loop over nodes to set normals
	for(ModelPart::NodesContainerType::iterator it =  rModelPart.NodesBegin();
                it !=rModelPart.NodesEnd(); it++)
        {
	  std::vector< array_1d<double,3> > N_Mat; 
	  N_Mat.reserve(10);
	  double nodal_area = 0.0;
	 
	  WeakPointerVector<Condition >& ng_cond = it->GetValue(NEIGHBOUR_CONDITIONS);
	  
	  if(ng_cond.size() != 0){	  
	    for(WeakPointerVector<Condition >::iterator ic = ng_cond.begin(); ic!=ng_cond.end(); ic++)
	    {
		Condition::GeometryType& pGeom = ic->GetGeometry();
		const array_1d<double,3>&  rNormal = ic->GetValue(NORMAL);
		const double Coef = 1.0 / pGeom.PointsNumber();
		double norm_normal = norm_2( rNormal );
		
		if(norm_normal != 0.0)
		{
		  nodal_area += Coef * norm_normal;
		  
		  if(N_Mat.size() == 0.0)
		      N_Mat.push_back( rNormal * Coef );
		  else{
			
			int added = 0;
			for(unsigned int ii=0; ii<N_Mat.size();++ii)
			  {
			  const array_1d<double,3>& temp_normal = N_Mat[ii];
			  double norm_temp = norm_2( temp_normal );
			  
			  double cos_alpha=temp_normal[0]*rNormal[0] + temp_normal[1]*rNormal[1] +temp_normal[2]*rNormal[2];
			  cos_alpha /= (norm_temp*norm_normal);
			  
			  if( cos_alpha > cos(0.017453293*rAlpha) ){
			    N_Mat[ii] += rNormal * Coef;
			    added = 1;}		      
			  }
			  
			if(!added)
			    N_Mat.push_back( rNormal*Coef );
			  
		    }
		  }
	      }
	  }
	  //compute NORMAL and mark 
	  array_1d<double,3> sum_Normal(3,0.0);	  
	  
	  for(unsigned int ii=0; ii<N_Mat.size(); ++ii){
	    sum_Normal += N_Mat[ii];
	  }
	  
	  noalias( it->FastGetSolutionStepValue(NORMAL) ) = sum_Normal;
	  it->FastGetSolutionStepValue(NODAL_PAUX) = nodal_area;
	  //assign IS_SLIP = 0 for vertices
	  if(N_Mat.size() == 2){
// 	    it->SetValue(IS_SLIP,0);	  
	    it->FastGetSolutionStepValue(IS_SLIP)=20.0;}	
	  else if(N_Mat.size() == 3)
	    it->FastGetSolutionStepValue(IS_SLIP)=30.0;	
	  else if(N_Mat.size() == 1)
	    it->FastGetSolutionStepValue(IS_SLIP)=10.0;	
	  
	   
	  
	}
               
        // For MPI: correct values on partition boundaries
        rModelPart.GetCommunicator().AssembleCurrentData(NORMAL);
        rModelPart.GetCommunicator().AssembleCurrentData(NODAL_PAUX);
	
        KRATOS_CATCH("");
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

