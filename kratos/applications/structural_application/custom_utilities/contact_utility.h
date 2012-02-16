/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
*   Last Modified by:    $Author: nagel $
*   Date:                $Date: 2009-03-25 08:14:58 $
*   Revision:            $Revision: 1.12 $
*
* ***********************************************************/

#if !defined(KRATOS_CONTACT_UTILITY_H_INCLUDED )
#define  KRATOS_CONTACT_UTILITY_H_INCLUDED

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"
// #include "solving_strategies/strategies/solving_strategy.h"
// #include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
// #include "solving_strategies/convergencecriterias/convergence_criteria.h"

#include "custom_conditions/contact_link_3D.h"
#include "custom_conditions/contact_link_3D_kinematic_linear.h"
#include "structural_application.h"

//default builder and solver
// #include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

namespace Kratos
{
    /**
     * auxiliary functions for contact problems
     */
    class ContactUtility
    {
        public:
            
            /**
             * Type Definitions 
             */
//             typedef ConvergenceCriteria<TSparseSpace,TDenseSpace> TConvergenceCriteriaType;
            
            /** 
             * Counted pointer of ContactUtility 
             */
            KRATOS_CLASS_POINTER_DEFINITION( ContactUtility );
            
            typedef std::size_t IndexType;
            
            typedef ModelPart::ConditionsContainerType ConditionsArrayType;
            
            typedef Geometry<Node<3> > GeometryType;
            
            typedef Properties PropertiesType;

            /**
             * Life Cycle 
             */
            
            /**
             * Constructor.
             */
            ContactUtility( int echo_level )
            {
                mEchoLevel = echo_level;
            }
            
            /**
             * Destructor.
             */
            virtual ~ContactUtility() {}
            
            /**
             * Operations
             */
            
            /**
             * Setting up the contact conditions.
             * This function performs the contact search and creates
             * temporary, virtual linking conditions for each integration point
             * on the slave surfaces.
             * @param mr_model_part the model part
             * @param initial_penalty the initial value for the penalty in normal direction
             * @param initial_penalty_t the initial value for the penalty in 
             *   transverse direction (for friction)
             * @param contact_double_check whether the problem should be solved using the master and the 
             *      slave surface as master and slave and additionally vice versa
             * @return the number of the last condition that is not a contact link
             */
            int SetUpContactConditions( ModelPart& mr_model_part, 
                                        double initial_penalty,
                                        double initial_penalty_t,
                                        bool contact_double_check
                                      )
            {
                k_contact = INT_MAX;
                k_contact_t = INT_MAX;
                //getting the array of the conditions
                ConditionsArrayType& ConditionsArray = mr_model_part.Conditions();
                //setting up candidate master conditions
                ConditionsArrayType MasterConditionsArray;
                //setting up an array of linking conditions
                ConditionsArrayType LinkingConditions;
                
                int lastRealCondition = mr_model_part.Conditions().size();

                if(lastRealCondition != 0 )
                {
                
                    for( ConditionsArrayType::ptr_iterator it = 
                         ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end();
                         ++it )
                    {
                        if( (*it)->GetValue( ACTIVATION_LEVEL ) == 0 )
                        {
                                if( (*it)->GetValue( IS_CONTACT_MASTER ) )
                                {
                                        MasterConditionsArray.push_back( *it );
                                }
                        }
                    }

                    GeometryType::Pointer tempGeometry =  GeometryType::Pointer( new Geometry<Node<3> >() );
                    
                    int properties_index = mr_model_part.NumberOfProperties();
                    PropertiesType::Pointer tempProperties( new PropertiesType(properties_index+1) );
                    mr_model_part.AddProperties( tempProperties );
                    
                    ConditionsArrayType::ptr_iterator it_begin=ConditionsArray.ptr_begin();
                    ConditionsArrayType::ptr_iterator it_end=ConditionsArray.ptr_end();

                    for (ConditionsArrayType::ptr_iterator it = it_begin; it!=it_end; ++it)
                    {
                        if( (*it)->GetValue( ACTIVATION_LEVEL ) == 0 )
                        {
                                if( (*it)->GetValue( IS_CONTACT_SLAVE ) )
                                {
                                        for( IndexType i = 0; i < (*it)->GetGeometry().IntegrationPoints().size(); i++ )
                                        {
                                                Point<3> MasterContactLocalPoint;
                                                Point<3> SlaveContactLocalPoint;
                                                (*it)->GetValue(PENALTY)[i]= initial_penalty;
                                                (*it)->GetValue(PENALTY_T)[i]= initial_penalty_t;
                                                Condition::Pointer CurrentMaster ;

                                                if( SearchPartner( 
                                                        mr_model_part,
                                                        (**it),
                                                        MasterConditionsArray, i,
                                                        MasterContactLocalPoint,
                                                        SlaveContactLocalPoint,
                                                        CurrentMaster ))
                                                {
                                                        IndexType newId = (mr_model_part.Conditions().end()-1)->Id()+LinkingConditions.size()+1;
                                                        //creating contact link element
                                                        Condition::Pointer newLink = Condition::Pointer( new ContactLink3D_Kinematic_Linear(newId,
                                                        tempGeometry,
                                                        tempProperties,
                                                        CurrentMaster, 
                                                        *it,
                                                        MasterContactLocalPoint,
                                                        SlaveContactLocalPoint, i) );

                                                        LinkingConditions.push_back( newLink );
                                                }
                                                else
                                                {
//                                                     std::cout << "no contact partner found for condition: " << (*it)->Id() <<", integration point " << i << std::endl;
                                                }
                                        }
                                }
                        }
                    }
                    
                    if( contact_double_check )
                    {
                        if( this->mEchoLevel > 1 )
                        {
                            std::cout << "double-check of contact problem activated: additional contact links are set up..." << std::endl;
                        }
                        MasterConditionsArray.clear();
                        for( ConditionsArrayType::ptr_iterator it = 
                             ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end();
                             ++it )
                        {
                            if( (*it)->GetValue( ACTIVATION_LEVEL ) == 0 )
                            {
                                if( (*it)->GetValue( IS_CONTACT_SLAVE ) )
                                {
                                    MasterConditionsArray.push_back( *it );
                                }
                            }
                        }
                        
                        GeometryType::Pointer tempGeometry =  GeometryType::Pointer( 
                                new Geometry<Node<3> >() );
                        
                        PropertiesType::Pointer tempProperties = ConditionsArray.begin()->pGetProperties();
                        
                        ConditionsArrayType::ptr_iterator it_begin=ConditionsArray.ptr_begin();
                        ConditionsArrayType::ptr_iterator it_end=ConditionsArray.ptr_end();
                        
                        for (ConditionsArrayType::ptr_iterator it = it_begin; it!=it_end; ++it)
                        {
                            if( (*it)->GetValue( ACTIVATION_LEVEL ) == 0 )
                            {
                                if( (*it)->GetValue( IS_CONTACT_MASTER ) )
                                {
                                    for( IndexType i = 0; 
                                         i < (*it)->GetGeometry().IntegrationPoints().size();
                                         i++ )
                                    {
                                        Point<3> MasterContactLocalPoint;
                                        Point<3> SlaveContactLocalPoint;
                                        Condition::Pointer CurrentMaster;
                                        if( SearchPartner(  
                                            mr_model_part,
                                            (**it),
                                               MasterConditionsArray, i,
                                               MasterContactLocalPoint,
                                               SlaveContactLocalPoint,
                                               CurrentMaster ))
                                        {
                                            IndexType newId = (
                                                    mr_model_part.Conditions().end()-1)->Id()
                                                    +LinkingConditions.size()+1;
                                            Condition::Pointer newLink = Condition::Pointer( 
                                                    new ContactLink3D_Kinematic_Linear( newId,
                                                    tempGeometry, tempProperties,
                                                            CurrentMaster, *it,
                                                                    MasterContactLocalPoint,
                                                                            SlaveContactLocalPoint, 
                                                                                i ) );
                                            LinkingConditions.push_back( newLink );
                                        }
                                        else
                                            std::cout<<"No Partner exists"<<std::endl;
                                    }
                                }
                            }
                        }
                    }
                    
                    //adding linking to model_part
                    KRATOS_WATCH(LinkingConditions.size());
                    if( LinkingConditions.size() == 0 )
                    {
                        std::cout << "!!!!!! NO LINKING CONDITIONS FOUND !!!!!!" << std::endl;
                    }
                    for(  ConditionsArrayType::ptr_iterator it=LinkingConditions.ptr_begin();
                          it != LinkingConditions.ptr_end(); ++it )
                    {
                        mr_model_part.Conditions().push_back( *it );
                    }
                    LinkingConditions.clear();
                }
                return lastRealCondition;
            }//SetUpContactConditions
            
//**********************************************************************
//**********************************************************************
            
            /**
             * This function updates the lagrangian multipliers. During the Uzawa loop
             * @param mr_model_part the model part
             * @param lastRealCondition number of the last condition that is not a contact linking
             * @param ... some parameters for an adaptation of the penalty values during the Uzawa loop
             */
            void Update( ModelPart& mr_model_part, int lastRealCondition, 
                         double friction_coefficient,
                         bool contact_ramp,
                         double ramp_criterion,
                         double ramp_criterion_t,
                         double ramp_factor,
                         double ramp_factor_t,
                         double max_penalty,
                         double max_penalty_t
                       )
            {
                
                ConditionsArrayType& ConditionsArray = mr_model_part.Conditions();
                
                /**
                 * Update of normal stress lagrange multipliers
                 */
                for (ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
                {
                    if( (*it)->GetValue( IS_CONTACT_SLAVE ) )
                    {
                        for( IndexType IntPoint = 0; 
                             IntPoint < (*it)->GetGeometry().IntegrationPoints().size();
                             ++IntPoint )
                        {

                            double newLambda = (*it)->GetValue( LAMBDAS )[IntPoint]
                                        +(*it)->GetValue( GAPS )[IntPoint]
                                        *(*it)->GetValue( PENALTY )[IntPoint];
                            if(newLambda < 0.0)
                            {
                                newLambda=0.0;
                            }

                            (*it)->GetValue( DELTA_LAMBDAS )[IntPoint] = newLambda
                                    - (*it)->GetValue( LAMBDAS )[IntPoint];
                            (*it)->GetValue( LAMBDAS )[IntPoint] = newLambda;
                        }
                    }
                }
                
                 /**
                 * Update of frictional stress lagrange multipliers
                 * 
                 */


//                 KRATOS_WATCH(BaseType::GetModelPart().GetMesh().GetProperties(2)[FRICTION_COEFFICIENT] );
                if( friction_coefficient > 0.0 )
                {
                    for (ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin()+lastRealCondition; 
                         it!=ConditionsArray.ptr_end(); ++it)
                    {
                        Vector newLambda_T_trial=ZeroVector(2);

                        Matrix m = (*it)->GetValue(CONTACT_LINK_M);

						Vector tangentialVelocity= GetRelativTangentialVelocity((*it)->GetValue(CONTACT_LINK_MASTER), (*it)->GetValue(CONTACT_LINK_SLAVE), (*it)->GetValue( SLAVE_CONTACT_LOCAL_POINT), (*it)->GetValue( MASTER_CONTACT_LOCAL_POINT));

                        newLambda_T_trial[0] =
                                (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(LAMBDAS_T)( (*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX),0) 
                                + (*it)->GetValue( CONTACT_LINK_SLAVE)->GetValue( PENALTY_T )[(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX)] 
                                * tangentialVelocity(0);
                             
                        newLambda_T_trial[1] =
                                (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(LAMBDAS_T)( (*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX),1)
                                + (*it)->GetValue( CONTACT_LINK_SLAVE)->GetValue( PENALTY_T )[(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX)]
                                * tangentialVelocity(1);

                        double NormLambda_T = sqrt( newLambda_T_trial[0]*m(0,0)*newLambda_T_trial[0]
                                    + newLambda_T_trial[0]*m(0,1)*newLambda_T_trial[1]
                                    + newLambda_T_trial[1]*m(1,0)*newLambda_T_trial[0]
                                    + newLambda_T_trial[1]*m(1,1)*newLambda_T_trial[1]);//;

                        if( ( NormLambda_T <= (((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS )[(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX)]) )*friction_coefficient))
                        {
                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( DELTA_LAMBDAS_T )( (*it)->GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX),0) = newLambda_T_trial[0]-(*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX),0);

                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( DELTA_LAMBDAS_T )( (*it)->GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX),1) = newLambda_T_trial[1]-(*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX),1);
                            
                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX),0) = newLambda_T_trial[0];
                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue( CONTACT_SLAVE_INTEGRATION_POINT_INDEX),1) = newLambda_T_trial[1];
                        }
                        else
                        {
                            newLambda_T_trial[0] = friction_coefficient
                                    * (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS )[(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX)]
                                    *newLambda_T_trial[0]/(NormLambda_T);

                            
                            newLambda_T_trial[1] = friction_coefficient
                                    * (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS )[(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX)]
                                    *newLambda_T_trial[1]/(NormLambda_T);
                            
                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( DELTA_LAMBDAS_T )( (*it)->GetValue(         
                                    CONTACT_SLAVE_INTEGRATION_POINT_INDEX),0) = newLambda_T_trial[0]-(*it)->GetValue( 
                                    CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue( 
                                    CONTACT_SLAVE_INTEGRATION_POINT_INDEX),0);
                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( DELTA_LAMBDAS_T )( (*it)->GetValue( 
                                    CONTACT_SLAVE_INTEGRATION_POINT_INDEX),1) = newLambda_T_trial[1]-(*it)->GetValue( 
                                    CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue( 
                                    CONTACT_SLAVE_INTEGRATION_POINT_INDEX),1);

                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue(
                                    CONTACT_SLAVE_INTEGRATION_POINT_INDEX),0) = newLambda_T_trial[0];
                            (*it)->GetValue( CONTACT_LINK_SLAVE )->GetValue( LAMBDAS_T )( (*it)->GetValue(
                                    CONTACT_SLAVE_INTEGRATION_POINT_INDEX),1) = newLambda_T_trial[1];
                        }
                    }
                }

                if(contact_ramp)
                {
                    std::cout << "##### ramping penalties ######" << std::endl;
                    
                    double alpha= ramp_criterion;
                    double alpha_T= ramp_criterion_t;
                    double beta= ramp_factor;
                    double beta_T= ramp_factor_t;
                    double Kmax=0.0;
                    double Kmax_T= 0.0;
//                     double max_penalty= max_penalty;
                    double max_penalty_T= max_penalty_t;
                    bool check_penalty = false;
                    bool check_penalty_T = false;
                    
                    for (ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
                    {
                        if( (*it)->GetValue( IS_CONTACT_SLAVE ) )
                        {
                            for( IndexType IntPoint = 0; 
                                 IntPoint < (*it)->GetGeometry().IntegrationPoints().size();
                                 ++IntPoint )
                            {
                                if(((*it)->GetValue( DELTA_LAMBDAS )[IntPoint]) > Kmax)
                                {
                                    Kmax = (*it)->GetValue( DELTA_LAMBDAS )[IntPoint];
                                }
                            }
                        }
                    }
                    if( friction_coefficient > 0.0 )
                    {
                        for (ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin()+lastRealCondition; 
                             it!=ConditionsArray.ptr_end(); ++it)
                        {
                            int i=(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX);
                            Matrix m = (*it)->GetValue(CONTACT_LINK_M);
                            double K_trial = sqrt((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                        *m(0,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                        +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                        *m(0,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                        +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                        *m(1,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                        +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                        *m(1,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(DELTA_LAMBDAS_T )(i,1));
                            if(K_trial > Kmax_T)
                            {
                                Kmax_T=K_trial;
                            }
                        }
                    }
                    KRATOS_WATCH( Kmax );
                    KRATOS_WATCH( k_contact );
                    if( Kmax >= k_contact)
                    {
                        std::cout << "###### KMAX > K_CONTACT #####" << std::endl;
                        for (ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
                        {
                            if( (*it)->GetValue( IS_CONTACT_SLAVE ) )
                            {
                                for( IndexType IntPoint = 0; 
                                     IntPoint < (*it)->GetGeometry().IntegrationPoints().size();
                                     ++IntPoint )
                                { 
                                    if( ((*it)->GetValue( DELTA_LAMBDAS )[IntPoint]) > 
                                           (k_contact/alpha))
                                    {
                                        (*it)->GetValue( PENALTY)[IntPoint] = beta*(*it)->GetValue( PENALTY )[IntPoint] ;
                                        std::cout << "new penalty: " << (*it)->GetValue( PENALTY)[IntPoint] << std::endl;
                                            
                                        if((*it)->GetValue( PENALTY)[IntPoint] > max_penalty)
                                        {
                                            (*it)->GetValue( LAMBDAS )[IntPoint]=((*it)->GetValue( LAMBDAS )[IntPoint]
                                                    -(*it)->GetValue( DELTA_LAMBDAS )[IntPoint])
                                                    +(*it)->GetValue( DELTA_LAMBDAS )[IntPoint]
                                                    /(max_penalty/
                                                    ((*it)->GetValue( PENALTY)[IntPoint])*beta);
                                            (*it)->GetValue( PENALTY)[IntPoint]=max_penalty;
                                        }
                                        else
                                        {
                                            (*it)->GetValue( LAMBDAS )[IntPoint]=((*it)->GetValue( LAMBDAS )[IntPoint]
                                                    -(*it)->GetValue( DELTA_LAMBDAS )[IntPoint])
                                                    +(*it)->GetValue( DELTA_LAMBDAS )[IntPoint]
                                            /beta;
                                        }
                                    }
                                }
                            } 
                        }
                        check_penalty =true;
                    }
                    
                    if( friction_coefficient > 0.0 )  
                    {
                        if(Kmax_T >= k_contact_t)
                        {
                            for (ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin()+lastRealCondition; 
                                 it!=ConditionsArray.ptr_end(); ++it)
                            {
                                int i=(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX);
                                Matrix m = (*it)->GetValue(CONTACT_LINK_M);

                                double K = sqrt((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            *m(0,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            *m(0,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                            +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                            *m(1,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                            *m(1,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(DELTA_LAMBDAS_T )(i,1));
                                
                                    
                                if( K > (k_contact_t/alpha_T))
                                {
                                    (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i]=beta_T* (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T )[i] ;

                                    if((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i]>max_penalty_T)
                                    {
                                        (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)
                                                -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0))
                                                +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0) /(max_penalty_T/((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i])*beta_T);
                                        (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)
                                                -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1))
                                                +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                        /(max_penalty_T/((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i])*beta_T);
                                        (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i]=max_penalty_T;
                                    }
                                    else
                                    {
                                        (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)
                                                -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0))
                                                +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                        /beta_T;
                                        (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)
                                                -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1))
                                                +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)/beta_T;
                                    }
                                }
                            }
                            check_penalty_T =true;
                        }
                    }

                    if(Kmax <= k_contact/alpha || !check_penalty)
                    {
                        k_contact = Kmax;
                        check_penalty =true;
                    }
                    
                    if( mr_model_part.GetProperties(1)[FRICTION_COEFFICIENT] > 0.0 )  
                    {
                        if(Kmax_T <= k_contact_t/alpha_T || !check_penalty_T)
                        {
                            k_contact_t = Kmax_T;
                            check_penalty_T =true;
                        }
                    }
                    
                    if(!check_penalty)
                    {
                        for (ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
                        {
                            if( (*it)->GetValue( IS_CONTACT_SLAVE ) )
                            {
                                for( IndexType IntPoint = 0; 
                                     IntPoint < (*it)->GetGeometry().IntegrationPoints().size();
                                     ++IntPoint )
                                { 
                                    if(((*it)->GetValue( DELTA_LAMBDAS )[IntPoint]) > (k_contact/alpha))
                                    {
                                        (*it)->GetValue( PENALTY )[IntPoint]=beta* ((*it)->GetValue( PENALTY )[IntPoint]);


                                        if((*it)->GetValue( PENALTY)[IntPoint]>max_penalty)
                                        {
                                            (*it)->GetValue( LAMBDAS )[IntPoint]=((*it)->GetValue( LAMBDAS )[IntPoint])
                                                    -(*it)->GetValue( DELTA_LAMBDAS )[IntPoint]
                                                    +(*it)->GetValue( DELTA_LAMBDAS )[IntPoint]
                                                    /(max_penalty/((*it)->GetValue( PENALTY)[IntPoint])*beta);
                                            (*it)->GetValue( PENALTY)[IntPoint]=max_penalty;
                                        }
                                        else
                                        {
                                            (*it)->GetValue( LAMBDAS )[IntPoint]=(*it)->GetValue( LAMBDAS )[IntPoint]
                                                    -(*it)->GetValue( DELTA_LAMBDAS )[IntPoint]
                                                    +(*it)->GetValue( DELTA_LAMBDAS )[IntPoint]
                                                    /beta;
                                        }
                                    } 
                                }
                            } 
                        }
                        k_contact = Kmax;
                        check_penalty =true;
                    }

                    if( friction_coefficient > 0.0 )  
                    {
                        if(!check_penalty_T)
                        {
                            for (ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin()+lastRealCondition; 
                                 it!=ConditionsArray.ptr_end(); ++it)
                            {
                                int i=(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX);
                                Matrix m = (*it)->GetValue(CONTACT_LINK_M);
                                double K = sqrt((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            *m(0,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            *m(0,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                            +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                            *m(1,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                            +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                            *m(1,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(DELTA_LAMBDAS_T )(i,1));
                                
                                if(K > (k_contact_t/alpha_T) )
                                {
                                    (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T )[i] = beta_T* ((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T )[i]);
                                    if((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i] > max_penalty_T)
                                    {
                                        (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0))
                                                -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                                +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                                        /(max_penalty_T/((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i])*beta_T);
                                        (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1))
                                                -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                                +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                                    /(max_penalty_T/((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i])*beta_T);
                                        (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( PENALTY_T)[i] = max_penalty_T;
                                    }
                                    else
                                    {
                                        (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0))
                                                -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                                +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                                   /beta_T;
                                        (*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)=((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1))
                                                -(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                                +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                                   /beta_T;
                                    }
                                } 
							
                            }
                            k_contact_t = Kmax_T;
                            check_penalty =true;
                        }
                    }
                }//if(BaseType::GetModelPart().GetMesh().GetProperties(1)[CONTACT_RAMP])
                return;
            }//Update
            
//**********************************************************************
//**********************************************************************
            
            /** 
             * Check whether the Uzawa loop is converged, check regarding the fulfillment of the 
             *constraints, the incremental update of the Lagrangian multipliers and the penalty energies
             * @param mr_model_part the model part
             * @param step SolutionStepNumber
             * @param lastRealCondition number of the last condition that is not a contact linking
             * @param friction_coefficient the fiction coefficient
             * @return is converged
             */
            bool IsConverged( ModelPart& mr_model_part, int step , int lastRealCondition, double friction_coefficient)
            {
                bool Converged = false;
                bool friction = false;
                double ratio = 0;
                double absolute = 0;
                double energy_contact = 0;
                double ratio_friction=0;
                double absolute_friction=0;
                double energy_friction = 0;
                double wriggers_crit= 0.0;
                double cumulative_penetration = 0.0;
//                 double gap = 0.0;
                double slip= 0.0;
                int Index = 1;
                int Index2 = 1;

                ConditionsArrayType& ConditionsArray = mr_model_part.Conditions();
                
                for( ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin();
                     it!=ConditionsArray.ptr_end(); ++it)
                {
                    if( (*it)->GetValue( IS_CONTACT_SLAVE ) )
                    {
                        for( IndexType i=0; i<(*it)->GetGeometry().IntegrationPoints().size();
                             i++ )
                        {
                            if((*it)->GetValue( LAMBDAS )[i]>0.0)
                            {
                                energy_contact += ((*it)->GetValue( GAPS )[i])
                                        * ((*it)->GetValue( GAPS )[i])
                                        *0.5*((*it)->GetValue(PENALTY)[i]);
                                absolute += (*it)->GetValue( LAMBDAS )[i];
                                ratio += (*it)->GetValue( DELTA_LAMBDAS )[i];

                                Index++;
                                cumulative_penetration += (*it)->GetValue( GAPS )[i]*(*it)->GetValue( GAPS )[i];

                                if(fabs(wriggers_crit) < fabs(((*it)->GetValue( GAPS )[i])))
                                        wriggers_crit= ((*it)->GetValue( GAPS )[i]);
                            }
                        }
                    }
                }
                
                if( friction_coefficient > 0.0 )
                {
                    friction=true;
                    
                    for (ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin() + lastRealCondition; 
                         it!=ConditionsArray.ptr_end(); ++it)
                    {
                        int i=(*it)->GetValue(CONTACT_SLAVE_INTEGRATION_POINT_INDEX);
                        
                        if((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS )[i]  > 0)
                        {
                            Matrix m = (*it)->GetValue(CONTACT_LINK_M);
                            absolute_friction +=
                                    sqrt((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)
                                    *m(0,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)
                                    +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)
                                    *m(0,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)
                                    +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)
                                    *m(1,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,0)
                                    +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1)
                                    *m(1,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( LAMBDAS_T )(i,1));
                            
                            ratio_friction +=  
                                    sqrt((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                    *m(0,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                    +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                    *m(0,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                    +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                    *m(1,0)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,0)
                                    +(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1)
                                    *m(1,1)*(*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue( DELTA_LAMBDAS_T )(i,1));
                            
                            if((*it)->GetValue(CONTACT_LINK_SLAVE)-> GetValue(STICK)(i) > 0.5)
                            {
                                Vector relVelo(2);
                                noalias(relVelo)=GetRelativTangentialVelocity((*it)->GetValue(CONTACT_LINK_MASTER), (*it)->GetValue(CONTACT_LINK_SLAVE), (*it)->GetValue( SLAVE_CONTACT_LOCAL_POINT), (*it)->GetValue( MASTER_CONTACT_LOCAL_POINT));
                                
                                energy_friction += 0.5*
                                        ((*it)->GetValue(CONTACT_LINK_SLAVE)->GetValue(PENALTY_T)[i])*(relVelo(0)*m(0,0)*relVelo(0)
                                        +relVelo(0)*m(0,1)*relVelo(1)
                                        +relVelo(1)*m(1,0)*relVelo(0)
                                        +relVelo(1)*m(1,1)*relVelo(1));
                                
                                slip+= (relVelo(0)*m(0,0)*relVelo(0)
                                        +relVelo(0)*m(0,1)*relVelo(1)
                                        +relVelo(1)*m(1,0)*relVelo(0)
                                        +relVelo(1)*m(1,1)*relVelo(1)); 
                                Index2++;
                            }
                        }
                    }
                }
                if( this->mEchoLevel > 1 )
                {
                    std::cout << "absolute Lambda: " << absolute/Index << std::endl;
                    std::cout << "relative Lambda: " << ratio/absolute << std::endl;
                    std::cout << "energy criterion: " << energy_contact/Index << std::endl;
                    std::cout << "normed gap: " << wriggers_crit<<std::endl;
                    std::cout << "weighted mean penetration: " << sqrt(cumulative_penetration/Index) << std::endl;
                    if( friction )
                    {
                        std::cout << "absolute Lambda friction: " << absolute_friction/Index << std::endl;
                        std::cout << "relative Lambda friction: " << ratio_friction/absolute_friction << std::endl;
                        std::cout << "energy criterion friction: " << energy_friction/Index2 << std::endl;
                        std::cout << "forbidden slip: " << sqrt(slip)/Index2  << std::endl;
                    }
                }
                else if( this->mEchoLevel > 0 )
                {
                    if( friction )
                    {
                        std::cout << "relative Lambda friction: " << ratio_friction/absolute_friction << std::endl;
                    }
                }
                if( friction_coefficient > 0.0 )
                {
                    if( (/*fabs(absolute/Index) < 1e-3 || fabs(ratio/absolute) < 1e-2) 
                          && */ sqrt(cumulative_penetration/Index) < 1e-5)
                          /*&& ((fabs(absolute_friction/Index) < 1e-3 ||
                          fabs(ratio_friction/absolute_friction) < 1e-2)*/
                          && sqrt(slip)/Index2< 1e-5 )
                    {
                        if( this->mEchoLevel > 0 )
                            std::cout << "Contact condition converged after " << step+1 << " steps" << std::endl;
                        Converged = true;
                    }
                    else
                    {
                        if( this->mEchoLevel > 0 )
                            std::cout << "Contact has been detected, next solution step required (UZAWA STEP: "<< step+1 << ")" << std::endl;
                    }
                }
                else
                {
                    if(/*(fabs(absolute/Index) < 1e-3 || fabs(ratio/absolute) < 1e-2) &&*/ /*fabs(wriggers_crit) < 1e-4)*/sqrt(cumulative_penetration/Index) < 1e-5 )
//                         sqrt(wriggers_crit)/Index < 3e-4)
                    {
                        if( this->mEchoLevel > 0 )
                            std::cout << "Contact condition converged after " << step+1 << " steps" << std::endl;
                        Converged = true;
                    }
                    else
                    {
                        if( this->mEchoLevel > 0 )
                            std::cout << "Contact has been detected, next solution step required (UZAWA STEP: "<< step+1 << ")" << std::endl;
                    }
                }
                return( Converged );
            }
            
            
            /**
             * This function cleans up the conditions list after uzawa iteration
             * @param mr_model_part the modelpart
             * @param lastRealCondition number of the last condition that is not a contact linking, 
             *   all conditions having a greater number/index will be erased
             */
            void Clean( ModelPart& mr_model_part, int lastRealCondition )
            {
                // cleaning up model part 
                mr_model_part.Conditions().erase(
                                         mr_model_part.Conditions().begin()+lastRealCondition,
                        mr_model_part.Conditions().end() );
            }
            
            
        private:
            /**
             * searches a contact partner for a given slave condition, 
             * by performing an outer and an inner search
             * @param mr_model_part the modelpart
             * @param Slave the current slave surfaces
             * @param AllMasterElements a container including all Master elements
             * @param IntegrationPointIndex Index of the current quadrature point on the slave surface
             * @param MasterContactLocalPoint local coordinates of the closest point projection
             *      on the master surface (output)
             * @param SlaveContactLocalPoint local coordinates of the quadrature point 
             *      on the slave surface
             * @param CurrentMaster pointer on the master facet where the closest point projection 
             *              is located
             */
            bool SearchPartner( ModelPart& mr_model_part, 
                                Condition& Slave,
                                ConditionsArrayType& AllMasterElements, 
                                const IndexType& IntegrationPointIndex,
                                Point<3>& MasterContactLocalPoint,
                                Point<3>& SlaveContactLocalPoint,
                                Condition::Pointer& CurrentMaster
                              )
            {
//                 std::cout << "searching partner" << std::endl;
                 KRATOS_TRY
                bool PartnerExists = false;
                
                //checking for existent master surfaces
                if( AllMasterElements.size() > 0 )
                {
                    PartnerExists = true; 
                    SlaveContactLocalPoint 
                            = Slave.GetGeometry().IntegrationPoints()[IntegrationPointIndex];
//                    KRATOS_WATCH(SlaveContactLocalPoint);
                    //calculating global coordinates of current integration point 
                    Point<3> SlaveContactGlobalPoint;

                    SlaveContactGlobalPoint = GlobalCoordinates( Slave, SlaveContactGlobalPoint, SlaveContactLocalPoint);

                    Point<3> GlobalCandidate;
                    //defining set of possible master surface elements
                    ConditionsArrayType::Pointer MasterSet( 
                            new ConditionsArrayType() );
                    double minDist = static_cast<double>(INT_MAX);
//                     KRATOS_WATCH(minDist);
//                     double minDist = INT_MAX;
                    //loop over all master surfaces (global search)
                    for( ConditionsArrayType::ptr_iterator it =
                         AllMasterElements.ptr_begin(); 
                         it != AllMasterElements.ptr_end();
                         ++it )
                    {
                        //loop over all nodes in current master surface
                        for( unsigned int n=0; n<(*it)->GetGeometry().PointsNumber(); n++ )
                        {
                            double dist = (((*it)->GetGeometry().GetPoint(n).X0()+
								(*it)->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
								-SlaveContactGlobalPoint[0])
                         		*
								(((*it)->GetGeometry().GetPoint(n).X0()+
								(*it)->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
								-SlaveContactGlobalPoint[0])
                         		+ 
								(((*it)->GetGeometry().GetPoint(n).Y0()+
								(*it)->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
								-SlaveContactGlobalPoint[1])
                         		* 
								(((*it)->GetGeometry().GetPoint(n).Y0()+
								(*it)->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
								-SlaveContactGlobalPoint[1])
                         		+ 
							(((*it)->GetGeometry().GetPoint(n).Z0()+
							(*it)->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
							-SlaveContactGlobalPoint[2])
                     	   	* 
							(((*it)->GetGeometry().GetPoint(n).Z0()+
							(*it)->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
							-SlaveContactGlobalPoint[2]);
//                             KRATOS_WATCH(dist);
                        	if( fabs(dist-minDist) < 1e-4 )
                       	 	{
                              	MasterSet->push_back(*it);
                        	} 
                        	else if( dist < minDist )
                        	{
                           		MasterSet->clear();

                           		GlobalCandidate = (*it)->GetGeometry().GetPoint(n);
//                                 KRATOS_WATCH( GlobalCandidate );
                           		minDist = dist;
                           
			   					MasterSet->push_back(*it);
                         	}
                      	}
                    }
                    //searching contact partner (local search)
                    Point<3> MasterContactGlobalPoint;
                    bool LocalPartnerExists = false;
//                    KRATOS_WATCH( GlobalCandidate );
                    for( ConditionsArrayType::ptr_iterator it = MasterSet->ptr_begin(); it != MasterSet->ptr_end(); ++it )
                    {
                        if( ClosestPoint( mr_model_part, *it, MasterContactGlobalPoint, 
                            MasterContactLocalPoint, SlaveContactGlobalPoint, GlobalCandidate ) )
                        {
                            CurrentMaster = *it;
                            LocalPartnerExists = true;
                            break;
                        }
                    }
                    PartnerExists = ( PartnerExists && LocalPartnerExists ); 
                }
                return PartnerExists;
                KRATOS_CATCH("")
            }//SearchPartner
            
            
            /**
             * This method searches for a given quadrature point the closest point 
             *   projection of a master surface. 
             * @param Surfrace given master surfaces
             * @param rResultGlobal global coordinates of the closest point projection
             * @param rResultLocal local coordinates of the closest point projection
             * @param rSlaveContactGlobalPoint global coordinates of the quadrature 
             *   point on the slave surface
             * @param rCandidateGlobal Global coordinates of the node that is closest to 
             *   the quadrature point on the slave surface
             */
            bool ClosestPoint( ModelPart& mr_model_part,
                               Condition::Pointer& Surface,
                               GeometryType::CoordinatesArrayType& rResultGlobal, 
                               GeometryType::CoordinatesArrayType& rResultLocal,
                               const GeometryType::CoordinatesArrayType& rSlaveContactGlobalPoint,
                               const GeometryType::CoordinatesArrayType& rCandidateGlobal
                             )
            {
               	double Xi1 = 0.0;
                double Xi2 = 0.0;
                double deltaXi1 = 0.0;
                double deltaXi2 = 0.0;
                Matrix localCoords;

                localCoords = (Surface->GetGeometry()).PointsLocalCoordinates( localCoords );
                //determining local coordinates for rResult
                for( unsigned int n=0; n<Surface->GetGeometry().PointsNumber(); n++ )
                {
                    if(    	fabs(rCandidateGlobal[0]-
							(Surface->GetGeometry().GetPoint(n).X0()
							+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
							) < 1e-7
                           	&& fabs(rCandidateGlobal[1]-
							(Surface->GetGeometry().GetPoint(n).Y0()
							+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
							) < 1e-7
                           	&& fabs(rCandidateGlobal[2]-
							(Surface->GetGeometry().GetPoint(n).Z0()
							+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
							) < 1e-7
                      )
                    {
                        Xi1 = localCoords(n,0);
                        Xi2 = localCoords(n,1);

                        break;
                    }
                }
                //setting up LocalPoint
                rResultLocal[0] = Xi1;
                rResultLocal[1] = Xi2;
                rResultLocal[2] = 0.0;
                //setting up rResult
                rResultGlobal = rCandidateGlobal;
//                 KRATOS_WATCH( rCandidateGlobal );
                //searching for orthogonal projection
                for( int k=0; k<1000; k++ )
				{
                    //setting up tangential vectors
                    Vector t1 = ZeroVector(3);//first tangential vector
                    Vector t2 = ZeroVector(3);//second tangential vector
                    //derivatives of tangential vectors
                    Vector dt11 = ZeroVector(3);
                    Vector dt12 = ZeroVector(3);
                    Vector dt21 = ZeroVector(3);
                    Vector dt22 = ZeroVector(3);
            
                    //retrieving first order derivatives in current solution point 
                    Matrix DN = ZeroMatrix(Surface->GetGeometry().PointsNumber(),2);
                    Surface->GetGeometry().ShapeFunctionsLocalGradients( DN, rResultLocal );
                    //retrieving second order derivatives in current solution point 
                    GeometryType::ShapeFunctionsSecondDerivativesType D2N;
                    Surface->GetGeometry().ShapeFunctionsSecondDerivatives( D2N, rResultLocal );
                    for( unsigned int n=0; n<Surface->GetGeometry().PointsNumber(); n++ )
                    {
                        //contribution to tangential vectors
                        t1[0] += (Surface->GetGeometry().GetPoint(n).X0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
								*DN(n,0);
                        t1[1] += (Surface->GetGeometry().GetPoint(n).Y0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
								*DN(n,0);
                        t1[2] += (Surface->GetGeometry().GetPoint(n).Z0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
								*DN(n,0);
                        t2[0] += (Surface->GetGeometry().GetPoint(n).X0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
								*DN(n,1);
                        t2[1] += (Surface->GetGeometry().GetPoint(n).Y0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
								*DN(n,1);
                        t2[2] += (Surface->GetGeometry().GetPoint(n).Z0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
								*DN(n,1);
                        //contribution to derivatives of tangential vectors
                        dt11[0] += (Surface->GetGeometry().GetPoint(n).X0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
								*D2N[n](0,0);
                        dt11[1] += (Surface->GetGeometry().GetPoint(n).Y0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
								*D2N[n](0,0);
                        dt11[2] += (Surface->GetGeometry().GetPoint(n).Z0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
								*D2N[n](0,0);
                        dt12[0] += (Surface->GetGeometry().GetPoint(n).X0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
								*D2N[n](0,1);
                        dt12[1] += (Surface->GetGeometry().GetPoint(n).Y0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
								*D2N[n](0,1);
                        dt12[2] += (Surface->GetGeometry().GetPoint(n).Z0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
								*D2N[n](0,1);
                        dt21[0] += (Surface->GetGeometry().GetPoint(n).X0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
								*D2N[n](1,0);
                        dt21[1] += (Surface->GetGeometry().GetPoint(n).Y0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
								*D2N[n](1,0);
                        dt21[2] += (Surface->GetGeometry().GetPoint(n).Z0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
								*D2N[n](1,0);
                        dt22[0] += (Surface->GetGeometry().GetPoint(n).X0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
								*D2N[n](1,1);
                        dt22[1] += (Surface->GetGeometry().GetPoint(n).Y0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
								*D2N[n](1,1);
                        dt22[2] += (Surface->GetGeometry().GetPoint(n).Z0()
								+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
								*D2N[n](1,1);
                    }
                    //defining auxiliary terms
                    double A1 = ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*t1[0])
                                +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*t1[1])
                                +((rSlaveContactGlobalPoint[2]-rResultGlobal[2])*t1[2]);
                    double A2 = ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*t2[0])
                                +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*t2[1])
                                +((rSlaveContactGlobalPoint[2]-rResultGlobal[2])*t2[2]);
                    double B11 = (-t1[0]*t1[0]-t1[1]*t1[1]-t1[2]*t1[2])
                                + ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*dt11[0])
                                +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt11[1])
                                +((rSlaveContactGlobalPoint[2]-rResultGlobal[2])*dt11[2]);
                    double B12 = (-t2[0]*t1[0]-t2[1]*t1[1]-t2[2]*t1[2])
                                + ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*dt12[0])
                                +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt12[1])
                                +((rSlaveContactGlobalPoint[2]-rResultGlobal[2])*dt12[2]);
                    double B21 = (-t1[0]*t2[0]-t1[1]*t2[1]-t1[2]*t2[2])
                                + ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*dt21[0])
                                +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt21[1])
                                +((rSlaveContactGlobalPoint[2]-rResultGlobal[2])*dt21[2]);
                    double B22 = (-t2[0]*t2[0]-t2[1]*t2[1]-t2[2]*t2[2])
                                + ((rSlaveContactGlobalPoint[0]-rResultGlobal[0])*dt22[0])
                                +((rSlaveContactGlobalPoint[1]-rResultGlobal[1])*dt22[1])
                                +((rSlaveContactGlobalPoint[2]-rResultGlobal[2])*dt22[2]);
                    //calculating update for Xi
                    deltaXi1 = -A1*B22/(B11*B22-B12*B21)+A2*B12/(B11*B22-B12*B21);
                    deltaXi2 =  A2*B21/(B11*B22-B12*B21)-A2*B11/(B11*B22-B12*B21);
                    //updating Xi
                    Xi1 += deltaXi1;
                    Xi2 += deltaXi2;
                    //updating LocalPoint
                    rResultLocal[0] = Xi1;
                    rResultLocal[1] = Xi2;
                    //updating rResult
                    rResultGlobal = ZeroVector( 3 );

                    rResultGlobal = GlobalCoordinates(*(Surface), rResultGlobal, rResultLocal );

		    		if( fabs(deltaXi1) < 1e-7 && fabs(deltaXi2) < 1e-7 )
                    {
                        //check whether contact point lies within elementary boundaries
                        if( (Surface->GetGeometry().size()==3 ||Surface->GetGeometry().size()==6) && (Xi1 <= 1.0 && Xi1 >= 0.0) && (Xi2 <= 1.0 && Xi2 >= 0.0) && ((Xi1+Xi2) <= 1.0))
                        {
                            return true;
                        }
						else if( (Surface->GetGeometry().size()==4 ||Surface->GetGeometry().size()==8 ||Surface->GetGeometry().size()==9) && (fabs(Xi1) <= 1.0 ) && (fabs(Xi2) <= 1.0) )
                        {
                            return true;
                        }
                        else
                        {
                            return false;
                        }
                    }
                }
                return false;
            }//ClosestPoint
            
    /**
    * Calculates for given Loacal coordinates the global coordinates 
    * @param Surface surface
    * @param rResult global coordinates
    * @param LocalCoordinates local coordinates
    * @return global coordinates
    */

     	GeometryType::CoordinatesArrayType& GlobalCoordinates(Condition& Surface, GeometryType::CoordinatesArrayType& rResult, GeometryType::CoordinatesArrayType const& LocalCoordinates)
    	{
		noalias(rResult)= ZeroVector(3);

		for(IndexType i = 0 ; i < Surface.GetGeometry().size() ; i++)
		{
			double shape_func= Surface.GetGeometry().ShapeFunctionValue(i,LocalCoordinates);

			rResult(0) += shape_func* 
				((Surface.GetGeometry()[i]).X0()
				+(Surface.GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_X));

			rResult(1) += shape_func* 
				((Surface.GetGeometry()[i]).Y0()
				+(Surface.GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_Y));

			rResult(2) += shape_func* 
				((Surface.GetGeometry()[i]).Z0()
				+(Surface.GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_Z));
		}
		return rResult;
    	}

    /**
    * Computes the relative tangential velocity between two points, one on the master surface 
    *    and one on the slave surface
    * @param Master surface
    * @param Slave surface
    * @param SlaveLocalCoordinates
    * @param MasterLocalCoordinates
    * @return relative tangential velocity between the point on the master surface, 
    *    and the point on the slave surface
    */

		Vector GetRelativTangentialVelocity(Condition::Pointer Master, Condition::Pointer Slave, Point<3> const& SlaveLocalCoordinates, Point<3> const& MasterLocalCoordinates)
		{

			Matrix T= TangentialVectors( Master, MasterLocalCoordinates );

			Vector result(2);

			Vector slave_velo(3);
			Vector master_velo(3);

			noalias(slave_velo)= ZeroVector(3);
			noalias(master_velo)= ZeroVector(3);

			for(IndexType i = 0 ; i < Slave->GetGeometry().size() ; i++)
			{
				double shape_func= Slave->GetGeometry().ShapeFunctionValue(i,SlaveLocalCoordinates );

				slave_velo+= ((Slave->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_DT))*
					shape_func;
			}
			for(IndexType i = 0 ; i < Master->GetGeometry().size() ; i++)
			{
				double shape_func= Master->GetGeometry().ShapeFunctionValue(i,MasterLocalCoordinates);
				master_velo+= ((Master->GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_DT))*
					shape_func;
			}
		
			Vector norm_T(2);

			norm_T(0)= sqrt(T(0,0)*T(0,0)+T(0,1)*T(0,1)+T(0,2)*T(0,2));

			norm_T(1)= sqrt(T(1,0)*T(1,0)+T(1,1)*T(1,1)+T(1,2)*T(1,2));

			result(0)= ((slave_velo(0)-master_velo(0))*T(0,0)+(slave_velo(1)-master_velo(1))*T(0,1)
				+(slave_velo(2)-master_velo(2))*T(0,2))/norm_T(0);

			result(1)= ((slave_velo(0)-master_velo(0))*T(1,0)+(slave_velo(1)-master_velo(1))*T(1,1)
				+(slave_velo(2)-master_velo(2))*T(1,2))/norm_T(1);

			return result;
		}

        /**
        * Computes the tangential vectors for a given point
        * @param Surface
        * @param rPoint Local Coordinates of the point on the surface
        * @return Tangential vectors in point rPoint on Surface
        */

   	 	Matrix TangentialVectors( Condition::Pointer Surface, 
                                             const GeometryType::CoordinatesArrayType& rPoint )
    	{
      	  	//setting up result matrix
        	Matrix T = ZeroMatrix( 2, 3 );
        //shape function gradients
       	 	Matrix DN = ZeroMatrix( Surface->GetGeometry().PointsNumber(),2);
        	Surface->GetGeometry().ShapeFunctionsLocalGradients( DN, rPoint );
        	//calculating tangential vectors
        	for( unsigned int n=0; n<Surface->GetGeometry().PointsNumber(); n++ )
        	{
            	T(0,0) += (Surface->GetGeometry().GetPoint(n).X0()
					+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
					*DN(n,0);
            	T(0,1) += (Surface->GetGeometry().GetPoint(n).Y0()
					+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
					*DN(n,0);
            	T(0,2) += (Surface->GetGeometry().GetPoint(n).Z0()
					+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
					*DN(n,0);
            	T(1,0) += (Surface->GetGeometry().GetPoint(n).X0()
					+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_X))
					*DN(n,1);
            	T(1,1) += (Surface->GetGeometry().GetPoint(n).Y0()
					+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Y))
					*DN(n,1);
            	T(1,2) += (Surface->GetGeometry().GetPoint(n).Z()
					+Surface->GetGeometry().GetPoint(n).GetSolutionStepValue(DISPLACEMENT_Z))
					*DN(n,1);
        }
        return( T );
    }


            int mEchoLevel;
            double k_contact;
            double k_contact_t;
            
    };//class ContactUtility
}  /* namespace Kratos.*/

#endif /* KRATOS_CONTACT_UTILITY_H_INCLUDED  defined */
