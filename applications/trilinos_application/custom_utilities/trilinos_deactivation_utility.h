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
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hurga $
//   Date:                $Date: 2009-02-26 13:49:42 $
//   Revision:            $Revision: 1.7 $
//
//


#if !defined(KRATOS_TRILINOS_DEACTIVATION_UTILITY_INCLUDED )
#define  KRATOS_TRILINOS_DEACTIVATION_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <algorithm>
#include <fstream>
#include <cmath>

#if !defined(isnan)
#define isnan(x) ((x)!=(x))
#endif
// External includes 

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "trilinos_application.h"

namespace Kratos
{
    /**
     * Deactivation and reactivation of elements
     * This process handles deactivation and reactivation
     * of elements and associated conditions.
     * In order to use this process, a variable calles
     * ACTIVATION_LEVEL of type int has to be defined.
     * The following values can be assigned to ACTIVATION_LEVEL:
     * ACTIVATION_LEVEL == 0: element is always active
     * ACTIVATION_LEVEL > 0: element will be deactivated on demand
     * ACTIVATION_LEVEL < 0: element is initially deactivated and
     *                       will be reactivated on demand
     */
    class TrilinosDeactivationUtility
    {
        public:
            
            typedef PointerVectorSet<Element> ElementsArrayType;
            typedef PointerVectorSet<Condition> ConditionsArrayType;
            
            /**
             * class pointer definition
             */
            KRATOS_CLASS_POINTER_DEFINITION( TrilinosDeactivationUtility );
            
            /**
             * Constructor.
             * The constructor takes the current model_part as argument.
             * Please note that reactivation of elements does only work
             * as long as the process that deactivated the elements before
             * is living.
                         */
            TrilinosDeactivationUtility()
            {
            }
            
            /**
             * Destructor.
             */
            virtual ~TrilinosDeactivationUtility()
            {
            }
            
            /**
             * Initializes all elements before performing any calculation.
             * This is done even for those elements that are deactivated
             * in the beginning of the simulation
             */
            void Initialize( ModelPart& model_part )
            {
                std::cout << "initializing deactivation utility" << std::endl;
                //initializing elements
                for ( ElementsArrayType::ptr_iterator it=model_part.Elements().ptr_begin();
                      it!=model_part.Elements().ptr_end(); ++it)
                {
                    if( (*it)->GetValue(ACTIVATION_LEVEL) < 0 ) (*it)->SetValue(IS_INACTIVE, true);
                    else (*it)->SetValue(IS_INACTIVE, false);
                    (*it)->Initialize();
                }
                for ( ConditionsArrayType::ptr_iterator it=model_part.Conditions().ptr_begin();
                      it != model_part.Conditions().ptr_end(); ++it)
                {
                    if( (*it)->GetValue(IS_CONTACT_MASTER) || (*it)->GetValue(IS_CONTACT_SLAVE) )
                    {
                        (*it)->SetValue(ACTIVATION_LEVEL, 0 );
                    }
                    if( (*it)->GetValue(ACTIVATION_LEVEL) < 0 ) (*it)->SetValue(IS_INACTIVE, true);
                    else (*it)->SetValue(IS_INACTIVE, false);
                    (*it)->Initialize();
                }
                std::cout << "deactivation utility initialized" << std::endl;
            }
            
            /**
             * Deactivates all elements and conditions marked with an
             * activation level in range (from_level, to_level).
             * Deactivated entities are stored in an intermediate
             * container and can be restored by calling Reactivate() or
             * ReactivateAll()
             */
            void Deactivate( ModelPart& model_part, int from_level, int to_level )
            {
                KRATOS_TRY;
                
                //first step: reactivate all elements and conditions
                ReactivateAll( model_part );
                //second step: deactivate elements and conditions to be deactivated currently
                // identify elements to be deactivated
                for ( ElementsArrayType::ptr_iterator it=model_part.Elements().ptr_begin();
                      it!=model_part.Elements().ptr_end(); ++it)
                {
                    //if element is to be deactivated
                    if( ( (*it)->GetValue( ACTIVATION_LEVEL ) >= from_level 
                            && (*it)->GetValue( ACTIVATION_LEVEL ) <= to_level
                            && (*it)->GetValue( ACTIVATION_LEVEL ) != 0 ) 
                            || ( (*it)->GetValue( ACTIVATION_LEVEL ) < 0 )
                      )
                    {
			(*it)->GetValue( IS_INACTIVE ) = true;
                    }
                }
                for( ConditionsArrayType::ptr_iterator it = model_part.Conditions().ptr_begin();
                     it != model_part.Conditions().ptr_end(); ++it )
                {
//                     std::cout << "condition: " << (*it)->Id() <<  "has activation level: " << (*it)->GetValue( ACTIVATION_LEVEL ) << std::endl;
                    if( ( (*it)->GetValue( ACTIVATION_LEVEL ) >= from_level
                            && (*it)->GetValue( ACTIVATION_LEVEL ) <= to_level
                            && (*it)->GetValue( ACTIVATION_LEVEL ) != 0 )
                            || ( (*it)->GetValue( ACTIVATION_LEVEL ) < 0 )
                      )
                    {
                        if( !( (*it)->GetValue( IS_CONTACT_MASTER ) || (*it)->GetValue( IS_CONTACT_SLAVE ) ) )
                        {
				(*it)->GetValue( IS_INACTIVE ) = true;
                        }
                    }
                }
				
                KRATOS_CATCH("")
            }
            
            /**
             * Reactivates all elements and conditions stored in the
             * intermediate containers
             */
            void ReactivateAll( ModelPart& model_part )
            {
                for ( ElementsArrayType::ptr_iterator it=model_part.Elements().ptr_begin();
                      it!=model_part.Elements().ptr_end(); ++it)
			(*it)->GetValue( IS_INACTIVE ) = false;
                for( ConditionsArrayType::ptr_iterator it = model_part.Conditions().ptr_begin();
                     it != model_part.Conditions().ptr_end(); ++it )
			(*it)->GetValue( IS_INACTIVE ) = false;
            }
            
            /**
             * reactivate all elements with activation level in range
             * (from_level, to_level)
             */
            void Reactivate( ModelPart& model_part, int from_level, int to_level )
            {
                KRATOS_TRY;
                for ( ElementsArrayType::ptr_iterator it=model_part.Elements().ptr_begin();
                      it!=model_part.Elements().ptr_end(); ++it)
                {
                    //if element is to be reactivated
                    if( (*it)->GetValue( ACTIVATION_LEVEL ) >= from_level 
                          && (*it)->GetValue( ACTIVATION_LEVEL ) <= to_level)
                    {
                        (*it)->SetValue(ACTIVATION_LEVEL, 0 );
                    }
                }
                for( ConditionsArrayType::ptr_iterator it = model_part.Conditions().ptr_begin();
                     it != model_part.Conditions().ptr_end(); ++it )
                {
                    if( ! ( (*it)->GetValue( IS_CONTACT_MASTER ) || (*it)->GetValue(IS_CONTACT_SLAVE) ) )
                    {
                        if( (*it)->GetValue( ACTIVATION_LEVEL ) >= from_level
                              && (*it)->GetValue(ACTIVATION_LEVEL) <= to_level)
                        {
                            (*it)->SetValue(ACTIVATION_LEVEL, 0);
                        }
                    }
                }
                KRATOS_CATCH("");
            }
            
            /**
             * reactivate all elements with activation level in range
             * (from_level, to_level)
             */
            void ReactivateStressFree( ModelPart& model_part, int from_level, int to_level )
            {
                KRATOS_TRY;
                
                for ( ElementsArrayType::ptr_iterator it=model_part.Elements().ptr_begin();
                      it!=model_part.Elements().ptr_end(); ++it)
                {
                    //if element is to be reactivated
                    if( (*it)->GetValue( ACTIVATION_LEVEL ) >= from_level 
                          && (*it)->GetValue( ACTIVATION_LEVEL ) <= to_level)
                    {
                        (*it)->Initialize();
                        (*it)->SetValue(ACTIVATION_LEVEL, 0 );
                    }
                }
                KRATOS_CATCH("");
            }
            
            /**
             * Turn back information as a string.
             */
            virtual std::string Info() const
            {
                return "TrilinosDeactivationUtility";
            }
            
            /**
             * Print information about this object.
             */
            virtual void PrintInfo(std::ostream& rOStream) const
            {
                rOStream << "TrilinosDeactivationUtility";
            }
            
            /**
             * Print object's data.
             */
            virtual void PrintData(std::ostream& rOStream) const
            {
            }
		
        private:
            
            /**
             * Containers for deactivated elements and conditions
             */
            PointerVector<Element> mDeactivatedElements;
            PointerVector<Condition> mDeactivatedConditions;
            
            /**
             * Assignment operator
             */
            //TrilinosDeactivationUtility& operator=(TrilinosDeactivationUtility const& rOther);
            
            /**
             * Copy constructor
             */
            //TrilinosDeactivationUtility(TrilinosDeactivationUtility const& rOther);
    
    };//class TrilinosDeactivationUtility

}  // namespace Kratos.

#endif // KRATOS_TRILINOS_DEACTIVATION_UTILITY_INCLUDED  defined 
