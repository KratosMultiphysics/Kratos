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


#if !defined(KRATOS_DEACTIVATION_UTILITY_INCLUDED )
#define  KRATOS_DEACTIVATION_UTILITY_INCLUDED

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
#include "includes/kratos_flags.h"
#include "utilities/openmp_utils.h"
#include "structural_application.h"

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
class DeactivationUtility
{
public:

    typedef PointerVectorSet<Element> ElementsArrayType;
    typedef PointerVectorSet<Condition> ConditionsArrayType;

    /**
     * class pointer definition
     */
    KRATOS_CLASS_POINTER_DEFINITION( DeactivationUtility );

    /**
     * Default Constructor.
     * Please note that reactivation of elements does only work
     * as long as the process that deactivated the elements before
     * is living.
     */
    DeactivationUtility()
    {
        mEchoLevel = 1;
    }

    /**
     * Constructor with echo level
     */
    DeactivationUtility(int EchoLevel)
    {
        mEchoLevel = EchoLevel;
    }

    /**
     * Destructor.
     */
    virtual ~DeactivationUtility()
    {
    }

    /**
     * Initializes all elements before performing any calculation.
     * This is done even for those elements that are deactivated
     * in the beginning of the simulation
     */
    void Initialize( ModelPart& model_part )
    {
        #ifndef _OPENMP
        if(mEchoLevel > 0)
            std::cout << "initializing deactivation utility" << std::endl;

        //initializing elements
        for ( ElementsArrayType::ptr_iterator it=model_part.Elements().ptr_begin();
                it!=model_part.Elements().ptr_end(); ++it )
        {
            if( (*it)->GetValue(ACTIVATION_LEVEL) < 0 )
            {
                (*it)->Set(ACTIVE, false);
            }
            else
            {
                (*it)->Set(ACTIVE, true);
            }
            (*it)->Initialize();
        }
        for ( ConditionsArrayType::ptr_iterator it=model_part.Conditions().ptr_begin();
                it != model_part.Conditions().ptr_end(); ++it )
        {
            if( (*it)->GetValue(IS_CONTACT_MASTER) || (*it)->GetValue(IS_CONTACT_SLAVE) )
            {
                (*it)->SetValue(ACTIVATION_LEVEL, 0 );
            }
            if( (*it)->GetValue(ACTIVATION_LEVEL) < 0 )
            {
                (*it)->Set(ACTIVE, false);
            }
            else
            {
                (*it)->Set(ACTIVE, true);
            }
            (*it)->Initialize();
        }
        #else
        if(mEchoLevel > 0)
            std::cout << "multithreaded initializing deactivation utility" << std::endl;

        int number_of_threads = omp_get_max_threads();

        vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, model_part.Elements().size(), element_partition);

        vector<unsigned int> condition_partition;
        OpenMPUtils::CreatePartition(number_of_threads, model_part.Conditions().size(), condition_partition);

        #pragma omp parallel for
        for(int k = 0; k < number_of_threads; ++k)
        {
            ElementsArrayType::ptr_iterator it_begin = model_part.Elements().ptr_begin() + element_partition[k];
            ElementsArrayType::ptr_iterator it_end = model_part.Elements().ptr_begin() + element_partition[k + 1];

            for (ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it)
            {
                if( (*it)->GetValue(ACTIVATION_LEVEL) < 0 )
                {
                    (*it)->Set(ACTIVE, false);
                }
                else
                {
                    (*it)->Set(ACTIVE, true);
                }
//                std::cout << "element of type " << typeid(*(*it)).name() << " Id = " << (*it)->Id() << " is going to be initialized" << std::endl;
                (*it)->Initialize();
            }
        }

        #pragma omp parallel for
        for(int k = 0; k < number_of_threads; ++k)
        {
            ConditionsArrayType::ptr_iterator it_begin = model_part.Conditions().ptr_begin() + condition_partition[k];
            ConditionsArrayType::ptr_iterator it_end = model_part.Conditions().ptr_begin() + condition_partition[k + 1];

            for (ConditionsArrayType::ptr_iterator it = it_begin; it != it_end; ++it)
            {
                if( (*it)->GetValue(IS_CONTACT_MASTER) || (*it)->GetValue(IS_CONTACT_SLAVE) )
                {
                    (*it)->SetValue(ACTIVATION_LEVEL, 0 );
                }
                if( (*it)->GetValue(ACTIVATION_LEVEL) < 0 )
                {
                    (*it)->Set(ACTIVE, false);
                }
                else
                {
                    (*it)->Set(ACTIVE, true);
                }
//                std::cout << "condition of type " << typeid(*(*it)).name() << " Id = " << (*it)->Id() << " is going to be initialized" << std::endl;
                (*it)->Initialize();
            }
        }
        #endif

        if(mEchoLevel > 0)
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
                (*it)->Set( ACTIVE, false );
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
                    (*it)->Set( ACTIVE, false );
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
        {
            (*it)->Set( ACTIVE, true );
        }
        for( ConditionsArrayType::ptr_iterator it = model_part.Conditions().ptr_begin();
                it != model_part.Conditions().ptr_end(); ++it )
        {
            (*it)->Set( ACTIVE, true );
        }
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

    /// Get name from an abject
    template<class Object>
    std::string GetName(const Object& o)
    {
        return std::string(typeid(o).name());
    }

    /**
     * Turn back information as a string.
     */
    virtual std::string Info() const
    {
        return "DeactivationUtility";
    }

    /**
     * Print information about this object.
     */
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DeactivationUtility";
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
    int mEchoLevel;

    /**
     * Assignment operator
     */
    //DeactivationUtility& operator=(DeactivationUtility const& rOther);

    /**
     * Copy constructor
     */
    //DeactivationUtility(DeactivationUtility const& rOther);

};//class DeactivationUtility

}  // namespace Kratos.

#endif // KRATOS_DEACTIVATION_UTILITY_INCLUDED  defined 
