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
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2008-01-24 16:55:10 $
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
             * Constructor.
             * The constructor takes the current model_part as argument.
             * Please note that reactivation of elements does only work
             * as long as the process that deactivated the elements before
             * is living.
                         */
            DeactivationUtility()
            {
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
                std::cout << "initializing deactivation utility" << std::endl;
                //initializing elements
                for ( ElementsArrayType::ptr_iterator it=model_part.Elements().ptr_begin();
                      it!=model_part.Elements().ptr_end(); ++it)
                {
                    (*it)->Initialize();
                }
                mDeactivatedElements.clear();
                mDeactivatedConditions.clear();
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
                std::vector<int> elements_to_be_deleted;
                for ( ElementsArrayType::ptr_iterator it=model_part.Elements().ptr_begin();
                      it!=model_part.Elements().ptr_end(); ++it)
                {
                    //if element is to be deactivated
                    if( ( (*it)->GetValue( ACTIVATION_LEVEL ) > from_level 
                            && (*it)->GetValue( ACTIVATION_LEVEL ) < to_level) 
                            || ( (*it)->GetValue( ACTIVATION_LEVEL ) < 0 )
                      )
                    {
                        elements_to_be_deleted.push_back((*it)->Id());
                        mDeactivatedElements.push_back( *(it) );
  //                      mDeactivatedElements.push_back( *(it.base()) );
                        //deactivate associated conditions (unless they are contact master)
                        Condition Cond = model_part.Conditions()[(*it)->Id()];
                        if( Cond.GetValue( IS_CONTACT_MASTER ) == false )
                        {
                            mDeactivatedConditions.push_back( Cond );
                            model_part.RemoveCondition( Cond );
                        }
                    }
                }
                // deactivate identified elements
                for( std::vector<int>::iterator it=elements_to_be_deleted.begin();
                     it != elements_to_be_deleted.end(); it++ )
                {
                    model_part.RemoveElement( *it );
                }
                //re-sort the elements
                model_part.Elements().Sort();
                model_part.Conditions().Sort();
				
                KRATOS_CATCH("")
            }
            
            /**
             * Reactivates all elements and conditions stored in the
             * intermediate containers
             */
            void ReactivateAll( ModelPart& model_part )
            {
                for( ElementsArrayType::iterator it = mDeactivatedElements.begin();
                     it != mDeactivatedElements.end(); it++ )
                {
                    model_part.AddElement( *(it.base()));
                }
                for( ConditionsArrayType::iterator jt = mDeactivatedConditions.begin();
                     jt != mDeactivatedConditions.end(); jt++ )
                {
                    model_part.AddCondition( *(jt.base()));
                }
                model_part.Elements().Sort();
                model_part.Conditions().Sort();

                mDeactivatedElements.clear();
                mDeactivatedConditions.clear();
            }
            
            /**
             * reactivate all elements with activation level in range
             * (from_level, to_level)
             */
            void Reactivate( ModelPart& model_part, int from_level, int to_level )
            {
                KRATOS_TRY;
                ElementsArrayType elements_to_be_restored;
                ConditionsArrayType conditions_to_be_restored;
                for( ElementsArrayType::iterator it = mDeactivatedElements.begin();
                     it != mDeactivatedElements.end(); it++ )
                {
					//if element is to be reactivated
                    if( (*it).GetValue( ACTIVATION_LEVEL ) > from_level 
                          && (*it).GetValue( ACTIVATION_LEVEL ) < to_level)
                    {
                        (*it).SetValue(ACTIVATION_LEVEL, 0 );
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
                
                /*for(ModelPart::NodeIterator i = model_part.NodesBegin() ; 
                i != model_part.NodesEnd() ; ++i)
                {
                (i)->X() = (i)->X0() + i->GetSolutionStepValue(DISPLACEMENT_NULL_X);
                (i)->Y() = (i)->Y0() + i->GetSolutionStepValue(DISPLACEMENT_NULL_Y);
                (i)->Z() = (i)->Z0() + i->GetSolutionStepValue(DISPLACEMENT_NULL_Z);

                outfile<<"NODE["<<(i)->Id()<<"]";
                outfile<<"["<<(i)->X()<<","<<(i)->Y()<<","<<(i)->Z()<<"]\n";
                }*/

                ElementsArrayType elements_to_be_restored;
                ConditionsArrayType conditions_to_be_restored;

                for( ElementsArrayType::iterator it = mDeactivatedElements.begin();
                     it != mDeactivatedElements.end(); it++ )
                {
                    //if element is to be reactivated
                    if( (*it).GetValue( ACTIVATION_LEVEL ) > from_level 
                          && (*it).GetValue( ACTIVATION_LEVEL ) < to_level)
                    {
                        (*it).Initialize();
                        (*it).SetValue(ACTIVATION_LEVEL, 0 );
                    }
                }
                //outfile.close();
                KRATOS_CATCH("");

            }
            
            /**
             * reactivate all elements with activation level in range
             * (from_level, to_level)
             */
/*            void ReactivateStressFree( ModelPart& model_part, int from_level, int to_level, const std::string& rDatafilename )
            {
                KRATOS_TRY;

                std::string nodesfile= rDatafilename+".deactivation.restart";
                std::ofstream outfile;
                outfile.open((char*)nodesfile.c_str());
                outfile.precision(12);

                for(ModelPart::NodeIterator i = model_part.NodesBegin() ; 
                    i != model_part.NodesEnd() ; ++i)
                {
                    (i)->X() = (i)->X0() + i->GetSolutionStepValue(DISPLACEMENT_NULL_X);
                    (i)->Y() = (i)->Y0() + i->GetSolutionStepValue(DISPLACEMENT_NULL_Y);
                    (i)->Z() = (i)->Z0() + i->GetSolutionStepValue(DISPLACEMENT_NULL_Z);

                    outfile<<"NODE["<<(i)->Id()<<"]";
                    outfile<<"["<<(i)->X()<<","<<(i)->Y()<<","<<(i)->Z()<<"]\n";
                }

                ElementsArrayType elements_to_be_restored;
                ConditionsArrayType conditions_to_be_restored;

                for( ElementsArrayType::iterator it = mDeactivatedElements.begin();
                     it != mDeactivatedElements.end(); it++ )
                {
                    //if element is to be reactivated
                    if( (*it).GetValue( ACTIVATION_LEVEL ) > from_level 
                          && (*it).GetValue( ACTIVATION_LEVEL ) < to_level)
                    {
                        (*it).Initialize();
                        (*it).SetValue(ACTIVATION_LEVEL, 0 );
                    }
                }

                outfile.close();
                KRATOS_CATCH("");

            }*/
            
            void ReactivateStressFreeFromFile( ModelPart& model_part, int from_level, int to_level, const std::string& rDatafilename )
            {
                KRATOS_TRY;
                unsigned int id;
                char c;
                std::string line;
                std::string out;
                std::string nodesfile= rDatafilename+".deactivation.restart";
                std::ifstream infile;
                infile.open((char*)nodesfile.c_str());
                infile.precision(20);
                infile.seekg(0,std::ios::beg);

                for(ModelPart::NodeIterator i = model_part.NodesBegin() ; 
                    i != model_part.NodesEnd() ; ++i)
                {
                    bool node_found= false;
                    std::streampos sp = infile.tellg(); 

                    for(int try_number=0; try_number<2; try_number++)//Loop searches first from last pointer position and thereafter from begining of file
                    {
                        while (!infile.eof()|| (try_number=1 && (infile.tellg()> sp) )) 
                        {
                            out="";
                            while(c!='\n')
                            {
                                infile.get(c);
                                if(c=='[') break;
                                out+=c;
                            }
                            if(out=="NODE")
                            {
                                out="";
                                while(c!='\n')
                                {
                                    infile.get(c);
                                    if(c==']') break;
                                    out+=c;
                                }
                                id= std::atoi((char*)out.c_str());
                                infile.get(c);
                                if(id== (i)->Id())
                                {
                                    node_found= true;
                                    out="";
                                    while(c!='\n')
                                    {
                                        infile.get(c);
                                        if(c== ',') break;
                                        out+=c;
                                    }
                                    double value= std::atof((char*)out.c_str());
									if(isnan(value))
                                    {
                                        std::cout<<"##### NaN FOUND ("<<i->Id()<<") ####"<<std::endl;
                                        continue;
                                    }

                                    (i)->X() = value;
	
                                    out="";
                                    while(c!='\n')
                                    {
                                        infile.get(c);
                                        if(c== ',') break;
                                        out+=c;
                                    }
                                    value= std::atof((char*)out.c_str());
									if(isnan(value))
                                    {
                                        std::cout<<"##### NaN FOUND ("<<i->Id()<<") ####"<<std::endl;
                                        continue;
                                    }
                                    (i)->Y() = value;
	
                                    out="";
                                    while(c!='\n')
                                    {
                                        infile.get(c);
                                        if(c== ']') break;
                                        out+=c;
                                    }
                                    value= std::atof((char*)out.c_str());
                                    if(isnan(value))
                                    {
                                        std::cout<<"##### NaN FOUND ("<<i->Id()<<") ####"<<std::endl;
                                        continue;
                                    }
	
                                    (i)->Z() = value;
	
                                    std::getline(infile, line);
                                    break;
                                }
                                else
                                {
                                    std::getline(infile, line);
                                    continue;
                                }
                            }
                        }
                        if(node_found)
                        {
                            break;
                        }
                        infile.seekg(0,std::ios::beg);
                    }
                    if(!node_found)
                        std::cout<<"##### NO CORRESPONDING NODE ("<<i->Id()<<") IN INPUT FILE ####"<<std::endl;
                }

                for( ElementsArrayType::iterator it = mDeactivatedElements.begin();
                     it != mDeactivatedElements.end(); it++ )
                {
					//if element is to be reactivated
                    if( (*it).GetValue( ACTIVATION_LEVEL ) > from_level 
                          && (*it).GetValue( ACTIVATION_LEVEL ) < to_level)
                    {
                        (*it).Initialize();
                        (*it).SetValue(ACTIVATION_LEVEL, 0 );
                    }
                }

                for( ElementsArrayType::iterator it = mDeactivatedElements.begin();
                     it != mDeactivatedElements.end(); it++ )
                {
					//if element is to be reactivated
                    if( (*it).GetValue( ACTIVATION_LEVEL ) > from_level 
                          && (*it).GetValue( ACTIVATION_LEVEL ) < to_level)
                    {
                        for(unsigned int node=0; node<(*it).GetGeometry().size(); node++)
                        {
                            (*it).GetGeometry()[node].X() = (*it).GetGeometry()[node].X0() +
                                    (*it).GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT_X);
                            (*it).GetGeometry()[node].Y() = (*it).GetGeometry()[node].Y0() +
                                    (*it).GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT_Y);
                            (*it).GetGeometry()[node].Z() = (*it).GetGeometry()[node].Z0() +
                                    (*it).GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT_Z);
                        }
                    }
                }
                KRATOS_CATCH("");

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
