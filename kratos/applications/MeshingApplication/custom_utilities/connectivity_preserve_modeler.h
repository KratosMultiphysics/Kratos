/*
==============================================================================
KratosPFEMApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


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
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-11-19 15:38:01 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_GENERATE_MODEL_PART_MODELER_INCLUDED )
#define  KRATOS_GENERATE_MODEL_PART_MODELER_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "modeler/modeler.h"

#include "meshing_application.h"
namespace Kratos {

    class ConnectivityPreserveModeler : public Modeler
    {
    public:

        ConnectivityPreserveModeler() {
        }

        /// Destructor.

        virtual ~ConnectivityPreserveModeler() {
        }

        //**********************************************************************************************
        //**********************************************************************************************
        ///This function fills the @param DestinationModelPart using the data obtained from @param  OriginModelPart
        ///the elements and conditions of the DestinationModelPart part use the same connectivity (and id) as the
        ///OriginModelPart but their type is determined by @param rReferenceElement and @param rReferenceBoundaryCondition

        void GenerateModelPart(
                ModelPart& OriginModelPart,
                ModelPart& DestinationModelPart,
                Element const& rReferenceElement,
                Condition const& rReferenceBoundaryCondition
                )
        {
            KRATOS_TRY;

            //assigning ProcessInfo
            DestinationModelPart.pGetProcessInfo() = OriginModelPart.pGetProcessInfo();

            //assigning Properties
            DestinationModelPart.pProperties() = OriginModelPart.pProperties();

            //assigning the nodes to the new model part
            DestinationModelPart.Nodes().clear();
            DestinationModelPart.Nodes() = OriginModelPart.Nodes();

            //generating the elements
            for (ModelPart::ElementsContainerType::iterator iii = OriginModelPart.ElementsBegin(); iii != OriginModelPart.ElementsEnd(); iii++) {
                Properties::Pointer properties = iii->pGetProperties();
                Element::Pointer p_element = rReferenceElement.Create(iii->Id(), iii->GetGeometry(), properties);
                DestinationModelPart.Elements().push_back(p_element);
            }
            std::cout << "Elements are generated" << std::endl;

            //generating the conditions
            for (ModelPart::ConditionsContainerType::iterator iii = OriginModelPart.ConditionsBegin(); iii != OriginModelPart.ConditionsEnd(); iii++) {
                Properties::Pointer properties = iii->pGetProperties();

                Condition::Pointer p_condition = rReferenceBoundaryCondition.Create(iii->Id(), iii->GetGeometry(), properties);
                DestinationModelPart.Conditions().push_back(p_condition);
            }
            std::cout << "Conditions are generated" << std::endl;

            KRATOS_CATCH("");
        }


        //**********************************************************************************************
        //**********************************************************************************************
    private:

        ConnectivityPreserveModeler & operator=(ConnectivityPreserveModeler const& rOther);

        /// Copy constructor.
        ConnectivityPreserveModeler(ConnectivityPreserveModeler const& rOther);


    };

} // namespace Kratos.

#endif //KRATOS_GENERATE_MODEL_PART_MODELER_INCLUDED  defined


