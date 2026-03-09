//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Hofer, https://github.com/PhiHo-eng
//                   Erich Wehrle, https://github.com/e-dub
//  based on original file from
//                   Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

#if !defined(KRATOS_TOPOLOGY_UPDATING_UTILITIES_H_INCLUDED)
#define  KRATOS_TOPOLOGY_UPDATING_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes


// Application includes
#include "topology_optimization_application.h"
#include "structure_response_function_utilities.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{


///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Solution utility that updates response values for next iteration.
/** Detail class definition.

 */

class TopologyUpdatingUtilities
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of TopologyUpdatingUtilities
    KRATOS_CLASS_POINTER_DEFINITION(TopologyUpdatingUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TopologyUpdatingUtilities( ModelPart& model_part )
    : mrModelPart(model_part)
    {
    }

    /// Destructor.
    virtual ~TopologyUpdatingUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------- UPDATE DENSITIES  -----------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------

    /// Finds the value of the X_PHYS (density) and updates it into the optimization problem
    void UpdateDensitiesUsingOCMethod( char update_type[], double volfrac, double greyscale , double OptItr , double qmax)
    {
        KRATOS_TRY;

        if ( strcmp( update_type , "oc_algorithm" ) == 0 ){
            BuiltinTimer timer;
            KRATOS_INFO("[TopOpt]") << "  Optimality Criterion Method (OC) chosen to solve the optimization problem" << std::endl;

            // Check if Grey Scale Filter should be used
            double q = 1;
            if (greyscale == 1)
            {
                if (OptItr < 10)
                    q = 1;
                else
                    q = std::min(qmax, 1.1*q);

                KRATOS_INFO("[TopOpt]") << "  Grey Scale Filter activated, q = " << q << std::endl;
            }
            else
            {
                KRATOS_INFO("[TopOpt]") << "  Grey Scale Filter deactivated, q = " << q << std::endl;
            }

            // Update Densities procedure
            double l1    = 0.0;
            double l2    = 1.0e12; ///10000000000.0
            const double move    = 0.2;
            double sum_X_Phys;
            int nele;
            double x_new = 0.0;
            double lmid = 0.0;
            double model_size;

            // Bisection algorithm to find Lagrange Multiplier so that volume constraint is satisfied (lmid)
            while ((l2-l1)/(l1+l2) > 0.001)
            {
                lmid = 0.5*(l2+l1);
                sum_X_Phys = 0.0;
                nele = 0;
                x_new = 0.0;
                model_size = 0.0;

                for( ModelPart::ElementIterator element_i = mrModelPart.ElementsBegin(); element_i!= mrModelPart.ElementsEnd(); element_i++ )
                {
                    const double x_old = element_i->GetValue(X_PHYS_OLD);
                    const int solid_void = element_i->GetValue(SOLID_VOID);
                    const double dcdx  = element_i->GetValue(DCDX);
                    const double dvdx  = element_i->GetValue(DVDX);
                    const double initial_element_size = element_i->GetValue(INITIAL_ELEMENT_SIZE);

                    // Update Density
                    // When q = 1, Grey Scale Filter is not activated, i.e., the results are in the classical OC update method
                    switch(solid_void)
                    {
                    // NORMAL elements
                    case 0:
                    {
                        x_new = std::max(0.0, std::max(x_old - move, std::min(1.0, pow(std::min(x_old + move, x_old * sqrt(-dcdx/dvdx/lmid)),q))));
                        break;
                    }
                    // ACTIVE elements (solid elements)
                    case 1:
                    {
                        x_new = 1;
                        break;
                    }
                    // PASSIVE elements (void elements)
                    case 2:
                    {
                        x_new = 0;
                        break;
                    }
                    default:
                    {
                        // If no element identification was found
                        KRATOS_INFO("[TopOpt]") << "This value for SOLID_VOID does not exist."<< std::endl;
                    }
                    }

                    // Update of the calculated X_PHYS for the next iteration
                    element_i->SetValue(X_PHYS, x_new);

                    // Updating additional quantities to determine the correct Lagrange Multiplier (lmid)
                    sum_X_Phys = sum_X_Phys + x_new*initial_element_size;
                    model_size += initial_element_size;
                    nele = nele + 1;
                }

                if( sum_X_Phys > (model_size*volfrac))
                    l1=lmid;
                else
                    l2=lmid;
            }

            // Printing of results
            KRATOS_INFO("[TopOpt]") << "  Updating of values performed               [ spent time =  " << timer.ElapsedSeconds() << " ] " << std::endl;
        } else {
            KRATOS_ERROR << "No valid optimization_algorithm selected for the simulation. Selected one: " << update_type << std::endl;
        }

        KRATOS_CATCH("");

    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "TopologyUpdatingUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TopologyUpdatingUtilities";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //TopologyUpdatingUtilities& operator=(TopologyUpdatingUtilities const& rOther);

    /// Copy constructor.
    //TopologyUpdatingUtilities(TopologyUpdatingUtilities const& rOther);


    ///@}

}; // Class TopologyUpdatingUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif	/* KRATOS_TOPOLOGY_UPDATING_UTILITIES_H_INCLUDED */
