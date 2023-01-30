//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//					 Philipp Hofer
//					 Erich Wehrle
//
// ==============================================================================

#if !defined(KRATOS_STRUCTURE_RESPONSE_FUNCTION_UTILITIES_H_INCLUDED)
#define  KRATOS_STRUCTURE_RESPONSE_FUNCTION_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes

// Application includes
#include "topology_optimization_application.h"
#include "utilities/builtin_timer.h"


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

/// Solution utility to compute structural analysis responses.
/** Detail class definition.

 */

class StructureResponseFunctionUtilities
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of StructureResponseFunctionUtilities
    KRATOS_CLASS_POINTER_DEFINITION(StructureResponseFunctionUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StructureResponseFunctionUtilities( ModelPart& model_part )
    : mr_structure_model_part(model_part)
    {
    }

    /// Destructor.
    virtual ~StructureResponseFunctionUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------- COMPUTE STRAIN ENERGY -------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------

    /// Computes the strain energy as the objective function of the optimization problem.
    double ComputeStrainEnergy()
    {
        KRATOS_TRY;

        BuiltinTimer timer;
        KRATOS_INFO("[TopOpt]") << "  Start calculating strain energy."<<std::endl;

        double Out = 0.0;
        double Global_Strain_Energy = 0.0;

        // Loop over all elements to calculate their local objective function and sum it into the global objective function (Global Strain Energy)
        for( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
                element_i++ )
        {

            element_i->Calculate(LOCAL_STRAIN_ENERGY, Out, mr_structure_model_part.GetProcessInfo());

            Global_Strain_Energy += element_i->GetValue(LOCAL_STRAIN_ENERGY);
            
        }

        KRATOS_INFO("[TopOpt]") <<  "  Strain energy calculated                  [ spent time =  " << timer.ElapsedSeconds() << " ] " << std::endl;

        // Return this obtained Global Strain Energy value as the objective function of the complete system
        return Global_Strain_Energy;

        KRATOS_CATCH("");
    }

    double ComputeVolumeFraction()
    {
        KRATOS_TRY;

        BuiltinTimer timer;
        KRATOS_INFO("[TopOpt]") <<"  Start calculating volume fraction."<<std::endl;

        double Global_Volume_Fraction = 0.0;
        double elemental_volume = 0.0;
        double design_variable = 0.0;
        double Total_volume = 0.0;


        // Loop over all elements to obtain their X_PHYS and know how many elements the model has
        for( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
                element_i++ )
        {

            elemental_volume = element_i->GetValue(INITIAL_ELEMENT_SIZE);
            design_variable = element_i->GetValue(X_PHYS);
            Global_Volume_Fraction += (elemental_volume*design_variable); //
            Total_volume += elemental_volume;
        }

        // Calculate and return the Global Volume Fraction by knowing how many elements the model has
        Global_Volume_Fraction = Global_Volume_Fraction/Total_volume;
        KRATOS_INFO("[TopOpt]") <<"  Global Volume Fraction: " << Global_Volume_Fraction << std::endl;
        KRATOS_INFO("[TopOpt]") <<"  Volume fraction calculated                [ spent time =  " << timer.ElapsedSeconds() << " ] " << std::endl;
        return Global_Volume_Fraction;

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
        return "StructureResponseFunctionUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "StructureResponseFunctionUtilities";
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

    ModelPart& mr_structure_model_part;

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
    //StructureResponseFunctionUtilities& operator=(StructureResponseFunctionUtilities const& rOther);

    /// Copy constructor.
    //StructureResponseFunctionUtilities(StructureResponseFunctionUtilities const& rOther);


    ///@}

}; // Class StructureResponseFunctionUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif	/* KRATOS_STRUCTURE_RESPONSE_FUNCTION_UTILITIES_H_INCLUDED */
