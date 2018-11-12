//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CASM_CEM_HARDENING_LAW_H_INCLUDED )
#define  KRATOS_CASM_CEM_HARDENING_LAW_H_INCLUDED



// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_hardening_laws/hardening_law.hpp"

namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

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

/// Short class definition.
/** Detail class definition.
*/
class CasmCemHardeningLaw 
        : public HardeningLaw 
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CasmCemHardeningLaw
    KRATOS_CLASS_POINTER_DEFINITION( CasmCemHardeningLaw );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CasmCemHardeningLaw();


    /// Copy constructor.
    CasmCemHardeningLaw(CasmCemHardeningLaw const& rOther);

    /// Assignment operator.
    CasmCemHardeningLaw& operator=(CasmCemHardeningLaw const& rOther);

    /// Destructor.
    ~CasmCemHardeningLaw();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    void CalculateHardening(PlasticVariablesType& rPlasticVariables, const double& rDeltaAlpha, const double& rDeltaBeta);
    
    Vector& CalculateHardening(Vector& rHardening, const double& rAlpha, const double& rBeta, const double& rAlphaCum, const double& rBetaCum, const double rTemperature = 0);
    
    double& CalculateHardening(double& rHardening, const double& rAlpha, const double& rBeta, const double& rAlphaCum, const double& rBetaCum, const double rTemperature = 0);
    
    double& CalculateBonding(double& rHardening, const double& rAlpha, const double& rBeta, const double& rAlphaCum, const double& rBetaCum, const double rTemperature = 0);


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
    // std::string Info() const;

    // /// Print information about this object.
    // void PrintInfo(std::ostream& rOStream) const;

    // /// Print object's data.
    // void PrintData(std::ostream& rOStream) const;


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
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class CasmCemHardeningLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


// /// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   LinearIsotropicHardeningLaw& rThis);

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const LinearIsotropicHardeningLaw& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CASM_CEM_HARDENING_LAW_H_INCLUDED defined

