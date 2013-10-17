/*
==============================================================================
KratosULFApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
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
//   Date:                $Date: 2007-11-06 12:28:50 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_ULF_FRAC2D_SWIMMING_H_INCLUDED )
#define  KRATOS_ULF_FRAC2D_SWIMMING_H_INCLUDED



// System includes


// External includes


// Project includes
#include "custom_elements/ulf_frac2d.h"

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

/// Short class definition.
/** Detail class definition.
*/
class UlfFrac2DSwimming
    : public UlfFrac2D
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of UlfFrac2DSwimming
    KRATOS_CLASS_POINTER_DEFINITION(UlfFrac2DSwimming);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UlfFrac2DSwimming(IndexType NewId, GeometryType::Pointer pGeometry);
    UlfFrac2DSwimming(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~UlfFrac2DSwimming();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);


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
        return "UlfFrac2DSwimming #" ;
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info() << Id();
    }
    /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;


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
    //static boost::numeric::ublas::bounded_matrix<double,3,3> msMassFactors;
    //static boost::numeric::ublas::bounded_matrix<double,3,2> msDN_Dx;
    //static array_1d<double,3> msN; //dimension = number of nodes
    //static Matrix msDN_DX;
    //static Matrix msMassFactors;
    //static array_1d<double,2> ms_vel_gauss; //dimesion coincides with space dimension
    //static array_1d<double,3> ms_temp_vec_np; //dimension = number of nodes
    //static array_1d<double,3> ms_u_DN;


    //static boost::numeric::ublas::bounded_matrix<double,3,6> msB;
    //static boost::numeric::ublas::bounded_matrix<double,3,3> ms_constitutive_matrix;
    //static boost::numeric::ublas::bounded_matrix<double,3,6> ms_temp;
    //static array_1d<double,6> ms_temp_vec;


    ///@}
    ///@name Member Variables
    ///@{
    double mA0; //stores area at the beginning of time-step

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    UlfFrac2DSwimming() : UlfFrac2D()
    {
    }

    virtual void save(Serializer& rSerializer) const
    {

        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, UlfFrac2D);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, UlfFrac2D);
    }

    ///@}
    ///@name Private Operators
    ///@{
    void CalculateLumpedMass();
    void PressureStep(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
    void VelocityStep(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
    //inline void CalculateGeometryData(Matrix& msDN_DX, Vector& N, double& Area)
    inline void CalculateGeometryData(boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, array_1d<double,3>& N, double& Area);

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
    //UlfFrac2DSwimming& operator=(const UlfFrac2DSwimming& rOther);

    /// Copy constructor.
    //UlfFrac2DSwimming(const UlfFrac2DSwimming& rOther);


    ///@}

}; // Class UlfFrac2DSwimming

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    UlfFrac2D& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const UlfFrac2D& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_ULF_FRAC2D_SWIMMING_H_INCLUDED   defined 


