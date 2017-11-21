//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#if !defined (KRATOS_EULER_LAW_2D_H_INCLUDED)
#define  KRATOS_EULER_LAW_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
    /**
     * Defines an Euler fluid constitutive law
     * This material law is defined by the parameters:
     * 1) DYNAMIC_VISCOSITY
     */
    class Euler2DLaw : public ConstitutiveLaw
    {
    public:
        /**
         * Type Definitions
         */
        typedef ProcessInfo      ProcessInfoType;
        typedef ConstitutiveLaw         BaseType;
        typedef std::size_t             SizeType;

        /**
         * Counted pointer of Euler2DLaw
         */
        KRATOS_CLASS_POINTER_DEFINITION(Euler2DLaw);

        /**
         * Life Cycle
         */

        /**
         * Default constructor.
         */
        Euler2DLaw();

        /**
         * Clone function (has to be implemented by any derived class)
         * @return a pointer to a new instance of this constitutive law
         */
        ConstitutiveLaw::Pointer Clone() const override;

        /**
         * Copy constructor.
         */
        Euler2DLaw (const Euler2DLaw& rOther);


        /**
         * Assignment operator.
         */

        //Newtonian3DLaw& operator=(const Newtonian3DLaw& rOther);


        /**
         * Destructor.
         */
        ~Euler2DLaw() override;

        /**
         * Operators
         */

        /**
         * Operations needed by the base class:
         */
        void CalculateMaterialResponseCauchy (Parameters& rValues) override;

        /**
         * This function is designed to be called once to check compatibility with element
         * @param rFeatures
         */
        void GetLawFeatures(Features& rFeatures) override;


        /**
         * This function is designed to be called once to perform all the checks needed
         * on the input provided. Checks can be "expensive" as the function is designed
         * to catch user's errors.
         * @param rMaterialProperties
         * @param rElementGeometry
         * @param rCurrentProcessInfo
         * @return
         */
        int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

        /**
         * Input and output
         */
        /**
         * Turn back information as a string.
         */
        //virtual String Info() const;
        /**
         * Print information about this object.
         */
        //virtual void PrintInfo(std::ostream& rOStream) const;
        /**
         * Print object's data.
         */
        //virtual void PrintData(std::ostream& rOStream) const;

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

        void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
        }

        void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
        }

    }; // Class Euler2DLaw
}  // namespace Kratos.
#endif // KRATOS_EULER_LAW_2D_H_INCLUDED  defined 
