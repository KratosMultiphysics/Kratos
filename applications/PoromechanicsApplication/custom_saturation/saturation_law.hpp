//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Danilo Cavalcanti
//

#pragma once

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/node.h"
#include "includes/properties.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos {

class KRATOS_API(POROMECHANICS_APPLICATION) SaturationLaw 
{

public:

    typedef ProcessInfo ProcessInfoType;
    typedef std::size_t SizeType;
    typedef Geometry<Node > GeometryType;

    ///------------------------------------------------------------------------------------------------

    KRATOS_CLASS_POINTER_DEFINITION(SaturationLaw);

    ///------------------------------------------------------------------------------------------------

    struct Parameters
    {

        KRATOS_CLASS_POINTER_DEFINITION(Parameters);

        ///------------------------------------------------------------------------------------------------

    private:

        double*             mpSl;
        double*             mpdSldPc;
        double*             mpkrl;
        double*             mpkrg;

        const Vector*       mpShapeFunctionsValues;
        const Matrix*       mpShapeFunctionsDerivatives;

        const ProcessInfo*  mpCurrentProcessInfo;
        const Properties*   mpMaterialProperties;
        const GeometryType* mpElementGeometry;

        ///------------------------------------------------------------------------------------------------

    public:

        Parameters()
        {
            mpSl=NULL;
            mpdSldPc=NULL;
            mpkrl=NULL;
            mpkrg=NULL;
            mpShapeFunctionsValues=NULL;
            mpShapeFunctionsDerivatives=NULL;
            mpCurrentProcessInfo=NULL;
            mpMaterialProperties=NULL;
            mpElementGeometry=NULL;
        };

        Parameters(const GeometryType& rElementGeometry,
                const Properties& rMaterialProperties,
                const ProcessInfo& rCurrentProcessInfo):
            mpCurrentProcessInfo(&rCurrentProcessInfo),
            mpMaterialProperties(&rMaterialProperties),
            mpElementGeometry(&rElementGeometry)
        {
            mpSl=NULL;
            mpdSldPc=NULL;
            mpkrl=NULL;
            mpkrg=NULL;
            mpShapeFunctionsValues=NULL;
            mpShapeFunctionsDerivatives=NULL;
        };

        Parameters(const Parameters & rNewParameters):
            mpSl(rNewParameters.mpSl),
            mpdSldPc(rNewParameters.mpdSldPc),
            mpkrl(rNewParameters.mpkrl),
            mpkrg(rNewParameters.mpkrg),
            mpShapeFunctionsValues(rNewParameters.mpShapeFunctionsValues),
            mpShapeFunctionsDerivatives(rNewParameters.mpShapeFunctionsDerivatives),
            mpCurrentProcessInfo(rNewParameters.mpCurrentProcessInfo),
            mpMaterialProperties(rNewParameters.mpMaterialProperties),
            mpElementGeometry(rNewParameters.mpElementGeometry)
        {

        };

        ~Parameters() = default;

        ///------------------------------------------------------------------------------------------------

        bool CheckAllParameters ()
        {
            if(CheckHydraulicVariables() &&  CheckShapeFunctions() && CheckInfoMaterialGeometry ())
                return 1;
            else
                return 0;
        }

        bool CheckHydraulicVariables ()
        {
            if(!mpSl)
                KRATOS_ERROR << "Sl NOT SET" << std::endl;

            if(!mpdSldPc)
                KRATOS_ERROR << "dSldPc NOT SET" << std::endl;

            if(!mpkrl)
                KRATOS_ERROR << "krl NOT SET" << std::endl;

            if(!mpkrg)
                KRATOS_ERROR << "krg NOT SET" << std::endl;

            return 1;
        }

        bool CheckShapeFunctions ()
        {
            if(!mpShapeFunctionsValues)
                KRATOS_ERROR << "ShapeFunctionsValues NOT SET" << std::endl;

            if(!mpShapeFunctionsDerivatives)
                KRATOS_ERROR << "ShapeFunctionsDerivatives NOT SET" << std::endl;

            return 1;
        }

        bool CheckInfoMaterialGeometry ()
        {
            if(!mpCurrentProcessInfo)
                KRATOS_ERROR << "rCurrentProcessInfo NOT SET" << std::endl;

            if(!mpMaterialProperties)
                KRATOS_ERROR << "MaterialProperties NOT SET" << std::endl;

            if(!mpElementGeometry)
                KRATOS_ERROR << "ElementGeometry NOT SET" << std::endl;

            return 1;
        }

        ///------------------------------------------------------------------------------------------------

        void SetSl (double& rSl)         {mpSl=&rSl;};
        void SetdSldPc (double& rdSldPc) {mpdSldPc=&rdSldPc;};
        void Setkrl (double& rkrl)       {mpkrl=&rkrl;};
        void Setkrg (double& rkrg)       {mpkrg=&rkrg;};

        void SetShapeFunctionsValues (const Vector& rShapeFunctionsValues)           {mpShapeFunctionsValues=&rShapeFunctionsValues;};
        void SetShapeFunctionsDerivatives (const Matrix& rShapeFunctionsDerivatives) {mpShapeFunctionsDerivatives=&rShapeFunctionsDerivatives;};

        void SetProcessInfo (const ProcessInfo& rProcessInfo)               {mpCurrentProcessInfo =&rProcessInfo;};
        void SetMaterialProperties (const Properties&  rMaterialProperties) {mpMaterialProperties =&rMaterialProperties;};
        void SetElementGeometry (const GeometryType& rElementGeometry)      {mpElementGeometry =&rElementGeometry;};

        ///------------------------------------------------------------------------------------------------

        double& GetSl()
        {
            KRATOS_DEBUG_ERROR_IF_NOT(IsSetSl()) << "Sl is not set!" << std::endl;
            return *mpSl;
        }
        double& GetdSldPc()
        {
            KRATOS_DEBUG_ERROR_IF_NOT(IsSetdSldPc()) << "dSldPc is not set!" << std::endl;
            return *mpdSldPc;
        }
        double& Getkrl()
        {
            KRATOS_DEBUG_ERROR_IF_NOT(IsSetkrl()) << "krl is not set!" << std::endl;
            return *mpkrl;
        }
        double& Getkrg()
        {
            KRATOS_DEBUG_ERROR_IF_NOT(IsSetkrg()) << "krg is not set!" << std::endl;
            return *mpkrg;
        }

        const Vector& GetShapeFunctionsValues()
        {
            KRATOS_DEBUG_ERROR_IF_NOT(IsSetShapeFunctionsValues()) << "ShapeFunctionsValues is not set!" << std::endl;
            return *mpShapeFunctionsValues;
        }
        const Matrix& GetShapeFunctionsDerivatives()
        {
            KRATOS_DEBUG_ERROR_IF_NOT(IsSetShapeFunctionsDerivatives()) << "ShapeFunctionsDerivatives is not set!" << std::endl;
            return *mpShapeFunctionsDerivatives;
        }

        const ProcessInfo& GetProcessInfo()
        {
            KRATOS_DEBUG_ERROR_IF_NOT(IsSetProcessInfo()) << "ProcessInfo is not set!" << std::endl;
            return *mpCurrentProcessInfo;
        }
        const Properties& GetMaterialProperties()
        {
            KRATOS_DEBUG_ERROR_IF_NOT(IsSetMaterialProperties()) << "MaterialProperties is not set!" << std::endl;
            return *mpMaterialProperties;
        }
        const GeometryType& GetElementGeometry()
        {
            KRATOS_DEBUG_ERROR_IF_NOT(IsSetElementGeometry()) << "ElementGeometry is not set!" << std::endl;
            return *mpElementGeometry;
        }

        ///------------------------------------------------------------------------------------------------

        double& GetSl (double& rSl)         {rSl=*mpSl; return rSl;};
        double& GetdSldPc (double& rdSldPc) {rdSldPc=*mpdSldPc; return rdSldPc;};
        double& Getkrl (double& rkrl)       {rkrl=*mpkrl; return rkrl;};
        double& Getkrg (double& rkrg)       {rkrg=*mpkrg; return rkrg;};

        ///------------------------------------------------------------------------------------------------

        bool IsSetSl ()     {return (mpSl != NULL);};
        bool IsSetdSldPc () {return (mpdSldPc != NULL);};
        bool IsSetkrl ()    {return (mpkrl != NULL);};
        bool IsSetkrg ()    {return (mpkrg != NULL);};

        bool IsSetShapeFunctionsValues ()      {return (mpShapeFunctionsValues != NULL);};
        bool IsSetShapeFunctionsDerivatives () {return (mpShapeFunctionsDerivatives != NULL);};

        bool IsSetProcessInfo ()        {return (mpCurrentProcessInfo != NULL);};
        bool IsSetMaterialProperties () {return (mpMaterialProperties != NULL);};
        bool IsSetElementGeometry ()    {return (mpElementGeometry != NULL);};

    };// struct Parameters end

    ///------------------------------------------------------------------------------------------------

    SaturationLaw() = default;

    SaturationLaw(const SaturationLaw& rOther)
    {
    }

    virtual ~SaturationLaw() = default;

    virtual SaturationLaw::Pointer Clone() const;

    ///------------------------------------------------------------------------------------------------
    
    // virtual SizeType WorkingSpaceDimension();

    ///------------------------------------------------------------------------------------------------

    virtual bool Has(const Variable<bool>& rThisVariable);
    virtual bool Has(const Variable<int>& rThisVariable);
    virtual bool Has(const Variable<double>& rThisVariable);
    virtual bool Has(const Variable<Vector>& rThisVariable);
    virtual bool Has(const Variable<Matrix>& rThisVariable);
    virtual bool Has(const Variable<array_1d<double, 3 > >& rThisVariable);
    virtual bool Has(const Variable<array_1d<double, 6 > >& rThisVariable);

    ///------------------------------------------------------------------------------------------------

    virtual bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue);
    virtual int& GetValue(const Variable<int>& rThisVariable, int& rValue);
    virtual double& GetValue(const Variable<double>& rThisVariable, double& rValue);
    virtual Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue);
    virtual Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue);
    virtual array_1d<double, 3>& GetValue(const Variable<array_1d<double, 3>>& rThisVariable, array_1d<double,3>& rValue);
    virtual array_1d<double, 6>& GetValue(const Variable<array_1d<double, 6>>& rThisVariable, array_1d<double, 6>& rValue);
    
    ///------------------------------------------------------------------------------------------------

    virtual void SetValue(const Variable<bool>& rVariable,
                          const bool& Value,
                          const ProcessInfo& rCurrentProcessInfo);
    virtual void SetValue(const Variable<int>& rVariable,
                          const int& Value,
                          const ProcessInfo& rCurrentProcessInfo);
    virtual void SetValue(const Variable<double>& rVariable,
                          const double& rValue,
                          const ProcessInfo& rCurrentProcessInfo);
    virtual void SetValue(const Variable<Vector >& rVariable,
                          const Vector& rValue,
                          const ProcessInfo& rCurrentProcessInfo);
    virtual void SetValue(const Variable<Matrix >& rVariable,
                          const Matrix& rValue,
                          const ProcessInfo& rCurrentProcessInfo);
    virtual void SetValue(const Variable<array_1d<double, 3>>& rVariable,
                          const array_1d<double, 3>& rValue,
                          const ProcessInfo& rCurrentProcessInfo);
    virtual void SetValue(const Variable<array_1d<double, 6>>& rVariable,
                          const array_1d<double, 6>& rValue,
                          const ProcessInfo& rCurrentProcessInfo);

    ///------------------------------------------------------------------------------------------------

    virtual bool& CalculateValue(Parameters& rParameterValues, const Variable<bool>& rThisVariable, bool& rValue);
    virtual int& CalculateValue(Parameters& rParameterValues, const Variable<int>& rThisVariable, int& rValue);
    virtual double& CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue);
    virtual Vector& CalculateValue(Parameters& rParameterValues, const Variable<Vector>& rThisVariable, Vector& rValue);
    virtual Matrix& CalculateValue(Parameters& rParameterValues, const Variable<Matrix>& rThisVariable, Matrix& rValue);
    virtual array_1d<double, 3>& CalculateValue(Parameters& rParameterValues, const Variable<array_1d<double, 3>>& rVariable, array_1d<double, 3>& rValue);
    virtual array_1d<double, 6>& CalculateValue(Parameters& rParameterValues, const Variable<array_1d<double, 6>>& rVariable, array_1d<double, 6>& rValue);

    ///------------------------------------------------------------------------------------------------

    virtual int Check(const Properties& rMaterialProperties,
                      const GeometryType& rElementGeometry,
                      const ProcessInfo& rCurrentProcessInfo) const;

    virtual bool ValidateInput(const Properties& rMaterialProperties);

    virtual void InitializeMaterial(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const Vector& rShapeFunctionsValues);

    virtual void CalculateMaterialResponse (Parameters& rValues);

    virtual bool RequiresInitializeMaterialResponse()
    {
        return false;
    }

    virtual void InitializeMaterialResponse (Parameters& rValues);

    virtual bool RequiresFinalizeMaterialResponse()
    {
        return false;
    }

    virtual void FinalizeMaterialResponse (Parameters& rValues);

    virtual void CalculateSaturation (Parameters& rValues);

    virtual void ResetMaterial(const Properties& rMaterialProperties,
                               const GeometryType& rElementGeometry,
                               const Vector& rShapeFunctionsValues);
    
    ///------------------------------------------------------------------------------------------------

    inline static bool HasSameType(const SaturationLaw& rLHS, const SaturationLaw& rRHS) {
        return (typeid(rLHS) == typeid(rRHS));
    }

    inline static bool HasSameType(const SaturationLaw* rLHS, const SaturationLaw* rRHS) {
        return SaturationLaw::HasSameType(*rLHS, *rRHS);
    }

    ///------------------------------------------------------------------------------------------------

    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "SaturationLaw";
        return buffer.str();
    }

    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "SaturationLaw";
    }

    virtual void PrintData(std::ostream& rOStream) const
    {
      rOStream << "SaturationLaw has no data";
    }

    ///------------------------------------------------------------------------------------------------

protected:

    struct SaturationLawVariables
    {
        double Slr; // Residual liquid saturation
        double Sgr; // Residual gas saturation
        double lambda; // Pore size factor
        double pb; // Gas entry pressure
        double pc; // Capillary pore pressure
        double Se; // Effective saturation
        double krmin; // Minimum relative permeability
    };

    /// Member Variables

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual void InitializeSaturationLawVariables(SaturationLawVariables& rVariables, Parameters& rValues);

    virtual void CalculateLiquidSaturationDegree(SaturationLawVariables& rVariables, Parameters& rValues);

    virtual void CalculateEffectiveSaturation(SaturationLawVariables& rVariables, Parameters& rValues);

    virtual void CalculateLiquidRelativePermeability(SaturationLawVariables& rVariables, Parameters& rValues);

    virtual void CalculateGasRelativePermeability(SaturationLawVariables& rVariables, Parameters& rValues);

    ///------------------------------------------------------------------------------------------------

private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

}; /* Class SaturationLaw */

inline std::istream& operator>>(std::istream& rIStream, SaturationLaw& rThis);

inline std::ostream& operator<<(std::ostream& rOStream, const SaturationLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/
