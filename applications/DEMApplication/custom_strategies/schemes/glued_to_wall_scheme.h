//
// Author: Miguel Angel Celigueta maceli@cimne.upc.edu
//

#if !defined(KRATOS_GLUED_TO_WALL_SCHEME_H_INCLUDED )
#define  KRATOS_GLUED_TO_WALL_SCHEME_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cfloat>

// Project includes
#include "dem_integration_scheme.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) GluedToWallScheme : public DEMIntegrationScheme {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;

        /// Pointer definition of GluedToWallScheme
        KRATOS_CLASS_POINTER_DEFINITION(GluedToWallScheme);

        /// Default constructor.
        GluedToWallScheme() {}

        GluedToWallScheme(Condition* p_wall, SphericParticle* p_sphere);

        GluedToWallScheme(Condition* p_wall, SphericParticle* p_sphere, bool& is_inside);

        /// Destructor.
        virtual ~GluedToWallScheme() {}

        DEMIntegrationScheme* CloneRaw() const override {
            DEMIntegrationScheme* cloned_scheme(new GluedToWallScheme(*this));
            return cloned_scheme;
        }

        DEMIntegrationScheme::Pointer CloneShared() const override {
            DEMIntegrationScheme::Pointer cloned_scheme(new GluedToWallScheme(*this));
            return cloned_scheme;
        }

        void SetTranslationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose = true) const override;
        void SetRotationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose = true) const override;

        void Move(Node & i, const double delta_t, const double force_reduction_factor, const int StepFlag) override;
        void Rotate(Node & i, const double delta_t, const double moment_reduction_factor, const int StepFlag) override;

        Condition* pGetCondition() {
            return mCondition;
        }

        Vector& GetShapeFunctionsValues() {
            return mShapeFunctionsValues;
        }

        double GetDistanceSignedWithNormal() {
            return mDistanceSignedWithNormal;
        }


        /// Turn back information as a string.
        virtual std::string Info() const override {
            std::stringstream buffer;
            buffer << "GluedToWallScheme";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const override {
            rOStream << "GluedToWallScheme";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const override {
        }


    protected:


    private:
        Condition* mCondition;
        Vector mShapeFunctionsValues;
        double mDistanceSignedWithNormal = 0.0;
        array_1d<double, 3> mInitialNormalToWall;
        array_1d<double, 3> mCurrentNormalToWall;

        /// Assignment operator.
        GluedToWallScheme& operator=(GluedToWallScheme const& rOther) {
            return *this;
        }

        /// Copy constructor.

        GluedToWallScheme(GluedToWallScheme const& rOther) {
            *this = rOther;
        }

        bool ShapeFunctionsValuesAreBetween0and1();
        ///@}

    }; // Class GluedToWallScheme

    inline std::istream& operator>>(std::istream& rIStream,
            GluedToWallScheme& rThis) {
        return rIStream;
    }

    inline std::ostream& operator<<(std::ostream& rOStream,
            const GluedToWallScheme& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_GLUED_TO_WALL_SCHEME_H_INCLUDED  defined
