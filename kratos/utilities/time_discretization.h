//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Ruben Zorrilla
//

#if !defined( KRATOS_TIME_DISCRETIZATION_H_INCLUDED )
#define  KRATOS_TIME_DISCRETIZATION_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes

namespace Kratos {

class ProcessInfo; // forward-declaring to not having to include it here

namespace TimeDiscretization {

class KRATOS_API(KRATOS_CORE) BDF
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(BDF);

    /**
     * @brief Construct a new BDF object
     * Constructor with time order
     * @param TimeOrder BDF time order
     */
    BDF(const unsigned int TimeOrder) : mTimeOrder(TimeOrder)
    {
        this->SetAuxBDFPointer(TimeOrder, mpAuxBDF);
    };

    /**
     * @brief Destroy the BDF object
     * Destructor of the BDF class
     */
    virtual ~BDF() = default;

    /**
     * @brief Return the BDF coefficients
     * This method computes the BDF coefficients for a provided time step
     * @param DeltaTime Time step to compute the BDF coefficients
     * @return std::vector<double> Vector containing the computed BDF coefficients
     */
    virtual std::vector<double> ComputeBDFCoefficients(double DeltaTime) const;

    /**
     * @brief Return the BDF coefficients
     * This method computes the BDF coefficients for a provided current and old time step
     * @param DeltaTime Time step to compute the BDF coefficients
     * @param PreviousDeltaTime Old time step to compute the BDF coefficients
     * @return std::vector<double> Vector containing the computed BDF coefficients
     */
    virtual std::vector<double> ComputeBDFCoefficients(double DeltaTime, double PreviousDeltaTime) const;

    /**
     * @brief Return the BDF coefficients
     * This method computes the BDF coefficients for the time step
     * stored in the variable DELTA_TIME of a ProcessInfo container
     * @param rProcessInfo ProcessInfo container with DELTA_TIME
     * @return std::vector<double> Vector containing the computed BDF coefficients
     */
    virtual std::vector<double> ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const;

    /**
     * @brief Computes the BDF coefficients
     * This method computes the BDF coefficients with the information
     * stored in the provided ProcessInfo. The computed coefficients
     * are saved in such ProcessInfo container.
     * @param rProcessInfo ProcessInfo container with DELTA_TIME
     */
    void ComputeAndSaveBDFCoefficients(ProcessInfo &rProcessInfo) const;

    /**
     * @brief Get the Time Order object
     * Auxiliary method to get the order of the BDF scheme
     * @return const std::size_t Order of the BDF scheme
     */
    std::size_t GetTimeOrder() const;

protected:

    /**
     * @brief Construct a new BDF object
     * Auxiliary constructor for derived classes
     */
    BDF(){};

private:

    const std::size_t mTimeOrder = 0; // Time order of the auxiliary BDF class
    Kratos::unique_ptr<BDF> mpAuxBDF = nullptr; // Pointer to an auxiliary BDF class with order

    /**
     * @brief Set the Aux BDF Pointer object
     * This method sets a pointer to one of the auxiliary BDF
     * classes according to the user defined time order
     * @param TimeOrder BDF scheme time order
     * @param rpAuxBDF Pointer to the auxiliary BDF scheme
     */
    static void SetAuxBDFPointer(
        const std::size_t TimeOrder,
        Kratos::unique_ptr<BDF> &rpAuxBDF);
};

class KRATOS_API(KRATOS_CORE) BDF1 : public BDF
{
public:
    std::vector<double> ComputeBDFCoefficients(double DeltaTime) const override;
    std::vector<double> ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const override;
};

class KRATOS_API(KRATOS_CORE) BDF2 : public BDF
{
public:
    std::vector<double> ComputeBDFCoefficients(double DeltaTime, double PreviousDeltaTime) const override;
    std::vector<double> ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const override;
};

class KRATOS_API(KRATOS_CORE) BDF3 : public BDF
{
public:
    std::vector<double> ComputeBDFCoefficients(double DeltaTime) const override;
    std::vector<double> ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const override;
};

class KRATOS_API(KRATOS_CORE) BDF4 : public BDF
{
public:
    std::vector<double> ComputeBDFCoefficients(double DeltaTime) const override;
    std::vector<double> ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const override;
};

class KRATOS_API(KRATOS_CORE) BDF5 : public BDF
{
public:
    std::vector<double> ComputeBDFCoefficients(double DeltaTime) const override;
    std::vector<double> ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const override;
};

class KRATOS_API(KRATOS_CORE) BDF6 : public BDF
{
public:
    std::vector<double> ComputeBDFCoefficients(double DeltaTime) const override;
    std::vector<double> ComputeBDFCoefficients(const ProcessInfo& rProcessInfo) const override;
};

class Newmark
{
public:
    Newmark() = default;

    explicit Newmark(const double NewmarkBeta,
                     const double NewmarkGamma)
        : mNewmarkBeta(NewmarkBeta),
          mNewmarkGamma(NewmarkGamma) {}

    double GetBeta()  const { return mNewmarkBeta; }
    double GetGamma() const { return mNewmarkGamma; }

private:
    double mNewmarkBeta=0.25;
    double mNewmarkGamma=0.5;
};

class Bossak
{
public:
    Bossak() = default;

    explicit Bossak(const double AlphaM)
        : mAlphaM(AlphaM) {}

    explicit Bossak(const double AlphaM,
                    const double NewmarkBeta,
                    const double NewmarkGamma)
        : mAlphaM(AlphaM),
          mNewmarkBeta(NewmarkBeta),
          mNewmarkGamma(NewmarkGamma) {}

    double GetAlphaM() const { return mAlphaM; }
    double GetBeta()   const { return mNewmarkBeta * (1-mAlphaM) * (1-mAlphaM); }
    double GetGamma()  const { return mNewmarkGamma - mAlphaM; }

private:
    double mAlphaM=-0.3;
    double mNewmarkBeta=0.25;
    double mNewmarkGamma=0.5;
};

class GeneralizedAlpha
{
public:
    GeneralizedAlpha() = default;

    explicit GeneralizedAlpha(const double AlphaM,
                              const double AlphaF)
        : mAlphaM(AlphaM),
          mAlphaF(AlphaF) {}

    explicit GeneralizedAlpha(const double AlphaM,
                              const double AlphaF,
                              const double NewmarkBeta,
                              const double NewmarkGamma)
        : mAlphaM(AlphaM),
          mAlphaF(AlphaF),
          mNewmarkBeta(NewmarkBeta),
          mNewmarkGamma(NewmarkGamma) {}

    double GetAlphaM() const { return mAlphaM; }
    double GetAlphaF() const { return mAlphaF; }
    double GetBeta()   const { return mNewmarkBeta * (1-mAlphaM+mAlphaF) * (1-mAlphaM+mAlphaF); }
    double GetGamma()  const { return mNewmarkGamma - mAlphaM + mAlphaF; }

private:
    double mAlphaM=-0.3;
    double mAlphaF=0.0;
    double mNewmarkBeta=0.25;
    double mNewmarkGamma=0.5;
};

inline std::size_t GetMinimumBufferSize(const BDF& rTimeDiscr) { return (rTimeDiscr.GetTimeOrder() + 1); }
inline std::size_t GetMinimumBufferSize(const BDF1& rTimeDiscr) { return 2;}
inline std::size_t GetMinimumBufferSize(const BDF2& rTimeDiscr) { return 3;}
inline std::size_t GetMinimumBufferSize(const BDF3& rTimeDiscr) { return 4;}
inline std::size_t GetMinimumBufferSize(const BDF4& rTimeDiscr) { return 5;}
inline std::size_t GetMinimumBufferSize(const BDF5& rTimeDiscr) { return 6;}
inline std::size_t GetMinimumBufferSize(const BDF6& rTimeDiscr) { return 7;}

inline std::size_t GetMinimumBufferSize(const Newmark& rTimeDiscr)          { return 2;}
inline std::size_t GetMinimumBufferSize(const Bossak& rTimeDiscr)           { return 2;}
inline std::size_t GetMinimumBufferSize(const GeneralizedAlpha& rTimeDiscr) { return 2;}

} // namespace TimeDiscretization.
}  // namespace Kratos.

#endif // KRATOS_TIME_DISCRETIZATION_H_INCLUDED  defined
