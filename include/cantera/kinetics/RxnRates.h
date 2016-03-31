/**
 *  @file RxnRates.h
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_RXNRATES_H
#define CT_RXNRATES_H

#include "cantera/kinetics/reaction_defs.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"

#include <iostream>

namespace Cantera
{

	//! Turbulent reaction rate correction coefficient

static doublereal Cc(doublereal m_b, doublereal m_E, doublereal recipT, doublereal TprimeOverT) {

		doublereal R_const = 1.9872041, t1, t2, t3, t4, t5, t6, t7, ta;
		doublereal recipT2 = pow(recipT, 2), m_E2 = pow(m_E, 2), TPOT2 = pow(TprimeOverT, 2), R2 = pow(1.9872041, 2);
		doublereal recipT3 = pow(recipT, 3), m_E3 = pow(m_E, 3), TPOT3 = pow(TprimeOverT, 3), R3 = pow(1.9872041, 3);
		doublereal recipT4 = pow(recipT, 4), m_E4 = pow(m_E, 4), TPOT4 = pow(TprimeOverT, 4), R4 = pow(1.9872041, 4);
		doublereal recipT5 = pow(recipT, 5), m_E5 = pow(m_E, 5), TPOT5 = pow(TprimeOverT, 5), R5 = pow(1.9872041, 5);
		doublereal recipT6 = pow(recipT, 6), m_E6 = pow(m_E, 6), TPOT6 = pow(TprimeOverT, 6), R6 = pow(1.9872041, 6);
		doublereal recipT7 = pow(recipT, 7), m_E7 = pow(m_E, 7), TPOT7 = pow(TprimeOverT, 7), R7 = pow(1.9872041, 6);

		t1 = (((((m_b*R_const) + m_E)*recipT)*TprimeOverT) / (R_const));
		t2 = (((R2*m_b*(m_b - 1)) + (2 * m_E*R_const*(m_b - 1)*recipT) + (m_E2*recipT2)) / (2 * R2))*TPOT2;
		t3 = (((R3*m_b*(m_b - 1)*(m_b - 2)) + (3 * m_E*R2*(m_b - 1)*(m_b - 2)*recipT) + (3 * m_E2*R_const*(m_b - 2)*recipT2) + (m_E3*recipT3)) / (6 * R3))*TPOT3;
		t4 = ((((((R4*m_b*(m_b - 1)*(m_b - 2))*(m_b - 3)) + (4 * m_E*R3*(m_b - 1)*(m_b - 2))*(m_b - 3)*recipT) + (6 * m_E2*R2*(m_b - 2))*(m_b - 3)*recipT2) + (4 * m_E3*R_const*(m_b - 3)*recipT3) + (m_E4*recipT4)) / (24 * R4))*TPOT4;
		t5 = ((((((R5*m_b*(m_b - 1)*(m_b - 2))*(m_b - 3)*(m_b - 4)) + (5 * m_E*R4*(m_b - 1)*(m_b - 2))*(m_b - 3)*(m_b - 4)*recipT) + (10 * m_E2*R3*(m_b - 2))*(m_b - 3)*(m_b - 4)*recipT2) + (10 * m_E3*R2*(m_b - 3)*(m_b - 4)*recipT3) + (5 * m_E4*R_const*(m_b - 4)*recipT4) + (m_E5*recipT5)) / (120 * R5))*TPOT5;
		t6 = ((((((R6*m_b*(m_b - 1)*(m_b - 2))*(m_b - 3)*(m_b - 4)*(m_b - 5)) + (6 * m_E*R5*(m_b - 1)*(m_b - 2))*(m_b - 3)*(m_b - 4)*(m_b - 5)*recipT) + (15 * m_E2*R4*(m_b - 2))*(m_b - 3)*(m_b - 4)*(m_b - 5)*recipT2) + (20 * m_E3*R3*(m_b - 3)*(m_b - 4)*(m_b - 5)*recipT3) + (15 * m_E4*R2*(m_b - 4)*(m_b - 5)*recipT4) + (6 * m_E5*R_const*(m_b - 5)*recipT5) + (m_E6*recipT6)) / (720 * R6))*TPOT6;
		t7 = ((((((R7*m_b*(m_b - 1)*(m_b - 2))*(m_b - 3)*(m_b - 4)*(m_b - 5)*(m_b - 6)) + (7 * m_E*R6*(m_b - 1)*(m_b - 2))*(m_b - 3)*(m_b - 4)*(m_b - 5)*(m_b - 6)*recipT) + (21 * m_E2*R5*(m_b - 2))*(m_b - 3)*(m_b - 4)*(m_b - 5)*(m_b - 6)*recipT2) + (35 * m_E3*R4*(m_b - 3)*(m_b - 4)*(m_b - 5)*(m_b - 6)*recipT3) + (35 * m_E4*R3*(m_b - 4)*(m_b - 5)*(m_b - 6)*recipT4) + (21 * m_E5*R2*(m_b - 5)*(m_b - 6)*recipT5)) + (7 * m_E6*R_const*(m_b - 6)*recipT6) + (m_E7*recipT7)) / (5040 * R7)*TPOT7;
		ta = (((R2*m_b*(m_b)) - (R2*m_b) + (2 * m_E*R_const*(m_b - 1)*recipT) + (m_E2*recipT2)) / (4 * R2))*TPOT2;;

		doublereal CorrectionCoefficient = 1 + t1 + t2 + t3 + t4 + t5 + t6 + t7;

		if (CorrectionCoefficient>1.e5){
			CorrectionCoefficient = 1.e5;
		}

		return CorrectionCoefficient;
	}


class Array2D;

//! Arrhenius reaction rate type depends only on temperature
/**
 * A reaction rate coefficient of the following form.
 *
 *   \f[
 *        k_f =  A T^b \exp (-E/RT)
 *   \f]
 */
class Arrhenius
{
public:
    //! return the rate coefficient type.
    static int type() {
        return ARRHENIUS_REACTION_RATECOEFF_TYPE;
    }

    //! Default constructor.
    Arrhenius();

    /// Constructor.
    /// @param A pre-exponential. The unit system is
    ///     (kmol, m, s). The actual units depend on the reaction
    ///     order and the dimensionality (surface or bulk).
    /// @param b Temperature exponent. Non-dimensional.
    /// @param E Activation energy in temperature units. Kelvin.
    Arrhenius(doublereal A, doublereal b, doublereal E);

    //! Update concentration-dependent parts of the rate coefficient.
    /*!
     *   For this class, there are no concentration-dependent parts, so this
     *   method does nothing.
     */
    void update_C(const doublereal* c) {
    }

    /**
     * Update the value of the natural logarithm of the rate constant.
     */
    doublereal updateLog(doublereal logT, doublereal recipT) const {
        return m_logA + m_b*logT - m_E*recipT;
    }

    /**
     * Update the value the rate constant.
     *
     * This function returns the actual value of the rate constant. It can be
     * safely called for negative values of the pre-exponential factor.
     */
    doublereal updateRC(doublereal logT, doublereal recipT) const {
        return m_A * std::exp(m_b*logT - m_E*recipT);
    }
	//Update the value of the turbulent rate constant.

	doublereal updateTurbulent(doublereal logT, doublereal recipT, doublereal TprimeOverT) const {
		return  (updateRC(logT, recipT))*(Cc_return(recipT, TprimeOverT));
	}

    //! Return the pre-exponential factor *A* (in m, kmol, s to powers depending
    //! on the reaction order)
    double preExponentialFactor() const {
        return m_A;
    }

	double Cc_return( doublereal recipT, doublereal TprimeOverT)const {
		return Cc(temperatureExponent(), activationEnergy_R(), recipT, TprimeOverT);
	}

	double Cc_Out(doublereal recipT, doublereal TprimeOverT)const {
		return Cc(0, 6500, recipT, TprimeOverT);
	}

    //! Return the temperature exponent *b*
    double temperatureExponent() const {
        return m_b;
    }

    //! Return the activation energy divided by the gas constant (i.e. the
    //! activation temperature) [K]
    doublereal activationEnergy_R() const {
        return m_E;
    }

protected:
    doublereal m_logA, m_b, m_E, m_A;
};


/**
 * An Arrhenius rate with coverage-dependent terms.
 *
 * The rate expression is given by [Kee, R. J., Coltrin, M. E., & Glarborg, P.
 * (2005). Chemically reacting flow: theory and practice. John Wiley & Sons.
 * Eq 11.113]:
 * \f[
 *     k_f = A T^b \exp \left(
 *             \ln 10 \sum a_k \theta_k
 *             - \frac{1}{RT} \left( E_a + \sum E_k\theta_k \right)
 *             + \sum m_k \ln \theta_k
 *             \right)
 *   \f]
 * or, equivalently, and as implemented in Cantera,
 * \f[
 *     k_f = A T^b \exp \left( - \frac{E_a}{RT} \right)
 *             \prod_k 10^{a_k \theta_k} \theta_k^{m_k}
 *             \exp \left( \frac{- E_k \theta_k}{RT} \right)
 *   \f]
 * where the parameters \f$ (a_k, E_k, m_k) \f$ describe the dependency on the
 * surface coverage of species \f$k, \theta_k \f$.
 */
class SurfaceArrhenius
{

public:
    static int type() {
        return SURF_ARRHENIUS_REACTION_RATECOEFF_TYPE;
    }

    SurfaceArrhenius();
    explicit SurfaceArrhenius(double A, double b, double Ta);

    //! Add a coverage dependency for species *k*, with pre-exponential
    //! dependence *a*, rate constant exponential dependency *m*, and activation
    //! energy dependence *e*, where *e* is in Kelvin, i.e. energy divided by
    //! the molar gas constant.
    void addCoverageDependence(size_t k, doublereal a,
                               doublereal m, doublereal e);

    void update_C(const doublereal* theta) {
        m_acov = 0.0;
        m_ecov = 0.0;
        m_mcov = 0.0;
        size_t k;
        doublereal th;
        for (size_t n = 0; n < m_ac.size(); n++) {
            k = m_sp[n];
            m_acov += m_ac[n] * theta[k];
            m_ecov += m_ec[n] * theta[k];
        }
        for (size_t n = 0; n < m_mc.size(); n++) {
            k = m_msp[n];
            th = std::max(theta[k], Tiny);
            m_mcov += m_mc[n]*std::log(th);
        }
    }

    /**
     * Update the value the rate constant.
     *
     * This function returns the actual value of the rate constant. It can be
     * safely called for negative values of the pre-exponential factor.
     */
    doublereal updateRC(doublereal logT, doublereal recipT) const {
        return m_A * std::exp(std::log(10.0)*m_acov + m_b*logT -
                              (m_E + m_ecov)*recipT + m_mcov);
    }

    //! Return the pre-exponential factor *A* (in m, kmol, s to powers depending
    //! on the reaction order) accounting coverage dependence.
    /*!
     *  Returns reaction prexponent accounting for both *a* and *m*.
     */
    doublereal preExponentialFactor() const {
        return m_A * std::exp(std::log(10.0)*m_acov + m_mcov);
    }

    //! Return effective temperature exponent
    doublereal temperatureExponent() const {
        return m_b;
    }

    //! Return the activation energy divided by the gas constant (i.e. the
    //! activation temperature) [K], accounting coverage dependence.
    doublereal activationEnergy_R() const {
        return m_E + m_ecov;
    }

protected:
    doublereal m_b, m_E, m_A;
    doublereal m_acov, m_ecov, m_mcov;
    std::vector<size_t> m_sp, m_msp;
    vector_fp m_ac, m_ec, m_mc;
};


//! Pressure-dependent reaction rate expressed by logarithmically interpolating
//! between Arrhenius rate expressions at various pressures.
class Plog
{
public:
    //! return the rate coefficient type.
    static int type() {
        return PLOG_REACTION_RATECOEFF_TYPE;
    }

    //! Default constructor.
    Plog() {}

    //! Constructor from Arrhenius rate expressions at a set of pressures
    explicit Plog(const std::multimap<double, Arrhenius>& rates);

    //! Update concentration-dependent parts of the rate coefficient.
    //! @param c natural log of the pressure in Pa
    void update_C(const doublereal* c) {
        logP_ = c[0];
        if (logP_ > logP1_ && logP_ < logP2_) {
            return;
        }

        auto iter = pressures_.upper_bound(c[0]);
        AssertThrowMsg(iter != pressures_.end(), "Plog::update_C",
                       "Pressure out of range: {}", logP_);
        AssertThrowMsg(iter != pressures_.begin(), "Plog::update_C",
                       "Pressure out of range: {}", logP_);

        // upper interpolation pressure
        logP2_ = iter->first;
        ihigh1_ = iter->second.first;
        ihigh2_ = iter->second.second;

        // lower interpolation pressure
        logP1_ = (--iter)->first;
        ilow1_ = iter->second.first;
        ilow2_ = iter->second.second;

        rDeltaP_ = 1.0 / (logP2_ - logP1_);
    }

    /**
     * Update the value the rate constant.
     *
     * This function returns the actual value of the rate constant.
     */
    doublereal updateRC(doublereal logT, doublereal recipT) const {
        double log_k1, log_k2;
        if (ilow1_ == ilow2_) {
            log_k1 = rates_[ilow1_].updateLog(logT, recipT);
        } else {
            double k = 1e-300; // non-zero to make log(k) finite
            for (size_t i = ilow1_; i < ilow2_; i++) {
                k += rates_[i].updateRC(logT, recipT);
            }
            log_k1 = std::log(k);
        }

        if (ihigh1_ == ihigh2_) {
            log_k2 = rates_[ihigh1_].updateLog(logT, recipT);
        } else {
            double k = 1e-300; // non-zero to make log(k) finite
            for (size_t i = ihigh1_; i < ihigh2_; i++) {
                k += rates_[i].updateRC(logT, recipT);
            }
            log_k2 = std::log(k);
        }

        return std::exp(log_k1 + (log_k2-log_k1) * (logP_-logP1_) * rDeltaP_);
    }

	/**
	* Update the value of the logarithm of the turbulent rate constant.
	*/
	doublereal updateTurbLog(doublereal logT, doublereal recipT, doublereal TprimeOverT) const {
		double log_k1, log_k2;
		if (ilow1_ == ilow2_) {
			log_k1 = rates_[ilow1_].updateLog(logT, recipT) * rates_[ilow1_].Cc_return(recipT, TprimeOverT);
		}
		else {
			double k = 1e-300; // non-zero to make log(k) finite
			double kTurb = 1e-300;
			for (size_t i = ilow1_; i < ilow2_; i++) {
				k += rates_[i].updateRC(logT, recipT);
				kTurb += k + (k * rates_[ilow1_].Cc_return(recipT, TprimeOverT));
				
			}
			log_k1 = std::log(kTurb);
		}

		if (ihigh1_ == ihigh2_) {
			log_k2 = rates_[ihigh1_].updateLog(logT, recipT)* rates_[ihigh1_].Cc_return(recipT, TprimeOverT);
		}
		else {
			double k = 1e-300; // non-zero to make log(k) finite
			double kTurb = 1e-300;
			for (size_t i = ihigh1_; i < ihigh2_; i++) {
				k += rates_[i].updateRC(logT, recipT);
				kTurb += k + k*rates_[ihigh1_].Cc_return(recipT, TprimeOverT);
				
			}
			log_k2 = std::log(kTurb);
		}

		return std::exp(log_k1 + (log_k2 - log_k1) * (logP_ - logP1_) * rDeltaP_);
	}

	/**
	* Update the value the turbulent rate constant.
	*
	* This function returns the actual value of the turbulent rate constant.
	*/
	doublereal updateTurbulent(doublereal logT, doublereal recipT, doublereal TprimeOverT) const {
		return std::exp(updateTurbLog(logT, recipT, TprimeOverT));
	}

    //! Check to make sure that the rate expression is finite over a range of
    //! temperatures at each interpolation pressure. This is potentially an
    //! issue when one of the Arrhenius expressions at a particular pressure
    //! has a negative pre-exponential factor.
    void validate(const std::string& equation);

    //! Return the pressures and Arrhenius expressions which comprise this
    //! reaction.
    std::vector<std::pair<double, Arrhenius> > rates() const;

protected:
    //! log(p) to (index range) in the rates_ vector
    std::map<double, std::pair<size_t, size_t> > pressures_;

    // Rate expressions which are referenced by the indices stored in pressures_
    std::vector<Arrhenius> rates_;

    double logP_; //!< log(p) at the current state
    double logP1_, logP2_; //!< log(p) at the lower / upper pressure reference

    //! Indices to the ranges within rates_ for the lower / upper pressure, such
    //! that rates_[ilow1_] through rates_[ilow2_] (inclusive) are the rates
    //! expressions which are combined to form the rate at the lower reference
    //! pressure.
    size_t ilow1_, ilow2_, ihigh1_, ihigh2_;

    double rDeltaP_; //!< reciprocal of (logP2 - logP1)
};

//! Pressure-dependent rate expression where the rate coefficient is expressed
//! as a bivariate Chebyshev polynomial in temperature and pressure.
class ChebyshevRate
{
public:
    //! return the rate coefficient type.
    static int type() {
        return CHEBYSHEV_REACTION_RATECOEFF_TYPE;
    }

    //! Default constructor.
    ChebyshevRate() {}

    //! Constructor directly from coefficient array
    /*
     *  @param Tmin    Minimum temperature [K]
     *  @param Tmax    Maximum temperature [K]
     *  @param Pmin    Minimum pressure [Pa]
     *  @param Pmax    Maximum pressure [Pa]
     *  @param coeffs  Coefficient array dimensioned `nT` by `nP` where `nT` and
     *      `nP` are the number of temperatures and pressures used in the fit,
     *      respectively.
     */
    ChebyshevRate(double Tmin, double Tmax, double Pmin, double Pmax,
                  const Array2D& coeffs);

    //! Update concentration-dependent parts of the rate coefficient.
    //! @param c base-10 logarithm of the pressure in Pa
    void update_C(const doublereal* c) {
        double Pr = (2 * c[0] + PrNum_) * PrDen_;
        double Cnm1 = 1;
        double Cn = Pr;
        double Cnp1;
        for (size_t j = 0; j < nT_; j++) {
            dotProd_[j] = chebCoeffs_[nP_*j] + Pr * chebCoeffs_[nP_*j+1];
        }
        for (size_t i = 2; i < nP_; i++) {
            Cnp1 = 2 * Pr * Cn - Cnm1;
            for (size_t j = 0; j < nT_; j++) {
                dotProd_[j] += Cnp1 * chebCoeffs_[nP_*j + i];
            }
            Cnm1 = Cn;
            Cn = Cnp1;
        }
    }

    /**
     * Update the value the rate constant.
     *
     * This function returns the actual value of the rate constant.
     */
    doublereal updateRC(doublereal logT, doublereal recipT) const {
        double Tr = (2 * recipT + TrNum_) * TrDen_;
        double Cnm1 = 1;
        double Cn = Tr;
        double Cnp1;
        double logk = dotProd_[0] + Tr * dotProd_[1];
        for (size_t i = 2; i < nT_; i++) {
            Cnp1 = 2 * Tr * Cn - Cnm1;
            logk += Cnp1 * dotProd_[i];
            Cnm1 = Cn;
            Cn = Cnp1;
        }
        return std::pow(10, logk);
    }

	doublereal updateTurbulent(doublereal logT, doublereal recipT, doublereal TprimeOverT) const {
		throw CanteraError("ChebyshevRate::updateTurbulent", "Not implemented");
	}

    //! Minimum valid temperature [K]
    double Tmin() const {
        return Tmin_;
    }

    //! Maximum valid temperature [K]
    double Tmax() const {
        return Tmax_;
    }

    //! Minimum valid pressure [Pa]
    double Pmin() const {
        return Pmin_;
    }

    //! Maximum valid pressure [Pa]
    double Pmax() const {
        return Pmax_;
    }

    //! Number of points in the pressure direction
    size_t nPressure() const {
        return nP_;
    }

    //! Number of points in the temperature direction
    size_t nTemperature() const {
        return nT_;
    }

    //! Access the Chebyshev coefficients.
    /*!
     *  \f$ \alpha_{t,p} = \mathrm{coeffs}[N_P*t + p] \f$ where
     *  \f$ 0 <= t < N_T \f$ and \f$ 0 <= p < N_P \f$.
     */
    const vector_fp& coeffs() const {
        return chebCoeffs_;
    }

protected:
    double Tmin_, Tmax_; //!< valid temperature range
    double Pmin_, Pmax_; //!< valid pressure range
    double TrNum_, TrDen_; //!< terms appearing in the reduced temperature
    double PrNum_, PrDen_; //!< terms appearing in the reduced pressure

    size_t nP_; //!< number of points in the pressure direction
    size_t nT_; //!< number of points in the temperature direction
    vector_fp chebCoeffs_; //!< Chebyshev coefficients, length nP * nT
    vector_fp dotProd_; //!< dot product of chebCoeffs with the reduced pressure polynomial
};

}

#endif
