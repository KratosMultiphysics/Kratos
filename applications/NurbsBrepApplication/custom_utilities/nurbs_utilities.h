#ifndef NURBS_UTILITIES
#define NURBS_UTILITIES

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <vector>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/math/special_functions/binomial.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/variables.h"
#include "includes/define.h"
//#include "nurbs_brep_application.h"
//#include "nurbs_brep_application_variables.h"

// ==============================================================================

namespace Kratos
{
  class NurbsUtilities
  {
  public:
    /// Destructor.
    virtual ~NurbsUtilities() { }

    /// Pointer definition of MapperUtilities
    KRATOS_CLASS_POINTER_DEFINITION(NurbsUtilities);

    // adapted from ALGORITHM A2.1, Page 68, "The NURBS Book" by Piegl, Tiller 
    // T.Oberbichler 02.2017
    //static int find_knot_span(const int &p, const Vector &knots, const double &u) {
    //  const auto upper = upper_bound(knots.begin() + p, knots.end() - p, u);

    //  const auto span = distance(knots.begin(), upper) - 1;

    //  return span;
    //}

    // #####################################################################################
    // #####################################################################################
    ///
    ///   \details   evaluates all non-zero B-Spline basis functions (recursive formula)
    ///              based on a specified u and its corresponding knot span
    ///              Piegl,L. and Tiller,W. "The NURBS Book - 2nd Edition", Springer Verlag
    ///              Algorithm A2.2
    ///
    /// ======================================================================================
    ///   \param[in]   _NBasisFct vector of non-zero B-Spline basis functions
    ///                           N = [Ni-p,p(u),...,Ni,p(u)]
    ///   \param[in]  _knotVec    knot vector in _par direction
    ///   \param[in]  _par        local parameter (curve -> 1D, isoline -> 2D)
    ///   \param[in]  _span       corresponding knot span
    ///   \param[in]  _pDeg       polynominal degree in _par direction
    ///
    /// ======================================================================================
    ///  \author     Daniel Baumgärtner (12/2016)
    //
    //########################################################################################
    static void eval_nonzero_basis_function(Vector& _NBasisFct, Vector& _knotVec, double _par, int _span, int _pDeg)
    {
      Vector left;
      Vector right;
      double saved;
      double temp;

      left.resize(_pDeg + 1);
      right.resize(_pDeg + 1);
      _NBasisFct.resize(_pDeg + 1);
      _NBasisFct(0) = 1.00;

      for (int i = 1; i <= _pDeg; i++)
      {
        left(i) = _par - _knotVec[_span + 1 - i];
        right(i) = _knotVec[_span + i] - _par;

        saved = 0.00;
        for (int j = 0; j < i; j++)
        {
          temp = _NBasisFct(j) / (right[j + 1] + left(i - j));
          _NBasisFct(j) = saved + right(j + 1)*temp;
          saved = left(i - j)*temp;
        }
        _NBasisFct(i) = saved;
      }
    }

    //  #####################################################################################
    // #######################################################################################
    ///
    ///   \details   evaluates all non-zero B-Spline basis functions (recursive formula)
    ///              and its n-th derivative based on a specified u and its corresponding
    ///              knot span
    ///              Piegl,L. and Tiller,W. "The NURBS Book - 2nd Edition", Springer Verlag
    ///              Algorithm A2.3
    ///
    /// ======================================================================================
    ///   \param[in] _dNBasisFct   matrix including the kth derivatives of the non-zero
    ///                            B-spline basis fuctions (reals)
    ///                            dNdu = |Ni-p,p(u)          ,...,  Ni,p(u)         |
    ///                                   |dNdu_i-p,p(u)      ,...,  dNdu_i,p(u)     |
    ///                                   |  :                  :       :            |
    ///                                   |(dNdu_i-p,p(u))^k  ,...,  (dNdu_i,p(u))^k |
    ///   \param[in]  _knotVec    knot vector in _par direction
    ///   \param[in]  _par        local parameter (curve -> 1D, isoline -> 2D) (knot of interest)
    ///   \param[in]  _span       corresponding knot span
    ///   \param[in]  _pDeg       polynominal degree in _par direction
    ///
    /// ======================================================================================
    ///  \author     Daniel Baumgärtner (12/2016)
    //
    //########################################################################################
    static void eval_nonzero_basis_function_with_derivatives(matrix<double>& _dNBasisFct, Vector _knotVec, double _par, int _span, int _pDeg, int _kth)
    {
      matrix<double> ndu;
      matrix<double> a;
      Vector left;
      Vector right;
      int s1, s2, j1, j2, jk, pk, jj, ll;
      double d, saved, temp;

      _dNBasisFct.resize(_kth + 1, _pDeg + 1);
      ndu.resize(_pDeg + 1, _pDeg + 1);
      a.resize(2, _pDeg + 1);
      left.resize(_pDeg + 1);
      right.resize(_pDeg + 1);
      ndu(0, 0) = 1.00;

      for (int i = 1; i <= _pDeg; i++) {
        left(i) = _par - _knotVec[_span + 1 - i];
        right(i) = _knotVec[_span + i] - _par;
        saved = 0.00;
        for (int j = 0; j < i; j++) {
          ndu(i, j) = right(j + 1) + left(i - j);
          temp = ndu(j, i - 1) / ndu(i, j);
          ndu(j, i) = saved + right(j + 1)*temp;
          saved = left(i - j)*temp;
        }
        ndu(i, i) = saved;
      }

      for (int i = 0; i <= _pDeg; i++) {
        _dNBasisFct(0, i) = ndu(i, _pDeg);
      }
      for (int j = 0; j <= _pDeg; j++) {
        s1 = 0;
        s2 = 1;
        a(0, 0) = 1.00;
        for (int k = 1; k <= _kth; k++) {
          d = 0.00;
          jk = j - k;
          pk = _pDeg - k;
          if (j >= k) {
            a(s2, 0) = a(s1, 0) / ndu(pk + 1, jk);
            d = a(s2, 0)*ndu(jk, pk);
          }
          if (jk >= -1) {
            j1 = 1;
          }
          else {
            j1 = -jk;
          }
          if (j - 1 <= pk) {
            j2 = k - 1;
          }
          else {
            j2 = _pDeg - j;
          }
          for (int l = j1; l <= j2; l++) {
            a(s2, l) = (a(s1, l) - a(s1, l - 1)) / ndu(pk + 1, jk + l);
            d += a(s2, l)*ndu(jk + l, pk);
          }
          if (j <= pk) {
            a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, j);
            d += a(s2, k)*ndu(j, pk);
          }
          _dNBasisFct(k, j) = d;
          ll = s1;
          s1 = s2;
          s2 = ll;
        }
      }
      jj = _pDeg;
      for (int k = 1; k <= _kth; k++) {
        for (int l = 0; l <= _pDeg; l++) {
          _dNBasisFct(k, l) *= jj;
        }
        jj *= (_pDeg - k);
      }
    }

    //  #####################################################################################
    // #######################################################################################
    //#
    //#                  ++++++++++++++++++++++++++++++++++++
    //#                  +++  NurbsBasis::find_Knot_Span  +++
    //#                  ++++++++++++++++++++++++++++++++++++
    //#
    ///   \details   returns the corresponding knot span based on a specified u
    ///              Piegl,L. and Tiller,W. "The NURBS Book - 2nd Edition", Springer Verlag
    ///              Algorithm A2.1
    ///
    /// ======================================================================================
    ///   \return     index of the first knot in the knot span in _par direction
    ///   \param[in]  _knotVec  knot vector in _par direction
    ///   \param[in]  _par      local parameter (curve -> 1D, isoline -> 2D)
    ///   \param[in]  _pDeg     polynominal degree in _par direction
    ///   \param[in]  _nCtrl    number of control points in _par direction
    ///
    /// ======================================================================================
    ///  \author     Daniel Baumgärtner (12/2016)
    //
    //########################################################################################
    static int find_knot_span(int _pDeg, Vector _knotVec, double _par)
    {
    	double epsilon = 1.0e-11;
    
    	if(_par<=_knotVec[0])
    		_par = _knotVec[0]+epsilon;
    
    	if(_par>=_knotVec[_knotVec.size()-1])
    		_par = _knotVec[_knotVec.size()-1]-epsilon;
    

      int _nCtrl = _knotVec.size() - _pDeg - 1;
    	int low = _pDeg;
    	int high = _nCtrl+1;
    	int span=(low+high)/2;
    
    	while (_par < _knotVec[span] || _par >= _knotVec[span+1])
    	{
    		if (_par < _knotVec[span])
    		{
    			high = span;
    		}
    		else
    		{
    			low = span;
    		}
    		span = (low+high)/2;
    	}
    
    	return span;
    }

  protected:

  private:
    /// Default constructor.
    NurbsUtilities() { }
    /// Assignment operator.
    NurbsUtilities& operator=(NurbsUtilities const& rOther);

    ///// input stream function
    //inline std::istream& operator >> (std::istream& rIStream,
    //  NurbsUtilities& rThis)
    //{
    //  return rIStream;
    //}

    ///// output stream function
    //inline std::ostream& operator << (std::ostream& rOStream,
    //  const NurbsUtilities& rThis)
    //{
    //  rThis.PrintInfo(rOStream);
    //  rOStream << std::endl;
    //  rThis.PrintData(rOStream);

    //  return rOStream;
    //}
  };
} // namespace Kratos.

#endif // NURBS_UTILITIES
