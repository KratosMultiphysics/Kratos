////////////////////////////////////////////////////////////////////////////////
// Copyright © 2018 Jérémie Dumas
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////

#include "GCMMASolver.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
// PUBLIC
////////////////////////////////////////////////////////////////////////////////

GCMMASolver::GCMMASolver(int nn, int mm, double ai, double ci, double di)
	: n(nn)
	, m(mm)
	, outeriter(0)
	, raa0eps(1e-6)
	, raaeps(raa0eps)
	, xmamieps(1e-5)
	//, epsimin(1e-7)
	, epsimin(std::sqrt(n + m) * 1e-9)
	, move(0.5)
	, albefa(0.1)
	, asyminit(0.5) // 0.2;
	, asymdec(0.7) // 0.65;
	, asyminc(1.2) // 1.08;
	, raa0(1.0)
	, raa(m, 1.0)
	, a(m, ai)
	, c(m, ci)
	, d(m, di)
	, y(m)
	, lam(m)
	, mu(m)
	, s(2 * m)
	, low(n)
	, upp(n)
	, alpha(n)
	, beta(n)
	, p0(n)
	, q0(n)
	, pij(n * m)
	, qij(n * m)
	, b(m)
	, grad(m)
	, hess(m * m)
	, r(m)
	, fapp(m)
	, xold1(n)
	, xold2(n)
{ }

void GCMMASolver::SetAsymptotes(double init, double decrease, double increase) {
	// Asymptotes initialization and increase/decrease
	asyminit = init;
	asymdec = decrease;
	asyminc = increase;
}

void GCMMASolver::OuterUpdate(double *xmma, const double *xval, double f0x, const double *df0dx, const double *fx,
                              const double *dfdx, const double *xmin, const double *xmax)
{
	// Compute asymptotes
	Asymp(xval, df0dx, dfdx, xmin, xmax);

	// Generate the subproblem
	GenSub(xval, f0x, df0dx, fx, dfdx, xmin, xmax);

	// Update xolds
	xold2 = xold1;
	std::copy_n(xval, n, xold1.data());

	// Solve the dual with an interior point method
	std::copy_n(xval, n, xmma);
	SolveDIP(xmma);

	// Solve the dual with a steepest ascent method
	// SolveDSA(xmma);

	// Compute approximation values
	ComputeApprox(xmma);
}

void GCMMASolver::InnerUpdate(double *xmma, double f0xnew, const double *fxnew, const double *xval, double f0x,
                              const double *df0dx, const double *fx, const double *dfdx, const double *xmin,
                              const double *xmax)
{
	// Update approximation factors
	RaaUpdate(xmma, xval, f0xnew, fxnew, xmin, xmax);

	// Generate the subproblem
	GenSub(xval, f0x, df0dx, fx, dfdx, xmin, xmax);

	// Solve the dual with an interior point method
	std::copy_n(xval, n, xmma);
	SolveDIP(xmma);

	// Solve the dual with a steepest ascent method
	// SolveDSA(xmma);

	// Compute approximation values
	ComputeApprox(xmma);
}

bool GCMMASolver::ConCheck(double f0xnew, const double *fxnew) const {
	if (f0app + epsimin < f0xnew) {
		return false;
	}
	for (int j = 0; j < m; ++j) {
		if (fapp[j] + epsimin < fxnew[j]) {
			return false;
		}
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// PRIVATE
////////////////////////////////////////////////////////////////////////////////

void GCMMASolver::SolveDIP(double *x) {

	for (int j = 0; j < m; ++j) {
		lam[j] = c[j] / 2.0;
		mu[j] = 1.0;
	}

	const double tol = epsimin; // 1.0e-9*sqrt(m+n);
	double epsi = 1.0;
	double err = 1.0;
	int loop;

	while (epsi > tol) {

		loop = 0;
		while (err > 0.9 * epsi && loop < 100) {
			loop++;

			// Set up Newton system
			XYZofLAMBDA(x);
			DualGrad(x);
			for (int j = 0; j < m; ++j) {
				grad[j] = -1.0 * grad[j] - epsi / lam[j];
			}
			DualHess(x);

			// Solve Newton system
			if (m > 1) {
				Factorize(hess.data(), m);
				Solve(hess.data(), grad.data(), m);
				for (int j = 0; j < m; ++j) {
					s[j] = grad[j];
				}
			} else if (m > 0) {
				s[0] = grad[0] / hess[0];
			}

			// Get the full search direction
			for (int i = 0; i < m; ++i) {
				s[m + i] = -mu[i] + epsi / lam[i] - s[i] * mu[i] / lam[i];
			}

			// Perform linesearch and update lam and mu
			DualLineSearch();

			XYZofLAMBDA(x);

			// Compute KKT res
			err = DualResidual(x, epsi);
		}
		epsi = epsi * 0.1;
	}
}

void GCMMASolver::SolveDSA(double *x) {

	for (int j = 0; j < m; ++j) {
		lam[j] = 1.0;
	}

	const double tol = epsimin; // 1.0e-9*sqrt(m+n);
	double err = 1.0;
	int loop = 0;

	while (err > tol && loop < 500) {
		loop++;
		XYZofLAMBDA(x);
		DualGrad(x);
		double theta = 1.0;
		err = 0.0;
		for (int j = 0; j < m; ++j) {
			lam[j] = std::max(0.0, lam[j] + theta * grad[j]);
			err += grad[j] * grad[j];
		}
		err = std::sqrt(err);
	}
}

double GCMMASolver::DualResidual(double *x, double epsi) {

	std::vector<double> res(2 * m);

	for (int j = 0; j < m; ++j) {
		res[j] = -b[j] - a[j] * z - y[j] + mu[j];
		res[j + m] = mu[j] * lam[j] - epsi;
		for (int i = 0; i < n; ++i) {
			res[j] += pij[i * m + j] / (upp[i] - x[i]) + qij[i * m + j] / (x[i] - low[i]);
		}
	}

	double nrI = 0.0;
	for (int i = 0; i < 2 * m; ++i) {
		if (nrI < std::abs(res[i])) {
			nrI = std::abs(res[i]);
		}
	}

	return nrI;
}

void GCMMASolver::DualLineSearch() {

	double theta = 1.005;
	for (int i = 0; i < m; ++i) {
		if (theta < -1.01 * s[i] / lam[i]) {
			theta = -1.01 * s[i] / lam[i];
		}
		if (theta < -1.01 * s[i + m] / mu[i]) {
			theta = -1.01 * s[i + m] / mu[i];
		}
	}
	theta = 1.0 / theta;

	for (int i = 0; i < m; ++i) {
		lam[i] = lam[i] + theta * s[i];
		mu[i] = mu[i] + theta * s[i + m];
	}
}

void GCMMASolver::DualHess(double *x) {

	std::vector<double> df2(n);
	std::vector<double> PQ(n * m);

	double pjlam, qjlam;
	for (int i = 0; i < n; ++i) {
		pjlam = p0[i];
		qjlam = q0[i];
		for (int j = 0; j < m; ++j) {
			pjlam += pij[i * m + j] * lam[j];
			qjlam += qij[i * m + j] * lam[j];
			PQ[i * m + j] = pij[i * m + j] / pow(upp[i] - x[i], 2.0) - qij[i * m + j] / pow(x[i] - low[i], 2.0);
		}
		df2[i] = -1.0 / (2.0 * pjlam / pow(upp[i] - x[i], 3.0) + 2.0 * qjlam / pow(x[i] - low[i], 3.0));
		double xp = (sqrt(pjlam) * low[i] + sqrt(qjlam) * upp[i]) / (sqrt(pjlam) + sqrt(qjlam));
		if (xp < alpha[i]) {
			df2[i] = 0.0;
		}
		if (xp > beta[i]) {
			df2[i] = 0.0;
		}
	}

	// Create the matrix/matrix/matrix product: PQ^T * diag(df2) * PQ
	std::vector<double> tmp(n * m);
	for (int j = 0; j < m; ++j) {
		for (int i = 0; i < n; ++i) {
			tmp[j * n + i] = 0.0;
			tmp[j * n + i] += PQ[i * m + j] * df2[i];
		}
	}

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			hess[i * m + j] = 0.0;
			for (int k = 0; k < n; k++) {
				hess[i * m + j] += tmp[i * n + k] * PQ[k * m + j];
			}
		}
	}

	double lamai = 0.0;
	for (int j = 0; j < m; ++j) {
		if (lam[j] < 0.0) {
			lam[j] = 0.0;
		}
		lamai += lam[j] * a[j];
		if (lam[j] > c[j]) {
			hess[j * m + j] += -1.0;
		}
		hess[j * m + j] += -mu[j] / lam[j];
	}

	if (lamai > 0.0) {
		for (int j = 0; j < m; ++j) {
			for (int k = 0; k < m; k++) {
				hess[j * m + k] += -10.0 * a[j] * a[k];
			}
		}
	}

	// pos def check
	double hessTrace = 0.0;
	for (int i = 0; i < m; ++i) {
		hessTrace += hess[i * m + i];
	}
	double hessCorr = 1e-4 * hessTrace / m;

	if (-1.0 * hessCorr < 1.0e-7) {
		hessCorr = -1.0e-7;
	}

	for (int i = 0; i < m; ++i) {
		hess[i * m + i] += hessCorr;
	}
}

void GCMMASolver::DualGrad(double *x) {
	for (int j = 0; j < m; ++j) {
		grad[j] = -b[j] - a[j] * z - y[j];
		for (int i = 0; i < n; ++i) {
			grad[j] += pij[i * m + j] / (upp[i] - x[i]) + qij[i * m + j] / (x[i] - low[i]);
		}
	}
}

void GCMMASolver::XYZofLAMBDA(double *x) {

	double lamai = 0.0;
	for (int i = 0; i < m; ++i) {
		if (lam[i] < 0.0) {
			lam[i] = 0;
		}
		y[i] = std::max(0.0, lam[i] - c[i]); // Note y=(lam-c)/d - however d is fixed at one !!
		lamai += lam[i] * a[i];
	}
	z = std::max(0.0, 10.0 * (lamai - 1.0)); // SINCE a0 = 1.0

	double pjlam, qjlam;
	for (int i = 0; i < n; ++i) {
		pjlam = p0[i];
		qjlam = q0[i];
		for (int j = 0; j < m; ++j) {
			pjlam += pij[i * m + j] * lam[j];
			qjlam += qij[i * m + j] * lam[j];
		}
		x[i] = (sqrt(pjlam) * low[i] + sqrt(qjlam) * upp[i]) / (sqrt(pjlam) + sqrt(qjlam));
		if (x[i] < alpha[i]) {
			x[i] = alpha[i];
		}
		if (x[i] > beta[i]) {
			x[i] = beta[i];
		}
	}
}

// Compute [low, upp, raa0, raa]
void GCMMASolver::Asymp(const double *xval, const double *df0dx, const double *dfdx, const double *xmin,
                        const double *xmax)
{
	// Forward the iterator
	++outeriter;

	// Set asymptotes
	if (outeriter < 3) {
		for (int i = 0; i < n; ++i) {
			low[i] = xval[i] - asyminit * (xmax[i] - xmin[i]);
			upp[i] = xval[i] + asyminit * (xmax[i] - xmin[i]);
		}
	} else {
		for (int i = 0; i < n; ++i) {
			double zzz = (xval[i] - xold1[i]) * (xold1[i] - xold2[i]);
			double gamma = 1.0;
			if (zzz < 0.0) {
				gamma = asymdec;
			} else if (zzz > 0.0) {
				gamma = asyminc;
			} else {
				gamma = 1.0;
			}
			low[i] = xval[i] - gamma * (xold1[i] - low[i]);
			upp[i] = xval[i] + gamma * (upp[i] - xold1[i]);

			double xmami = std::max(xmamieps, xmax[i] - xmin[i]);
			// double xmami = xmax[i] - xmin[i];
			low[i] = std::max(low[i], xval[i] - 100.0 * xmami);
			low[i] = std::min(low[i], xval[i] - 1e-5 * xmami);
			upp[i] = std::max(upp[i], xval[i] + 1e-5 * xmami);
			upp[i] = std::min(upp[i], xval[i] + 100.0 * xmami);

			double xmi = xmin[i] - 1.0e-6;
			double xma = xmax[i] + 1.0e-6;
			if (xval[i] < xmi) {
				low[i] = xval[i] - (xma - xval[i]) / 0.9;
				upp[i] = xval[i] + (xma - xval[i]) / 0.9;
			}
			if (xval[i] > xma) {
				low[i] = xval[i] - (xval[i] - xmi) / 0.9;
				upp[i] = xval[i] + (xval[i] - xmi) / 0.9;
			}
		}
	}

	// Set raa0, raa
	// raa0 = 0;
	// std::fill(raa.begin(), raa.end(), 0);
	// for (int i = 0; i < n; ++i) {
	// 	double xmami = std::max(xmamieps, xmax[i] - xmin[i]);
	// 	raa0 += std::abs(df0dx[i]) * xmami;
	// 	for (int j = 0; j < m; ++j) {
	// 		raa[j] += std::abs(dfdx[i*m+j])*xmami;
	// 	}
	// }
	raa0 = std::max(raa0eps, 0.1 * raa0); // std::max(raa0eps, (0.1/n)*raa0);
	for (int j = 0; j < m; ++j) {
		raa[j] = std::max(raaeps, 0.1 * raa[j]); // std::max(raaeps, (0.1/n)*raa[j]);
	}
}

// Update [raa0, raa]
void GCMMASolver::RaaUpdate(const double *xmma, const double *xval, double f0xnew, const double *fxnew,
                            const double *xmin, const double *xmax)
{
	const double raacofmin = 1e-12;
	double raacof = 0.0;
	for (int i = 0; i < n; ++i) {
		double xmami = std::max(xmamieps, xmax[i] - xmin[i]);
		double xxux = (xmma[i] - xval[i]) / (upp[i] - xmma[i]);
		double xxxl = (xmma[i] - xval[i]) / (xmma[i] - low[i]);
		double xxul = xxux * xxxl;
		double ulxx = (upp[i] - low[i]) / xmami;
		raacof += xxul * ulxx;
	}
	raacof = std::max(raacofmin, raacof);
	// std::cout << "raacof: " << raacof << std::endl;

	if (f0xnew > f0app + 0.5 * epsimin) {
		double deltaraa0 = (1.0 / raacof) * (f0xnew - f0app);
		raa0 = std::min(1.1 * (raa0 + deltaraa0), 10.0 * raa0);
	}
	// std::cout << "raa0: " << raa0 << std::endl << "raaj: ";
	for (int j = 0; j < m; ++j) {
		if (fxnew[j] > fapp[j] + 0.5 * epsimin) {
			double deltaraa = (1.0 / raacof) * (fxnew[j] - fapp[j]);
			raa[j] = std::min(1.1 * (raa[j] + deltaraa), 10.0 * raa[j]);
		}
		// std::cout << raa[j] << ' ';
	}
	// std::cout << std::endl;
}

void GCMMASolver::GenSub(const double *xval, double f0x, const double *df0dx, const double *fx, const double *dfdx,
                         const double *xmin, const double *xmax)
{
	// Set bounds and the coefficients for the approximation
	r0 = 0;
	std::fill(r.begin(), r.end(), 0);
	for (int i = 0; i < n; ++i) {
		// Compute bounds alpha and beta
		alpha[i] = std::max(xmin[i], low[i] + albefa * (xval[i] - low[i]));
		alpha[i] = std::max(alpha[i], xval[i] - move * (xmax[i] - xmin[i]));
		// alpha[i] = std::min(alpha[i], xmax[i]);
		beta[i] = std::min(xmax[i], upp[i] - albefa * (upp[i] - xval[i]));
		beta[i] = std::min(beta[i], xval[i] + move * (xmax[i] - xmin[i]));
		// beta[i]  = std::max(beta[i], xmin[i]);

		double uxinv = 1.0 / (upp[i] - xval[i]);
		double xlinv = 1.0 / (xval[i] - low[i]);

		// Objective function
		{
			double df0dxp = std::max(0.0, df0dx[i]);
			double df0dxm = std::max(0.0, -1.0 * df0dx[i]);
			double xmamiinv = 1.0 / std::max(xmamieps, xmax[i] - xmin[i]);
			double pq = 0.001 * std::abs(df0dx[i]) + raa0 * xmamiinv;
			p0[i] = std::pow(upp[i] - xval[i], 2.0) * (df0dxp + pq);
			q0[i] = std::pow(xval[i] - low[i], 2.0) * (df0dxm + pq);
			r0 += p0[i] * uxinv + q0[i] * xlinv;
		}

		// Constraints
		for (int j = 0; j < m; ++j) {
			double dfdxp = std::max(0.0, dfdx[i * m + j]);
			double dfdxm = std::max(0.0, -1.0 * dfdx[i * m + j]);
			double xmamiinv = 1.0 / std::max(xmamieps, xmax[i] - xmin[i]);
			double pq = 0.001 * std::abs(dfdx[i * m + j]) + raa[j] * xmamiinv;
			pij[i * m + j] = std::pow(upp[i] - xval[i], 2.0) * (dfdxp + pq);
			qij[i * m + j] = std::pow(xval[i] - low[i], 2.0) * (dfdxm + pq);
			r[j] += pij[i * m + j] * uxinv + qij[i * m + j] * xlinv;
		}
	}

	r0 = f0x - r0;
	for (int j = 0; j < m; ++j) {
		r[j] = fx[j] - r[j];
		b[j] = -r[j]; // The constant for the constraints
	}
}

void GCMMASolver::ComputeApprox(const double *xmma) {
	f0app = 0;
	std::fill(fapp.begin(), fapp.end(), 0);
	for (int i = 0; i < n; ++i) {
		double uxinv = 1.0 / (upp[i] - xmma[i]);
		double xlinv = 1.0 / (xmma[i] - low[i]);
		f0app += p0[i] * uxinv + q0[i] * xlinv;
		for (int j = 0; j < m; ++j) {
			fapp[j] += pij[i * m + j] * uxinv + qij[i * m + j] * xlinv;
		}
	}
	f0app += r0;
	// std::cout << "f0: " << f0app << std::endl;
	for (int j = 0; j < m; ++j) {
		fapp[j] += r[j];
		// std::cout << "fj: " << fapp[j] << ' ';
	}
	// std::cout << std::endl;
}

void GCMMASolver::Factorize(double *K, int n) {

	for (int s = 0; s < n - 1; ++s) {
		for (int i = s + 1; i < n; ++i) {
			K[i * n + s] = K[i * n + s] / K[s * n + s];
			for (int j = s + 1; j < n; ++j) {
				K[i * n + j] = K[i * n + j] - K[i * n + s] * K[s * n + j];
			}
		}
	}
}

void GCMMASolver::Solve(double *K, double *x, int n) {

	for (int i = 1; i < n; ++i) {
		double a = 0.0;
		for (int j = 0; j < i; ++j) {
			a = a - K[i * n + j] * x[j];
		}
		x[i] = x[i] + a;
	}

	x[n - 1] = x[n - 1] / K[(n - 1) * n + (n - 1)];
	for (int i = n - 2; i >= 0; --i) {
		double a = x[i];
		for (int j = i + 1; j < n; ++j) {
			a = a - K[i * n + j] * x[j];
		}
		x[i] = a / K[i * n + i];
	}
}
