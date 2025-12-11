#include <iostream>
#include <fstream>
#include <cmath>
#include "nr.h"

using namespace std;

// Global variables for odeint
int kmax, kount;
DP dxsav;
Vec_DP* xp_p;
Mat_DP* yp_p;

// derivs function for the system dx/dt = v, dv/dt = -0.2v + x - 0.1x^3 + cos(t)
void derivs(const DP x, Vec_I_DP& y, Vec_O_DP& dydx)
{
    // y[0] = x (position), y[1] = v (velocity)
    // x is the independent variable (time t)
    dydx[0] = y[1]; // dx/dt = v
    dydx[1] = -0.2 * y[1] + y[0] - 0.1 * y[0] * y[0] * y[0] + cos(x); // dv/dt = -0.2v + x - 0.1x^3 + cos(t)
}

// rkck function (Runge-Kutta-Cash-Karp stepper)
void NR::rkck(Vec_I_DP& y, Vec_I_DP& dydx, const DP x,
    const DP h, Vec_O_DP& yout, Vec_O_DP& yerr,
    void derivs(const DP, Vec_I_DP&, Vec_O_DP&))
{
    static const DP a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875,
        b21 = 0.2, b31 = 3.0 / 40.0, b32 = 9.0 / 40.0, b41 = 0.3, b42 = -0.9,
        b43 = 1.2, b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0,
        b54 = 35.0 / 27.0, b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0,
        b63 = 575.0 / 13824.0, b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0,
        c1 = 37.0 / 378.0, c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0,
        dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0,
        dc4 = c4 - 13525.0 / 55296.0, dc5 = -277.00 / 14336.0, dc6 = c6 - 0.25;
    int i;

    int n = y.size();
    Vec_DP ak2(n), ak3(n), ak4(n), ak5(n), ak6(n), ytemp(n);
    for (i = 0; i < n; i++)
        ytemp[i] = y[i] + b21 * h * dydx[i];
    derivs(x + a2 * h, ytemp, ak2);
    for (i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i]);
    derivs(x + a3 * h, ytemp, ak3);
    for (i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
    derivs(x + a4 * h, ytemp, ak4);
    for (i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
    derivs(x + a5 * h, ytemp, ak5);
    for (i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);
    derivs(x + a6 * h, ytemp, ak6);
    for (i = 0; i < n; i++)
        yout[i] = y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
    for (i = 0; i < n; i++)
        yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);
}


// rkqs function (adaptive stepsize control)
void NR::rkqs(Vec_IO_DP& y, Vec_IO_DP& dydx, DP& x, const DP htry,
    const DP eps, Vec_I_DP& yscal, DP& hdid, DP& hnext,
    void derivs(const DP, Vec_I_DP&, Vec_O_DP&))
{
    const DP SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89e-4;
    int i;
    DP errmax, h, htemp, xnew;

    int n = y.size();
    h = htry;
    Vec_DP yerr(n), ytemp(n);
    for (;;) {
        rkck(y, dydx, x, h, ytemp, yerr, derivs);
        errmax = 0.0;
        for (i = 0; i < n; i++) errmax = MAX(errmax, fabs(yerr[i] / yscal[i]));
        errmax /= eps;
        if (errmax <= 1.0) break;
        htemp = SAFETY * h * pow(errmax, PSHRNK);
        h = (h >= 0.0 ? MAX(htemp, 0.1 * h) : MIN(htemp, 0.1 * h));
        xnew = x + h;
        if (xnew == x) nrerror("stepsize underflow in rkqs");
    }
    if (errmax > ERRCON) hnext = SAFETY * h * pow(errmax, PGROW);
    else hnext = 5.0 * h;
    x += (hdid = h);
    for (i = 0; i < n; i++) y[i] = ytemp[i];
}



// odeint function (Runge-Kutta driver)
void NR::odeint(Vec_IO_DP& ystart, const DP x1, const DP x2, const DP eps,
    const DP h1, const DP hmin, int& nok, int& nbad,
    void derivs(const DP, Vec_I_DP&, Vec_O_DP&),
    void rkqs(Vec_IO_DP&, Vec_IO_DP&, DP&, const DP, const DP,
        Vec_I_DP&, DP&, DP&, void (*)(const DP, Vec_I_DP&, Vec_O_DP&)))
{
    const int MAXSTP = 10000;
    const DP TINY = 1.0e-30;
    int i, nstp;
    DP xsav, x, hnext, hdid, h;

    int nvar = ystart.size();
    Vec_DP yscal(nvar), y(nvar), dydx(nvar);
    Vec_DP& xp = *xp_p;
    Mat_DP& yp = *yp_p;
    x = x1;
    h = SIGN(h1, x2 - x1);
    nok = nbad = kount = 0;
    for (i = 0; i < nvar; i++) y[i] = ystart[i];
    if (kmax > 0) xsav = x - dxsav * 2.0;
    for (nstp = 0; nstp < MAXSTP; nstp++) {
        derivs(x, y, dydx);
        for (i = 0; i < nvar; i++)
            yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + TINY;
        if (kmax > 0 && kount < kmax - 1 && fabs(x - xsav) > fabs(dxsav)) {
            for (i = 0; i < nvar; i++) yp[i][kount] = y[i];
            xp[kount++] = x;
            xsav = x;
        }
        if ((x + h - x2) * (x + h - x1) > 0.0) h = x2 - x;
        rkqs(y, dydx, x, h, eps, yscal, hdid, hnext, derivs);
        if (hdid == h) ++nok; else ++nbad;
        if ((x - x2) * (x2 - x1) >= 0.0) {
            for (i = 0; i < nvar; i++) ystart[i] = y[i];
            if (kmax != 0) {
                for (i = 0; i < nvar; i++) yp[i][kount] = y[i];
                xp[kount++] = x;
            }
            return;
        }
        if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
        h = hnext;
    }
    nrerror("Too many steps in routine odeint");
}





int main()
{
    const int nvar = 2; // Two variables: x and v
    Vec_DP ystart(nvar);
    ystart[0] = 0.0; // Initial x(0) = 0
    ystart[1] = 0.0; // Initial v(0) = 0
    DP x1 = 0.0, x2 = 400.0; // Integrate from t=0 to t=400
    DP eps = 1.0e-4, h1 = 0.1, hmin = 0.0; // Integration parameters
    int nok, nbad;

    // Set up storage for intermediate results
    kmax = 1000; // Maximum number of stored points
    dxsav = 0.1; // Save points every 0.1 time units
    Vec_DP xp_data(kmax);
    Mat_DP yp_data(nvar, kmax);
    xp_p = &xp_data;
    yp_p = &yp_data;

    // Integrate the system
    try {
        NR::odeint(ystart, x1, x2, eps, h1, hmin, nok, nbad, derivs, NR::rkqs);
    }
    catch (const std::exception& e) {
        std::cerr << "Error during integration: " << e.what() << std::endl;
        return 1;
    }

    // Output results to a CSV file
    std::ofstream out_file("results.csv");
    if (!out_file) {
        std::cerr << "Error opening output file" << std::endl;
        return 1;
    }
    out_file << "t,x,v\n"; // CSV header
    for (int i = 0; i < kount; i++) {
        out_file << xp_data[i] << "," << yp_data[0][i] << "," << yp_data[1][i] << "\n";
    }
    out_file.close();

    // Print integration summary
    std::cout << "Integration complete. Good steps: " << nok << ", Bad steps: " << nbad << std::endl;
    std::cout << "Final x: " << ystart[0] << ", v: " << ystart[1] << std::endl;
    std::cout << "Results saved to results.csv" << std::endl;
    

    return 0;
}