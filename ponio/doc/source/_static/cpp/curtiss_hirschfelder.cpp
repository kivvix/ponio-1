#include <cmath>
#include <filesystem>
#include <string>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/splitting.hpp>
#include <ponio/time_span.hpp>

/*
solve Curtiss and Hirschfelder problem:

$$
    \begin{aligned}
        \dot{y} =  k(\cos(t) - y) \\
        y(0) = y_0
    \end{aligned}
$$

$y_0 = 2$
*/

int
main()
{
    using namespace ponio::observer;

    double const k   = 50.;
    double const y_0 = 2.;

    ponio::time_span<double> t_span = { 0., 4.0 };
    double dt                       = 0.05;

    {
    // explicit Runge-Kutta method
#include "curtiss_hirschfelder_all/erk.cxx"
    }

    {
    // embedded Runge-Kutta method
#include "curtiss_hirschfelder_all/dp.cxx"
    }

    {
    // diagonal implicit Runge-Kutta method
#include "curtiss_hirschfelder_all/dirk.cxx"
    }

    {
    // Lawson method
#include "curtiss_hirschfelder_all/lrk.cxx"
    }

    {
    // exponential Runge-Kutta method
#include "curtiss_hirschfelder_all/exprk.cxx"
    }

    {
    // Runge-Kutta Chebyshev method
#include "curtiss_hirschfelder_all/rkc.cxx"
    }

    {
    // ROCK2 method
#include "curtiss_hirschfelder_all/rock2.cxx"
    }

    {
    // ROCK4 method
#include "curtiss_hirschfelder_all/rock4.cxx"
    }

    {
    // Runge-Kutta Legendre first-order method
#include "curtiss_hirschfelder_all/rkl1.cxx"
    }

    {
    // Runge-Kutta Legendre second-order method
#include "curtiss_hirschfelder_all/rkl2.cxx"
    }

    {
    // Lie splitting method
#include "curtiss_hirschfelder_all/lie.cxx"
    }

    {
    // Strang splitting method
#include "curtiss_hirschfelder_all/strang.cxx"
    }

    {
        // adaptive Strang splitting method
        // #include "curtiss_hirschfelder_all/adaptive_strang.cxx"
    }

    {
    // pirock
#include "curtiss_hirschfelder_all/pirock.cxx"
    }

    { // additive Runge-Kutta
#include "curtiss_hirschfelder_all/ark.cxx"
    }
}
