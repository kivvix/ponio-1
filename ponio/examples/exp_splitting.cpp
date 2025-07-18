// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <cmath>
#include <filesystem>
#include <string>
#include <tuple>
#include <utility>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>
#include <ponio/splitting.hpp>
#include <ponio/user_defined_method.hpp>

// solve $\dot{y} = y$ with $y(t=0) = 1$, and $t\in[0,2]$.

int
main()
{
    std::string const dirname = "exp_splitting_data";

    double const lambda = 0.3;

    auto f1 = [=]( double, double y )
    {
        return lambda * y;
    };

    auto f2 = [=]( double, double y )
    {
        return ( 1.0 - lambda ) * y;
    };

    auto pb = ponio::make_problem( ponio::make_simple_problem( f1 ), ponio::make_simple_problem( f2 ) );

    auto exact_solver_f1 = [=]( auto /* f */, double tn, double yn, double dt ) -> std::tuple<double, double, double>
    {
        return { tn + dt, std::exp( lambda * dt ) * yn, dt };
    };

    auto exact_solver = [=]( auto /* f */, double tn, double yn, double dt ) -> std::tuple<double, double, double>
    {
        return { tn + dt, std::exp( dt ) * yn, dt };
    };

    double const y0 = 1.0;
    double const dt = 0.5;

    // 3 examples
    {
        auto strang = ponio::splitting::make_strang_tuple( std::make_pair( ponio::make_user_defined_method( exact_solver_f1 ), dt ),
            std::make_pair( ponio::runge_kutta::rk_22_ralston(), 0.5 * dt ) );

        auto filename = std::filesystem::path( dirname ) / "exp_strang.dat";
        ponio::observer::file_observer fobs( filename );
        ponio::solve( pb, strang, y0, { 0., 5.0 }, dt, fobs );
    }
    {
        auto exact = ponio::make_user_defined_method( exact_solver );

        auto filename = std::filesystem::path( dirname ) / "exp_exact.dat";
        ponio::observer::file_observer fobs( filename );
        ponio::solve( pb, exact, y0, { 0., 5.0 }, dt, fobs );
    }
    {
        auto rk2 = ponio::runge_kutta::rk_22_ralston();

        auto filename = std::filesystem::path( dirname ) / "exp_rk2.dat";
        ponio::observer::file_observer fobs( filename );
        ponio::solve( pb, rk2, y0, { 0., 5.0 }, dt, fobs );
    }

    return 0;
}
