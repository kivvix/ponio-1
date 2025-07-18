// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// NOLINTBEGIN(misc-include-cleaner)

#include <filesystem>
#include <string>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/solver.hpp>

// NOLINTEND(misc-include-cleaner)

// solve $\dot{u} = u$ with $u(t=0) = 1$, and $t\in[0,2]$.

int
main()
{
    std::string const dirname = "exp_data";
    auto filename             = std::filesystem::path( dirname ) / "exp.dat";
    ponio::observer::file_observer fobs( filename );

    auto identity = ponio::make_simple_problem(
        []( double, double u )
        {
            return u;
        } );
    double const x0 = 1.0;
    double const dt = 0.1;

    ponio::solve( identity, ponio::runge_kutta::rk_nssp_21(), x0, { 0., 2.0 }, dt, fobs );

    return 0;
}
