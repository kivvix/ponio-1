// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <cmath>
#include <iostream>
#include <numbers>
#include <numeric>

#include <algorithm>
#include <array>
#include <concepts>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <tuple>
#include <valarray>

namespace fs = std::filesystem;

#include <solver/butcher_methods.hpp>
#include <solver/detail.hpp>
#include <solver/observer.hpp>
#include <solver/problem.hpp>
#include <solver/solver.hpp>

#include <samurai/algorithm.hpp>
#include <samurai/bc.hpp>
#include <samurai/field.hpp>
#include <samurai/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/stencil_field.hpp>
#include <samurai/subset/subset_op.hpp>

#include <cmath>

/*
solve advection problem:

$$
    \begin{aligned}
        u_t + a u_x = 0 \\
        y(0) = u_0
    \end{aligned}
$$
*/

auto init = []( auto& mesh )
{
    auto u = samurai::make_field<double, 1>( "u", mesh );
    u.fill( 0. );

    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            auto center         = cell.center();
            const double radius = .5;

            const double x_center = 0;
            u[cell]               = std::cos( center[0] * 2. * 3.1415926535 / 4.0 );

            // if ( std::abs( center[0] - x_center ) <= radius )
            // {
            //     u[cell] = 1;
            // }

            // if ( -1 < center[0] && center[0] < 0 )
            // {
            //     u[cell] = center[0] + 1;
            // }
            // else if ( 0 < center[0] && center[0] < 1 )
            // {
            //     u[cell] = -center[0] + 1;
            // }
        } );

    return u;
};

// int
// main()
// {
//     constexpr std::size_t dim = 1;

//     using Config = samurai::MRConfig<dim>;

//     // Multiresolution parameters
//     std::size_t min_level = 4;
//     std::size_t max_level = 10;

//     const samurai::Box<double, dim> box( { -1 }, { 1 } );
//     samurai::MRMesh<Config> mesh( box, min_level, max_level, { true } );

//     auto u = init( mesh );

//     auto a = detail::init_fill_array<5>( init, mesh );

//     return 0;
// }

template <class Field>
void
save( fs::path const& path, std::string const& filename, Field const& u, std::string const& suffix = "" )
{
    auto& mesh  = u.mesh();
    auto level_ = samurai::make_field<std::size_t, 1>( "level", mesh );

    if ( !fs::exists( path ) )
    {
        fs::create_directory( path );
    }

    samurai::for_each_cell( mesh,
        [&]( auto const& cell )
        {
            level_[cell] = cell.level;
        } );

    samurai::save( path, fmt::format( "{}{}", filename, suffix ), mesh, u, level_ );
}

int
main()
{
    constexpr std::size_t dim = 1;
    using Config              = samurai::MRConfig<dim>;

    // Simulation parameters
    double left_box  = -2;
    double right_box = 2;
    bool is_periodic = true;

    // Multiresolution parameters
    std::size_t min_level = 2;
    std::size_t max_level = 6;
    double mr_epsilon     = 2.e-4; // Threshold used by multiresolution
    double mr_regularity  = 1.;    // Regularity guess for multiresolution

    const samurai::Box<double, dim> box( { left_box }, { right_box } );
    samurai::MRMesh<Config> mesh( box, min_level, max_level, { is_periodic } );

    auto u0 = init( mesh );

    std::cout << std::boolalpha;
    std::cout << "std::ranges::range<field>: " << std::ranges::range<decltype( u0 )> << "\n";

    return 0;
}

/*
int
main()
{
    std::string dirname = "trp_data";

    constexpr std::size_t dim = 1;
    using Config              = samurai::MRConfig<dim>;

    // Simulation parameters
    double left_box  = -2;
    double right_box = 2;
    bool is_periodic = true;

    // Multiresolution parameters
    std::size_t min_level = 2;
    std::size_t max_level = 6;
    double mr_epsilon     = 2.e-4; // Threshold used by multiresolution
    double mr_regularity  = 1.;    // Regularity guess for multiresolution

    // Output parameters
    fs::path path = fs::current_path() / "transport_data";
    fs::create_directories( path );
    std::string filename = "transport1d_";

    const samurai::Box<double, dim> box( { left_box }, { right_box } );
    samurai::MRMesh<Config> mesh( box, min_level, max_level, { is_periodic } );

    double a   = 1.;
    double Tf  = 4.0;
    double cfl = 0.5;
    double dt  = cfl * samurai::cell_length( max_level );

    auto u0 = init( mesh );
    //  auto unp1 = samurai::make_field<double, 1>( "unp1", mesh );

    using state_t  = decltype( u0 );
    auto advection = [=]( double, auto&& u ) -> auto
    {
        // auto tmp = u;
        // tmp      = -samurai::upwind( a, u );
        // return tmp;

        samurai::update_ghost_mr( u );
        return -samurai::upwind( a, u );
    };

    auto sol_range = ode::make_solver_range( advection,
        ode::butcher::rk_44(),
        u0,
        { 0., Tf },
        dt,
        [&]()
        {
            return init( mesh );
        } );

    auto it           = sol_range.begin();
    auto MRadaptation = samurai::make_MRAdapt( it->state );
    MRadaptation( mr_epsilon, mr_regularity );
    for ( auto& ki : it.meth.kis )
    {
        ki.resize();
        ki.fill( 0 );
    }

    // save( path, filename, it->state, "_init" );

    // std::size_t count = 0;
    // for ( auto& ui : sol_range )
    // // while ( it->time <= Tf )
    // {
    //     std::cout << "t : " << it->time << " " << it->time_step << std::flush << "\r";
    //     // MRadaptation( mr_epsilon, mr_regularity );

    //     // std::cout << "\n" << it->state << std::endl;

    //     // if ( count % 10 == 0 )
    //     {
    //         std::string suffix = fmt::format( "_ite_{}", count );
    //         ui.state.name()    = "u";
    //         save( path, filename, ui.state, suffix );
    //     }
    //     // std::cout << "it->state: " << it->state << "\n";
    //     samurai::update_ghost_mr( ui.state );

    //     ++count;
    //     ++it;
    // }
    // std::cout << std::endl;

    std::size_t count = 0;
    while ( it->time <= Tf )
    {
        // std::cout << "\n" << it->state << std::endl;

        // if ( count % 10 == 0 )
        {
            std::string suffix = fmt::format( "_ite_{}", count );
            it->state.name()   = "u";
            save( path, filename, it->state, suffix );
        }
        // std::cout << "it->state: " << it->state << "\n";

        std::cout << "t : " << it->time << " " << it->time_step << std::flush << "\r";
        ++count;
        ++it;

        MRadaptation( mr_epsilon, mr_regularity );
        samurai::update_ghost_mr( it->state );
        for ( auto& ki : it.meth.kis )
        {
            ki.resize();
            ki.fill( 0 );
        }
    }
    std::cout << std::endl;

    return 0;
}
*/
