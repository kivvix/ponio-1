// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <iostream>
#include <valarray>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <numbers>
#include <sstream>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/runge_kutta/pirock.hpp>
#include <ponio/samurai_linear_algebra.hpp>
#include <ponio/solver.hpp>
#include <ponio/time_span.hpp>

#include <samurai/field.hpp>
#include <samurai/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>
#include <samurai/schemes/fv.hpp>

#include <filesystem>
namespace fs = std::filesystem;

template <class field_t>
void
save( fs::path const& path, std::string const& filename, field_t& u, std::string const& suffix = "" )
{
    auto mesh   = u.mesh();
    auto level_ = samurai::make_field<std::size_t, 1>( "level", mesh );
    u.name()    = "u";

    if ( !fs::exists( path ) )
    {
        fs::create_directory( path );
    }

    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            level_[cell] = cell.level;
        } );

    samurai::save( path, fmt::format( "{}{}", filename, suffix ), mesh, u, level_ );
}

int
main( int argc, char** argv )
{
    PetscInitialize( &argc, &argv, 0, nullptr );

    constexpr std::size_t dim = 2; // cppcheck-suppress unreadVariable
    using config_t            = samurai::MRConfig<dim>;
    using box_t               = samurai::Box<double, dim>;
    using point_t             = typename box_t::point_t;

    // simulation parameters --------------------------------------------------
    constexpr double Da = 2.5e-3;
    constexpr double Db = 2.5e-3;
    constexpr double Dc = 1.5e-3;

    constexpr double epsilon = 1e-2;
    constexpr double mu      = 1e-5;
    constexpr double f       = 1.6;
    constexpr double q       = 2e-3;

    constexpr double left_box  = 0.;
    constexpr double right_box = 1.;
    constexpr double t_ini     = 0.;
    constexpr double t_end     = 2.;

    // multiresolution parameters
    std::size_t min_level = 2;
    std::size_t max_level = 6;
    double mr_epsilon     = 1e-2; // Threshold used by multiresolution
    double mr_regularity  = 0.;   // Regularity guess for multiresolution

    // output parameters
    std::string const dirname = "belousov_zhabotinsky_2d_3eq_pirock_data";
    fs::path path             = std::filesystem::path( dirname );
    std::string filename      = "u";
    fs::create_directories( path );

    // define mesh
    point_t box_corner1, box_corner2;
    box_corner1.fill( left_box );
    box_corner2.fill( right_box );
    box_t box( box_corner1, box_corner2 );
    samurai::MRMesh<config_t> mesh{ box, min_level, max_level };

    // init solution ----------------------------------------------------------
    // auto u_ini = init( mesh );
    auto u_ini = samurai::make_field<3>( "u", mesh );

    u_ini.fill( 0 );
    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            auto& a = u_ini[cell]( 0 );
            auto& b = u_ini[cell]( 1 );
            auto& c = u_ini[cell]( 2 );

            double x = cell.center()[0] - 0.5;
            double y = cell.center()[1] - 0.5;

            double theta = std::atan( y / x );

            if ( x < 0. )
            {
                theta = std::numbers::pi + theta;
            }
            else if ( y < 0. )
            {
                theta = 2 * std::numbers::pi + theta;
            }

            double bss = q * ( ( f + 1. ) / ( f - 1. ) );
            double css = bss;

            if ( theta < 0.5 )
            {
                b = 0.8;
            }
            else
            {
                b = bss;
            }
            c = css + theta / ( 8. * std::numbers::pi * f );

            a = ( f * c ) / ( q + b );
        } );
    samurai::make_bc<samurai::Neumann<1>>( u_ini, 0., 0., 0. );

    // define problem ---------------------------------------------------------

    // diffusion terme
    auto d    = samurai::DiffCoeff<3>( { Da, Db, Dc } );
    auto diff = samurai::make_multi_diffusion_order2<decltype( u_ini )>( d );
    auto fd   = [&]( double /* t */, auto&& u )
    {
        samurai::make_bc<samurai::Neumann<1>>( u, 0., 0., 0. );
        samurai::update_ghost_mr( u );
        return -diff( u );
    };

    // reaction terme
    using cfg  = samurai::LocalCellSchemeConfig<samurai::SchemeType::NonLinear, decltype( u_ini )::size, decltype( u_ini )>;
    auto react = samurai::make_cell_based_scheme<cfg>();
    react.set_name( "Reaction" );
    react.set_scheme_function(
        [&]( auto const& cell, auto const& field ) -> samurai::SchemeValue<cfg>
        {
            auto u  = field[cell];
            auto& a = u[0];
            auto& b = u[1];
            auto& c = u[2];

            // clang-format off
            return {
                1. / mu * ( -q * a - a * b + f * c ),
                1. / epsilon * ( q * a - a * b + b * ( 1. - b ) ),
                b - c
            };
            // clang-format on
        } );
    // or set option in command line with : -snes_fd -pc_type none
    react.set_jacobian_function(
        [&]( auto const& cell, auto const& field ) -> samurai::JacobianMatrix<cfg>
        {
            auto u  = field[cell];
            auto& a = u[0];
            auto& b = u[1];
            // auto& c = u[2];

            return {
                {( -q - b ) / mu,      -a / mu,                           f / mu},
                { ( q - b ) / epsilon, 1. / epsilon * ( -a - 2 * b + 1 ), 0.    },
                { 0.,                  1.,                                -1.   }
            };
        } );
    auto fr_t = [&]( double /* t */ )
    {
        return react;
    };
    auto fr = [&]( double t, auto&& u )
    {
        u.get_bc().clear();
        samurai::make_bc<samurai::Neumann<1>>( u, 0., 0., 0. );
        samurai::update_ghost_mr( u );
        return fr_t( t )( u );
    };

    ponio::time_span<double> const t_span = { t_ini, t_end };
    double dt                             = ( t_end - t_ini ) / 4000;

    auto eigmax_computer = [=]( auto&, double, auto&, double )
    {
        double dx = samurai::cell_length( max_level );
        // return 10. * Da * 8. / ( dx * dx );
        return 4. / ( dx * dx );
    };

    auto pb = ponio::make_imex_operator_problem( fd, fr, fr_t );

    // time loop  -------------------------------------------------------------
    static constexpr bool is_embedded = false;

    auto pirock_b0    = ponio::runge_kutta::pirock::pirock<1>( ponio::runge_kutta::pirock::beta_0<double>(), eigmax_computer );
    auto pirock_a1    = ponio::runge_kutta::pirock::pirock<2>( ponio::runge_kutta::pirock::alpha_fixed<double>( 1.0 ), eigmax_computer );
    auto pirock_b0_st = ponio::runge_kutta::pirock::pirock<1, is_embedded>( ponio::runge_kutta::pirock::beta_0<double>(),
        eigmax_computer,
        ponio::shampine_trick::shampine_trick<decltype( u_ini )>(),
        1e-3 );
    auto pirock_a1_st = ponio::runge_kutta::pirock::pirock<2, is_embedded>( ponio::runge_kutta::pirock::alpha_fixed<double>( 1.0 ),
        eigmax_computer,
        ponio::shampine_trick::shampine_trick<decltype( u_ini )>(),
        1e-3 );
    auto rock4        = ponio::runge_kutta::rock::rock4<is_embedded>( eigmax_computer );

    // auto sol_range = ponio::make_solver_range( pb, rock4, u_ini, t_span, dt );
    auto sol_range = ponio::make_solver_range( pb, pirock_b0_st, u_ini, t_span, dt );

    auto it_sol = sol_range.begin();

    // preapre MR for solution on iterator
    auto mr_adaptation = samurai::make_MRAdapt( it_sol->state );
    samurai::make_bc<samurai::Neumann<1>>( it_sol->state, 0., 0., 0. );
    mr_adaptation( mr_epsilon, mr_regularity );
    samurai::update_ghost_mr( it_sol->state );

    std::size_t n_save = 0;
    save( path, filename, it_sol->state, "_initial_condition" );
    save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save++ ) );

    while ( it_sol->time < t_end )
    {
        it_sol->state.get_bc().clear();
        samurai::make_bc<samurai::Neumann<1>>( it_sol->state, 0., 0., 0. );
        //  TODO: add a callback function to make this before each iteration
        for ( auto& ki : it_sol.meth.kis )
        {
            ki.resize();
            ki.fill( 0. );
        }

        ++it_sol;
        std::cout << "tⁿ: " << std::setw( 8 ) << it_sol->time << " (Δt: " << it_sol->time_step << ") " << n_save << std::endl;

        mr_adaptation( mr_epsilon, mr_regularity );
        samurai::update_ghost_mr( it_sol->state );

        save( path, filename, it_sol->state, fmt::format( "_ite_{}", n_save++ ) );
    }
    std::cout << std::endl;

    PetscFinalize();

    return 0;
}
