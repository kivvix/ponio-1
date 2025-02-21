// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <samurai/field.hpp>
#include <samurai/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/schemes/fv.hpp>

#include <ponio/observer.hpp>
#include <ponio/problem.hpp>
#include <ponio/runge_kutta.hpp>
#include <ponio/samurai_linear_algebra.hpp>
#include <ponio/solver.hpp>

#include <filesystem>

template <class field1_t, class field2_t>
void
save( std::filesystem::path const& path, std::string const& filename, field1_t& u, field2_t& v, std::string const& suffix = "" )
{
    auto mesh   = u.mesh();
    auto level_ = samurai::make_field<std::size_t, 1>( "level", mesh );
    u.name()    = "c";

    if ( !std::filesystem::exists( path ) )
    {
        std::filesystem::create_directory( path );
    }

    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            level_[cell] = cell.level;
        } );

    samurai::save( path, fmt::format( "{}{}", filename, suffix ), mesh, u, v, level_ );
}

int
main( int argc, char** argv )
{
    PetscInitialize( &argc, &argv, 0, nullptr );

    constexpr std::size_t dim = 2;
    using config_t            = samurai::MRConfig<dim, 3>;
    using box_t               = samurai::Box<double, dim>;
    using point_t             = typename box_t::point_t;

    // simulation parameters --------------------------------------------------
    double const L     = 1.0;
    double const D     = 1.0; // ???????????????????
    double const tau_d = L * L / D;

    double const Da      = 2.5e9;
    double const Ta      = 20000;
    double const To      = 300;
    double const tau     = 6.72;
    double const Re      = 1000;
    double const Sc      = 1.0;
    double const Delta_s = 0.02;

    double const left_box  = -L;
    double const right_box = L;

    auto x_s = [=]( double x )
    {
        return x / L;
    };
    auto y_s = [=]( double y )
    {
        return y / L;
    };
    auto t_s = [=]( double t )
    {
        return t / tau_d;
    };

    double y_0s  = -0.5;
    double x_0sp = 0.25;
    double x_0sm = -x_0sp;

    double const t_ini = 0.;
    double const t_end = 0.4e-3 * tau_d;

    // multiresolution parameters
    std::size_t const min_level = 2;
    std::size_t const max_level = 5;
    double const mr_epsilon     = 1e-3; // Threshold used by multiresolution
    double const mr_regularity  = 1.;   // Regularity guess for multiresolution

    // output parameters
    std::string const dirname  = std::string( argv[0] ) + std::string( "_data" );
    std::filesystem::path path = std::filesystem::path( dirname );
    std::string filename       = "sol";

    // define mesh
    point_t box_corner1, box_corner2;
    box_corner1.fill( left_box );
    box_corner2.fill( right_box );
    box_t box( box_corner1, box_corner2 );
    std::array<bool, dim> periodic = { false, false };
    samurai::MRMesh<config_t> mesh{ box, min_level, max_level, periodic };

    // init solution ----------------------------------------------------------

    auto c_init = samurai::make_field<double, 1>( "c", mesh );
    c_init.fill( 0. );
    samurai::for_each_cell( mesh,
        [&]( auto& cell )
        {
            // double x = cell.center( 0 );
            double y = cell.center( 1 );

            if ( y_s( y ) <= y_0s )
            {
                c_init[cell] = std::exp( ( y_s( y ) - y_0s ) / Delta_s );
            }
            else
            {
                c_init[cell] = 1.;
            }
        } );

    using field_t = decltype( c_init );

    // define problem ---------------------------------------------------------

    // diffusion terme: (d_xx ∙ + d_yy ∙)c
    samurai::DiffCoeff<dim> diff_coeff;
    diff_coeff.fill( 1.0 );
    auto diff = samurai::make_diffusion_order2<field_t>( diff_coeff );
    auto fd   = [&]( double /* t */, auto&& c_ )
    {
        samurai::update_ghost_mr( c_ );
        return -diff( c_ );
    };

    // reaction terme: Da (1-c) exp( - Ta/(To(1+tau c)) )
    using cfg_react = samurai::LocalCellSchemeConfig<samurai::SchemeType::NonLinear, field_t::size, field_t>;
    auto react      = samurai::make_cell_based_scheme<cfg_react>();
    react.set_name( "source" );
    react.set_scheme_function(
        [&]( auto const& cell, auto const& c_ ) -> samurai::SchemeValue<cfg_react>
        {
            return Da * ( 1 - c_[cell] ) * std::exp( -Ta / ( To * ( 1 + tau * c_[cell] ) ) );
        } );
    // or set option in command line with : -snes_fd -pc_type none
    react.set_jacobian_function(
        [&]( auto const& cell, auto const& c_ ) -> samurai::JacobianMatrix<cfg_react>
        {
            return Da * ( -Ta * tau * ( c_[cell] - 1 ) - To * ( c_[cell] * tau + 1 ) * ( c_[cell] * tau + 1 ) )
                 * exp( -Ta / ( To * ( c_[cell] * tau + 1 ) ) ) / ( ( To * ( c_[cell] * tau + 1 ) ) * ( To * ( c_[cell] * tau + 1 ) ) );
        } );
    auto fr_t = [&]( double /* t */ )
    {
        return react;
    };
    auto fr = [&]( double t, auto&& c_ )
    {
        // samurai::update_ghost_mr( c_ );
        return fr_t( t )( c_ );
    };

    // advection terme: (v_x d_x ∙ + v_y d_y ∙)c
    auto velocity        = samurai::make_field<double, dim>( "velocity", mesh );
    auto update_velocity = [&]( auto& vel, double t )
    {
        samurai::for_each_cell( mesh[decltype( mesh )::mesh_id_t::reference],
            [&]( auto& cell )
            {
                double x = cell.center( 0 );
                double y = cell.center( 1 );

                double x_0s = ( x < 0. ) ? x_0sm : x_0sp;
                double sign = ( x < 0. ) ? 1.0 : -1.0;

                double r_s2     = ( x_s( x ) - x_0s ) * ( x_s( x ) - x_0s ) + ( y_s( y ) - y_0s ) * ( y_s( y ) - y_0s ); // compute r_s^2
                double v_thetas = Re * Sc * ( 1.0 - std::exp( -r_s2 / ( 4. * Sc * t_s( t ) ) ) ); // remove r_s from this expression

                vel[cell][0] = ( y_s( y ) - y_0s ) * v_thetas / r_s2;        // divide by r_s^2
                vel[cell][1] = sign * ( x_s( x ) - x_0s ) * v_thetas / r_s2; // divide by r_s^2
            } );
    };
    update_velocity( velocity, t_ini );
    auto conv = samurai::make_convection_weno5<field_t>( velocity );
    auto fa   = [&]( double /* t */, auto&& c_ )
    {
        samurai::update_ghost_mr( c_ );
        return -conv( c_ );
    };

    // make ponio problems

    auto diffusion_reaction_pb = ponio::make_imex_operator_problem( fd, fr, fr_t );
    auto advection_pb          = ponio::make_problem( fa );

    auto pb = ponio::make_problem( diffusion_reaction_pb, advection_pb );

    // time loop  -------------------------------------------------------------
    auto eigmax_computer = [&]( auto&, double, auto&, double )
    {
        double dx = mesh.cell_length( mesh.max_level() );
        return 4. / ( dx * dx );
    };

    auto pirock = ponio::runge_kutta::pirock::pirock<1>( ponio::runge_kutta::pirock::beta_0<double>(),
        eigmax_computer,
        ponio::shampine_trick::shampine_trick<field_t>() );
    auto rk3    = ponio::runge_kutta::rk_ssp_33();

    auto strang = ponio::splitting::make_strang_tuple( std::make_pair( pirock, 1e-6 ), std::make_pair( rk3, 1e-8 ) );

    double dt = 1e-7;

    // range to iterate over solution
    auto sol_range = ponio::make_solver_range( pb, strang, c_init, { t_ini, t_end }, dt );

    auto it_sol = sol_range.begin();
    // std::cout << it_sol.tolerance << "\n";

    samurai::make_bc<samurai::Neumann<1>>( it_sol->state, 0. );

    // Prepare MR for solution on iterator
    auto MRadaptation = samurai::make_MRAdapt( it_sol->state );
    MRadaptation( mr_epsilon, mr_regularity );
    samurai::update_ghost_mr( it_sol->state );

    std::size_t n_iteration = 0;
    std::size_t n_save      = 0;
    save( path, filename, it_sol->state, velocity, fmt::format( "_ite_{}", n_save++ ) );

    std::cerr << "> time loop" << std::endl;
    while ( it_sol->time < t_end )
    {
        // samurai::make_bc<samurai::Neumann<1>>( it_sol->state, 0. );

        for ( auto& ki : it_sol.meth.stages( std::integral_constant<std::size_t, 0>{} ) )
        {
            ki.resize();
            ki.fill( 0. );
        }
        update_velocity( velocity, it_sol->time );

        ++it_sol;
        ++n_iteration;

        std::cout << "tⁿ: " << std::setw( 8 ) << it_sol->time << " (Δt: " << it_sol->time_step << ")  "
                  << "N stages:" << it_sol.info().get( std::integral_constant<std::size_t, 0>() ).number_of_stages << " it:" << n_iteration
                  << "   \r";

        MRadaptation( mr_epsilon, mr_regularity );
        samurai::update_ghost_mr( it_sol->state );

        if ( n_iteration % 10 == 0 )
        {
            save( path, filename, it_sol->state, velocity, fmt::format( "_ite_{}", n_save++ ) );
        }
    }
    std::cout << std::endl;

    PetscFinalize();

    return 0;
}
