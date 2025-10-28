// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private, include "../runge_kutta.hpp"

#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <string_view> // NOLINT(misc-include-cleaner)

#include "../butcher_tableau.hpp"
#include "../detail.hpp"
#include "../iteration_info.hpp"
#include "../ponio_config.hpp"
#include "../stage.hpp" // NOLINT(misc-include-cleaner)

namespace ponio::runge_kutta::additive_runge_kutta
{

    template <typename tableau_im_t, typename tableau_ex_t>
    struct additive_runge_kutta
    {
        static_assert( tableau_im_t::N_stages == tableau_ex_t::N_stages );

        tableau_im_t butcher_im;
        tableau_ex_t butcher_ex;
        static constexpr std::size_t N_stages = tableau_im_t::N_stages;
        static constexpr bool is_embedded     = butcher::is_embedded_tableau<tableau_im_t> && butcher::is_embedded_tableau<tableau_ex_t>;
        static constexpr std::size_t order    = tableau_im_t::order;
        static constexpr std::string_view id  = detail::join_id_v<tableau_im_t, detail::separator, tableau_ex_t>;

        using value_t = typename tableau_im_t::value_t;

        additive_runge_kutta( double tolerance = default_config::tol )
            : butcher_im()
            , butcher_ex()
            , _info( tolerance )
        {
            _info.number_of_eval = N_stages;
        }

        template <typename problem_t, typename state_t, typename array_kj_t, std::size_t I>
        void
        stage( Stage<I>, problem_t& pb, value_t tn, state_t& un, array_kj_t const& Kj, value_t dt, state_t& Ki )
        {
            yi = detail::tpl_inner_product<I>( butcher_ex.A[I], K_ex_j, un, dt );
            yi = detail::tpl_inner_product<I>( butcher_im.A[I], K_im_j, yi, dt ) + ... implicit_part;
            // solve yi

            Ki_ex_j[I + 1] = pb.explicit_part( tn + butcher_ex.c[I] * dt, yi );
            Ki_im_j[I + 1] = pb.implicit_part( tn + butcher_im.c[I] * dt, yi );
        }

        template <typename problem_t, typename state_t, typename array_kj_t>
        void
        stage( Stage<N_stages>, problem_t&, value_t, state_t& un, array_kj_t const& Kj, value_t dt, state_t& Ki )
        {
            Ki = detail::tpl_inner_product<N_stages>( butcher_ex.b, K_ex_j, un, dt );
            Ki = detail::tpl_inner_product<N_stages>( butcher_im.b, K_im_j, Ki, dt );
        }

        template <typename problem_t, typename state_t, typename array_kj_t, typename tab_t = tableau_t>
            requires std::same_as<tab_t, tableau_t> && is_embedded
        void
        stage( Stage<N_stages + 1>, problem_t&, value_t, state_t& un, array_kj_t const& Kj, value_t dt, state_t& Ki )
        {
            Ki = detail::tpl_inner_product<N_stages>( butcher_ex.b2, K_ex_j, un, dt );
            Ki = detail::tpl_inner_product<N_stages>( butcher_im.b2, K_im_j, Ki, dt );
        }

        /**
         * @brief gets `iteration_info` object
         */
        auto&
        info()
        {
            return _info;
        }

        /**
         * @brief gets `iteration_info` object (constant version)
         */
        auto const&
        info() const
        {
            return _info;
        }

        iteration_info<tableau_t> _info;
    };

} // namespace ponio::runge_kutta::additive_runge_kutta
