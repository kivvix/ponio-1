// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// IWYU pragma: private, include "../runge_kutta.hpp"

#pragma once

#include <concepts>
#include <cstddef>
#include <string_view> // NOLINT(misc-include-cleaner)

#include "../butcher_tableau.hpp"
#include "../detail.hpp"
#include "../iteration_info.hpp"
#include "../ponio_config.hpp"
#include "../stage.hpp" // NOLINT(misc-include-cleaner)

namespace ponio::runge_kutta::explicit_runge_kutta
{

    template <typename tableau_t>
    struct explicit_runge_kutta
    {
        tableau_t butcher;
        static constexpr std::size_t N_stages = tableau_t::N_stages;
        static constexpr bool is_embedded     = butcher::is_embedded_tableau<tableau_t>;
        static constexpr std::size_t order    = tableau_t::order;
        static constexpr std::string_view id  = tableau_t::id;

        using value_t = typename tableau_t::value_t;

        explicit_runge_kutta( double tolerance = default_config::tol )
            : butcher()
            , _info( tolerance )
        {
            _info.number_of_eval = N_stages;
        }

        template <typename problem_t, typename state_t, typename array_ki_t, std::size_t I>
        state_t
        stage( Stage<I>, problem_t& f, value_t tn, state_t& un, array_ki_t const& Ki, value_t dt )
        {
            return f( tn + butcher.c[I] * dt, detail::tpl_inner_product<I>( butcher.A[I], Ki, un, dt ) );
        }

        template <typename problem_t, typename state_t, typename array_ki_t>
        state_t
        stage( Stage<N_stages>, problem_t&, value_t, state_t& un, array_ki_t const& Ki, value_t dt )
        {
            return detail::tpl_inner_product<N_stages>( butcher.b, Ki, un, dt );
        }

        template <typename problem_t, typename state_t, typename array_ki_t, typename tab_t = tableau_t>
            requires std::same_as<tab_t, tableau_t> && is_embedded
        state_t
        stage( Stage<N_stages + 1>, problem_t&, value_t, state_t& un, array_ki_t const& Ki, value_t dt )
        {
            return detail::tpl_inner_product<N_stages>( butcher.b2, Ki, un, dt );
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

} // namespace ponio::runge_kutta::explicit_runge_kutta
