// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <ranges>
#include <tuple>
#include <type_traits>

#include "detail.hpp"
#include "splitting.hpp"
#include "stage.hpp"

namespace ponio
{

    template <typename state_t>
    auto
    error_estimate( state_t const& un, state_t const& unp1, state_t const& unp1bis )
    {
        return std::abs( ( unp1 - unp1bis ) / ( 1.0 + std::max( std::abs( un ), std::abs( unp1 ) ) ) );
    }

    template <typename state_t>
        requires std::ranges::range<state_t>
    auto
    error_estimate( state_t const& un, state_t const& unp1, state_t const& unp1bis )
    {
        auto it_unp1    = std::ranges::cbegin( unp1 );
        auto it_unp1bis = std::ranges::cbegin( unp1bis );
        auto last       = std::ranges::cend( un );

        auto n_elm = std::distance( std::ranges::cbegin( un ), last );

        using value_t = std::remove_cvref_t<decltype( *it_unp1 )>;
        auto r        = static_cast<value_t>( 0. );

        for ( auto it_un = std::ranges::cbegin( un ); it_un != last; ++it_un, ++it_unp1, ++it_unp1bis )
        {
            auto tmp = ( *it_unp1 - *it_unp1bis ) / ( 1.0 + std::max( std::abs( *it_un ), std::abs( *it_unp1 ) ) );
            r += tmp * tmp;
        }
        return std::sqrt( ( 1. / static_cast<double>( n_elm ) ) * r );

        /*
        return std::sqrt(
          std::accumulate(
            std::cbegin(un), std::cend(un),
            [&it_unp1,&it_unp1bis]( auto r , auto uni ) mutable {
              return std::pow(
                  (*it_unp1 - *it_unp1bis++)/(1.0 + std::max(uni,*it_unp1++))
                , 2u );
            }
          )
        );
        */
    }

    /** @class method
     *  @brief define a time method
     *
     *  A method is define by its algorithm (how compute \f$u^{n+1}\f$ from
     *  \f$(t^n,u^n,\Delta t)\f$) and store possible substeps. Actually this class
     *  avoid the user to specify a template parameter, while making it easy to
     *  add new methods.
     *
     *  @tparam Algorithm_t type of the algorithm which define properly the method
     *  @tparam state_t type of \f$u^n\f$
     */
    template <typename Algorithm_t, typename state_t>
    struct method
    {
        static constexpr bool is_embedded = Algorithm_t::is_embedded;
        using step_storage_t              = typename std::
            conditional<is_embedded, std::array<state_t, Algorithm_t::N_stages + 2>, std::array<state_t, Algorithm_t::N_stages + 1>>::type;

        Algorithm_t alg;
        step_storage_t kis;

        method( Algorithm_t const& alg_, state_t const& shadow_of_u0 );

        template <typename Function_t>
            requires std::invocable<Function_t>
        method( Algorithm_t const& alg_, Function_t&& ki_generator );

        // method( method&& rhs );

        template <typename Problem_t, typename value_t>
        inline std::tuple<value_t, state_t, value_t>
        operator()( Problem_t& f, value_t tn, state_t& un, value_t dt );

        template <std::size_t I = 0, typename Problem_t, typename value_t, typename Algo_t = Algorithm_t>
            requires std::same_as<Algo_t, Algorithm_t> && Algorithm_t::is_embedded
        typename std::enable_if<( I == Algorithm_t::N_stages + 1 ), void>::type
        _call_stage( Problem_t& f, value_t tn, state_t& un, value_t dt );

        template <std::size_t I = 0, typename Problem_t, typename value_t, typename Algo_t = Algorithm_t>
            requires std::same_as<Algo_t, Algorithm_t>
        typename std::enable_if<( I == Algorithm_t::N_stages + 1 ), void>::type
        _call_stage( Problem_t& f, value_t tn, state_t& un, value_t dt );

        template <std::size_t I = 0, typename Problem_t, typename value_t, typename Algo_t = Algorithm_t>
            requires std::same_as<Algo_t, Algorithm_t>
        typename std::enable_if<( I < Algorithm_t::N_stages + 1 ), void>::type
        _call_stage( Problem_t& f, value_t tn, state_t& un, value_t dt );

        template <typename value_t, typename Algo_t = Algorithm_t>
        std::tuple<value_t, state_t, value_t>
        _return( value_t tn, state_t const& un, value_t dt );

        template <typename value_t, typename Algo_t = Algorithm_t>
            requires std::same_as<Algo_t, Algorithm_t> && Algorithm_t::is_embedded
        std::tuple<value_t, state_t, value_t>
        _return( value_t tn, state_t const& un, value_t dt );
    };

    /**
     * constructor of \ref method from its stages and a \f$u_0\f$ (only for preallocation)
     * of temporary substeps
     * @param algo         a `Algorithm_t` objet with predifined stages of the method
     * @param shadow_of_u0 an object with the same size of computed value for allocation
     */
    template <typename Algorithm_t, typename state_t>
    method<Algorithm_t, state_t>::method( Algorithm_t const& alg_, state_t const& shadow_of_u0 )
        : alg( alg_ )
        , kis( ::detail::init_fill_array<std::tuple_size<step_storage_t>::value>( shadow_of_u0 ) )
    {
    }

    /**
     * constructor of \ref method from its stages and a \f$u_0\f$ (only for preallocation)
     * of temporary substeps
     * @param algo         a `Algorithm_t` objet with predifined stages of the method
     * @param ki_generator an invocable object returns an object of the same shape of \f$u^n\f$ or \f$k_i\f$
     */
    template <typename Algorithm_t, typename state_t>
    template <typename Function_t>
        requires std::invocable<Function_t>
    method<Algorithm_t, state_t>::method( Algorithm_t const& alg_, Function_t&& ki_generator )
        : alg( alg_ )
        , kis( ::detail::init_fill_array<std::tuple_size<step_storage_t>::value>( ki_generator ) )
    {
    }

    /**
     * constructor of \ref method from its stages and a \f$u_0\f$ (only for preallocation)
     * of temporary substeps
     * @param rhs a move constructor from this value
     */
    // template <typename Algorithm_t, typename state_t>
    // method<Algorithm_t, state_t>::method( method<Algorithm_t, state_t>&& rhs )
    //     : alg( std::move(rhs.alg) )
    //     , kis( std::move(rhs.kis) )
    // {
    // }

    /**
     * call operator which process all stages of underlying algorithm
     * @param f  callable obect which represents the problem to solve
     * @param tn time \f$t^n\f$ last time where solution is computed
     * @param un computed solution \f$u^n\f$ à time \f$t^n\f$
     * @param dt time step
     * @return tuple of result of iteration \f$(t^{n+1},u^{n+1},\Delta t_{opt})\f$ if
     * iteration is accepted, \f$(t^{n},u^{n},\Delta t_{opt})\f$ otherwise. If `Algorithm_t`
     * is a constant time step method, so \f$\Delta t_{opt}\f$ is alway equal to the same
     * initial value.
     */
    template <typename Algorithm_t, typename state_t>
    template <typename Problem_t, typename value_t>
    inline std::tuple<value_t, state_t, value_t>
    method<Algorithm_t, state_t>::operator()( Problem_t& f, value_t tn, state_t& un, value_t dt )
    {
        _call_stage( f, tn, un, dt );

        return _return( tn, un, dt );
    }

    template <typename Algorithm_t, typename state_t>
    template <std::size_t I, typename Problem_t, typename value_t, typename Algo_t>
        requires std::same_as<Algo_t, Algorithm_t> && Algorithm_t::is_embedded
    typename std::enable_if<( I == Algorithm_t::N_stages + 1 ), void>::type
    method<Algorithm_t, state_t>::_call_stage( Problem_t& f, value_t tn, state_t& un, value_t dt )
    {
        kis[I] = alg.stage( Stage<I>{}, f, tn, un, kis, dt );
    }

    template <typename Algorithm_t, typename state_t>
    template <std::size_t I, typename Problem_t, typename value_t, typename Algo_t>
        requires std::same_as<Algo_t, Algorithm_t>
    typename std::enable_if<( I == Algorithm_t::N_stages + 1 ), void>::type
    method<Algorithm_t, state_t>::_call_stage( Problem_t&, value_t, state_t&, value_t )
    {
    }

    /**
     * unroll all stages of `Algorithm_t` with templated recursion
     * @tparam I stage of `Algorithm_t` to compute
     * @param f  problem whose solution must be computed
     * @param tn time \f$t^n\f$ last time where solution is computed
     * @param un computed solution \f$u^n\f$ à time \f$t^n\f$
     * @param dt time step
     * @return this function store its result in specific attribut of \ref method
     */
    template <typename Algorithm_t, typename state_t>
    template <std::size_t I, typename Problem_t, typename value_t, typename Algo_t>
        requires std::same_as<Algo_t, Algorithm_t>
    typename std::enable_if<( I < Algorithm_t::N_stages + 1 ), void>::type
    method<Algorithm_t, state_t>::_call_stage( Problem_t& f, value_t tn, state_t& un, value_t dt )
    {
        kis[I] = alg.stage( Stage<I>{}, f, tn, un, kis, dt );
        _call_stage<I + 1>( f, tn, un, dt );
    }

    /**
     * return values \f$(t^n,u^n,\Delta t)\f$ after call of all stages
     * @param tn time at the begining of the step
     * @param un state at the begining of the step
     * @param dt time step of the step
     * @return return \f$(t^{n+1},u^{n+1},\Delta t_{opt})\f$ if iteration is accepted,
     * and return \f$(t^{n},u^{n},\Delta t_{opt})\f$ otherwise.
     * @details This member function differs if the algorithm is adaptive time stepping or not.
     */
    template <typename Algorithm_t, typename state_t>
    template <typename value_t, typename Algo_t>
    std::tuple<value_t, state_t, value_t>
    method<Algorithm_t, state_t>::_return( value_t tn, state_t const&, value_t dt )
    {
        return std::make_tuple( tn + dt, kis.back(), dt );
    }

    template <typename Algorithm_t, typename state_t>
    template <typename value_t, typename Algo_t>
        requires std::same_as<Algo_t, Algorithm_t> && Algorithm_t::is_embedded
    std::tuple<value_t, state_t, value_t>
    method<Algorithm_t, state_t>::_return( value_t tn, state_t const& un, value_t dt )
    {
        auto error = error_estimate( un, kis[Algorithm_t::N_stages], kis[Algorithm_t::N_stages + 1] );

        value_t new_dt = 0.9 * std::pow( alg.tol / error, 1. / static_cast<value_t>( Algorithm_t::order ) ) * dt;
        new_dt         = std::min( std::max( 0.2 * dt, new_dt ), 5. * dt );

        if ( error > alg.tol )
        {
            return std::make_tuple( tn, un, new_dt );
        }
        return std::make_tuple( tn + dt, kis[Algorithm_t::N_stages], new_dt );
    }

    /**
     *  generic factory to build a method from an algoritm, it only reuses `method`
     *  constructor
     *  @param algo         a `Algorithm_t` objet with predifined stages of the method
     *  @param shadow_of_u0 an object with the same size of computed value for allocation
     */
    template <typename state_t, typename Algorithm_t>
    auto
    make_method( Algorithm_t const& algo, state_t const& shadow_of_u0 )
    {
        return method<Algorithm_t, state_t>( algo, shadow_of_u0 );
    }

    /**
     *  generic factory to build a method from an algoritm, it only reuses `method`
     *  constructor
     *  @param algo         a `Algorithm_t` objet with predifined stages of the method
     *  @param ki_generator a generator for \f$k_i\f$ in substeps of algorithm
     */
    template <typename state_t, typename Algorithm_t, typename function_t>
        requires std::invocable<function_t>
    auto
    make_method( Algorithm_t const& algo, function_t&& ki_generator )
    {
        return method<Algorithm_t, state_t>( algo, std::forward<function_t>( ki_generator ) );
    }

    /**
     * @brief factory of tuple of methods from a tuple `tuple_algorithms_t`
     *
     * @tparam state_t            type of state variable \f$u^n\f$
     * @tparam tuple_algorithms_t type of tuple of `Algorithm_t`
     * @tparam Arg_t              type of argument
     *
     * @param algos a tuple of `Algorithm_t` objets with predifined stages
     * @param arg   only one `state_t` value if it is from a value (which is copy in each state) or a generator
     */
    template <typename state_t, typename tuple_algorithms_t, typename Arg_t>
    auto
    make_tuple_methods( tuple_algorithms_t const& algos, Arg_t&& arg )
    {
        return std::apply(
            [&]( auto const&... algos_pack )
            {
                auto helper_make_method = [&]( auto const& algo )
                {
                    return make_method( algo, std::forward<Arg_t>( arg ) );
                };
                return std::make_tuple( helper_make_method( algos_pack )... );
            },
            algos );
    }

    /**
     * factory of tuple of methods from a tuple of `Algorithm_t`
     * @param algos        a tuple of `Algorithm_t` objets with predifined stages
     * @param shadow_of_u0 an object with the same size of computed value for allocation
     * @details this factory is to prevent duplucation of code in factory of methods for
     * splitting methods (Lie or Strang method).
     */
    template <typename state_t, typename... Algorithms_t>
    auto
    make_tuple_methods( std::tuple<Algorithms_t...> const& algos, state_t const& shadow_of_u0 )
    {
        return make_tuple_methods_impl<state_t>( algos, shadow_of_u0 );
        // return std::apply(
        //     [&]( auto const&... algos_pack )
        //     {
        //         auto helper_make_method = [&]( auto const& algo )
        //         {
        //             return make_method( algo, shadow_of_u0 );
        //         };
        //         return std::make_tuple( helper_make_method( algos_pack )... );
        //     },
        //     algos );
    }

    template <typename state_t, typename... Algorithms_t, typename function_t>
        requires std::invocable<function_t>
    auto
    make_tuple_methods( std::tuple<Algorithms_t...> const& algos, function_t&& ki_generator )
    {
        return make_tuple_methods_impl<state_t>( algos, std::forward<function_t>( ki_generator ) );
    }

    /**
     *  generic factory to build a method from an algoritm, it only reuses `method`
     *  constructor
     *  @param algos        a variadic `splitting::lie_tuple` of `Algorithms_t`
     *  @param shadow_of_u0 an object with the same size of computed value for allocation
     */
    template <typename state_t, typename... Algorithms_t>
    auto
    make_method( splitting::lie::lie_tuple<Algorithms_t...> const& algos, state_t const& shadow_of_u0 )
    {
        auto methods = make_tuple_methods( algos.algos, shadow_of_u0 );
        return splitting::lie::make_lie_from_tuple( methods, algos.time_steps );
    }

    template <typename state_t, typename... Algorithms_t, typename function_t>
        requires std::invocable<function_t>
    auto
    make_method( splitting::lie_tuple<Algorithms_t...> const& algos, function_t&& ki_generator )
    {
        auto methods = make_tuple_methods<state_t>( algos.algos, std::forward<function_t>( ki_generator ) );
        return splitting::make_lie_from_tuple( methods, algos.time_steps );
    }

    /**
     *  generic factory to build a method from an algoritm, it only reuses `method`
     *  constructor
     *  @param algos        a variadic `splitting::strang_tuple` of `Algorithms_t`
     *  @param shadow_of_u0 an object with the same size of computed value for allocation
     */
    template <typename state_t, typename... Algorithms_t>
    auto
    make_method( splitting::strang::strang_tuple<Algorithms_t...> const& algos, state_t const& shadow_of_u0 )
    {
        auto methods = make_tuple_methods( algos.algos, shadow_of_u0 );
        return splitting::strang::make_strang_from_tuple( methods, algos.time_steps );
    }

    template <typename state_t, typename... Algorithms_t, typename function_t>
        requires std::invocable<function_t>
    auto
    make_method( splitting::strang_tuple<Algorithms_t...> const& algos, function_t&& ki_generator )
    {
        auto methods = make_tuple_methods<state_t>( algos.algos, std::forward<function_t>( ki_generator ) );
        return splitting::make_strang_from_tuple( methods, algos.time_steps );
    }

} // namespace ponio
