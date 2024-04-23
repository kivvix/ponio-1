// Copyright 2022 PONIO TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#pragma once
// autogenerated file

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string_view>
#include <tuple>

#include "../butcher_tableau.hpp"
#include "../ponio_config.hpp"
#include "../runge_kutta/dirk.hpp"
#include "../runge_kutta/erk.hpp"
#include "../runge_kutta/exprk.hpp"
#include "../runge_kutta/lrk.hpp"
#include "../runge_kutta/rkc.hpp"

namespace ponio::runge_kutta
{
    // clang-format off

// explicit Runge-Kutta methods
{% for rk in list_erk %}
/**
 * @brief Butcher tableau of {{ rk.label }} method
 * @tparam value_t type of coefficient (``double``by default)
 */
template <typename value_t=double>
struct butcher_{{ rk.id }} : public butcher::{{ "adaptive_" if 'b2' in rk else "" }}butcher_tableau<{{ rk.A|length }}, value_t>
{
  using base_t = butcher::{{ "adaptive_" if 'b2' in rk else "" }}butcher_tableau<{{ rk.A|length }}, value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  static constexpr std::size_t order    = {{ rk.order }};
  static constexpr std::string_view id  = "{{ rk.id }}";

  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_{{ rk.id }}()
  : base_t(
    {{ '{{' }}
    {% for ai in rk.A -%}
      { {{ ai }} }{{ ",\n      " if not loop.last else "" }}
    {%- endfor %}
    {{ '}}' }}, // A
    { {{ rk.b }} }, // b
    {% if 'b2' in rk -%}{ {{ rk.b2 }} }, // b2 {%- endif %}
    { {{ rk.c }} }  // c
  )
  {}
};

/**
 * @brief {{ rk.label }} method
 * @tparam value_t type of coefficient (``double``by default)
 *
 * @details see more on [ponio](https://josselin.massot.gitlab.labos.polytechnique.fr/ponio/viewer.html#{{ rk.id }})
 *
 * This method is based on the following Butcher tableau:
 *
 * \f[
 *  \begin{array}{c|{%- for ci in rk.butcher.c -%}c{%- endfor -%}}
      {%- for ai in rk.butcher.A %}
 *      {{ rk.butcher.c[loop.index0] }} & {{ ai|join(' & ') }} \\
 {%- endfor %}
 *    \hline
 *      & {{ rk.butcher.b|join(' & ') }} {% if 'b2' in rk %} \\
 *    \hline
 *      & {{ rk.butcher.b2|join(' & ') }}
{%- endif %}
 *  \end{array}
 * \f]
 *
 * + **stages:** {{ rk.A|length }}
 * + **order:** {{ rk.order }}
 * + **stages order:** {{ rk.stage_order }}
 * + **stability function:** \f[ {{ rk.stability_function }} \f]
 */
template <typename value_t>
using {{ rk.id }}_t = explicit_runge_kutta::explicit_runge_kutta<butcher_{{ rk.id }}<value_t>>;

using {{ rk.id }} = explicit_runge_kutta::explicit_runge_kutta<butcher_{{ rk.id }}<double>>;

/**
 * @brief Lawson l{{ rk.label }} method
 *
 * @tparam value_t type of coefficient (``double``by default)
 * @tparam Exp_t   type of exponential function passing in argument
 *
 * @param exp_ exponential function for the Lawson method
 * @param tol  tolerance (only for adaptive time step methods)
 *
 * @details see more on [ponio](https://josselin.massot.gitlab.labos.polytechnique.fr/ponio/viewer.html#{{ rk.id }})
 */
template <typename value_t, typename Exp_t>
constexpr auto l{{ rk.id }}_t = []( Exp_t exp_ , double tol=ponio::default_config::tol )
{
  return lawson_runge_kutta::make_lawson<butcher_{{ rk.id }}<value_t>,Exp_t>(exp_,tol);
};

template <typename Exp_t>
auto
l{{ rk.id }}( Exp_t exp_ , double tol=ponio::default_config::tol )
{
  return lawson_runge_kutta::make_lawson<butcher_{{ rk.id }}<double>,Exp_t>(exp_,tol);
}

{% endfor %}

/**
 * @brief Type of tuple that contains all classical explicit Runge-Kutta methods of ponio
 */
template <typename value_t>
using erk_tuple = std::tuple< {{ list_erk | sformat("{}_t<value_t>", attribute="id") | join(", ") }} >;

/**
 * @brief Type of tuple that contains all Lawson methods of ponio
*/
template <typename value_t, typename Exp_t>
using lrk_tuple = std::tuple< {{ list_erk | sformat("decltype(l{}_t<value_t, Exp_t>)", attribute="id") | join(", ") }} >;

// exponential Runge-Kutta methods
{% for rk in list_exprk %}
/**
 * @brief Butcher tableau of {{ rk.label }} method
 * @tparam _value_t  type of coefficient (``double`` by default)
 * @tparam _linear_t type of linear part (``double`` by default)
 */
template <typename _value_t=double, typename _linear_t=double>
struct butcher_{{ rk.id }}
{
  using value_t  = _value_t;
  using linear_t = _linear_t;
  using func_t   = std::function<linear_t(linear_t &&)>;
  static constexpr std::size_t N_stages = {{ rk.b.type|length }};
  static constexpr std::size_t order    = 1;
  static constexpr std::string_view id  = "{{ rk.id }}";

  std::tuple<{% for t in rk.A.type -%}{{ t }}{{ " , " if not loop.last else "" }}{%- endfor %}> a;
  std::tuple<{% for t in rk.b.type -%}{{ t }}{{ " , " if not loop.last else "" }}{%- endfor %}> b;
  std::array<value_t, N_stages> c;

  butcher_{{ rk.id }}()
  : a( {% for aij in rk.A.code %}{{ aij }}{{ " , " if not loop.last else "" }}{%- endfor %} )
  , b( {% for bi  in rk.b.code %}{{ bi  }}{{ " , " if not loop.last else "" }}{%- endfor %} )
  , c({ {{ rk.c }} })
  {}
};

/**
 * @brief {{ rk.label }} method
 * @tparam value_t  type of coefficient (``double``by default)
 * @tparam linear_t type of coefficient (``double``by default)
 */
template <typename value_t, typename linear_t>
using {{ rk.id }}_t = exponential_runge_kutta::explicit_exp_rk_butcher<butcher_{{ rk.id }}<value_t, linear_t>>;

using {{ rk.id }} = exponential_runge_kutta::explicit_exp_rk_butcher<butcher_{{ rk.id }}<double, double>>;

{% endfor %}

/**
 * @brief Type of tuple that contains all exponential Runge-Kutta methods of ponio
*/
template <typename value_t, typename linear_t>
using exprk_tuple = std::tuple< {{ list_exprk | sformat("{}_t<value_t, linear_t>", attribute="id") | join(", ") }} >;

// diagonal implicit Runge-Kutta methods
{% for rk in list_dirk %}
/**
 * @brief Butcher tableau of {{ rk.label }} method
 * @tparam value_t type of coefficient (``double``by default)
 */
template <typename value_t=double>
struct butcher_{{ rk.id }} : public butcher::{{ "adaptive_" if 'b2' in rk else "" }}butcher_tableau<{{ rk.A|length }},value_t>
{
  using base_t = butcher::{{ "adaptive_" if 'b2' in rk else "" }}butcher_tableau<{{ rk.A|length }},value_t>;
  static constexpr std::size_t N_stages = base_t::N_stages;
  static constexpr std::size_t order    = {{ rk.order }};
  static constexpr std::string_view id  = "{{ rk.id }}";

  using base_t::A;
  using base_t::b;
  using base_t::c;

  butcher_{{ rk.id }}()
  : base_t(
    {{ '{{' }}
    {% for ai in rk.A -%}
      { {{ ai }} }{{ ",\n      " if not loop.last else "" }}
    {%- endfor %}
    {{ '}}' }}, // A
    { {{ rk.b }} }, // b
    {% if 'b2' in rk -%}{ {{ rk.b2 }} }, // b2 {%- endif %}
    { {{ rk.c }} }  // c
  )
  {}
};

/**
 * @brief {{ rk.label }} method
 *
 * @tparam value_t          type of coefficient (``double``by default)
 * @tparam linear_algebra_t type that provides linear algebra if it is undefined for state_t (``void`` by default)
 * @tparam Args optional    arguments to build linear_algebra_t object
 */
template <typename value_t, typename linear_algebra_t=void, typename ... Args>
auto
{{ rk.id }}_t ( Args ... args )
{
  return diagonal_implicit_runge_kutta::make_dirk<butcher_{{ rk.id }}<value_t>, linear_algebra_t>(args...);
}

template <typename linear_algebra_t=void, typename ... Args>
auto
{{ rk.id }} ( Args ... args )
{
  return diagonal_implicit_runge_kutta::make_dirk<butcher_{{ rk.id }}<double>, linear_algebra_t>(args...);
}

{% endfor %}

    // clang-format on

} // namespace ponio::runge_kutta
