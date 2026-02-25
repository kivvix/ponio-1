List of additive Runge-Kutta methods
====================================

{% for rk in list_ark %}
.. doxygenvariable:: ponio::runge_kutta::{{ rk.id }}_t
  :project: ponio

{% endfor %}
