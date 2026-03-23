List of additive Runge-Kutta methods
====================================

{% for rk in list_ark %}{% if rk.explicit.b2 is undefined %}
.. doxygenvariable:: ponio::runge_kutta::{{ rk.id }}_t
  :project: ponio

{% endif %}{% endfor %}


Embedded methods
~~~~~~~~~~~~~~~~

{% for rk in list_ark %}{% if rk.explicit.b2 is defined %}
.. doxygenvariable:: ponio::runge_kutta::{{ rk.id }}_t
  :project: ponio

{% endif %}{% endfor %}
