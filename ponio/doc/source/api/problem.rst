Problems
========

Simple problem
--------------

.. doxygenclass:: ponio::simple_problem
   :project: ponio
   :members:

.. doxygenfunction:: ponio::make_simple_problem
   :project: ponio


Lawson problem
--------------

.. doxygenclass:: ponio::lawson_problem
   :project: ponio
   :members:

.. doxygenfunction:: ponio::make_lawson_problem
   :project: ponio


Implicit problem
----------------

This kind of problems are defined by a function and its Jacobian.

.. doxygenclass:: ponio::implicit_problem
   :project: ponio
   :members:

.. doxygenfunction:: ponio::make_implicit_problem
   :project: ponio


Implicit problem with operator
------------------------------

This kind of problems are defined by a function and an operator, for example to solve :

.. math::

  \dot{u} = f(t, u)

- the function is :math:`f` ;
- the operator is :math:`t \mapsto f(t, \cdot)` which is a function.

.. doxygenclass:: ponio::implicit_operator_problem
   :project: ponio
   :members:

.. doxygenfunction:: ponio::make_implicit_operator_problem
   :project: ponio

Coupled Implicit-Explicit problem
---------------------------------

This kind of problems are defined by a function and an operator, for example to solve :

.. math::

  \dot{u} = f(t, u) + g(t,u)

where we would like to solve explicitly :math:`f` and implicitly :math:`g`.

.. doxygenclass:: ponio::imex_problem
   :project: ponio
   :members:

Like for implicit problems we can define the implicit part with its Jacobian matrix or with a samurai operator.

.. doxygenfunction:: ponio::make_imex_jacobian_problem
   :project: ponio

.. doxygenfunction:: ponio::make_imex_operator_problem
   :project: ponio

Problem
-------

A problem is a collection of sub-problems defined by previous problems.

.. doxygenclass:: ponio::problem
   :project: ponio
   :members:
