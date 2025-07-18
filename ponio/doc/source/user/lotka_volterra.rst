Understand observers with Lotka-Volterra equations
==================================================

In this example, we present what can we do with :doc:`observers <../api/observer>` in ponio, with Lotka-Voletrra equations. The equations are the following

.. math::

  \begin{cases}
    \frac{\mathrm{d}x}{\mathrm{d}t} = \alpha x - \beta xy\\
    \frac{\mathrm{d}y}{\mathrm{d}t} = \delta xy - \gamma y\\
  \end{cases}

with :math:`\alpha = \frac{2}{3}`, :math:`\beta = \frac{4}{3}` and :math:`\delta = \gamma = 1`. Initial condition :math:`(x, y)(t=0) = (x_0, x_0)`, and :math:`x_0 = 1`.

We write the problem like

.. literalinclude:: ../_static/cpp/lotka_volterra_fobs.cpp
  :language: cpp
  :lines: 19-29
  :lineno-start: 19
  :linenos:

And launch the simulation between :math:`0` and :math:`15` with the time step :math:`\Delta t = 0.1`

.. literalinclude:: ../_static/cpp/lotka_volterra_fobs.cpp
  :language: cpp
  :lines: 34-39
  :lineno-start: 34
  :linenos:


The next sections specify how to define the observer :code:`obs`, to write into a file, in standard output, a stream, or a user-defined observer. The observer is called at each iteration in time and take the current time :math:`t^n`, the current state :math:`u^n` and the current time step :math:`\Delta t^n`.

The file observer
-----------------

In this example we export data into a file with :cpp:class:`ponio::observer::file_observer` class. The ponio library provides this output for simple types and containers (thanks concepts). The output is formate as following

.. literalinclude:: ../_static/cpp/lotka_volterra_fobs.txt
  :language: text
  :lines: 1-5

First column is the current time, the last one is the current time step, and between we get the multiple values of current state. Other observers provide by ponio have the same format.

The :cpp:class:`ponio::observer::file_observer` class can be build with a `std::string <https://en.cppreference.com/w/cpp/string/basic_string>`_ or a `std::filesystem::path <https://en.cppreference.com/w/cpp/filesystem/path>`_. The constructor makes necessary directory if needed.

.. literalinclude:: ../_static/cpp/lotka_volterra_fobs.cpp
  :language: cpp
  :lines: 31-32
  :lineno-start: 31
  :linenos:

.. note::

  If name is fixed at compile time you can use literal :code:`_fobs` by

  .. code-block:: cpp

    using namespace ponio::observer;

    auto obs = "output.txt"_fobs;

.. note::

  Data are not flushed into the output file, you have to wait the flush of the buffer.

The full example can be found in :download:`lotka_volterra_fobs.cpp <../_static/cpp/lotka_volterra_fobs.cpp>`.


.. figure:: ../_static/cpp/lotka_volterra_fobs.png
    :width: 500 px
    :alt: Solution in phase space (x, y)

    Solution in phase space :math:`(x, y)`


The :code:`cout` observer
-------------------------

In this example we export data into the standard output with :cpp:class:`ponio::observer::cout_observer` class.

.. literalinclude:: ../_static/cpp/lotka_volterra_cobs.cpp
  :language: cpp
  :lines: 29
  :lineno-start: 29
  :linenos:

.. note::

  Data are not flushed into the standard output, you have to wait the flush of the buffer.

The full example can be found in :download:`lotka_volterra_cobs.cpp <../_static/cpp/lotka_volterra_cobs.cpp>`.


.. figure:: ../_static/cpp/lotka_volterra_cobs.png
    :width: 500 px
    :alt: Solution in phase space (x, y)

    Solution in phase space :math:`(x, y)`


The stream observer
-------------------

In this example we export data into the standard output with :cpp:class:`ponio::observer::stream_observer` class. This observer is the more generally one and can be plug with any other type of stream, for example a `std::stringstream <https://en.cppreference.com/w/cpp/io/basic_stringstream>`_. The :cpp:class:`ponio::observer::file_observer` is a specialization with a `std::ofstream <https://en.cppreference.com/w/cpp/io/basic_ofstream>`_ and :cpp:class:`ponio::observer::cout_observer` with `std::cout <https://en.cppreference.com/w/cpp/io/cout>`_.


.. literalinclude:: ../_static/cpp/lotka_volterra_sobs.cpp
  :language: cpp
  :lines: 31-32
  :lineno-start: 31
  :linenos:

Next we get all informations into the observer or our buffer.

.. literalinclude:: ../_static/cpp/lotka_volterra_sobs.cpp
  :language: cpp
  :lines: 41
  :lineno-start: 41
  :linenos:


The full example can be found in :download:`lotka_volterra_sobs.cpp <../_static/cpp/lotka_volterra_sobs.cpp>`.


.. figure:: ../_static/cpp/lotka_volterra_sobs.png
    :width: 500 px
    :alt: Solution in phase space (x, y)

    Solution in phase space :math:`(x, y)`


The user-defined observer
-------------------------

An observer is a invocable object (function, functor, lambda) with the following signature

.. code-block:: cpp

  void operator() ( value_t current_time, state_t & current_state, value_t current_time_step );

In our case, :code:`value_t` is :code:`double` and :code:`state_t` is :code:`std::valarray<double>`. We propose to export also the invariant of Lotka-Voletrra equations given by

.. math::

  V = \delta x - \ln(x) + \beta y - \alpha \ln(y)

We implement the observer into a class

.. literalinclude:: ../_static/cpp/lotka_volterra_uobs.cpp
  :language: cpp
  :lines: 13-50
  :lineno-start: 13
  :linenos:

Now we just have to build an instance of this class

.. literalinclude:: ../_static/cpp/lotka_volterra_uobs.cpp
  :language: cpp
  :lines: 69
  :lineno-start: 69
  :linenos:

The output of this observer looks like this


.. literalinclude:: ../_static/cpp/lotka_volterra_uobs.txt
  :language: text
  :lines: 1-5

The full example can be found in :download:`lotka_volterra_uobs.cpp <../_static/cpp/lotka_volterra_uobs.cpp>`.

.. figure:: ../_static/cpp/lotka_volterra_uobs.png
    :width: 500 px
    :alt: Solution in phase space (x, y)

    Solution in phase space :math:`(x, y)`


We can compute from previous simulation the invariant :math:`V`, or display it with out user-defined observer. We don't change the default precision in C++ which is set to 6 by default (see `std::ios_base::precision <https://en.cppreference.com/w/cpp/io/ios_base/precision>`_), in ponio observers the precision is set to the maximum of precision of given type.

.. figure:: ../_static/cpp/lotka_volterra.png
    :width: 500 px
    :alt: Relative error of invariant V

    Relative error of invariant :math:`V`
