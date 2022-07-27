epydemic-signals: Epidemics as signals
======================================

.. image:: https://badge.fury.io/py/epydemic-signals.svg
    :target: https://badge.fury.io/py/epydemic-signals

.. image:: https://github.com/simoninireland/epydemic-signals/actions/workflows/ci.yaml/badge.svg
    :target: https://github.com/simoninireland/epydemic-signals/actions/workflows/ci.yaml

.. image:: https://zenodo.org/badge/428717218.svg
   :target: https://zenodo.org/badge/latestdoi/428717218

.. image:: https://www.gnu.org/graphics/gplv3-88x31.png
    :target: https://www.gnu.org/licenses/gpl-3.0.en.html

Overview
--------

``epydemic-signals`` is an experiment in treating epidemic processes
on networks as "signals" that can be analysed using the tools of
signal processing, discrete algebra, and algebraic topology. It is
closely integrated with the ``epydemic`` epidemic process simulator.

A graph signal associates a value with each node in a network at each
point in time during the execution of a process. The signal can change
at each event in the process' evolution.

``epydemic`` already provides some monitoring of processes, and
further monitoring can be added either by sub-classing the different
processes or by adding additional process
instances. ``epydemic-signals`` provides a third approach: tapping the
raw event streams of simulations into one or more signal generators
that can be used to construct time-varying signals derived from the
process. The result is one or more time series, which are saved to
``epyc`` notebooks alongside other experimental results.


Installation
------------

You can install ``epydemic-signals`` directly from PyPi using ``pip``:

::

   pip install epydemic-signals

The master distribution of ``epydemic-signals`` is hosted on GitHub. To obtain a
copy, just clone the repo:

::

    git clone git@github.com:simoninireland/epydemic-signals.git
    cd epydemic
    pip install .



Author and license
------------------

Copyright (c) 2021--2022, Simon Dobson <simoninireland@gmail.com>

Licensed under the `GNU General Public Licence v3 <https://www.gnu.org/licenses/gpl-3.0.en.html>`_.
