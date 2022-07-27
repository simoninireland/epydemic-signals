# Setup for epydemic-signals
#
# Copyright (C) 2021--2022 Simon Dobson
#
# This file is part of epydemic-signals, an experiment in epidemic processes.
#
# epydemic-signals is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# epydemic-signals is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with epydemic-signals. If not, see <http://www.gnu.org/licenses/gpl.html>.

from setuptools import setup

with open('README.rst') as f:
    longDescription = f.read()

setup(name='epydemic-signals',
      version='0.1.1',
      description='An experiment in spidemics as signals',
      long_description=longDescription,
      url='http://github.com/simoninireland/epydemic-signals',
      author='Simon Dobson',
      author_email='simoninireland@gmail.com',
      license='License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      classifiers=['Development Status :: 2 - Pre-Alpha',
                   'Intended Audience :: Science/Research',
                   'Intended Audience :: Developers',
                   'Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7',
                   'Programming Language :: Python :: 3.8',
                   'Programming Language :: Python :: 3.9',
                   'Programming Language :: Python :: 3.10',
                   'Topic :: Scientific/Engineering'],
      python_requires='>=3.6',
      packages=['epydemic_signals',
                'epydemic_signals.plot'],
      package_data={'epydemic_signals': ['py.typed']},
      zip_safe=False,
      install_requires=["epydemic >= 1.11.1", "networkx", "pandas", "pygsp", "pyunlocbox", "matplotlib", "mypy", "jedi", "jedi-language-server", "black", ],
      extra_requires={':python_version < 3.8': ['typing_extensions']},
)
