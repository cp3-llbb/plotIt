""" Setuptools-based setup module for the plotIt python library

derived from the pypa example, see https://github.com/pypa/sampleproject
"""

from setuptools import setup
import os, os.path

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the relevant file
from io import open
with open(os.path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="plotIt",

    version="0.1.0",

    description="A utility to plot ROOT histograms",
    long_description=long_description,
    long_description_content_type="text/markdown",

    url="https://github.com/cp3-llbb/plotIt",

    author="Sébastien Brochet, IPNL and CP3-llbb teams",

    license="unknown",

    classifiers=[
        'Development Status :: 5 - Production/Stable',

        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Software Development :: Libraries :: Python Modules',

        'License :: Other/Proprietary License',

        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
    ],

    install_requires=["mplbplot @ git+https://github.com/pieterdavid/mplbplot.git@master"],

    packages=["plotIt"],
    package_dir={"plotIt": "python"},
)
