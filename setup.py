“””
Setup configuration for geometric-optimization package.
“””

from setuptools import setup, find_packages

with open(“README.md”, “r”, encoding=“utf-8”) as fh:
long_description = fh.read()

setup(
name=“geometric-optimization”,
version=“0.1.0”,
author=“JinnZ2 Collaborative Research”,
author_email=“anonymous@research.local”,
description=“Geometric optimization via E₈ lattice structure”,
long_description=long_description,
long_description_content_type=“text/markdown”,
url=“https://github.com/JinnZ2/geometric-optimization”,
packages=find_packages(),
classifiers=[
“Development Status :: 3 - Alpha”,
“Intended Audience :: Science/Research”,
“Topic :: Scientific/Engineering :: Mathematics”,
“Topic :: Scientific/Engineering :: Physics”,
“License :: OSI Approved :: Apache Software License”,
“Programming Language :: Python :: 3”,
“Programming Language :: Python :: 3.8”,
“Programming Language :: Python :: 3.9”,
“Programming Language :: Python :: 3.10”,
“Programming Language :: Python :: 3.11”,
],
python_requires=”>=3.8”,
install_requires=[
“numpy>=1.20.0”,
“scipy>=1.7.0”,
],
extras_require={
“dev”: [
“pytest>=6.0”,
“pytest-cov>=2.12”,
“black>=21.0”,
“flake8>=3.9”,
“mypy>=0.900”,
],
“gpu”: [
“jax[cuda]>=0.3.0”,
“jaxlib>=0.3.0”,
],
“viz”: [
“matplotlib>=3.3.0”,
“seaborn>=0.11.0”,
“plotly>=5.0.0”,
],
“notebooks”: [
“jupyter>=1.0.0”,
“ipywidgets>=7.6.0”,
],
},
entry_points={
“console_scripts”: [
“gopt-optimize=gas.cli:main”,
],
},
)
