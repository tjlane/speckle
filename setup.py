#coding: utf8

"""
Setup script for speckle.
"""

from glob import glob

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


setup(name='speckle',
      version='0.0.1',
      author="TJ Lane",
      author_email="tjlane@stanford.edu",
      description='XPCS code',
      packages=["speckle", "speckle"],
      package_dir={"speckle": "speckle"},
      scripts=[s for s in glob('scripts/*') if not s.endswith('__.py')],
      test_suite="test")
