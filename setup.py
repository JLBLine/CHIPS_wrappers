import setuptools
from setuptools import setup
# from setuptools.command.build_py import build_py as _build_py
# import os
# from wodenpy.wodenpy_setup.git_helper import make_gitdict
# import numpy as np

# class GitInfo(setuptools.Command):
#   '''A custom command to create a json file containing wodenpy git information.'''

#   description = 'Create the file "wodenpy/wodenpy_gitinfo.json" containing git information '
#   user_options = []

#   def initialize_options(self):
#     '''Set default values for options (this has to be included for
#     setuptools.Command to work)'''
#     # Each user option must be listed here with their default value.
#     self.git_info = True

#   def finalize_options(self):
#     '''Post-process options (this has to be included for
#     setuptools.Command to work)'''
#     if self.git_info:
#         print('Creating file wodenpy/wodenpy_gitinfo.npz')

#   def run(self):
#     '''Write the wodenpy git npz file.'''

#     ##Find where we are running the pip install from, and add in a sensible
#     ##place to save the git dictionary
#     save_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),'wodenpy', 'wodenpy_gitinfo.npz')

#     git_dict = make_gitdict()
#     np.savez(save_path, **git_dict)


# class BuildPyCommand(_build_py):
#   '''Custom build command to run the gitinfo command during build'''

#   def run(self):
#     self.run_command('gitinfo')
#     _build_py.run(self)

setup(
    name = "chips_wrappers",
    version = '1.0.0',
    author = "Jack L. B. Line, Aman Chokshi, Jaiden Cook, Dev Null",
    url = "https://github.com/JLBLine/CHIPS_wrappers.git",
    python_requires=">=3.7",
    description = 'Python wrappers to run and plot CHIPS outputs',
    long_description = open("README.md").read(),
    long_description_content_type = 'text/markdown',
    classifiers = [
        "Programming Language :: Python :: 3",
        "License ::  GNU GENERAL PUBLIC LICENSE Version 3.0",
        "Operating System :: Linux",
    ],
    packages = ['chips_wrappers',
                'chips_wrappers.plotting',
    'chips_wrappers.ps_methods',
    'chips_wrappers.run_chips',
    'chips_wrappers.setup'],

    scripts=["scripts/plotchips_all.py",
             "scripts/chips1D_tsv.py",
             "scripts/run_CHIPS.py",],
    install_requires=[
        "numpy",
        "matplotlib",
        "pandas",
        "scipy",
        "astropy",
    ],
)
