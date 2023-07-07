from setuptools import setup, find_packages
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.rst").read_text()

setup(
    name='pycrispr',
    version='1.0.1',
    py_modules=['pycrispr'], 
    description='A package for CRISPR-Cas9 screen analysis',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    project_urls={
        'Documentation': 'https://pycrispr.readthedocs.io/',
        'Source': 'https://github.com/niekwit/pycrispr',
    },
    author='Niek Wit',
    author_email='nw416@cam.ac.uk',
    license='GPL-3.0 license',
    packages=find_packages(),
    install_requires=['pyyaml','Click','sphinx-click'
                      ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)', 
        'Operating System :: POSIX :: Linux', 
        'Programming Language :: Python :: 3',
    ],
    entry_points={
        'console_scripts': [
            'pycrispr = pycrispr.scripts.pycrispr:cli',
        ],
    },
    include_package_data=True,
)
