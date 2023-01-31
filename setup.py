from setuptools import setup, find_packages

setup(
    name='pycrispr',
    version='0.1.0',
    py_modules=['pycrispr'], 
    description='A package for CRISPR-Cas9 screen analysis',
    url='https://github.com/niekwit/pycrispr',
    author='Niek Wit',
    author_email='nw416@cam.ac.uk',
    license='GPL-3.0 license',
    packages=find_packages(),
    install_requires=['seaborn','matplotlib',
                      'numpy','pandas','pyyaml',
                      'tqdm','gseapy','cutadapt', 'multiqc','Click'
                      ],
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GPL-3.0 License', 
        'Operating System :: POSIX :: Linux', 
        'Programming Language :: Python :: 3',
    ],
    entry_points={
        'console_scripts': [
            'pycrispr = pycrispr.scripts.pycrispr:cli',
        ],
    },
)
