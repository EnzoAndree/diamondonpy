from setuptools import setup, find_packages

setup(
    name="diamondonpy",
    version="0.1.1",
    description="A Python wrapper for the DIAMOND bioinformatics tool",
    author="Enzo Guerrero-Araya",
    packages=find_packages(),
    python_requires=">=3.6",
    install_requires=[
        'pandas>=1.0.0',
        'numpy>=1.18.0'
    ],
    extras_require={
        'dev': [
            'pytest>=6.0',
            'pytest-cov>=2.0',
            'pytest-mock>=3.0'
        ]
    },
    entry_points={
        "console_scripts": [
            "diamondonpy=diamondonpy.__main__:main"
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
) 