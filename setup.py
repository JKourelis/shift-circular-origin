"""
Setup configuration for shift-circular-origin package.
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="shift-circular-origin",
    version="1.0.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="Rotate circular DNA sequences to start at specified origin sequences",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/shift-circular-origin",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.7",
    install_requires=[
        "biopython>=1.70",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "flake8",
            "black",
        ],
    },
    entry_points={
        "console_scripts": [
            "shift-circular-origin=shift_circular_origin.cli:main",
        ],
    },
    keywords="bioinformatics dna sequences plasmids origin rotation circular",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/shift-circular-origin/issues",
        "Source": "https://github.com/yourusername/shift-circular-origin",
    },
)