#!/usr/bin/env python3
"""Setup script for retina-spatial-scrna-analysis package."""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="retina-spatial-scrna",
    version="0.1.0",
    author="ChangHee Lee Lab",
    description="Spatial gene expression analysis toolkit for retinal single-cell RNA-seq data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/chlee-tabin/retina-spatial-scrna-analysis",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.10",
    install_requires=requirements,
    include_package_data=True,
    package_data={
        "retina_spatial_scrna": ["control_genes.yaml"],
    },
)
