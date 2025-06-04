from setuptools import setup, find_packages


setup(
    name="gpu-coloc",
    version="0.1",
    packages=find_packages(),
    license="MIT",
    
    description="Ultra-fast GPU-enabled Bayesian colocalisation",
    url="https://github.com/mjesse-github/gpu-coloc",
    author="Mihkel Jesse",

    install_requires=[
        "filelock>=3.17.0",
        "fsspec>=2025.2.0",
        "Jinja2>=3.1.5",
        "MarkupSafe>=3.0.2",
        "mpmath>=1.3.0",
        "networkx>=3.4.2",
        "numpy>=2.2.3",
        "pandas>=2.2.3",
        "pyarrow>=19.0.0",
        "python-dateutil>=2.9.0.post0",
        "pytz>=2025.1",
        "six>=1.17.0",
        "sympy>=1.13.1",
        "torch>=2.6.0",
        "tqdm>=4.67.1",
        "typing_extensions>=4.12.2",
        "tzdata>=2025.1"
    ],
    entry_points={
        "console_scripts": [
            "gpu-coloc = gpu_coloc.cli:main",
        ],
    },
    python_requires=">=3.12",
)
