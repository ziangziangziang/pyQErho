import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="qe-rho",
    version="1.0.0",
    author="Ziang",
    author_email="ziang.zhang@kaust.edu.sa",
    description="A package to generate charge-density.dat for Quantum ESPRESSO",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ziangziangziang/pyQErho",
    project_urls={
        "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)
