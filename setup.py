import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="qe-rho-ziang",
    version="0.0.2",
    author="Ziang",
    author_email="ziang.zhang@kaust.com",
    description="A package to generate charge density initial guess for Quantum ESPRESSO",
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
