conda activate deepchem
python setup.py sdist bdist_wheel
twine upload --repository testpypi dist/* --verbose
pip install --index-url https://test.pypi.org/simple/ qe_rho_ziang
