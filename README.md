# Tournser
Computing homology of tournaplexes (based on https://github.com/luetge/flagser)

An online version is available [here](https://homepages.abdn.ac.uk/neurotopology/tournser.html).

To install, first download repository with:
```sh
git clone --recursive https://github.com/JasonPSmith/tournser.git
```
Next, compile flagser count with
```
make
```

To install pyflagsercount run:
```sh
pip install .
```
Requirements: For pyflagsercount the python packages numpy and pybind11 are required, and cmake version ≥ 2.8.12.
