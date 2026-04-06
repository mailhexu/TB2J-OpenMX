from setuptools import setup

# cffi_modules is not expressible in pyproject.toml, so a minimal setup.py
# is kept solely for this purpose.  All other metadata lives in pyproject.toml.
setup(
    cffi_modules=["TB2J_OpenMX/ffimod.py:ffi"],
    scripts=["scripts/openmx2J.py"],
)
