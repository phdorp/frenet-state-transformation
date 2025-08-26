import os

import pybind11_stubgen as stubgen

# import enables process by stubgen
from frenetTransform import _core

# generate stubs after building the package
stubgen.main(
    ["-o", os.environ["PATH_PYBINDS"], "--root-suffix", "Stubs", "frenetTransform"]
)

# mypy does not recognize .pyi files without a py.typed marker file
open(os.environ["PATH_PYSTUBS"] + "/py.typed", "w", encoding="utf-8")
