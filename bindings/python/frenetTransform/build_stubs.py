import os
import shutil
import sys

import pybind11_stubgen as stubgen
from frenetTransform import _core as transform

# stub out directory defaults to build if not provided
pathBuild = "build"
if len(sys.argv) == 2:
    pathBuild = sys.argv[1]
pathPybinds = pathBuild + "/bindings/python"

# generate stubs after building the package
stubgen.main(["-o", pathPybinds, "--root-suffix", "Stubs", "frenetTransform"])

# copy stubs to the package directory and preserve files in build
shutil.copytree(
    pathPybinds + "/frenetTransformStubs",
    os.path.dirname(transform.__file__),
    dirs_exist_ok=True,
)

# mypy does not recognize .pyi files without a py.typed marker file
open(os.path.dirname(transform.__file__) + "/py.typed", "w", encoding="utf-8")
