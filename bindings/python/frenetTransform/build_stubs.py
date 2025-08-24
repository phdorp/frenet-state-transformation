import os

import pybind11_stubgen as stubgen
from frenetTransform import _core as transform

# generate stubs after building the package
stubgen.main(["-o", os.path.dirname(transform.__file__) + "/..", "frenetTransform"])

# mypy does not recognize .pyi files without a py.typed marker file
open(os.path.dirname(transform.__file__) + "/py.typed", 'w', encoding='utf-8')
