FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11
  GIT_TAG v3.0.0
  DOWNLOAD_EXTRACT_TIMESTAMP NEW
)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(_core frenetTransform.cpp)
target_link_libraries(
  _core
  PUBLIC
  transform
)
install(TARGETS _core DESTINATION frenetTransform)
install(DIRECTORY stubs/ DESTINATION frenetTransform)