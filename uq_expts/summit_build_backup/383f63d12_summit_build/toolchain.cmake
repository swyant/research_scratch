include(summit_init)
set(TOOLS_PREFIX "$ENV{HOME}/tools")
set(HDF5_PREFER_PARALLEL ON)
set(PYBIND11_FINDPYTHON ON)

# Python configuration
set(Python3_ROOT_DIR "$ENV{CONDA_PREFIX}")
set(Python3_EXECUTABLE "$ENV{CONDA_PREFIX}/bin/python")

# Dynamic Python path detection
execute_process(
    COMMAND python -c "import sysconfig; print(sysconfig.get_path('include'))"
    OUTPUT_VARIABLE PYTHON_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(PYTHON_INCLUDE_DIR)
    set(Python3_INCLUDE_DIR "${PYTHON_INCLUDE_DIR}")
    message(STATUS "Found Python include directory: ${Python3_INCLUDE_DIR}")
else()
    set(Python3_INCLUDE_DIR "$ENV{CONDA_PREFIX}/include/python3.11")
endif()

execute_process(
    COMMAND python -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))"
    OUTPUT_VARIABLE PYTHON_LIB_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(PYTHON_LIB_DIR)
    set(Python3_LIBRARY "${PYTHON_LIB_DIR}/libpython3.12.so")
    message(STATUS "Found Python library directory: ${PYTHON_LIB_DIR}")
else()
    set(Python3_LIBRARY "$ENV{CONDA_PREFIX}/lib/libpython3.12.so")
endif()

set(Python3_NumPy_INCLUDE_DIRS "$ENV{CONDA_PREFIX}/lib/python3.11/site-packages/numpy/core/include")
set(Python3_FIND_VERSION "3.12")
set(Python3_FIND_VERSION_MAJOR "3")
set(Python3_FIND_VERSION_MINOR "12")

# Dependency configuration
summit_set_dependency(Eigen3     "$ENV{CONDA_PREFIX}"                         "Eigen3::Eigen")
summit_set_dependency(GTEST      "$ENV{CONDA_PREFIX}"                         "GTest::gtest")
summit_set_dependency(GSL        "$ENV{CONDA_PREFIX}"                         "GSL::gsl")
summit_set_dependency(HDF5       "$ENV{CONDA_PREFIX}"                         "HDF5::HDF5")
summit_set_dependency(METIS      "$ENV{CONDA_PREFIX}"                         "METIS::METIS")
summit_set_dependency(MPI        "$ENV{CONDA_PREFIX}"                         "MPI::MPI_CXX")
summit_set_dependency(ParMETIS   "$ENV{CONDA_PREFIX}"                         "ParMETIS::ParMETIS")
summit_set_dependency(PETSc      "$ENV{HOME}/local_software/petsc/arch-linux-cxx-debug"                         "PETSc::PETSc")
#summit_set_dependency(PETSc      "$ENV{CONDA_PREFIX}"                         "PETSc::PETSc")
summit_set_dependency(pybind11   "$ENV{CONDA_PREFIX}"                         "pybind11::module")
summit_set_dependency(PYRE       "$ENV{CONDA_PREFIX}"                       "pyre::pyre")
summit_set_dependency(Python     "$ENV{CONDA_PREFIX}"                         "Python3::Python")
#summit_set_dependency(SLEPc      "$ENV{CONDA_PREFIX}"                         "SLEPc::SLEPc")
summit_set_dependency(VTK        "$ENV{CONDA_PREFIX}"                         "VTK::CommonCore")
summit_set_dependency(VTK        "$ENV{CONDA_PREFIX}"                         "VTK::IOXML")
summit_set_dependency(VTK        "$ENV{CONDA_PREFIX}"                         "VTK::CommonDataModel")
summit_set_dependency(yaml-cpp   "$ENV{CONDA_PREFIX}"                         "yaml-cpp")
#summit_set_dependency(CGAL       "$ENV{CONDA_PREFIX}"                         "CGAL::CGAL")

# Optional dependencies (comment out if not needed)
# summit_set_dependency(CANTERA    "$ENV{CONDA_PREFIX}"                              "Cantera::Cantera")
# summit_set_dependency(CUDA       ""                                           "CUDA::CUDA")
# summit_set_dependency(GMSH       "$ENV{CONDA_PREFIX}"                              "GMSH::GMSH")
# summit_set_dependency(KOKKOS     "${TOOLS_PREFIX}/kokkos/cuda"                "Kokkos::kokkos")
# summit_set_dependency(MKL        "/opt/intel/oneapi/mkl/2025.1"               "MKL::MKL")
# summit_set_dependency(SUNDIALS   "$ENV{CONDA_PREFIX}"                              "SUNDIALS::sundials_cvode")

list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" CACHE STRING "" FORCE)
summit_generate_features_file()
