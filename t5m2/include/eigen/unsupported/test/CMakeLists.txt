
include_directories(../../test ../../unsupported ../../Eigen)

set(SPARSE_LIBS "")

find_package(Cholmod)
if(CHOLMOD_FOUND)
  add_definitions("-DEIGEN_CHOLMOD_SUPPORT")
  include_directories(${CHOLMOD_INCLUDES})
  set(SPARSE_LIBS ${SPARSE_LIBS} ${CHOLMOD_LIBRARIES})
  ei_add_property(EIGEN_TESTED_BACKENDS "Cholmod, ")
else(CHOLMOD_FOUND)
  ei_add_property(EIGEN_MISSING_BACKENDS "Cholmod, ")
endif(CHOLMOD_FOUND)

find_package(Umfpack)
if(UMFPACK_FOUND)
  add_definitions("-DEIGEN_UMFPACK_SUPPORT")
  include_directories(${UMFPACK_INCLUDES})
  set(SPARSE_LIBS ${SPARSE_LIBS} ${UMFPACK_LIBRARIES})
  ei_add_property(EIGEN_TESTED_BACKENDS "UmfPack, ")
else(UMFPACK_FOUND)
  ei_add_property(EIGEN_MISSING_BACKENDS "UmfPack, ")
endif(UMFPACK_FOUND)

find_package(SuperLU)
if(SUPERLU_FOUND)
  add_definitions("-DEIGEN_SUPERLU_SUPPORT")
  include_directories(${SUPERLU_INCLUDES})
  set(SPARSE_LIBS ${SPARSE_LIBS} ${SUPERLU_LIBRARIES})
  ei_add_property(EIGEN_TESTED_BACKENDS  "SuperLU, ")
else(SUPERLU_FOUND)
  ei_add_property(EIGEN_MISSING_BACKENDS  "SuperLU, ")
endif(SUPERLU_FOUND)

find_package(GoogleHash)
if(GOOGLEHASH_FOUND)
  add_definitions("-DEIGEN_GOOGLEHASH_SUPPORT")
  include_directories(${GOOGLEHASH_INCLUDES})
  ei_add_property(EIGEN_TESTED_BACKENDS  "GoogleHash, ")
else(GOOGLEHASH_FOUND)
  ei_add_property(EIGEN_MISSING_BACKENDS  "GoogleHash, ")
endif(GOOGLEHASH_FOUND)

find_package(Adolc)
if(ADOLC_FOUND)
  include_directories(${ADOLC_INCLUDES})
  ei_add_property(EIGEN_TESTED_BACKENDS "Adolc, ")
  ei_add_test(forward_adolc " " ${ADOLC_LIBRARIES})
else(ADOLC_FOUND)
  ei_add_property(EIGEN_MISSING_BACKENDS "Adolc, ")
endif(ADOLC_FOUND)

# this test seems to never have been successful on x87, so is considered to contain a FP-related bug.
# see thread: "non-linear optimization test summary"
#ei_add_test(NonLinearOptimization)

ei_add_test(NumericalDiff)
ei_add_test(autodiff)
ei_add_test(BVH)
ei_add_test(matrix_exponential)
ei_add_test(matrix_function)
ei_add_test(alignedvector3)
ei_add_test(FFT)

find_package(MPFR 2.3.0)
if(MPFR_FOUND)
  include_directories(${MPFR_INCLUDES})
  ei_add_property(EIGEN_TESTED_BACKENDS "MPFR C++, ")
  ei_add_test(mpreal_support " " ${MPFR_LIBRARIES} )
else()
  ei_add_property(EIGEN_MISSING_BACKENDS "MPFR C++, ")
endif()

ei_add_test(sparse_llt  " " "${SPARSE_LIBS}")
ei_add_test(sparse_ldlt " " "${SPARSE_LIBS}")
ei_add_test(sparse_lu   " " "${SPARSE_LIBS}")
ei_add_test(sparse_extra   " " " ")

find_package(FFTW)
if(FFTW_FOUND)
  ei_add_property(EIGEN_TESTED_BACKENDS "fftw, ")
  ei_add_test(FFTW  "-DEIGEN_FFTW_DEFAULT " "-lfftw3 -lfftw3f -lfftw3l" )
else()
  ei_add_property(EIGEN_MISSING_BACKENDS "fftw, ")
endif()

option(EIGEN_TEST_NO_OPENGL "Disable OpenGL support in unit tests" OFF)
if(NOT EIGEN_TEST_NO_OPENGL)
find_package(OpenGL)
find_package(GLUT)
find_package(GLEW)
if(OPENGL_FOUND AND GLUT_FOUND AND GLEW_FOUND)
  ei_add_property(EIGEN_TESTED_BACKENDS "OpenGL, ")
  set(EIGEN_GL_LIB ${GLUT_LIBRARIES} ${GLEW_LIBRARIES})
  ei_add_test(openglsupport  "" "${EIGEN_GL_LIB}" )
else()
  ei_add_property(EIGEN_MISSING_BACKENDS "OpenGL, ")
endif()
endif(NOT EIGEN_TEST_NO_OPENGL)

find_package(GSL)
if(GSL_FOUND AND GSL_VERSION_MINOR LESS 9)
  set(GSL_FOUND "")
endif(GSL_FOUND AND GSL_VERSION_MINOR LESS 9)
if(GSL_FOUND)
  add_definitions("-DHAS_GSL" ${GSL_DEFINITIONS})
  include_directories(${GSL_INCLUDE_DIR})
  ei_add_property(EIGEN_TESTED_BACKENDS "GSL, ")
else(GSL_FOUND)
  ei_add_property(EIGEN_MISSING_BACKENDS "GSL, ")
  set(GSL_LIBRARIES " ")
endif(GSL_FOUND)

ei_add_test(polynomialutils)
ei_add_test(polynomialsolver " " "${GSL_LIBRARIES}" )

