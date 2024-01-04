IF (DEFINED VCPKG_INSTALLED_DIR)
  # List `vcpkg` dependencies from VCPKG_INSTALLED_DIR so we can install them
  # together with the python API.
  FILE(
    GLOB VCPKG_DEPENDENCIES
    RELATIVE "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin/"
    "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin/*"
  )
ENDIF()

IF ( CGAL_FOUND )
    SET(
      pcsr_cgal_SRCS
      pcsr_cgal.cpp
      pcsr_pcl.cpp
      pcsr_cgal_py.cpp
    )

    ADD_LIBRARY(_pcsr_cgal MODULE ${pcsr_cgal_SRCS})
    TARGET_LINK_LIBRARIES(
      _pcsr_cgal
      ${CGAL_LIBRARIES}
      ${BOOST_LIBRARIES}
      ${PYTHON_LIBRARIES}
      ${MPFR_LIBRARIES}
      ${GMP_LIBRARIES}
      ${TBB_LIBRARIES}
    )
    PYTHON_EXTENSION_MODULE(_pcsr_cgal)

    # Install the `vcpkg` dependencies that are used in our ${PROJECT_NAME}
    # library in the same DESTINATION as the python API. Note that it
    # only supports collecting the runtime dependencies for Windows, Linux and
    # macOS platforms.
    IF (DEFINED SKBUILD)
      INSTALL(
        TARGETS _pcsr_cgal
        RUNTIME_DEPENDENCIES
        PRE_INCLUDE_REGEXES ${VCPKG_DEPENDENCIES}
        PRE_EXCLUDE_REGEXES ".*"
        DESTINATION pcsr
    )
    ELSE()
      MESSAGE(STATUS "PCSR_PYTHON_SITE_PACKAGES_INSTALL_DIR=${PCSR_PYTHON_SITE_PACKAGES_INSTALL_DIR}")
      INSTALL(
        TARGETS _pcsr_cgal
        RUNTIME_DEPENDENCIES
        PRE_INCLUDE_REGEXES ${VCPKG_DEPENDENCIES}
        PRE_EXCLUDE_REGEXES ".*"
        DESTINATION ${PCSR_PYTHON_SITE_PACKAGES_INSTALL_DIR}
      )
    ENDIF()
ENDIF()

IF ( PCL_FOUND )
    SET(
      pcsr_pcl_SRCS
      pcsr_pcl.cpp
    )

    ADD_LIBRARY(_pcsr_pcl MODULE ${pcsr_pcl_SRCS})
    TARGET_LINK_LIBRARIES(
      _pcsr_pcl
      ${PCSR_PCL_LIBRARIES}
      ${BOOST_LIBRARIES}
      ${PYTHON_LIBRARIES}
    )
    PYTHON_EXTENSION_MODULE(_pcsr_pcl)

    IF (DEFINED SKBUILD)
      INSTALL(
        TARGETS _pcsr_pcl
        RUNTIME_DEPENDENCIES
        PRE_INCLUDE_REGEXES ${VCPKG_DEPENDENCIES}
        PRE_EXCLUDE_REGEXES ".*"
        DESTINATION pcsr
    )
    ELSE()
      MESSAGE(STATUS "PCSR_PYTHON_SITE_PACKAGES_INSTALL_DIR=${PCSR_PYTHON_SITE_PACKAGES_INSTALL_DIR}")
      INSTALL(
        TARGETS _pcsr_pcl
        RUNTIME_DEPENDENCIES
        PRE_INCLUDE_REGEXES ${VCPKG_DEPENDENCIES}
        PRE_EXCLUDE_REGEXES ".*"
        DESTINATION ${PCSR_PYTHON_SITE_PACKAGES_INSTALL_DIR}
      )
    ENDIF()
ENDIF()
