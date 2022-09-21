
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
      ${PCSR_PCL_LIBRARIES}
      ${BOOST_LIBRARIES}
      ${PCSR_BOOST_PYTHON_LIBRARIES}
      ${PYTHON_LIBRARIES}
      ${MPFR_LIBRARIES}
      ${GMP_LIBRARIES}
      ${TBB_LIBRARIES}
    )
    PYTHON_EXTENSION_MODULE(_pcsr_cgal)

    INSTALL(TARGETS _pcsr_cgal LIBRARY DESTINATION pcsr)
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
      ${PCSR_BOOST_PYTHON_LIBRARIES}
      ${PYTHON_LIBRARIES}
    )
    PYTHON_EXTENSION_MODULE(_pcsr_pcl)
    INSTALL(TARGETS _pcsr_pcl LIBRARY DESTINATION pcsr)
ENDIF()
