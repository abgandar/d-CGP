# Keplerian_toolbox lib source files.
SET(dCGP_LIB_SRC_LIST
	${CMAKE_CURRENT_SOURCE_DIR}/expression.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/function_set.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/rng.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/wrapped_functions.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/fitness_functions.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/basis_function.cpp
)

#Build Static Library
ADD_LIBRARY(dcgp_s STATIC ${dCGP_LIB_SRC_LIST})

#Build Dynamic Library (only if needed by PyKEP)
#SET(LIB_INSTALL_PATH "lib")
#ADD_LIBRARY(dcgp SHARED ${dCGP_LIB_SRC_LIST})
#TARGET_LINK_LIBRARIES(dcgp ${LIBRARIES_4_DYNAMIC})
#INSTALL(TARGETS dcgp
#		RUNTIME DESTINATION ${LIB_INSTALL_PATH}
#		LIBRARY DESTINATION ${LIB_INSTALL_PATH}
#)



