cmake_minimum_required (VERSION 3.13)
project(topoprocessor)

include(CTest)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if ( CYGWIN )
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++1y -fdiagnostics-color=always \
-DENABLE_CYGWIN_MATH_FIX -fpermissive -Wwrite-strings")
endif ( CYGWIN )

# -Wsign-conversion

if ( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -std=c++1y -fdiagnostics-color=always \
	-fpermissive -Wc++11-extensions")
	set(LINK_LIBS ${LINK_LIBS} -lpthread) # hack, must be fixed
endif( )

if ( ${CMAKE_SYSTEM_NAME} MATCHES "Apple" )
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fdiagnostics-color=always \
	-fpermissive")
	set(LINK_LIBS ${LINK_LIBS})
endif( )

set(SOURCE lean_cohomology.cpp generic_two_manifold.cpp mapped_file.cpp timecounter.cpp)

include_directories(./)

add_library(topoprocessor_mesh STATIC mesh.cpp)
target_include_directories(topoprocessor_mesh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(topoprocessor_mesh PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_executable(topoprocessor ./main/main.cpp ${SOURCE})
target_link_libraries(topoprocessor PRIVATE topoprocessor_mesh ${LINK_LIBS})
install (TARGETS topoprocessor DESTINATION .)

option(TOPOPROCESSOR_BUILD_PYTHON "Build the Python mesh bindings" OFF)

if (TOPOPROCESSOR_BUILD_PYTHON)
	find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
	find_package(pybind11 CONFIG REQUIRED)
	pybind11_add_module(topoprocessor_python MODULE python_bindings.cpp)
	target_link_libraries(topoprocessor_python PRIVATE topoprocessor_mesh)
	set_target_properties(topoprocessor_python PROPERTIES OUTPUT_NAME topoprocessor)
	install(TARGETS topoprocessor_python LIBRARY DESTINATION . RUNTIME DESTINATION .)
endif()

if (BUILD_TESTING)
	find_package(Python3 REQUIRED COMPONENTS Interpreter)
	add_executable(test_mesh ../tests/test_mesh.cpp)
	target_link_libraries(test_mesh PRIVATE topoprocessor_mesh)
	add_test(NAME mesh_core COMMAND test_mesh)
	add_test(
		NAME recovered_mesh_invariants
		COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_recovered_meshes.py
	)
	if (TOPOPROCESSOR_BUILD_PYTHON)
		add_test(
			NAME python_mesh_bindings
			COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_python_bindings.py
		)
		set_tests_properties(
			python_mesh_bindings
			PROPERTIES ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:topoprocessor_python>"
		)
	endif()
endif()
