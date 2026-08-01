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

add_library(topoprocessor_mesh STATIC mesh.cpp mesh_adapters.cpp)
target_include_directories(topoprocessor_mesh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(
	topoprocessor_mesh PUBLIC $<$<CXX_COMPILER_ID:MSVC>:/EHsc>)
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
	set(TOPOPROCESSOR_PYTHON_PACKAGE_DIR
		${CMAKE_CURRENT_BINARY_DIR}/python/topoprocessor)
	file(MAKE_DIRECTORY ${TOPOPROCESSOR_PYTHON_PACKAGE_DIR})
	configure_file(
		${CMAKE_CURRENT_SOURCE_DIR}/../python/topoprocessor/__init__.py
		${TOPOPROCESSOR_PYTHON_PACKAGE_DIR}/__init__.py COPYONLY)
	configure_file(
		${CMAKE_CURRENT_SOURCE_DIR}/../python/topoprocessor/adapters.py
		${TOPOPROCESSOR_PYTHON_PACKAGE_DIR}/adapters.py COPYONLY)
	set_target_properties(
		topoprocessor_python
		PROPERTIES
			OUTPUT_NAME _core
			LIBRARY_OUTPUT_DIRECTORY ${TOPOPROCESSOR_PYTHON_PACKAGE_DIR}
			RUNTIME_OUTPUT_DIRECTORY ${TOPOPROCESSOR_PYTHON_PACKAGE_DIR})
	foreach(configuration DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
		set_target_properties(
			topoprocessor_python
			PROPERTIES
				LIBRARY_OUTPUT_DIRECTORY_${configuration}
					${TOPOPROCESSOR_PYTHON_PACKAGE_DIR}
				RUNTIME_OUTPUT_DIRECTORY_${configuration}
					${TOPOPROCESSOR_PYTHON_PACKAGE_DIR})
	endforeach()
	install(
		TARGETS topoprocessor_python
		LIBRARY DESTINATION topoprocessor
		RUNTIME DESTINATION topoprocessor)
	install(
		FILES
			${CMAKE_CURRENT_SOURCE_DIR}/../python/topoprocessor/__init__.py
			${CMAKE_CURRENT_SOURCE_DIR}/../python/topoprocessor/adapters.py
		DESTINATION topoprocessor)
endif()

if (BUILD_TESTING)
	find_package(Python3 REQUIRED COMPONENTS Interpreter)
	add_executable(test_mesh ../tests/test_mesh.cpp)
	target_link_libraries(test_mesh PRIVATE topoprocessor_mesh)
	add_test(NAME mesh_core COMMAND test_mesh)
	add_executable(test_mesh_adapters ../tests/test_mesh_adapters.cpp)
	target_link_libraries(test_mesh_adapters PRIVATE topoprocessor_mesh)
	add_test(
		NAME mesh_adapters
		COMMAND test_mesh_adapters ${CMAKE_CURRENT_SOURCE_DIR}/../tests/fixtures
	)
	add_test(
		NAME recovered_mesh_invariants
		COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_recovered_meshes.py
	)
	add_test(
		NAME python_in_memory_adapters
		COMMAND ${Python3_EXECUTABLE}
			${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_in_memory_adapters.py
	)
	if (TOPOPROCESSOR_BUILD_PYTHON)
		add_test(
			NAME python_mesh_bindings
			COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_python_bindings.py
		)
		set_tests_properties(
			python_in_memory_adapters python_mesh_bindings
			PROPERTIES ENVIRONMENT
				"PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/python"
		)
	endif()
endif()
