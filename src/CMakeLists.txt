cmake_minimum_required(VERSION 3.15)
project(topoprocessor LANGUAGES CXX)

include(CTest)
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(topoprocessor_mesh STATIC mesh.cpp mesh_adapters.cpp topology.cpp)
add_library(topoprocessor::core ALIAS topoprocessor_mesh)
target_include_directories(topoprocessor_mesh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(
	topoprocessor_mesh PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/EHsc>)
set_target_properties(topoprocessor_mesh PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(topoprocessor_io STATIC mapped_file.cpp)
target_include_directories(topoprocessor_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(
	topoprocessor_io PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/EHsc>)

option(TOPOPROCESSOR_BUILD_LEGACY_CLI "Build the historical CLI" ON)
if (TOPOPROCESSOR_BUILD_LEGACY_CLI)
	find_package(Threads REQUIRED)
	add_executable(
		topoprocessor
		./main/main.cpp lean_cohomology.cpp generic_two_manifold.cpp
		timecounter.cpp)
	target_link_libraries(
		topoprocessor
		PRIVATE topoprocessor_mesh topoprocessor_io Threads::Threads)
	target_include_directories(topoprocessor PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
	target_compile_options(
		topoprocessor PRIVATE
			$<$<CXX_COMPILER_ID:GNU>:-fpermissive>)
	target_compile_definitions(
		topoprocessor PRIVATE $<$<BOOL:${CYGWIN}>:ENABLE_CYGWIN_MATH_FIX>)
	install(TARGETS topoprocessor RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
install(
	TARGETS topoprocessor_mesh
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
	FILES mesh.hpp mesh_adapters.hpp topology.hpp
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/topoprocessor)

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
	configure_file(
		${CMAKE_CURRENT_SOURCE_DIR}/../python/topoprocessor/loaders.py
		${TOPOPROCESSOR_PYTHON_PACKAGE_DIR}/loaders.py COPYONLY)
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
			${CMAKE_CURRENT_SOURCE_DIR}/../python/topoprocessor/loaders.py
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
	add_executable(test_topology ../tests/test_topology.cpp)
	target_link_libraries(test_topology PRIVATE topoprocessor_mesh)
	add_test(NAME topology_core COMMAND test_topology)
	add_executable(
		test_recovered_topology ../tests/test_recovered_topology.cpp)
	target_link_libraries(test_recovered_topology PRIVATE topoprocessor_mesh)
	add_test(
		NAME recovered_topology
		COMMAND test_recovered_topology
			${CMAKE_CURRENT_SOURCE_DIR}/../tests/fixtures
	)
	add_executable(test_mapped_file ../tests/test_mapped_file.cpp)
	target_link_libraries(test_mapped_file PRIVATE topoprocessor_io)
	add_test(
		NAME mapped_file
		COMMAND test_mapped_file ${CMAKE_CURRENT_BINARY_DIR}
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
			NAME python_current_mesher_files
			COMMAND ${Python3_EXECUTABLE}
				${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_current_mesher_files.py
		)
		add_test(
			NAME python_mesh_bindings
			COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_python_bindings.py
		)
		set_tests_properties(
			python_in_memory_adapters python_current_mesher_files
			python_mesh_bindings
			PROPERTIES ENVIRONMENT
				"PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/python"
		)
	endif()
endif()
