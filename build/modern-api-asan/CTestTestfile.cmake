# CMake generated Testfile for 
# Source directory: C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src
# Build directory: C:/Users/Bernard Kapidani/repos/topoprocessor/build/modern-api-asan
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(mesh_core "C:/Users/Bernard Kapidani/repos/topoprocessor/build/modern-api-asan/test_mesh.exe")
set_tests_properties(mesh_core PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/CMakeLists.txt;98;add_test;C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/CMakeLists.txt;0;")
add_test(mesh_adapters "C:/Users/Bernard Kapidani/repos/topoprocessor/build/modern-api-asan/test_mesh_adapters.exe" "C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/../tests/fixtures")
set_tests_properties(mesh_adapters PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/CMakeLists.txt;101;add_test;C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/CMakeLists.txt;0;")
add_test(topology_core "C:/Users/Bernard Kapidani/repos/topoprocessor/build/modern-api-asan/test_topology.exe")
set_tests_properties(topology_core PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/CMakeLists.txt;107;add_test;C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/CMakeLists.txt;0;")
add_test(recovered_topology "C:/Users/Bernard Kapidani/repos/topoprocessor/build/modern-api-asan/test_recovered_topology.exe" "C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/../tests/fixtures")
set_tests_properties(recovered_topology PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/CMakeLists.txt;111;add_test;C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/CMakeLists.txt;0;")
add_test(mapped_file "C:/Users/Bernard Kapidani/repos/topoprocessor/build/modern-api-asan/test_mapped_file.exe" "C:/Users/Bernard Kapidani/repos/topoprocessor/build/modern-api-asan")
set_tests_properties(mapped_file PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/CMakeLists.txt;118;add_test;C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/CMakeLists.txt;0;")
add_test(recovered_mesh_invariants "C:/Users/Bernard Kapidani/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe" "C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/../tests/test_recovered_meshes.py")
set_tests_properties(recovered_mesh_invariants PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/CMakeLists.txt;122;add_test;C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/CMakeLists.txt;0;")
add_test(python_in_memory_adapters "C:/Users/Bernard Kapidani/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe" "C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/../tests/test_in_memory_adapters.py")
set_tests_properties(python_in_memory_adapters PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/CMakeLists.txt;126;add_test;C:/Users/Bernard Kapidani/AppData/Local/Temp/topoprocessor-publish-in-memory-portable-io/src/CMakeLists.txt;0;")
