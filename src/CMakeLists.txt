cmake_minimum_required (VERSION 3.13)
project(topoprocessor)

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

add_executable(topoprocessor ./main/main.cpp ${SOURCE})
target_link_libraries(topoprocessor ${LINK_LIBS})
install (TARGETS topoprocessor)
#install (TARGETS topoprocessor DESTINATION ${CMAKE_INSTALL_PREFIX})
