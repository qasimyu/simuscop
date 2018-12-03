cmake_minimum_required(VERSION 2.8)
project(SimuSCoP)

SET(CMAKE_BUILD_TYPE Release)

add_definitions(-std=c++11)

add_subdirectory(lib)

include_directories (
	${SimuSCoP_SOURCE_DIR}/lib/config
	${SimuSCoP_SOURCE_DIR}/lib/fastahack
	${SimuSCoP_SOURCE_DIR}/lib/genome
	${SimuSCoP_SOURCE_DIR}/lib/matrix
	${SimuSCoP_SOURCE_DIR}/lib/profile
	${SimuSCoP_SOURCE_DIR}/lib/mydefine
	${SimuSCoP_SOURCE_DIR}/lib/segment
	${SimuSCoP_SOURCE_DIR}/lib/seqwriter
	${SimuSCoP_SOURCE_DIR}/lib/snp
	${SimuSCoP_SOURCE_DIR}/lib/split
	${SimuSCoP_SOURCE_DIR}/lib/threadpool
	${SimuSCoP_SOURCE_DIR}/lib/psifunc
)

add_executable(seqToProfile src/seqToProfile.cpp)
add_executable(simuReads src/simuReads.cpp)

target_link_libraries(seqToProfile genome profile mydefine)
target_link_libraries(simuReads config genome profile mydefine split threadpool)

SET_TARGET_PROPERTIES(seqToProfile PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${SimuSCoP_BINARY_DIR}/bin)
SET_TARGET_PROPERTIES(simuReads PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${SimuSCoP_BINARY_DIR}/bin)

find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
	set_property(TARGET seqToProfile PROPERTY COMPILE_OPTIONS "-pthread")
	set_property(TARGET seqToProfile PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
	set_property(TARGET simuReads PROPERTY COMPILE_OPTIONS "-pthread")
	set_property(TARGET simuReads PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
	target_link_libraries(seqToProfile "${CMAKE_THREAD_LIBS_INIT}")
	target_link_libraries(simuReads "${CMAKE_THREAD_LIBS_INIT}")
endif()
