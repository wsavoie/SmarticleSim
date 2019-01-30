#=============================================================================
# CMake configuration file for Chrono-Projects
# 
#=============================================================================

cmake_minimum_required(VERSION 2.8)

project(Smarticles)

# Set location of additional CMAKE modules
#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Set location of executable
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})


#-----------------------------------------------------------------------------
# Always use full RPATH (differentiating between the build and install trees)
#-----------------------------------------------------------------------------

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

#-----------------------------------------------------------------------------
# Force C++11
#-----------------------------------------------------------------------------

#if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#    set(CH_LINKERFLAG_EXE  "${CH_LINKERFLAG_EXE} -framework IOKit -framework Cocoa -framework OpenGL")
#endif()

#-----------------------------------------------------------------------------
# Disable some warnings
#-----------------------------`------------------------------------------------

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    if(MSVC)
        add_definitions( "-D_CRT_SECURE_NO_DEPRECATE" )  # avoids deprecation warnings
        add_definitions( "-D_SCL_SECURE_NO_DEPRECATE" )  # avoids deprecation warnings
        add_definitions( "-DNOMINMAX" ) # do not use MSVC's min/max macros
        add_definitions( "-MP" ) # perform parallel builds
    endif(MSVC)
endif()
#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(DEMOS
	#SmarticlesSystem
	SmarticlesSystem_Dynamic
)

#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode.


find_package(Chrono
             # COMPONENTS Parallel
             # OPTIONAL_COMPONENTS OpenGL
             COMPONENTS Irrlicht
             CONFIG
		
)

# If Chrono and/or the required component(s) were not found, return now.

if(NOT Chrono_FOUND)
  message("Could not find requirements for PARALLEL projects")
  return()
endif()

#--------------------------------------------------------------
# Include paths and libraries
#--------------------------------------------------------------

# (A) Path to the Chrono include headers
# - If using an installed version of Chrono, this will be the path 
#   to the installed headers (the configuration headers are also
#   available there)
# - If using a build version of Chrono, this will contain both the
#   path to the Chrono sources and the path to the chrono BUILD tree
#   (the latter for the configuration headers)
# 
# (B) Path to the top of the source tree for this project
# - for access to utility headers
set(JSON_LIB_PATH "" CACHE FILEPATH "json lib path")
include_directories(

    ${CHRONO_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}
	${JSON_LIB_PATH}
	
)

#--------------------------------------------------------------
# Append to the parent's list of DLLs (and make it visible up)
#--------------------------------------------------------------

list(APPEND ALL_DLLS "${CHRONO_DLLS}")
set(ALL_DLLS "${ALL_DLLS}" PARENT_SCOPE)
set(CHRONO_BIN_PATH "D:/ChronoCode/chronoPkgs/ChronoBuild/bin/" CACHE FILEPATH "chrono build path")

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

message(STATUS "Parallel programs...")

foreach(PROGRAM ${DEMOS})

  message(STATUS "...add ${PROGRAM}")
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.h*")
  add_executable(${PROGRAM}  "${PROGRAM}.cpp" ${HEADERS} ${SOURCES})
  source_group(""  FILES "${PROGRAM}.cpp")

  set_target_properties(${PROGRAM} PROPERTIES
    FOLDER demos
    COMPILE_FLAGS "${CHRONO_CXX_FLAGS}"
    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
    LINK_FLAGS "${CHRONO_CXX_FLAGS} ${CHRONO_LINKER_FLAGS}"
  )

  target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES})
  
  #install(TARGETS ${PROGRAM} DESTINATION bin)

endforeach(PROGRAM)

#--------------------------------------------------------------
# Copy DLLs
#--------------------------------------------------------------

# On Windows, traverse the list of DLLs and copy them next to
# the binaries.

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# set specifically only release and debug build types
# set(CMAKE_CONFIGURATION_TYPES Release Debug CACHE TYPE INTERNAL FORCE )


  # Remove duplicates in the list of DLLs
  list(REMOVE_DUPLICATES ALL_DLLS)

  message(STATUS "==== Add custom commands for copying DLLs ====")
  message(STATUS "DLL list...")

  # Create custom target for copying DLLs; add it to the default build target
  add_custom_target(COPY_DLLS ALL)

  # Loop over the list of all DLLs and create a custom command
  
	
  	file(GLOB ReleaseDLL "${Chrono_DIR}/../bin/Release/*.dll")
	file(GLOB DebugDLL "${Chrono_DIR}/../bin/Debug/*.dll")
	file(GLOB RelWithDebInfoDLL "${Chrono_DIR}/../bin/RelWithDebInfo/*.dll")
	list(REMOVE_DUPLICATES ReleaseDLL)
	list(REMOVE_DUPLICATES DebugDLL)
	list(REMOVE_DUPLICATES RelWithDebInfoDLL)
	message(STATUS ${DebugDLL})

	
	# add_custom_command(
        # TARGET COPY_DLLS POST_BUILD
        # COMMAND ${CMAKE_COMMAND} -E copy_if_different
            # "${sDLL}"
            # "${CMAKE_BINARY_DIR}/$<CONFIGURATION>"			
	# )
	#have to use one a specific build here since I don't know how to add this to the post build command below
	# file(GLOB DLLS "${CHRONO_BIN_PATH}/Release/*.dll")
	
	foreach(DLL ${CHRONO_DLL_NAMES})
	add_custom_command(
        TARGET COPY_DLLS POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CHRONO_BIN_PATH}/$<CONFIGURATION>/${DLL}"
            "${CMAKE_BINARY_DIR}/$<CONFIGURATION>/${DLL}"			
	)
	endforeach()
# message(STATUS "\nParallel programs... " ${CMAKE_BINARY_DIR}/$<CONFIGURATION>/${DLL}".dll"	)
endif()
