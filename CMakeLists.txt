cmake_minimum_required (VERSION 3.19 FATAL_ERROR)

find_package(cetmodules)

project(event-display VERSION 7.0.1)

include(CetCMakeEnv)
cet_cmake_env()
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

cet_set_compiler_flags(DIAGS VIGILANT
 WERROR
 NO_UNDEFINED
 EXTRA_FLAGS -pedantic
 -Wno-unused-parameter -Wno-type-limits -Wno-ignored-qualifiers -Wno-unused-but-set-parameter -Wno-extra -Wno-non-virtual-dtor -Wno-implicit-fallthrough -Wno-deprecated-copy
)

file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/EventDisplay SYMBOLIC)


find_package(Offline REQUIRED EXPORT)
find_package(art_root_io REQUIRED EXPORT)
find_package(ROOT REQUIRED EXPORT)
find_package(nlohmann_json REQUIRED EXPORT)
get_property(importTargets DIRECTORY "${CMAKE_SOURCE_DIR}" PROPERTY IMPORTED_TARGETS)
message("---> Found importtargets ${importTargets}")

include(ArtDictionary)
include(BuildPlugins)


add_subdirectory(src)

install_headers(SUBDIRNAME EventDisplay SUBDIRS inc)

#install_headers(SUBDIRNAME EventDisplay SUBDIRS config)
#install_headers(SUBDIRNAME EventDisplay SUBDIRS gdml)
#install_headers(SUBDIRNAME EventDisplay SUBDIRS tutorial)
#install(DIRECTORY utils DESTINATION python)
#install(DIRECTORY tutorial doc example-analysis-scripts validation DESTINATION ${CMAKE_INSTALL_DATAROOTDIR})
install(DIRECTORY config gdml tutorial DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/EventDisplay)

install_source(SUBDIRS src)
install_headers(NO_RELATIVE_SUBDIR SUBDIRS inc SUBDIRNAME EventDisplay)
install_fhicl(SUBDIRS fcl SUBDIRNAME EventDisplay/fcl)
install_fhicl(SUBDIRS examples SUBDIRNAME EventDisplay/examples)

cet_cmake_config()
