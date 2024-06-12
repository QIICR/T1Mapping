cmake_minimum_required(VERSION 2.8.9)

project(T1Mapping)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/T1Mapping")
set(EXTENSION_CATEGORY "Quantification")
set(EXTENSION_CONTRIBUTORS "Xiao Da(MGH)")
set(EXTENSION_DESCRIPTION "T1 mapping estimates effective tissue parameter maps (T1) from multi-spectral FLASH MRI scans with different flip angles.")
set(EXTENSION_ICONURL "https://slicer.org/slicerWiki/images/3/32/T1_Mapping_Logo_Resized.png")
set(EXTENSION_SCREENSHOTURLS "https://slicer.org/slicerWiki/images/6/68/T1_Mapping_CPP_GUI.png")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(T1Mapping)
## NEXT_MODULE

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()



