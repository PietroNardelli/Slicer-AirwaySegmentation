cmake_minimum_required(VERSION 2.8.9)

project(Slicer-AirwaySegmentation)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/AirwaySegmentation")
set(EXTENSION_CATEGORY "Segmentation")
set(EXTENSION_CONTRIBUTORS "Pietro Nardelli (University College Cork)")
set(EXTENSION_DESCRIPTION "CLI module for airway segmentation of chest CT images")
set(EXTENSION_ICONURL "http://www.slicer.org/slicerWiki/images/7/7d/AirwaySegmentation-ICON.png")
set(EXTENSION_SCREENSHOTURLS "http://www.slicer.org/slicerWiki/images/2/24/Slicer4-AirwaySegmentation-GUI.png")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(AirwaySegmentationCLI)
add_subdirectory(AirwaySegmentation)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
