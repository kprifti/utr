# Register the CMake cache variable TARGET with type string and default value 154Sm if no value is already in the cache
set(TARGET "Sm_No8" CACHE STRING "Set the activation target to be used in the DetectorConstruction")

# Register the valid choices (PROPERTY STRINGS) in ccmake for the CMake cache variable TARGET
set_property(CACHE TARGET PROPERTY STRINGS "MixedSourceAtDet3" "MixedSourceAtDet4" "PointSourceAtDet3" "PointSourceAtDet4" "Ce_No1" "Ce_No2" "Ce_No3" "Ce_No4" "Ce_No5" "Ce_No6" "Ce_No7" "Sm_No8" "Sm_No9" "Sm_No10" "Sm_No11" "Sm_No12" "Sm_No13" "Sm_No14" "Au_No1" "Au_No2" "Au_No3" "Au_No4" "Au_No5" "Au_No6" "Au_No7" "Au_No8" "Au_No9" "Au_No10" "Au_No11" "Au_No12" )

# Define the preprocessor definition TARGET with the value being the CMake variable TARGET's value but with quotes added
# add_compile_definitions(TARGET="${TARGET}")

# Alternative to the preprocessor definition:
# Use configure_file to create the header DetectorConstructionConfig.hh from the template DetectorConstructionConfig.configure_file.hh, which then contains the definition of an proper c++ string variable TARGET
# Advantage: Is only in scope if this header is included and doesn't require a full project recompilation if the CMake variable TARGET is changed (unlike when using a global preprocessor definition).
# Disadvantage: The automatically generated header file pollutes the source directory (one could copy it to some directory in the build directory instead though).
configure_file(
  "${CMAKE_CURRENT_LIST_DIR}/DetectorConstructionConfig.configure_file.hh"
  "${CMAKE_CURRENT_LIST_DIR}/DetectorConstructionConfig.hh"
)
