set(Launcher_Command ${Slicer_LAUNCH_COMMAND})

set(TEMP ${T1Mapping_BINARY_DIR}/Testing/Temporary)

set(CLP ${MODULE_NAME})

#-----------------------------------------------------------------------------
add_executable(${CLP}Test ${CLP}Test.cxx)
target_link_libraries(${CLP}Test ${CLP}Lib ${SlicerExecutionModel_EXTRA_EXECUTABLE_TARGET_LIBRARIES})
set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})

#-----------------------------------------------------------------------------
set(testName ${CLP}_TestHelp)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPoint --help
  )
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
# Regression Tests QIBA001
#-----------------------------------------------------------------------------
set(inputDataBaseName ${CMAKE_SOURCE_DIR}/T1Mapping/Data/RegressionTests/QIBA001/Input/QIBA001)
set(referenceDataBaseDir ${CMAKE_SOURCE_DIR}/T1Mapping/Data/RegressionTests/QIBA001/Reference/)

set(testName QIBA001_AllFlipAngles)
set(tempOutDataBaseName ${TEMP}/${testName})
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  --compareIntensityTolerance 1e-6
  --compare ${referenceDataBaseDir}${testName}-T1VFA.nrrd ${tempOutDataBaseName}-T1VFA.nrrd
  --compare ${referenceDataBaseDir}${testName}-rsqrVolume.nrrd ${tempOutDataBaseName}-rsqrVolume.nrrd
  --compare ${referenceDataBaseDir}${testName}-fittedVolume.nrrd ${tempOutDataBaseName}-fittedVolume.nrrd
  ModuleEntryPoint
    --T1VFA ${tempOutDataBaseName}-T1VFA.nrrd
    --rsqrVolume ${tempOutDataBaseName}-rsqrVolume.nrrd
    --fittedVolume ${tempOutDataBaseName}-fittedVolume.nrrd
    ${inputDataBaseName}-phantom.nrrd
  )
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

set(testName QIBA001_InclFlipAngles)
set(tempOutDataBaseName ${TEMP}/${testName})
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  --compareIntensityTolerance 1e-6
  --compare ${referenceDataBaseDir}${testName}-T1VFA.nrrd ${tempOutDataBaseName}-T1VFA.nrrd
  --compare ${referenceDataBaseDir}${testName}-rsqrVolume.nrrd ${tempOutDataBaseName}-rsqrVolume.nrrd
  --compare ${referenceDataBaseDir}${testName}-fittedVolume.nrrd ${tempOutDataBaseName}-fittedVolume.nrrd
  ModuleEntryPoint
    --FAInclude "3, 6, 9"
    --T1VFA ${tempOutDataBaseName}-T1VFA.nrrd
    --rsqrVolume ${tempOutDataBaseName}-rsqrVolume.nrrd
    --fittedVolume ${tempOutDataBaseName}-fittedVolume.nrrd
    ${inputDataBaseName}-phantom.nrrd
  )
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

set(testName QIBA001_ExclFlipAngles)
set(tempOutDataBaseName ${TEMP}/${testName})
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  --compareIntensityTolerance 1e-6
  --compare ${referenceDataBaseDir}${testName}-T1VFA.nrrd ${tempOutDataBaseName}-T1VFA.nrrd
  --compare ${referenceDataBaseDir}${testName}-rsqrVolume.nrrd ${tempOutDataBaseName}-rsqrVolume.nrrd
  --compare ${referenceDataBaseDir}${testName}-fittedVolume.nrrd ${tempOutDataBaseName}-fittedVolume.nrrd
  ModuleEntryPoint
    --FAExclude "3, 6, 9"
    --T1VFA ${tempOutDataBaseName}-T1VFA.nrrd
    --rsqrVolume ${tempOutDataBaseName}-rsqrVolume.nrrd
    --fittedVolume ${tempOutDataBaseName}-fittedVolume.nrrd
    ${inputDataBaseName}-phantom.nrrd
  )
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

set(testName QIBA001_MaskedInput)
set(tempOutDataBaseName ${TEMP}/${testName})
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  --compareIntensityTolerance 1e-6
  --compare ${referenceDataBaseDir}${testName}-T1VFA.nrrd ${tempOutDataBaseName}-T1VFA.nrrd
  --compare ${referenceDataBaseDir}${testName}-rsqrVolume.nrrd ${tempOutDataBaseName}-rsqrVolume.nrrd
  --compare ${referenceDataBaseDir}${testName}-fittedVolume.nrrd ${tempOutDataBaseName}-fittedVolume.nrrd
  ModuleEntryPoint
    --T1VFA ${tempOutDataBaseName}-T1VFA.nrrd
    --rsqrVolume ${tempOutDataBaseName}-rsqrVolume.nrrd
    --fittedVolume ${tempOutDataBaseName}-fittedVolume.nrrd
    --mask ${inputDataBaseName}-mask-label.nrrd
    ${inputDataBaseName}-phantom.nrrd
  )
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

set(testName QIBA001_FailOnInclAndExclFlipAngles)
set(tempOutDataBaseName ${TEMP}/${testName})
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPointExpectFail
    --FAInclude "3, 6"
    --FAExclude "9, 15"
    --T1VFA ${tempOutDataBaseName}-T1VFA.nrrd
    --rsqrVolume ${tempOutDataBaseName}-rsqrVolume.nrrd
    --fittedVolume ${tempOutDataBaseName}-fittedVolume.nrrd
    ${inputDataBaseName}-phantom.nrrd
  )
set_property(TEST ${testName} PROPERTY LABELS ${CLP})


#-----------------------------------------------------------------------------











