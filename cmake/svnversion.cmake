if (ENABLE_SVN)
	if (OPENMS_HAS_SVNVERSION)
		execute_process(COMMAND ${SVNVERSION_EXECUTABLE} -n ${SOURCE_DIR} OUTPUT_VARIABLE SVNVERSION)
	endif()
	if(NOT SVNVERSION)
		set(SVNVERSION "unknown")
	endif()
	file(WRITE ${SVN_REVISION_FILE}.tmp "// generated by the build system! do not edit this manually\n#define OPENMS_SVN_REVISION \"${SVNVERSION}\"\n")
	execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SVN_REVISION_FILE}.tmp ${SVN_REVISION_FILE})
	file(REMOVE ${SVN_REVISION_FILE}.tmp)
else()
	# write the file once:
	if (NOT EXISTS ${SVN_REVISION_FILE})
		set(SVNVERSION "unknown")
		file(WRITE ${SVN_REVISION_FILE} "// generated by the build system! do not edit this manually\n#define OPENMS_SVN_REVISION \"${SVNVERSION}\"\n")
	endif()
endif()
