### the directory name
set(directory include/OpenMS/ANALYSIS/OPENSWATH)

### list all header files of the directory here
set(sources_list_h
  CachedmzML.h
  DIAHelper.h
  MRMDecoy.h
  MRMTransitionGroupPicker.h
  SpectrumAddition.h
  ChromatogramExtractor.h
  DIAPrescoring.h
  MRMFeatureFinderScoring.h
  TransitionTSVReader.h
  DIAScoring.h
  MRMRTNormalizer.h
  OpenSwathHelper.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\ANALYSIS\\OPENSWATH" FILES ${sources_h})
set_source_files_properties(${directory}/sources.cmake PROPERTIES HEADER_FILE_ONLY TRUE)

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

