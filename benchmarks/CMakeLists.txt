add_executable(LArNESTMeanYieldsBenchmarks LArNESTMeanYieldsBenchmarks.cxx)
target_link_libraries(LArNESTMeanYieldsBenchmarks PUBLIC LArNEST::Core)

add_executable(LArNESTFluctuationBenchmarks LArNESTFluctuationBenchmarks.cxx)
target_link_libraries(LArNESTFluctuationBenchmarks PUBLIC LArNEST::Core)