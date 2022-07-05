#compile the project as a standalone application
add_executable(opacity main.cxx Conductivity.hxx FileIO.hxx Opacity.hxx roots.hxx)
target_compile_features(opacity PRIVATE cxx_std_17)