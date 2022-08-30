#compile the project as a standalone application
add_executable(OpacityTests CodeTests.cxx Conductivity.hxx FileIO.hxx Opacity.hxx roots.hxx allocators.hxx)
#add_library(OpacityTools Conductivity.hxx FileIO.hxx Opacity.hxx roots.hxx)
