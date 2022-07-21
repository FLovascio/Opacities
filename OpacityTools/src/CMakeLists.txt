#compile the project as a standalone application
add_executable(OpacityTests CodeTests.cxx Conductivity.hxx FileIO.hxx Opacity.hxx roots.hxx)
#add_library(OpacityTools Conductivity.hxx FileIO.hxx Opacity.hxx roots.hxx)
