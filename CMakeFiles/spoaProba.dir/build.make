# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/mlipovac/Matea/zavrsniRad

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/mlipovac/Matea/zavrsniRad

# Include any dependencies generated for this target.
include CMakeFiles/spoaProba.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/spoaProba.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/spoaProba.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/spoaProba.dir/flags.make

CMakeFiles/spoaProba.dir/main.cpp.o: CMakeFiles/spoaProba.dir/flags.make
CMakeFiles/spoaProba.dir/main.cpp.o: main.cpp
CMakeFiles/spoaProba.dir/main.cpp.o: CMakeFiles/spoaProba.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/mlipovac/Matea/zavrsniRad/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/spoaProba.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/spoaProba.dir/main.cpp.o -MF CMakeFiles/spoaProba.dir/main.cpp.o.d -o CMakeFiles/spoaProba.dir/main.cpp.o -c /Users/mlipovac/Matea/zavrsniRad/main.cpp

CMakeFiles/spoaProba.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/spoaProba.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mlipovac/Matea/zavrsniRad/main.cpp > CMakeFiles/spoaProba.dir/main.cpp.i

CMakeFiles/spoaProba.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/spoaProba.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mlipovac/Matea/zavrsniRad/main.cpp -o CMakeFiles/spoaProba.dir/main.cpp.s

# Object files for target spoaProba
spoaProba_OBJECTS = \
"CMakeFiles/spoaProba.dir/main.cpp.o"

# External object files for target spoaProba
spoaProba_EXTERNAL_OBJECTS =

spoaProba: CMakeFiles/spoaProba.dir/main.cpp.o
spoaProba: CMakeFiles/spoaProba.dir/build.make
spoaProba: spoa/lib/libspoa.a
spoaProba: alignment/libalignment.a
spoaProba: /Library/Developer/CommandLineTools/SDKs/MacOSX11.3.sdk/usr/lib/libz.tbd
spoaProba: CMakeFiles/spoaProba.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/mlipovac/Matea/zavrsniRad/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable spoaProba"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/spoaProba.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/spoaProba.dir/build: spoaProba
.PHONY : CMakeFiles/spoaProba.dir/build

CMakeFiles/spoaProba.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/spoaProba.dir/cmake_clean.cmake
.PHONY : CMakeFiles/spoaProba.dir/clean

CMakeFiles/spoaProba.dir/depend:
	cd /Users/mlipovac/Matea/zavrsniRad && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mlipovac/Matea/zavrsniRad /Users/mlipovac/Matea/zavrsniRad /Users/mlipovac/Matea/zavrsniRad /Users/mlipovac/Matea/zavrsniRad /Users/mlipovac/Matea/zavrsniRad/CMakeFiles/spoaProba.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/spoaProba.dir/depend

