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
include alignment/CMakeFiles/alignment.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include alignment/CMakeFiles/alignment.dir/compiler_depend.make

# Include the progress variables for this target.
include alignment/CMakeFiles/alignment.dir/progress.make

# Include the compile flags for this target's objects.
include alignment/CMakeFiles/alignment.dir/flags.make

alignment/CMakeFiles/alignment.dir/alignment.cpp.o: alignment/CMakeFiles/alignment.dir/flags.make
alignment/CMakeFiles/alignment.dir/alignment.cpp.o: alignment/alignment.cpp
alignment/CMakeFiles/alignment.dir/alignment.cpp.o: alignment/CMakeFiles/alignment.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/mlipovac/Matea/zavrsniRad/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object alignment/CMakeFiles/alignment.dir/alignment.cpp.o"
	cd /Users/mlipovac/Matea/zavrsniRad/alignment && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT alignment/CMakeFiles/alignment.dir/alignment.cpp.o -MF CMakeFiles/alignment.dir/alignment.cpp.o.d -o CMakeFiles/alignment.dir/alignment.cpp.o -c /Users/mlipovac/Matea/zavrsniRad/alignment/alignment.cpp

alignment/CMakeFiles/alignment.dir/alignment.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/alignment.dir/alignment.cpp.i"
	cd /Users/mlipovac/Matea/zavrsniRad/alignment && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mlipovac/Matea/zavrsniRad/alignment/alignment.cpp > CMakeFiles/alignment.dir/alignment.cpp.i

alignment/CMakeFiles/alignment.dir/alignment.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/alignment.dir/alignment.cpp.s"
	cd /Users/mlipovac/Matea/zavrsniRad/alignment && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mlipovac/Matea/zavrsniRad/alignment/alignment.cpp -o CMakeFiles/alignment.dir/alignment.cpp.s

# Object files for target alignment
alignment_OBJECTS = \
"CMakeFiles/alignment.dir/alignment.cpp.o"

# External object files for target alignment
alignment_EXTERNAL_OBJECTS =

alignment/libalignment.a: alignment/CMakeFiles/alignment.dir/alignment.cpp.o
alignment/libalignment.a: alignment/CMakeFiles/alignment.dir/build.make
alignment/libalignment.a: alignment/CMakeFiles/alignment.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/mlipovac/Matea/zavrsniRad/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libalignment.a"
	cd /Users/mlipovac/Matea/zavrsniRad/alignment && $(CMAKE_COMMAND) -P CMakeFiles/alignment.dir/cmake_clean_target.cmake
	cd /Users/mlipovac/Matea/zavrsniRad/alignment && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/alignment.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
alignment/CMakeFiles/alignment.dir/build: alignment/libalignment.a
.PHONY : alignment/CMakeFiles/alignment.dir/build

alignment/CMakeFiles/alignment.dir/clean:
	cd /Users/mlipovac/Matea/zavrsniRad/alignment && $(CMAKE_COMMAND) -P CMakeFiles/alignment.dir/cmake_clean.cmake
.PHONY : alignment/CMakeFiles/alignment.dir/clean

alignment/CMakeFiles/alignment.dir/depend:
	cd /Users/mlipovac/Matea/zavrsniRad && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mlipovac/Matea/zavrsniRad /Users/mlipovac/Matea/zavrsniRad/alignment /Users/mlipovac/Matea/zavrsniRad /Users/mlipovac/Matea/zavrsniRad/alignment /Users/mlipovac/Matea/zavrsniRad/alignment/CMakeFiles/alignment.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : alignment/CMakeFiles/alignment.dir/depend

