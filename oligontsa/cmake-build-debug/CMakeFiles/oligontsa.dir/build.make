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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/oligontsa.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/oligontsa.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/oligontsa.dir/flags.make

CMakeFiles/oligontsa.dir/main.cpp.o: CMakeFiles/oligontsa.dir/flags.make
CMakeFiles/oligontsa.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/oligontsa.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/oligontsa.dir/main.cpp.o -c /Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa/main.cpp

CMakeFiles/oligontsa.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/oligontsa.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa/main.cpp > CMakeFiles/oligontsa.dir/main.cpp.i

CMakeFiles/oligontsa.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/oligontsa.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa/main.cpp -o CMakeFiles/oligontsa.dir/main.cpp.s

CMakeFiles/oligontsa.dir/osa.cpp.o: CMakeFiles/oligontsa.dir/flags.make
CMakeFiles/oligontsa.dir/osa.cpp.o: ../osa.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/oligontsa.dir/osa.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/oligontsa.dir/osa.cpp.o -c /Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa/osa.cpp

CMakeFiles/oligontsa.dir/osa.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/oligontsa.dir/osa.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa/osa.cpp > CMakeFiles/oligontsa.dir/osa.cpp.i

CMakeFiles/oligontsa.dir/osa.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/oligontsa.dir/osa.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa/osa.cpp -o CMakeFiles/oligontsa.dir/osa.cpp.s

# Object files for target oligontsa
oligontsa_OBJECTS = \
"CMakeFiles/oligontsa.dir/main.cpp.o" \
"CMakeFiles/oligontsa.dir/osa.cpp.o"

# External object files for target oligontsa
oligontsa_EXTERNAL_OBJECTS =

oligontsa: CMakeFiles/oligontsa.dir/main.cpp.o
oligontsa: CMakeFiles/oligontsa.dir/osa.cpp.o
oligontsa: CMakeFiles/oligontsa.dir/build.make
oligontsa: CMakeFiles/oligontsa.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable oligontsa"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/oligontsa.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/oligontsa.dir/build: oligontsa
.PHONY : CMakeFiles/oligontsa.dir/build

CMakeFiles/oligontsa.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/oligontsa.dir/cmake_clean.cmake
.PHONY : CMakeFiles/oligontsa.dir/clean

CMakeFiles/oligontsa.dir/depend:
	cd /Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa /Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa /Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa/cmake-build-debug /Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa/cmake-build-debug /Users/kmive/Desktop/GWU/Summer-2023/HSCI-6273/WK05/oligontsa/cmake-build-debug/CMakeFiles/oligontsa.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/oligontsa.dir/depend

