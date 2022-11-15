# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/nems/parallel-physical-partical-stimulation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/nems/parallel-physical-partical-stimulation

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/nems/parallel-physical-partical-stimulation/CMakeFiles /home/nems/parallel-physical-partical-stimulation/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/nems/parallel-physical-partical-stimulation/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named spatial-benchmark

# Build rule for target.
spatial-benchmark: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 spatial-benchmark
.PHONY : spatial-benchmark

# fast build rule for target.
spatial-benchmark/fast:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/build
.PHONY : spatial-benchmark/fast

BruteForce.o: BruteForce.cpp.o

.PHONY : BruteForce.o

# target to build an object file
BruteForce.cpp.o:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/BruteForce.cpp.o
.PHONY : BruteForce.cpp.o

BruteForce.i: BruteForce.cpp.i

.PHONY : BruteForce.i

# target to preprocess a source file
BruteForce.cpp.i:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/BruteForce.cpp.i
.PHONY : BruteForce.cpp.i

BruteForce.s: BruteForce.cpp.s

.PHONY : BruteForce.s

# target to generate assembly for a file
BruteForce.cpp.s:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/BruteForce.cpp.s
.PHONY : BruteForce.cpp.s

Core.o: Core.cpp.o

.PHONY : Core.o

# target to build an object file
Core.cpp.o:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/Core.cpp.o
.PHONY : Core.cpp.o

Core.i: Core.cpp.i

.PHONY : Core.i

# target to preprocess a source file
Core.cpp.i:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/Core.cpp.i
.PHONY : Core.cpp.i

Core.s: Core.cpp.s

.PHONY : Core.s

# target to generate assembly for a file
Core.cpp.s:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/Core.cpp.s
.PHONY : Core.cpp.s

HierarchicalGrid.o: HierarchicalGrid.cpp.o

.PHONY : HierarchicalGrid.o

# target to build an object file
HierarchicalGrid.cpp.o:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/HierarchicalGrid.cpp.o
.PHONY : HierarchicalGrid.cpp.o

HierarchicalGrid.i: HierarchicalGrid.cpp.i

.PHONY : HierarchicalGrid.i

# target to preprocess a source file
HierarchicalGrid.cpp.i:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/HierarchicalGrid.cpp.i
.PHONY : HierarchicalGrid.cpp.i

HierarchicalGrid.s: HierarchicalGrid.cpp.s

.PHONY : HierarchicalGrid.s

# target to generate assembly for a file
HierarchicalGrid.cpp.s:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/HierarchicalGrid.cpp.s
.PHONY : HierarchicalGrid.cpp.s

Kdtree.o: Kdtree.cpp.o

.PHONY : Kdtree.o

# target to build an object file
Kdtree.cpp.o:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/Kdtree.cpp.o
.PHONY : Kdtree.cpp.o

Kdtree.i: Kdtree.cpp.i

.PHONY : Kdtree.i

# target to preprocess a source file
Kdtree.cpp.i:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/Kdtree.cpp.i
.PHONY : Kdtree.cpp.i

Kdtree.s: Kdtree.cpp.s

.PHONY : Kdtree.s

# target to generate assembly for a file
Kdtree.cpp.s:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/Kdtree.cpp.s
.PHONY : Kdtree.cpp.s

LooseOctree.o: LooseOctree.cpp.o

.PHONY : LooseOctree.o

# target to build an object file
LooseOctree.cpp.o:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/LooseOctree.cpp.o
.PHONY : LooseOctree.cpp.o

LooseOctree.i: LooseOctree.cpp.i

.PHONY : LooseOctree.i

# target to preprocess a source file
LooseOctree.cpp.i:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/LooseOctree.cpp.i
.PHONY : LooseOctree.cpp.i

LooseOctree.s: LooseOctree.cpp.s

.PHONY : LooseOctree.s

# target to generate assembly for a file
LooseOctree.cpp.s:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/LooseOctree.cpp.s
.PHONY : LooseOctree.cpp.s

Octree.o: Octree.cpp.o

.PHONY : Octree.o

# target to build an object file
Octree.cpp.o:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/Octree.cpp.o
.PHONY : Octree.cpp.o

Octree.i: Octree.cpp.i

.PHONY : Octree.i

# target to preprocess a source file
Octree.cpp.i:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/Octree.cpp.i
.PHONY : Octree.cpp.i

Octree.s: Octree.cpp.s

.PHONY : Octree.s

# target to generate assembly for a file
Octree.cpp.s:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/Octree.cpp.s
.PHONY : Octree.cpp.s

SortAndSweep.o: SortAndSweep.cpp.o

.PHONY : SortAndSweep.o

# target to build an object file
SortAndSweep.cpp.o:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/SortAndSweep.cpp.o
.PHONY : SortAndSweep.cpp.o

SortAndSweep.i: SortAndSweep.cpp.i

.PHONY : SortAndSweep.i

# target to preprocess a source file
SortAndSweep.cpp.i:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/SortAndSweep.cpp.i
.PHONY : SortAndSweep.cpp.i

SortAndSweep.s: SortAndSweep.cpp.s

.PHONY : SortAndSweep.s

# target to generate assembly for a file
SortAndSweep.cpp.s:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/SortAndSweep.cpp.s
.PHONY : SortAndSweep.cpp.s

SphereObject.o: SphereObject.cpp.o

.PHONY : SphereObject.o

# target to build an object file
SphereObject.cpp.o:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/SphereObject.cpp.o
.PHONY : SphereObject.cpp.o

SphereObject.i: SphereObject.cpp.i

.PHONY : SphereObject.i

# target to preprocess a source file
SphereObject.cpp.i:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/SphereObject.cpp.i
.PHONY : SphereObject.cpp.i

SphereObject.s: SphereObject.cpp.s

.PHONY : SphereObject.s

# target to generate assembly for a file
SphereObject.cpp.s:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/SphereObject.cpp.s
.PHONY : SphereObject.cpp.s

UniformGrid.o: UniformGrid.cpp.o

.PHONY : UniformGrid.o

# target to build an object file
UniformGrid.cpp.o:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/UniformGrid.cpp.o
.PHONY : UniformGrid.cpp.o

UniformGrid.i: UniformGrid.cpp.i

.PHONY : UniformGrid.i

# target to preprocess a source file
UniformGrid.cpp.i:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/UniformGrid.cpp.i
.PHONY : UniformGrid.cpp.i

UniformGrid.s: UniformGrid.cpp.s

.PHONY : UniformGrid.s

# target to generate assembly for a file
UniformGrid.cpp.s:
	$(MAKE) -f CMakeFiles/spatial-benchmark.dir/build.make CMakeFiles/spatial-benchmark.dir/UniformGrid.cpp.s
.PHONY : UniformGrid.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... spatial-benchmark"
	@echo "... BruteForce.o"
	@echo "... BruteForce.i"
	@echo "... BruteForce.s"
	@echo "... Core.o"
	@echo "... Core.i"
	@echo "... Core.s"
	@echo "... HierarchicalGrid.o"
	@echo "... HierarchicalGrid.i"
	@echo "... HierarchicalGrid.s"
	@echo "... Kdtree.o"
	@echo "... Kdtree.i"
	@echo "... Kdtree.s"
	@echo "... LooseOctree.o"
	@echo "... LooseOctree.i"
	@echo "... LooseOctree.s"
	@echo "... Octree.o"
	@echo "... Octree.i"
	@echo "... Octree.s"
	@echo "... SortAndSweep.o"
	@echo "... SortAndSweep.i"
	@echo "... SortAndSweep.s"
	@echo "... SphereObject.o"
	@echo "... SphereObject.i"
	@echo "... SphereObject.s"
	@echo "... UniformGrid.o"
	@echo "... UniformGrid.i"
	@echo "... UniformGrid.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

