CMAKE_MINIMUM_REQUIRED(VERSION 2.8...3.13)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
CMAKE_POLICY(SET CMP0003 NEW)
CMAKE_POLICY(SET CMP0015 NEW)

PROJECT(spatial-benchmark)

if (${CMAKE_SYSTEM_NAME} MATCHES Windows)
    set (SB_FREEGLUT_DIR "${CMAKE_SOURCE_DIR}/freeglut")
    include_directories (${SB_FREEGLUT_DIR}/include)
    link_directories (${SB_FREEGLUT_DIR}/lib)
endif ()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set (sb_sources
    BruteForce.cpp
    BruteForce_omp.cpp
    HierarchicalGrid.cpp
    HierarchicalGrid_omp.cpp
    Kdtree.cpp
    Kdtree_omp.cpp
    LooseOctree.cpp
    LooseOctree_omp.cpp
    Octree.cpp
    Octree_omp.cpp
    SortAndSweep.cpp
    SortAndSweep_omp.cpp
    SphereObject.cpp
    UniformGrid.cpp
    UniformGrid_omp.cpp
    Core.cpp
)

set (sb_headers
    Base.h
    BruteForce.h
    BruteForce_omp.h
    HierarchicalGrid.h
    HierarchicalGrid_omp.h
    ISpatialCell.h
    ISpatialObject.h
    ISpatialStructure.h
    Kdtree.h
    Kdtree_omp.h
    LooseOctree.h
    LooseOctree_omp.h
    Octree.h
    Octree_omp.h
    SortAndSweep.h
    SortAndSweep_omp.h
    SphereObject.h
    UniformGrid.h
    UniformGrid_omp.h
    Vector3.h
)

add_executable (spatial-benchmark ${sb_sources} ${sb_headers})

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(spatial-benchmark Threads::Threads)

if (${CMAKE_SYSTEM_NAME} MATCHES Darwin)
    target_link_libraries (spatial-benchmark "-framework OpenGL -framework GLUT -framework ApplicationServices")
elseif (${CMAKE_SYSTEM_NAME} MATCHES Linux)
    target_link_libraries (spatial-benchmark glut GLU GL)
elseif (${CMAKE_SYSTEM_NAME} MATCHES Windows)
    target_link_libraries (spatial-benchmark freeglut opengl32)
endif()
