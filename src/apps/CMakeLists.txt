SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCMAKE_BUILD_TYPE=Release")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPERF_PROF")

add_app_exec(deepwalk)
add_app_exec(ppr)
add_app_exec(metapath)
add_app_exec(node2vec)
