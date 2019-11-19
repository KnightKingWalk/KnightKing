function(add_app_exec EXEC_NAME)
    add_executable(${EXEC_NAME} ${EXEC_NAME}.cpp)
    target_link_libraries(${EXEC_NAME} PUBLIC ${MPI_LIBRARIES})
endfunction(add_app_exec)

function(add_test_exec EXEC_NAME)
    add_executable(${EXEC_NAME} ${EXEC_NAME}.cpp)
    target_link_libraries(${EXEC_NAME} PUBLIC ${GTEST_LIBRARIES} ${MPI_LIBRARIES})
endfunction(add_test_exec)

function(add_tool_exec EXEC_NAME)
    add_executable(${EXEC_NAME} ${EXEC_NAME}.cpp)
    target_link_libraries(${EXEC_NAME} PUBLIC ${MPI_LIBRARIES})
endfunction(add_tool_exec)

foreach(prog "test_storage")
    add_test("${prog}" "${KTK_RUNTIME_OUTPUT_DIRECTORY}/${prog}")
endforeach(prog)

foreach(prog "test_graph" "test_path" "test_walker" "test_bound" "test_outlier" "test_deepwalk" "test_ppr" "test_metapath" "test_node2vec")
    add_test("${prog}" "./bin/${prog}")
    add_test("distributed_${prog}" "mpirun" "-n" "2" "${KTK_RUNTIME_OUTPUT_DIRECTORY}/${prog}")
endforeach(prog)
