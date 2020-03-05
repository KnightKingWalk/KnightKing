/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2019 Ke Yang, Tsinghua University 
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <vector>

#include <gtest/gtest.h>

#include "storage.hpp"
#include "graph.hpp"
#include "walk.hpp"
#include "path.hpp"
#include "util.hpp"
#include "test.hpp"
#include "../apps/ppr.hpp"
#include "../apps/node2vec.hpp"

void check_path_data(PathSet *ps, std::vector<std::vector<vertex_id_t> > &std_ps)
{
    for (int s_i = 0; s_i < ps->seg_num; s_i++)
    {
        for (walker_id_t p_i = 0; p_i < ps->path_num[s_i]; p_i++)
        {
            walker_id_t walker = ps->walker_id[s_i][p_i];
            ASSERT_TRUE(0 <= walker && walker < std_ps.size());
            ASSERT_EQ(ps->path_length[s_i][p_i], std_ps[walker].size());
            ASSERT_EQ(ps->path_length[s_i][p_i], ps->path_end[s_i][p_i] - ps->path_begin[s_i][p_i]);
            step_t step = 0;
            for (auto p = ps->path_begin[s_i][p_i]; p != ps->path_end[s_i][p_i]; p++)
            {
                EXPECT_EQ(*p, std_ps[walker][step]);
                step++;
            }
        }
    }

    walker_id_t local_path_num = 0;
    for (int s_i = 0; s_i < ps->seg_num; s_i++)
    {
        local_path_num += ps->path_num[s_i];
    }
    walker_id_t tot_path_num;
    MPI_Allreduce(&local_path_num, &tot_path_num, 1, get_mpi_data_type<walker_id_t>(), MPI_SUM, MPI_COMM_WORLD);
    ASSERT_EQ(tot_path_num, std_ps.size());

    std::thread send_thread([&]() {
        for (int s_i = 0; s_i < ps->seg_num; s_i++)
        {
            MPI_Send(ps->walker_id[s_i], ps->path_num[s_i], get_mpi_data_type<walker_id_t>(), 0, 0, MPI_COMM_WORLD);
        }
    });
    if (get_mpi_rank() == 0)
    {
        std::vector<bool> vis(tot_path_num, false);
        for (partition_id_t p_i = 0; p_i < get_mpi_size(); p_i++)
        {
            for (int s_i = 0; s_i < ps->seg_num; s_i++)
            {
                int recv_size = 0;
                MPI_Status recv_status;
                MPI_Probe(p_i, 0, MPI_COMM_WORLD, &recv_status);
                MPI_Get_count(&recv_status, get_mpi_data_type<walker_id_t>(), &recv_size);
                std::vector<walker_id_t> recv_data(recv_size);
                MPI_Recv(recv_data.data(), recv_size, get_mpi_data_type<walker_id_t>(), p_i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for (int w_i = 0; w_i < recv_size; w_i++)
                {
                    EXPECT_EQ(vis[recv_data[w_i]], false);
                    vis[recv_data[w_i]] = true;
                }
            }
        }

    }
    send_thread.join();
}

void broadcast_rw_sequences(std::vector<std::vector<vertex_id_t> > &seq)
{
    size_t seq_size = seq.size();
    MPI_Bcast(&seq_size, 1, get_mpi_data_type<size_t>(), 0, MPI_COMM_WORLD);
    if (get_mpi_rank() != 0)
    {
        seq.resize(seq_size);
    }
    for (size_t s_i = 0; s_i < seq_size; s_i++)
    {
        size_t seq_length = seq[s_i].size();
        MPI_Bcast(&seq_length, 1, get_mpi_data_type<size_t>(), 0, MPI_COMM_WORLD);
        if (get_mpi_rank() != 0)
        {
            seq[s_i].resize(seq_length);
        }
        MPI_Bcast(seq[s_i].data(), seq_length, get_mpi_data_type<vertex_id_t>(), 0, MPI_COMM_WORLD);
    }
}

void test_first_order_path(vertex_id_t v_num, int worker_number)
{
    std::vector<std::vector<vertex_id_t> > rw_sequences;
    WalkEngine<EmptyData, EmptyData> graph;
    graph.set_concurrency(worker_number);
    graph.load_graph(v_num, test_data_file);

    PathSet *ps = nullptr;
    WalkConfig walk_conf;
    walk_conf.set_output_consumer(
        [&] (PathSet* ps_param) {
            // Assume only has one iteration
            assert(ps == nullptr);
            ps = ps_param;
        }
    );

    real_t terminate_prob = 1.0 / 80;
    walker_id_t walker_num = graph.get_vertex_num() * 50 + graph.get_edge_num() * 10 + rand() % 100;
    MPI_Bcast(&walker_num, 1, get_mpi_data_type<walker_id_t>(), 0, MPI_COMM_WORLD);

    ppr(&graph, walker_num, terminate_prob, nullptr, &walk_conf);
    graph.collect_walk_sequence(rw_sequences, walker_num);
    broadcast_rw_sequences(rw_sequences);
    check_path_data(ps, rw_sequences);
}

void test_second_order_path(vertex_id_t v_num, int worker_number)
{
    std::vector<std::vector<vertex_id_t> > rw_sequences;
    WalkEngine<EmptyData, Node2vecState> graph;
    graph.set_concurrency(worker_number);
    graph.load_graph(v_num, test_data_file);

    PathSet *ps;
    WalkConfig walk_conf;
    walk_conf.set_output_consumer(
        [&] (PathSet* ps_param) {
            ps = ps_param;
        }
    );

    Node2vecConf n2v_conf;
    n2v_conf.walker_num = graph.get_vertex_num() * 50 + graph.get_edge_num() * 10 + rand() % 100;
    n2v_conf.walk_length = rand() % 20 + 60;
    n2v_conf.p = rand() % 4 + 1;
    n2v_conf.q = rand() % 4 + 1;
    if (rand() % 2 == 0)
    {
        n2v_conf.p = 1.0 / n2v_conf.p;
    }
    if (rand() % 2 == 0)
    {
        n2v_conf.q = 1.0 / n2v_conf.q;
    }
    MPI_Bcast(&n2v_conf, sizeof(n2v_conf), get_mpi_data_type<char>(), 0, MPI_COMM_WORLD);

    node2vec(&graph, n2v_conf, &walk_conf);
    graph.collect_walk_sequence(rw_sequences, n2v_conf.walker_num);
    broadcast_rw_sequences(rw_sequences);
    check_path_data(ps, rw_sequences);
    delete ps;
}

void test_path(int order)
{
    edge_id_t e_nums_arr[] = {1000, 2000, 4000, 5556, 8888, 10000, 20000, 100000};
    vertex_id_t v_num = 1000 + rand() % 500;
    std::vector<edge_id_t> e_nums(e_nums_arr, e_nums_arr + 8);
    /*
    size_t e_nums_arr[] = {30};
    vertex_id_t v_num = 10;
    std::vector<size_t> e_nums(e_nums_arr, e_nums_arr + 1);
    */

    MPI_Bcast(&v_num, 1, get_mpi_data_type<vertex_id_t>(), 0, MPI_COMM_WORLD);

    for (auto &e_num : e_nums_arr)
    {
        if (get_mpi_rank() == 0)
        {
            gen_undirected_graph_file<EmptyData>(v_num, e_num);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        int worker_number = rand() % 8 + 1;
        MPI_Bcast(&worker_number, 1, get_mpi_data_type<int>(), 0, MPI_COMM_WORLD);
        if (order == 1)
        {
            test_first_order_path(v_num, worker_number);
        } else if (order == 2)
        {
            test_second_order_path(v_num, worker_number);
        } else {
            exit(1);
        }
    }
    if (get_mpi_rank() == 0)
    {
        rm_test_graph_temp_file();
    }
}

TEST(PathOutput, FirstOrderRandomWalk)
{
    test_path(1);
}

TEST(PathOutput, SecondOrderRandomWalk)
{
    test_path(2);
}

GTEST_API_ int main(int argc, char *argv[])
{
    MPI_Instance mpi_instance(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    mute_nonroot_gtest_events();
    int result = RUN_ALL_TESTS();
    return result;
}
