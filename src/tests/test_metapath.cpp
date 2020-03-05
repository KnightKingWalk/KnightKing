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
#include <utility>
#include <queue>
#include <type_traits>

#include <gtest/gtest.h>

#include "storage.hpp"
#include "graph.hpp"
#include "walk.hpp"
#include "util.hpp"
#include "test.hpp"
#include "test_walk.hpp"
#include "../apps/metapath.hpp"

const int edge_type_num = 5;

template<typename edge_data_t>
void check_metapath_random_walk(vertex_id_t v_num, Edge<edge_data_t> *edges, edge_id_t e_num, walker_id_t walker_num, step_t walk_lenght, std::vector<std::vector<std::vector<bool> > > schemes, std::vector<std::vector<vertex_id_t> > &seq, std::vector<Walker<MetapathState> > &walker_init_state)
{
    size_t max_state_num = 0;
    for (auto &sch : schemes)
    {
        max_state_num += sch.size();
    }
    auto get_state_id = [&] (scheme_id_t scheme_id, meta_state_t meta_state)
    {
        assert(scheme_id < schemes.size());
        size_t state_id = 0;
        for (scheme_id_t s_i = 0; s_i < scheme_id; s_i++)
        {
            state_id += schemes[s_i].size();
        }
        return state_id + meta_state;
    };
    std::vector<std::vector<Edge<edge_data_t> > > graph(v_num);
    for (edge_id_t e_i = 0; e_i < e_num; e_i++)
    {
        graph[edges[e_i].src].push_back(edges[e_i]);
    }
    for (auto &adj : graph)
    {
        std::sort(adj.begin(), adj.end(), [](const Edge<edge_data_t> a, const Edge<edge_data_t> b){return a.dst < b.dst;});
    }
    auto get_edge_idx = [&] (vertex_id_t src, vertex_id_t dst)
    {
        Edge<edge_data_t> target;
        target.dst = dst;
        auto p = std::lower_bound(graph[src].begin(), graph[src].end(), target, [](const Edge<edge_data_t> &a, const Edge<edge_data_t> &b) { return a.dst < b.dst; });
        assert(p != graph[src].end());
        return p - graph[src].begin();
    };
    std::vector<std::vector<bool> > vis(v_num);
    for (auto &x : vis)
    {
        x.resize(max_state_num, false);
    }
    struct QueueItem
    {
        scheme_id_t sch_id;
        meta_state_t state;
        vertex_id_t vertex;
    };
    for (walker_id_t v_i = 0; v_i < v_num; v_i ++)
    {
        for (scheme_id_t sch_i = 0; sch_i < schemes.size(); sch_i++)
        {
            std::queue<QueueItem> q;
            auto expand_func = [&] (QueueItem current)
            {
                size_t state_id = get_state_id(current.sch_id, current.state);
                if (!vis[current.vertex][state_id])
                {
                    vis[current.vertex][state_id] = true;
                    for (auto edge : graph[current.vertex])
                    {
                        if (!schemes[current.sch_id][current.state][edge.data.get_meta()])
                        {
                            continue;
                        }
                        QueueItem next;
                        next.sch_id = current.sch_id;
                        next.state = (current.state + 1) % schemes[current.sch_id].size();
                        next.vertex = edge.dst;
                        q.push(next);
                    }
                }

            };
            QueueItem start;
            start.vertex = v_i;
            start.sch_id = sch_i;
            start.state = 0;
            expand_func(start);
            while (q.empty() == false)
            {
                QueueItem current = q.front();
                q.pop();
                expand_func(current);
            }
        }
    }
    std::vector<std::vector<std::vector<double> > > std_trans_mat(v_num);
    for (vertex_id_t v_i = 0; v_i < v_num; v_i++)
    {
        std_trans_mat[v_i].resize(max_state_num);
        for (scheme_id_t sch_i = 0; sch_i < schemes.size(); sch_i++)
        {
            for (meta_state_t ms_i = 0; ms_i < schemes[sch_i].size(); ms_i ++)
            {
                size_t state_id = get_state_id(sch_i, ms_i);
                auto &dist = std_trans_mat[v_i][state_id];
                dist.resize(graph[v_i].size(), 0.0);
                if (!vis[v_i][state_id])
                {
                    continue;
                }
                for (size_t e_i = 0; e_i < graph[v_i].size(); e_i++)
                {
                    auto &edge = graph[v_i][e_i];
                    double val = schemes[sch_i][ms_i][edge.data.get_meta()]? 1.0 : 0.0;
                    val *= edge.data.get_weight();
                    dist[e_i] = val;
                }
            }
        }
    }
    std::vector<std::vector<std::vector<double> > > real_trans_mat(v_num);
    for (vertex_id_t v_i = 0; v_i < v_num; v_i++)
    {
        real_trans_mat[v_i].resize(max_state_num);
        for (size_t s_i = 0; s_i < max_state_num; s_i++)
        {
            real_trans_mat[v_i][s_i].resize(graph[v_i].size(), 0.0);
        }
    }
    for (walker_id_t w_i = 0; w_i < seq.size(); w_i++)
    {
        meta_state_t init_meta = walker_init_state[w_i].data.state;
        meta_state_t scheme_id = walker_init_state[w_i].data.scheme_id;
        for (step_t s_i = 0; s_i + 1 < seq[w_i].size(); s_i++)
        {
            vertex_id_t current_v = seq[w_i][s_i];
            meta_state_t ms = (init_meta + s_i) % schemes[scheme_id].size();
            size_t state_id = get_state_id(scheme_id, ms);
            size_t edge_idx = get_edge_idx(seq[w_i][s_i], seq[w_i][s_i + 1]);
            real_trans_mat[current_v][state_id][edge_idx] += 1.0;
            ASSERT_TRUE(schemes[scheme_id][ms][graph[current_v][edge_idx].data.get_meta()]);
        }
    }
    auto get_flat_mat = [] (std::vector<std::vector<std::vector<double> > > &three_d_mat)
    {
        std::vector<std::vector<double> > two_d_mat;
        for (auto &x : three_d_mat)
        {
            for (auto &y : x)
            {
                two_d_mat.push_back(y);
            }
        }
        return two_d_mat;
    };
    auto std_flat_mat = get_flat_mat(std_trans_mat);
    auto real_flat_mat = get_flat_mat(real_trans_mat);
    mat_normalization(std_flat_mat);
    mat_normalization(real_flat_mat);
    cmp_trans_matrix(real_flat_mat, std_flat_mat);
}

template<typename edge_data_t>
void test_metapath(vertex_id_t v_num, int worker_number)
{
    WalkEngine<edge_data_t, MetapathState> graph;
    graph.set_concurrency(worker_number);
    graph.load_graph(v_num, test_data_file);

    step_t walk_length = 80 + rand() % 20;
    walker_id_t walker_num = graph.get_vertex_num() * 500 + graph.get_edge_num() * 500 + rand() % 100;
    MPI_Bcast(&walk_length, sizeof(walk_length), get_mpi_data_type<char>(), 0, MPI_COMM_WORLD);
    MPI_Bcast(&walker_num, sizeof(walker_num), get_mpi_data_type<char>(), 0, MPI_COMM_WORLD);
    std::vector<std::vector<std::vector<bool> > > schemes;
    if (get_mpi_rank() == 0)
    {
        int scheme_num = 3 + rand() % 2;
        for (int s_i = 0; s_i < scheme_num; s_i++)
        {
            int scheme_length = 3 + rand() % 2;
            std::vector<std::vector<bool> > sch;
            for (int l_i = 0; l_i < scheme_length; l_i++)
            {
                std::vector<bool> val;
                for (int v_i = 0; v_i < edge_type_num; v_i++)
                {
                    val.push_back((bool)(rand() % 2));
                }
                sch.push_back(val);
            }
            schemes.push_back(sch);
        }
    }
    size_t scheme_num = schemes.size();
    MPI_Bcast(&scheme_num, sizeof(scheme_num), get_mpi_data_type<char>(), 0, MPI_COMM_WORLD);
    schemes.resize(scheme_num);
    for (int s_i = 0; s_i < scheme_num; s_i++)
    {
        size_t scheme_length = schemes[s_i].size();
        MPI_Bcast(&scheme_length, sizeof(scheme_length), get_mpi_data_type<char>(), 0, MPI_COMM_WORLD);
        schemes[s_i].resize(scheme_length);
        for (int l_i = 0; l_i < scheme_length; l_i++)
        {
            schemes[s_i][l_i].resize(edge_type_num);
            for (int v_i = 0; v_i < edge_type_num; v_i++)
            {
                bool val = schemes[s_i][l_i][v_i];
                MPI_Bcast(&val, sizeof(val), get_mpi_data_type<char>(), 0, MPI_COMM_WORLD);
                schemes[s_i][l_i][v_i] = val;
            }
        }
    }

    metapath(&graph, schemes, walker_num, walk_length);

    std::vector<std::vector<vertex_id_t> > rw_sequences;
    std::vector<Walker<MetapathState> > walker_init_state;
    graph.collect_walk_sequence(rw_sequences, walker_num);
    graph.collect_walker_init_state(walker_init_state);

    if (get_mpi_rank() == 0)
    {
        Edge<edge_data_t> *std_edges;
        edge_id_t std_edge_num;
        read_graph(test_data_file, 0, 1, std_edges, std_edge_num);
        check_metapath_random_walk(v_num, std_edges, std_edge_num, walker_num, walk_length, schemes, rw_sequences, walker_init_state);
    }
}

template<typename T>
std::function<void(T&)> get_edge_data_gen_func()
{
    printf("[error] Undefined edge data generator\n");
    exit(1);
}

template<>
std::function<void(UnweightedMetaData&)> get_edge_data_gen_func<UnweightedMetaData>()
{
    auto func = [] (UnweightedMetaData &data)
    {
        data.meta_info = rand() % edge_type_num;
    };
    return func;
}

template<>
std::function<void(WeightedMetaData&)> get_edge_data_gen_func<WeightedMetaData>()
{
    auto func = [] (WeightedMetaData &data)
    {
        data.meta_info = rand() % edge_type_num;
        gen_rand_edge_data<real_t>(data.weight);
    };
    return func;
}

template<typename edge_data_t>
void test_metapath()
{
    edge_id_t e_nums_arr[] = {100, 200, 300, 400, 500, 600};
    vertex_id_t v_num = 50 + rand() % 20;
    std::vector<edge_id_t> e_nums(e_nums_arr, e_nums_arr + 6);
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
            gen_undirected_graph_file<edge_data_t>(v_num, e_num, get_edge_data_gen_func<edge_data_t>());
        }
        MPI_Barrier(MPI_COMM_WORLD);
        int worker_number = rand() % 8 + 1;
        MPI_Bcast(&worker_number, 1, get_mpi_data_type<int>(), 0, MPI_COMM_WORLD);
        test_metapath<edge_data_t>(v_num, worker_number);
    }
    if (get_mpi_rank() == 0)
    {
        rm_test_graph_temp_file();
    }
}

TEST(Metapath, Unbiased)
{
    test_metapath<UnweightedMetaData>();
}

TEST(Outlier, Biased)
{
    test_metapath<WeightedMetaData>();
}

GTEST_API_ int main(int argc, char *argv[])
{
    MPI_Instance mpi_instance(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    mute_nonroot_gtest_events();
    int result = RUN_ALL_TESTS();
    return result;
}
