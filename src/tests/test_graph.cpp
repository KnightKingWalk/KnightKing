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
#include <map>
#include <set>
#include <type_traits>

#include <gtest/gtest.h>

#include "storage.hpp"
#include "graph.hpp"
#include "util.hpp"
#include "test.hpp"

template<typename edge_data_t>
class GraphTester: public GraphEngine<edge_data_t>
{
public:
    void get_edges(EdgeContainer<edge_data_t> *ec, std::vector<Edge<edge_data_t>> &ret)
    {
        ret.clear();
        for (vertex_id_t v_i = this->vertex_partition_begin[this->local_partition_id]; v_i < this->vertex_partition_end[this->local_partition_id]; v_i++)
        {
            for (auto p = ec->adj_lists[v_i].begin; p != ec->adj_lists[v_i].end; p++)
            {
                Edge<edge_data_t> e;
                e.src = v_i;
                e.dst = p->neighbour;
				if (!std::is_same<edge_data_t, EmptyData>::value)
				{
					e.data = p->data;
				}
                ret.push_back(e);
            }
        }
    }
    void set_concurrency(int worker_number)
    {
        this->set_graph_engine_concurrency(worker_number);
    }
};

template<typename edge_data_t>
void check_edges(GraphTester<edge_data_t>* graph, bool load_as_undirected = false)
{
    std::vector<Edge<edge_data_t> > local_graph_edges;
    graph->get_edges(graph->csr, local_graph_edges);
    if (get_mpi_rank() == 0)
    {
        Edge<edge_data_t> *std_edges;
        edge_id_t std_edge_num;
        read_graph(test_data_file, 0, 1, std_edges, std_edge_num);
        if (load_as_undirected)
        {
            std::vector<Edge<edge_data_t> > temp;
            for (edge_id_t e_i = 0; e_i < std_edge_num; e_i++)
            {
                temp.push_back(std_edges[e_i]);
                std::swap(std_edges[e_i].src, std_edges[e_i].dst);
                temp.push_back(std_edges[e_i]);
            }
            delete []std_edges;
            std_edge_num *= 2;
            std_edges = new Edge<edge_data_t>[std_edge_num];
            memcpy(std_edges, temp.data(), sizeof(Edge<edge_data_t>) * std_edge_num);
        }
        auto graph_edges = local_graph_edges;
        for (partition_id_t p_i = 1; p_i < get_mpi_size(); p_i++)
        {
            int recv_size = 0;
            MPI_Status recv_status;
            MPI_Probe(p_i, Tag_ShuffleGraph, MPI_COMM_WORLD, &recv_status);
            MPI_Get_count(&recv_status, get_mpi_data_type<char>(), &recv_size);
            std::vector<Edge<edge_data_t> > remote_edges(recv_size / sizeof(Edge<edge_data_t>));
            MPI_Recv(remote_edges.data(), recv_size, get_mpi_data_type<char>(), p_i, Tag_ShuffleGraph, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (auto e : remote_edges)
            {
                graph_edges.push_back(e);
            }
        }
        cmp_edges(graph_edges.data(), graph_edges.size(), std_edges, std_edge_num); 
        delete []std_edges;
    } else
    {
        MPI_Send(local_graph_edges.data(), local_graph_edges.size() * sizeof(Edge<edge_data_t>), get_mpi_data_type<char>(), 0, Tag_ShuffleGraph, MPI_COMM_WORLD);
    }
}

template<typename edge_data_t>
void test_static_edge(vertex_id_t v_num, bool load_as_undirected = false)
{
    GraphTester<edge_data_t> graph;
    int worker_number = rand() % 8 + 1;
    graph.set_concurrency(worker_number);
    graph.load_graph(v_num, test_data_file, load_as_undirected);
    check_edges(&graph, load_as_undirected);
}


template<typename edge_data_t>
void test_edges(bool load_as_undirected = false)
{
    edge_id_t e_nums_arr[] = {0, 2, 6, 16, 8888, 10000, 20000, 100000};
    vertex_id_t v_num = 1000 + rand() % 1000;
    std::vector<edge_id_t> e_nums(e_nums_arr, e_nums_arr + 8);
    /*
    size_t e_nums_arr[] = {20};
    vertex_id_t v_num = 20;
    std::vector<size_t> e_nums(e_nums_arr, e_nums_arr + 1);
    */

    for (auto &e_num : e_nums_arr)
    {
        if (get_mpi_rank() == 0)
        {
            if (load_as_undirected)
            {
                gen_directed_graph_file<edge_data_t>(v_num, e_num);
            } else
            {
                gen_undirected_graph_file<edge_data_t>(v_num, e_num);
            }
        }
        MPI_Bcast(&v_num, 1, get_mpi_data_type<vertex_id_t>(), 0, MPI_COMM_WORLD);
        test_static_edge<edge_data_t>(v_num, load_as_undirected);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (get_mpi_rank() == 0)
    {
        rm_test_graph_temp_file();
    }
}

TEST(GraphEngine, DefaultLoad)
{
    test_edges<EmptyData>();
    test_edges<real_t>();
}

TEST(GraphEngine, LoadAsUndirected)
{
    test_edges<EmptyData>(true);
    test_edges<real_t>(true);
}

GTEST_API_ int main(int argc, char *argv[])
{
    MPI_Instance mpi_instance(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    mute_nonroot_gtest_events();
    int result = RUN_ALL_TESTS();
    return result;
}

