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

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <fstream>
#include <vector>
#include <utility>
#include <map>
#include <set>
#include <type_traits>

#include <gtest/gtest.h>

#include "storage.hpp"
#include "generator_helper.hpp"
#include "mpi_helper.hpp"
#include "graph.hpp"

//temporary file for holding randomly generated graph
//which is used for unit test
const char *test_data_file = "746123_embedding_test_temp_data";

template<typename edge_data_t>
void gen_undirected_graph_file(vertex_id_t v_num, edge_id_t e_num, std::vector<Edge<edge_data_t> > &glb_edges, std::function<void(edge_data_t&)> edge_data_gen_func = nullptr, GraphFormat gf = GF_Binary)
{
    //bi-directional edge
    assert(e_num % 2 == 0);
    std::set< std::pair<vertex_id_t, vertex_id_t> > filter;
    glb_edges.clear();
    for (edge_id_t i = 0; i < e_num / 2; i++)
    {
        bool ok = false;
        vertex_id_t s;
        vertex_id_t t;
        while (!ok)
        {
            s = rand() % v_num;
            t = rand() % v_num;
            if (s == t) continue;
            if (s > t) std::swap(s, t);
            if (filter.find(std::make_pair(s, t)) == filter.end())
            {
                filter.insert(std::make_pair(s, t));
                ok = true;
            }
        }
        Edge<edge_data_t> e;
        e.src = s;
        e.dst = t;
        if (edge_data_gen_func != nullptr)
        {
            edge_data_gen_func(e.data);
        } else
        {
            gen_rand_edge_data<edge_data_t>(e.data);
        }
        glb_edges.push_back(e);
        std::swap(e.src, e.dst);
        glb_edges.push_back(e);
    }
    assert(e_num == glb_edges.size());
    if (gf == GF_Binary)
    {
        write_graph(test_data_file, glb_edges.data(), e_num);
    } else if (gf == GF_Edgelist)
    {
        write_edgelist(test_data_file, glb_edges.data(), e_num);
    } else
    {
        fprintf(stderr, "Unsupported graph format");
        exit(1);
    }
}

template<typename edge_data_t>
void gen_undirected_graph_file(vertex_id_t v_num, edge_id_t e_num, std::function<void(edge_data_t&)> edge_data_gen_func = nullptr, GraphFormat gf = GF_Binary)
{
    std::vector<Edge<edge_data_t> > es;
    gen_undirected_graph_file(v_num, e_num, es, edge_data_gen_func, gf);
}

template<typename edge_data_t>
void gen_directed_graph_file(vertex_id_t v_num, edge_id_t e_num, std::vector<Edge<edge_data_t> > &glb_edges, std::function<void(edge_data_t&)> edge_data_gen_func = nullptr, GraphFormat gf = GF_Binary)
{
    std::set< std::pair<vertex_id_t, vertex_id_t> > filter;
    glb_edges.clear();
    for (edge_id_t i = 0; i < e_num; i++)
    {
        bool ok = false;
        vertex_id_t s;
        vertex_id_t t;
        while (!ok)
        {
            s = rand() % v_num;
            t = rand() % v_num;
            if (s == t) continue;
            if (filter.find(std::make_pair(s, t)) == filter.end())
            {
                filter.insert(std::make_pair(s, t));
                ok = true;
            }
        }
        Edge<edge_data_t> e;
        e.src = s;
        e.dst = t;
        if (edge_data_gen_func != nullptr)
        {
            edge_data_gen_func(e.data);
        } else
        {
            gen_rand_edge_data<edge_data_t>(e.data);
        }
        glb_edges.push_back(e);
    }
    assert(e_num == glb_edges.size());
    if (gf == GF_Binary)
    {
        write_graph(test_data_file, glb_edges.data(), e_num);
    } else if (gf == GF_Edgelist)
    {
        write_edgelist(test_data_file, glb_edges.data(), e_num);
    } else
    {
        fprintf(stderr, "Unsupported graph format\n");
        exit(1);
    }
}

template<typename edge_data_t>
void gen_directed_graph_file(vertex_id_t v_num, edge_id_t e_num, std::function<void(edge_data_t&)> edge_data_gen_func = nullptr, GraphFormat gf = GF_Binary)
{
    std::vector<Edge<edge_data_t> > es;
    gen_directed_graph_file(v_num, e_num, es, edge_data_gen_func, gf);
}

void rm_test_graph_temp_file()
{
    std::remove(test_data_file);
}

template <typename edge_data_t>
void edge_data_expect_equal(typename std::enable_if<std::is_same<EmptyData, edge_data_t>::value, edge_data_t>::type &a, edge_data_t &b)
{
}

template <typename edge_data_t>
void edge_data_expect_equal(typename std::enable_if<std::is_integral<edge_data_t>::value, edge_data_t>::type &a, edge_data_t &b)
{
    EXPECT_EQ(a, b);
}

template <typename edge_data_t>
void edge_data_expect_equal(typename std::enable_if<std::is_same<float, edge_data_t>::value, edge_data_t>::type &a, edge_data_t &b)
{
    EXPECT_FLOAT_EQ(a, b);
}

template <typename edge_data_t>
struct CmpEdgeLess
{
private:
    template<typename T>
    bool cmp_edge(const typename std::enable_if<std::is_same<EmptyData, T>::value, Edge<T> >::type &a, const Edge<T> &b) const
    {
        if (a.src != b.src)
        {
            return a.src < b.src;
        } else
        {
            return a.dst < b.dst;
        }
    }
    template<typename T>
    bool cmp_edge(const typename std::enable_if<std::is_arithmetic<T>::value, Edge<T> >::type &a, const Edge<T> &b) const
    {
        if (a.src != b.src)
        {
            return a.src < b.src;
        } else if (a.dst != b.dst)
        {
            return a.dst < b.dst;
        } else
        {
            return a.data < b.data;
        }
    }
public:
    bool operator()(const Edge<edge_data_t> &a, const Edge<edge_data_t> &b) const
    {
        return cmp_edge<edge_data_t>(a, b);
    }
};

template<typename edge_data_t>
void cmp_edges(Edge<edge_data_t> *edges, edge_id_t e_num, Edge<edge_data_t> *std_edges, edge_id_t std_e_num)
{
    CmpEdgeLess<edge_data_t> cmp_edge_less;
    ASSERT_EQ(e_num, std_e_num);
    std::sort(edges, edges + e_num, cmp_edge_less);
    std::sort(std_edges, std_edges + std_e_num, cmp_edge_less);
    for (edge_id_t i = 0; i < std_e_num; i++)
    {
//        printf("%u->%u %u->%u\n", edges[i].src, edges[i].dst, std_edges[i].src, std_edges[i].dst);
        EXPECT_EQ(edges[i].src, std_edges[i].src);
        EXPECT_EQ(edges[i].dst, std_edges[i].dst);
        edge_data_expect_equal(edges[i].data, std_edges[i].data);
    }
}

void mute_nonroot_gtest_events()
{
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (get_mpi_rank() != 0)
    {
        delete listeners.Release(listeners.default_result_printer());
    }
}
