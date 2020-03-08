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
#include "util.hpp"
#include "test.hpp"

template<typename edge_data_t>
void test_load_graph(partition_id_t partition_num, edge_id_t e_num, std::vector<Edge<edge_data_t> > &glb_edges, GraphFormat gf)
{
    std::vector<Edge<edge_data_t>*> loaded_edges(partition_num);
    std::vector<size_t> loaded_e_nums(partition_num);
    for (int i = 0; i < partition_num; i++)
    {
        if (gf == GF_Binary)
        {
            read_graph(test_data_file, i, partition_num, loaded_edges[i], loaded_e_nums[i]);
        } else if (gf == GF_Edgelist)
        {
            read_edgelist(test_data_file, i, partition_num, loaded_edges[i], loaded_e_nums[i]);
        } else
        {
            fprintf(stderr, "Unsupported graph format\n");
            exit(1);
        }
    }
    size_t max_e_num = 0;
    size_t min_e_num = glb_edges.size();
    for (auto e_num : loaded_e_nums)
    {
        max_e_num = std::max(max_e_num, e_num);
        min_e_num = std::min(min_e_num, e_num);
    }
    //for load balance
    EXPECT_TRUE(max_e_num - min_e_num <= 10);

    size_t tot_loaded_e_num = 0;
    for (auto &en : loaded_e_nums)
    {
        tot_loaded_e_num += en;
    }
    EXPECT_EQ(tot_loaded_e_num, e_num);
    size_t p = 0;
    for (int i = 0; i < partition_num; i++)
    {
        for(size_t j = 0; j < loaded_e_nums[i]; j++)
        {
            auto &le = loaded_edges[i][j];
            auto &e = glb_edges[p++];
            EXPECT_EQ(e.src, le.src);
            EXPECT_EQ(e.dst, le.dst);
            edge_data_expect_equal(e.data, le.data);
        }
    }
    for (int i = 0; i < partition_num; i++)
    {
        delete []loaded_edges[i];
    }
}

void test(GraphFormat gf)
{
    //edge_id_t e_nums_arr[] = {0, 2, 4, 6, 10, 16, 6536, 65538};
    edge_id_t e_nums_arr[] = {0, 2, 4, 6, 10, 16, 64, 644};
    std::vector<edge_id_t> e_nums(e_nums_arr, e_nums_arr + 8);
    for (auto &e_num : e_nums_arr)
    {
        for (int partition_num = 1; partition_num <= 8; partition_num ++)
        {
            std::vector<Edge<EmptyData> > glb_edges;
            gen_undirected_graph_file<EmptyData>(500 + rand() % 10000, e_num, glb_edges, nullptr, gf);
            test_load_graph(partition_num, e_num, glb_edges, gf);
            std::vector<Edge<real_t> > glb_weighted_edges;
            gen_undirected_graph_file<real_t>(500 + rand() % 10000, e_num, glb_weighted_edges, nullptr, gf);
            test_load_graph(partition_num, e_num, glb_weighted_edges, gf);
        }
    }
    rm_test_graph_temp_file();
}

TEST(Storage, ReadWriteBinaryGraph)
{
    test(GF_Binary);
}

TEST(Storage, ReadWriteEdgeList)
{
    test(GF_Edgelist);
}

GTEST_API_ int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    return result;
}

