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

void print_mat(std::vector<std::vector<double> > &mat)
{
    for (auto &vec : mat)
    {
        for (auto v : vec)
        {
            printf("%lf ", v);
        }
        printf("\n");
    }
    printf("\n");
}

void mat_normalization(std::vector<std::vector<double> > &mat)
{
    for (auto &vec : mat)
    {
        double sum = 0;
        for (auto val : vec)
        {
            sum += val;
        }
        if (sum != 0)
        {
            for (auto &val : vec)
            {
                val = val / sum;
            }
        }
    }
}

void cmp_trans_matrix(std::vector<std::vector<double> > &a, std::vector<std::vector<double> > &b, double variance_upper_bound = 1.0)
{
    ASSERT_EQ(a.size(), b.size());
    double max_row_variance = 0;
    double variance = 0;
    for (size_t r_i = 0; r_i < a.size(); r_i++)
    {
        double row_variance = 0;
        ASSERT_EQ(a[r_i].size(), b[r_i].size());
        bool err_in_row = false;
        for (size_t c_i = 0; c_i < a[r_i].size(); c_i++)
        {
            if (!((a[r_i][c_i] == 0 && b[r_i][c_i] == 0) || (a[r_i][c_i] != 0 && b[r_i][c_i] != 0)))
            {
                printf("(%lf %lf)\n", a[r_i][c_i], b[r_i][c_i]);
                err_in_row = true;
            }
            EXPECT_TRUE((a[r_i][c_i] == 0 && b[r_i][c_i] == 0) || (a[r_i][c_i] != 0 && b[r_i][c_i] != 0));
            double diff = a[r_i][c_i] - b[r_i][c_i];
            row_variance += diff * diff;
        }
        if (err_in_row)
        {
            for (size_t c_i = 0; c_i < a[r_i].size(); c_i++)
            {
                if (a[r_i][c_i] != 0 || b[r_i][c_i] != 0)
                {
                    printf("%.3lf ", a[r_i][c_i]);
                }
            }
            printf("\n");
            for (size_t c_i = 0; c_i < a[r_i].size(); c_i++)
            {
                if (a[r_i][c_i] != 0 || b[r_i][c_i] != 0)
                {
                    printf("%.3lf ", b[r_i][c_i]);
                }
            }
            printf("\n");
        }
        max_row_variance = std::max(max_row_variance, row_variance);
        variance += row_variance;
    }
    printf("max_row_variance : %lf, variance : %lf\n", max_row_variance, variance);
    ASSERT_TRUE(variance < variance_upper_bound);
}

double get_edge_trans_weight(Edge<EmptyData> &e)
{
    return 1;
}

double get_edge_trans_weight(Edge<float> &e)
{
    return e.data;
}

template<typename edge_data_t>
void get_static_walk_trans_matrix(vertex_id_t v_num, Edge<edge_data_t> *edges, edge_id_t e_num, std::vector<std::vector<double> > &trans_mat)
{
    std::vector<double> weight_sum(v_num, 0.0);
    for (edge_id_t e_i = 0; e_i < e_num; e_i++)
    {
        weight_sum[edges[e_i].src] += get_edge_trans_weight(edges[e_i]);
    }
    for (vertex_id_t v_i = 0; v_i < v_num; v_i++)
    {
        if (weight_sum[v_i] == 0)
        {
            trans_mat[v_i][v_i] = 0;
        }
    }
    for (edge_id_t e_i = 0; e_i < e_num; e_i++)
    {
        auto &e = edges[e_i];
        trans_mat[e.src][e.dst] += get_edge_trans_weight(e) / weight_sum[e.src];
    }
}

template<typename edge_data_t>
void check_static_random_walk(vertex_id_t v_num, Edge<edge_data_t> *edges, edge_id_t e_num, std::vector<std::vector<vertex_id_t> > rw_sequences)
{
    std::vector<std::vector<double> > trans_mat(v_num);
    for (auto &vec : trans_mat)
    {
        vec.resize(v_num, 0.0);
    }
    get_static_walk_trans_matrix(v_num, edges, e_num, trans_mat);
    
    //check if sequences are legal
    for (auto &s : rw_sequences)
    {
        for (size_t v_i = 0; v_i + 1 < s.size(); v_i++)
        {
            ASSERT_TRUE(trans_mat[s[v_i]][s[v_i + 1]] != 0);
        }
    }

    std::vector<std::vector<double> > real_trans_mat(v_num);
    for (auto &vec : real_trans_mat)
    {
        vec.resize(v_num, 0.0);
    }
    for (auto &s : rw_sequences)
    {
        for (size_t v_i = 0; v_i + 1 < s.size(); v_i++)
        {
            real_trans_mat[s[v_i]][s[v_i + 1]] += 1;
        }
    }
    mat_normalization(real_trans_mat);
    //check if trans_mat is obeyed during random walk
    cmp_trans_matrix(real_trans_mat, trans_mat);
}
