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

#include "walk.hpp"
#include "static_comp.hpp"

struct Node2vecState
{
    vertex_id_t previous_vertex;
};

struct Node2vecConf
{
    real_t p;
    real_t q;
    walker_id_t walker_num;
    step_t walk_length;
};

template<typename edge_data_t>
void node2vec(WalkEngine<edge_data_t, Node2vecState> *graph, Node2vecConf conf)
{
    MPI_Barrier(MPI_COMM_WORLD);
    Timer timer;

    real_t p = conf.p;
    real_t q = conf.q;
    step_t walk_length = conf.walk_length;
    walker_id_t walker_num = conf.walker_num;

    vertex_id_t local_vertex_begin = graph->get_local_vertex_begin();
    vertex_id_t local_vertex_end = graph->get_local_vertex_end();

    graph->template process_vertices<vertex_id_t>([&](vertex_id_t v_i) {
        std::sort(graph->csr->adj_lists[v_i].begin, graph->csr->adj_lists[v_i].end, [](const AdjUnit<edge_data_t> a, const AdjUnit<edge_data_t> b){return a.neighbour < b.neighbour;});
        return 0;
    });
    real_t upperbound = std::max(1.0 / p, std::max(1.0, 1.0 / q));
    real_t lowerbound = std::min(1.0 / p, std::min(1.0, 1.0 / q));
    graph->set_walkers(
        walker_num,
        nullptr,
        [&] (Walker<Node2vecState> &walker, vertex_id_t current_v, AdjUnit<edge_data_t> *edge)
        {
            walker.data.previous_vertex = current_v;
        }
    );
    graph->template second_order_random_walk<vertex_id_t, bool> (
        [&] (Walker<Node2vecState> &walker, vertex_id_t current_v)
        {
            return walker.step >= walk_length ? 0.0 : 1.0;
        },
        get_trivial_static_comp(graph),
        [&] (Walker<Node2vecState> &walker, walker_id_t walker_idx, vertex_id_t current_v, AdjUnit<edge_data_t> *edge)
        {
            if (walker.step != 0)
            {
                stateQuery<vertex_id_t> query;
                query.src_v = current_v;
                query.walker_idx = walker_idx;
                query.data = edge->neighbour;
                graph->emit(walker.data.previous_vertex, query);
            }
        },
        [&] (vertex_id_t vtx, stateQuery<vertex_id_t> query, AdjList<edge_data_t>* adj_list)
        {
            stateResponse<bool> response;
            response.walker_idx = query.walker_idx;
            AdjUnit<edge_data_t> target;
            target.neighbour = query.data;
            response.data = std::binary_search(adj_list->begin, adj_list->end, target, [](const AdjUnit<edge_data_t> &a, const AdjUnit<edge_data_t> &b) { return a.neighbour < b.neighbour; });
            graph->emit(query.src_v, response);
        },
        [&] (Walker<Node2vecState> &walker, stateResponse<bool> &response, vertex_id_t current_v, AdjUnit<edge_data_t> *edge)
        {
            if (walker.step == 0)
            {
                return upperbound;
            } else
            {
                if (walker.data.previous_vertex == edge->neighbour)
                {
                    return 1 / p;
                } else if (response.data)
                {
                    return (real_t) 1;
                } else
                {
                    return 1 / q;
                }
            }
        },
        [&] (vertex_id_t v_id, AdjList<edge_data_t> *adj_lists)
        {
            return upperbound;
        },
        [&] (vertex_id_t v_id, AdjList<edge_data_t> *adj_lists)
        {
            return lowerbound;
        }
    );

#ifndef UNIT_TEST
    printf("total time %lfs\n", timer.duration());
#endif
}
