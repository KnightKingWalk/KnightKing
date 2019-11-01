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
#include "metascheme.hpp"

typedef uint8_t scheme_mask_t;

std::vector<std::vector<scheme_mask_t> > get_scheme_mask(std::vector<std::vector<std::vector<bool> > > &schemes)
{
    std::vector<std::vector<scheme_mask_t> > ret(schemes.size());

    for (size_t i = 0; i < schemes.size(); i++)
    {
        ret[i].resize(schemes[i].size());
        for (size_t j = 0; j < schemes[i].size(); j++)
        {
            scheme_mask_t val = 0;
            for (size_t k = 0; k < schemes[i][j].size(); k++)
            {
                if (schemes[i][j][k])
                {
                    val |= (1 << k);
                }
            }
            ret[i][j] = val;
        }
    }
    return ret;
}

void metapath(WalkEngine<int, MetapathState> *graph, std::vector<std::vector<std::vector<bool> > > schemes, walker_id_t walker_num, step_t walk_length)
{
    MPI_Barrier(MPI_COMM_WORLD);
    Timer timer;

    auto scheme_masks = get_scheme_mask(schemes);
    scheme_mask_t* vertex_masks = graph->alloc_vertex_array<scheme_mask_t>();
    graph->template process_vertices<vertex_id_t>(
        [&] (vertex_id_t v_id)
        {
            vertex_masks[v_id] = 0;
            for (auto p = graph->csr->adj_lists[v_id].begin; p != graph->csr->adj_lists[v_id].end; p++)
            {
                vertex_masks[v_id] |= (1 << p->data);
            }
            return 0;
        }
    );
    graph->set_walkers(
        walker_num,
        [&] (Walker<MetapathState> &walker, vertex_id_t start_vertex)
        {
            walker.data.scheme_id = graph->get_thread_local_rand_gen()->gen(schemes.size());
            walker.data.state = 0;
        },
        [&] (Walker<MetapathState> &walker, vertex_id_t current_v, AdjUnit<int> *edge)
        {
            walker.data.state = (walker.data.state + 1) % schemes[walker.data.scheme_id].size();
        }
    );
    graph->random_walk(
        [&] (Walker<MetapathState> &walker, vertex_id_t current_v)
        {
            return (walker.step >= walk_length || !(vertex_masks[current_v] & scheme_masks[walker.data.scheme_id][walker.data.state])) ? 0.0 : 1.0;
        },
        nullptr,
        [&] (Walker<MetapathState> &walker, vertex_id_t current_v, AdjUnit<int> *edge)
        {
            if (schemes[walker.data.scheme_id][walker.data.state][edge->data])
            //if (sh[walker.data.scheme_id][walker.data.state] & (1 << edge->data))
            {
                return 1.0;
            } else
            {
                return 0.0;
            }
        },
        [&] (vertex_id_t v_id, AdjList<int> *adj_lists)
        {
            return 1.0;
        }
    );
    graph->dealloc_vertex_array(vertex_masks);

#ifndef UNIT_TEST
    printf("total time %lfs\n", timer.duration());
#endif
}

void metapath(WalkEngine<WeightedMetaData, MetapathState> *graph, std::vector<std::vector<std::vector<bool> > > schemes, walker_id_t walker_num, step_t walk_length)
{
    MPI_Barrier(MPI_COMM_WORLD);
    Timer timer;

    auto scheme_masks = get_scheme_mask(schemes);
    scheme_mask_t* vertex_masks = graph->alloc_vertex_array<scheme_mask_t>();
    graph->template process_vertices<vertex_id_t>(
        [&] (vertex_id_t v_id)
        {
            vertex_masks[v_id] = 0;
            for (auto p = graph->csr->adj_lists[v_id].begin; p != graph->csr->adj_lists[v_id].end; p++)
            {
                vertex_masks[v_id] |= (1 << p->data.meta_info);
            }
            return 0;
        }
    );
    graph->set_walkers(
        walker_num,
        [&] (Walker<MetapathState> &walker, vertex_id_t start_vertex)
        {
            walker.data.scheme_id = graph->get_thread_local_rand_gen()->gen(schemes.size());
            walker.data.state = 0;
        },
        [&] (Walker<MetapathState> &walker, vertex_id_t current_v, AdjUnit<WeightedMetaData> *edge)
        {
            walker.data.state = (walker.data.state + 1) % schemes[walker.data.scheme_id].size();
        }
    );
    graph->random_walk(
        [&] (Walker<MetapathState> &walker, vertex_id_t current_v)
        {
            return (walker.step >= walk_length || !(vertex_masks[current_v] & scheme_masks[walker.data.scheme_id][walker.data.state])) ? 0.0 : 1.0;
        },
        [&] (vertex_id_t v, AdjUnit<WeightedMetaData> *edge)
        {
            return edge->data.weight;
        },
        [&] (Walker<MetapathState> &walker, vertex_id_t current_v, AdjUnit<WeightedMetaData> *edge)
        {
            if (schemes[walker.data.scheme_id][walker.data.state][edge->data.meta_info])
            {
                return 1.0;
            } else
            {
                return 0.0;
            }
        },
        [&] (vertex_id_t v_id, AdjList<WeightedMetaData> *adj_lists)
        {
            return 1.0;
        }
    );
    graph->dealloc_vertex_array(vertex_masks);

#ifndef UNIT_TEST
    printf("total time %lfs\n", timer.duration());
#endif
}
