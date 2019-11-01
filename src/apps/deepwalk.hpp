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

void deepwalk(WalkEngine<EmptyData, EmptyData> *graph, walker_id_t walker_num, step_t walk_length)
{
    MPI_Barrier(MPI_COMM_WORLD);
    Timer timer;

    graph->set_walkers(walker_num);
    graph->random_walk(
        [&] (Walker<EmptyData>& walker, vertex_id_t current_v)
        {
            return walker.step >= walk_length ? 0.0 : 1.0;
        }
    );

#ifndef UNIT_TEST
    printf("total time %lfs\n", timer.duration());
#endif
}

void deepwalk(WalkEngine<real_t, EmptyData> *graph, walker_id_t walker_num, step_t walk_length)
{
    MPI_Barrier(MPI_COMM_WORLD);
    Timer timer;

    graph->set_walkers(walker_num);
    graph->random_walk(
        [&] (Walker<EmptyData>& walker, vertex_id_t current_v)
        {
            return walker.step >= walk_length ? 0.0 : 1.0;
        },
        [&] (vertex_id_t v, AdjUnit<real_t> *edge)
        {
            return edge->data;
        }
    );

#ifndef UNIT_TEST
    printf("total time %lfs\n", timer.duration());
#endif
}
