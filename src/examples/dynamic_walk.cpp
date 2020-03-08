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

#include "walk.hpp"
#include "option_helper.hpp"

struct WalkState
{
    vertex_id_t last_vertex;
};

int main(int argc, char** argv)
{
    MPI_Instance mpi_instance(&argc, &argv);

    TruncatedRandomWalkOptionHelper opt;
    opt.parse(argc, argv);

    WalkEngine<real_t, WalkState> graph;
    graph.load_graph(opt.v_num, opt.graph_path.c_str(), opt.make_undirected);
    WalkConfig walk_conf;
    if (!opt.output_path.empty())
    {
        walk_conf.set_output_file(opt.output_path.c_str());
    }
    if (opt.set_rate)
    {
        walk_conf.set_walk_rate(opt.rate);
    }

    auto init_walker_func = [&] (Walker<WalkState> &walker, vertex_id_t start_vertex)
    {
        /*At first, the last vertex is not defined*/
        walker.data.last_vertex = UINT_MAX;
    };
    auto update_walker_func = [&] (Walker<WalkState> &walker, vertex_id_t current_v, AdjUnit<real_t> *edge)
    {
        walker.data.last_vertex = current_v;
    };
    WalkerConfig<real_t, WalkState> walker_conf(34, init_walker_func, update_walker_func);

    auto extension_comp = [&] (Walker<WalkState>& walker, vertex_id_t current_v)
    {
        return walker.step >= opt.walk_length ? 0.0 : 1.0; /*walk walk_length steps then terminate*/
    };
    auto static_comp = [&] (vertex_id_t v, AdjUnit<real_t> *edge)
    {
        return edge->data; /*edge->data is a real number denoting edge weight*/
    };
    auto dynamic_comp = [&] (Walker<WalkState> &walker, vertex_id_t current_v, AdjUnit<real_t> *edge)
    {
        if (walker.step == 0)
        {
            /*No return edge for the first step*/
            return 1.0;
        } else if (edge->neighbour == walker.data.last_vertex)
        {
            /*if return edge, double the un-normalized transition probability*/
            return 2.0;
        } else
        {
            /*if not return edge*/
            return 1.0;
        }
    };
    auto dynamic_comp_upperbound = [&] (vertex_id_t v_id, AdjList<real_t> *adj_lists)
    {
        return 2.0;
    };

    TransitionConfig<real_t, WalkState> tr_conf(extension_comp, static_comp, dynamic_comp, dynamic_comp_upperbound);
    graph.random_walk(&walker_conf, &tr_conf, &walk_conf);

    return 0;
}
