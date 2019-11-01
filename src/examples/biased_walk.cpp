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

int main(int argc, char** argv)
{
    MPI_Instance mpi_instance(&argc, &argv);

    TruncatedRandomWalkOptionHelper opt;
    opt.parse(argc, argv);

    WalkEngine<real_t, EmptyData> graph;
    graph.load_graph(opt.v_num, opt.graph_path.c_str());
    if (!opt.output_path.empty())
    {
        graph.set_output();
    }
    graph.set_walkers(opt.walker_num);
    auto extension_comp = [&] (Walker<EmptyData>& walker, vertex_id_t current_v)
    {
        return walker.step >= opt.walk_length ? 0.0 : 1.0; /*walk opt.walk_length steps then terminate*/
    };
    auto static_comp = [&] (vertex_id_t v, AdjUnit<real_t> *edge)
    {
        return edge->data; /*edge->data is a real number denoting edge weight*/
    };
    graph.random_walk(extension_comp, static_comp);
    if (!opt.output_path.empty())
    {
        PathSet path_data = graph.get_path_data();
        graph.dump_path_data(path_data, opt.output_path.c_str());
        graph.free_path_data(path_data);
    }
    return 0;
}
