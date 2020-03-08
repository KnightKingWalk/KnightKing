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

    RandomWalkOptionHelper opt;
    opt.parse(argc, argv);

    WalkEngine<real_t, EmptyData> graph;
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
    WalkerConfig<real_t, EmptyData> walker_conf(opt.walker_num);
    auto extension_comp = [&] (Walker<EmptyData>& walker, vertex_id_t current_v)
    {
        return 0.875; /*the probability to continue the walk*/
    };
    TransitionConfig<real_t, EmptyData> tr_conf(extension_comp);
    graph.random_walk(&walker_conf, &tr_conf, &walk_conf);
    return 0;
}
