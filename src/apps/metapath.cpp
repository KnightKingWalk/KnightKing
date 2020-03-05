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
#include "metapath.hpp"

class MetapathOptionHelper : public STruncatedRandomWalkOptionHelper
{
private:
    args::ValueFlag<std::string> schemes_path_flag;
public:
    std::string schemes_path;
    MetapathOptionHelper():
        schemes_path_flag(parser, "schemes", "schemes file path", {"schemes"})
    {}
    virtual void parse(int argc, char** argv)
    {
        STruncatedRandomWalkOptionHelper::parse(argc, argv);

        assert(schemes_path_flag);
        schemes_path = args::get(schemes_path_flag);
    }
};

template<typename edge_data_t>
void run(WalkEngine<edge_data_t,MetapathState > *graph, MetapathOptionHelper *opt)
{
    graph->load_graph(opt->v_num, opt->graph_path.c_str(), opt->make_undirected);
    WalkConfig walk_conf;
    if (!opt->output_path.empty())
    {
        walk_conf.set_output_file(opt->output_path.c_str());
    }
    if (opt->set_rate)
    {
        walk_conf.set_walk_rate(opt->rate);
    }
    metapath(graph, read_metapath_schemes(opt->schemes_path.c_str()), opt->walker_num, opt->walk_length, &walk_conf);
}

int main(int argc, char** argv)
{
    MPI_Instance mpi_instance(&argc, &argv);

    MetapathOptionHelper opt;
    opt.parse(argc, argv);

    if (opt.static_comp.compare("weighted") == 0)
    {
        WalkEngine<WeightedMetaData, MetapathState> graph;
        run(&graph, &opt);
    } else if(opt.static_comp.compare("unweighted") == 0)
    {
        WalkEngine<UnweightedMetaData, MetapathState> graph;
        run(&graph, &opt);
    } else
    {
        exit(1);
    }
    return 0;
}
