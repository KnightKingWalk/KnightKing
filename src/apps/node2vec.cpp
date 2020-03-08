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
#include "node2vec.hpp"

class Node2vecOptionHelper : public STruncatedRandomWalkOptionHelper
{
private:
    args::ValueFlag<float> p_flag;
    args::ValueFlag<float> q_flag;
public:
    float p;
    float q;
    Node2vecOptionHelper():
        p_flag(parser, "p", "hyperparameter p", {'p'}),
        q_flag(parser, "q", "hyperparameter q", {'q'})
    {}
    virtual void parse(int argc, char** argv)
    {
        STruncatedRandomWalkOptionHelper::parse(argc, argv);

        assert(p_flag);
        p = args::get(p_flag);

        assert(q_flag);
        q = args::get(q_flag);
    }
    virtual Node2vecConf get_n2v_conf()
    {
        Node2vecConf n2v_conf;
        n2v_conf.p = p;
        n2v_conf.q = q;
        n2v_conf.walker_num = walker_num;
        n2v_conf.walk_length = walk_length;
        return n2v_conf;
    }
};

template<typename edge_data_t>
void run(WalkEngine<edge_data_t, Node2vecState> *graph, Node2vecOptionHelper *opt)
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
    Node2vecConf n2v_conf = opt->get_n2v_conf();
    node2vec(graph, n2v_conf, &walk_conf);
}

int main(int argc, char** argv)
{
    MPI_Instance mpi_instance(&argc, &argv);

    Node2vecOptionHelper opt;
    opt.parse(argc, argv);

    puts(opt.static_comp.c_str());
    if (opt.static_comp.compare("weighted") == 0)
    {
        WalkEngine<real_t, Node2vecState> graph;
        run(&graph, &opt);
    } else if(opt.static_comp.compare("unweighted") == 0)
    {
        WalkEngine<EmptyData, Node2vecState> graph;
        run(&graph, &opt);
    } else
    {
        exit(1);
    }
    return 0;
}
