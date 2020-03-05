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

#include "type.hpp"
#include "graph.hpp"
#include "path.hpp"

template<typename walker_data_t>
struct Walker
{
public:
    walker_id_t id;
    step_t step;
    walker_data_t data;
};

template<>
struct Walker<EmptyData>
{
public:
    uint32_t id;
    union
    {
        step_t step;
        EmptyData data;
    };
};

template<typename edge_data_t>
struct AliasBucket
{
    real_t p;
    AdjUnit<edge_data_t> *p_ptr, *q_ptr;
};

template<typename edge_data_t>
struct AliasTableContainer
{
    AliasBucket<edge_data_t> *buckets;
    AliasBucket<edge_data_t> **index;
    AliasTableContainer() : buckets(nullptr), index(nullptr) {}
    ~AliasTableContainer()
    {
    }
};

template <typename query_data_t>
struct stateQuery
{
    vertex_id_t src_v;
    walker_id_t walker_idx;
    query_data_t data;
};

template<>
struct stateQuery <EmptyData>
{
    vertex_id_t src_v;
    union
    {
        walker_id_t walker_idx;
        EmptyData data;
    };
};

template<typename response_data_t>
struct stateResponse
{
    walker_id_t walker_idx;
    response_data_t data;
};

template<typename edge_data_t>
struct SecondOrderCandidate
{
    AdjUnit<edge_data_t> *candidate;
    real_t randval;
    bool accepted;
};


class WalkConfig
{
public:
    // output setting
    bool output_file_flag;
    std::string output_path_prefix;
    bool print_with_head_info;
    bool output_consumer_flag;
    std::function<void (PathSet*)> output_consumer_func;
    double rate;

    WalkConfig()
    {
        output_file_flag = false;
        print_with_head_info = false;
        output_consumer_flag = false;
        output_consumer_func = nullptr;
        rate = 1;
    }

    void set_output_file(const char* output_path_prefix_param, bool with_head_info = false)
    {
        assert(output_path_prefix_param != nullptr);
        this->output_file_flag = true;
        this->output_path_prefix = std::string(output_path_prefix_param);
        this->print_with_head_info = with_head_info;
    }

    void set_output_consumer(std::function<void (PathSet*)> output_consumer_func_param)
    {
        assert(output_consumer_func_param != nullptr);
        this->output_consumer_flag = true;
        this->output_consumer_func = output_consumer_func_param;
    }

    void set_walk_rate(double rate_param)
    {
        assert(rate_param > 0 && rate_param <= 1);
        this->rate = rate_param;
    }
};

template<typename edge_data_t, typename walker_data_t>
class WalkerConfig
{
public:
    // walker setting
    walker_id_t walker_num;
    std::function<vertex_id_t (walker_id_t)> walker_init_dist_func;
    std::function<void (Walker<walker_data_t>&, vertex_id_t)> walker_init_state_func;
    std::function<void (Walker<walker_data_t>&, vertex_id_t, AdjUnit<edge_data_t> *)> walker_update_state_func;

    WalkerConfig()
    {
        walker_num = 0;
        walker_init_dist_func = nullptr;
        walker_init_state_func = nullptr;
        walker_update_state_func = nullptr;
    }

    WalkerConfig(
        walker_id_t walker_num_param,
        std::function<void (Walker<walker_data_t>&, vertex_id_t)> walker_init_state_func_param = nullptr,
        std::function<void (Walker<walker_data_t>&, vertex_id_t, AdjUnit<edge_data_t> *)> walker_update_state_func_param = nullptr,
        std::function<vertex_id_t (walker_id_t)> walker_init_dist_func_param = nullptr
    )
    {
        WalkerConfig();
        set_walkers(
            walker_num_param,
            walker_init_state_func_param,
            walker_update_state_func_param,
            walker_init_dist_func_param
        );
    }

    void set_walkers(
        walker_id_t walker_num_param,
        std::function<void (Walker<walker_data_t>&, vertex_id_t)> walker_init_state_func_param = nullptr,
        std::function<void (Walker<walker_data_t>&, vertex_id_t, AdjUnit<edge_data_t> *)> walker_update_state_func_param = nullptr,
        std::function<vertex_id_t (walker_id_t)> walker_init_dist_func_param = nullptr
    )
    {
        this->walker_num = walker_num_param;
        this->walker_init_state_func = walker_init_state_func_param;
        this->walker_update_state_func = walker_update_state_func_param;
        this->walker_init_dist_func = walker_init_dist_func_param;
    }
};

template<typename edge_data_t, typename walker_data_t>
class TransitionConfig
{
public:
	// random walk setting
	std::function<real_t (Walker<walker_data_t>&, vertex_id_t)> extension_comp_func;
	std::function<real_t(vertex_id_t, AdjUnit<edge_data_t>*)> static_comp_func;
	std::function<real_t (Walker<walker_data_t>&, vertex_id_t, AdjUnit<edge_data_t> *)> dynamic_comp_func;
	std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_upperbound_func;
	std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_lowerbound_func;
	std::function<void(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, real_t&, vertex_id_t&)> outlier_upperbound_func;
	std::function<AdjUnit<edge_data_t>*(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, vertex_id_t)> outlier_search_func;

    TransitionConfig()
    {
        extension_comp_func = nullptr;
        static_comp_func = nullptr;
        dynamic_comp_func = nullptr;
        dcomp_upperbound_func = nullptr;
        dcomp_lowerbound_func = nullptr;
        outlier_upperbound_func = nullptr;
        outlier_search_func = nullptr;
    }

    TransitionConfig(
        std::function<real_t (Walker<walker_data_t>&, vertex_id_t)> extension_comp_func_param,
        std::function<real_t(vertex_id_t, AdjUnit<edge_data_t>*)> static_comp_func_param = nullptr,
        std::function<real_t (Walker<walker_data_t>&, vertex_id_t, AdjUnit<edge_data_t> *)> dynamic_comp_func_param = nullptr,
        std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_upperbound_func_param = nullptr,
        std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_lowerbound_func_param = nullptr,
        std::function<void(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, real_t&, vertex_id_t&)> outlier_upperbound_func_param = nullptr,
        std::function<AdjUnit<edge_data_t>*(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, vertex_id_t)> outlier_search_func_param = nullptr
    )
    {
        TransitionConfig();
        set_transition(
            extension_comp_func_param,
            static_comp_func_param,
            dynamic_comp_func_param,
            dcomp_upperbound_func_param,
            dcomp_lowerbound_func_param,
            outlier_upperbound_func_param,
            outlier_search_func_param
        );
    }

    void set_transition(
        std::function<real_t (Walker<walker_data_t>&, vertex_id_t)> extension_comp_func_param,
        std::function<real_t(vertex_id_t, AdjUnit<edge_data_t>*)> static_comp_func_param = nullptr,
        std::function<real_t (Walker<walker_data_t>&, vertex_id_t, AdjUnit<edge_data_t> *)> dynamic_comp_func_param = nullptr,
        std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_upperbound_func_param = nullptr,
        std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_lowerbound_func_param = nullptr,
        std::function<void(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, real_t&, vertex_id_t&)> outlier_upperbound_func_param = nullptr,
        std::function<AdjUnit<edge_data_t>*(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, vertex_id_t)> outlier_search_func_param = nullptr
    )
    {
        if (dynamic_comp_func_param != nullptr)
        {
            assert(dcomp_upperbound_func_param != nullptr);
        } else
        {
            assert(dcomp_upperbound_func_param == nullptr);
            assert(dcomp_lowerbound_func_param == nullptr);
            assert(outlier_upperbound_func_param == nullptr);
            assert(outlier_search_func_param == nullptr);
        }
        assert(outlier_upperbound_func_param == nullptr && outlier_search_func_param == nullptr || outlier_upperbound_func_param != nullptr &&  outlier_search_func_param != nullptr);

        this->extension_comp_func = extension_comp_func_param;
        this->static_comp_func = static_comp_func_param;
        this->dynamic_comp_func = dynamic_comp_func_param;
        this->dcomp_upperbound_func = dcomp_upperbound_func_param;
        this->dcomp_lowerbound_func = dcomp_lowerbound_func_param;
        this->outlier_upperbound_func = outlier_upperbound_func_param;
        this->outlier_search_func = outlier_search_func_param;
    }
};

template<typename edge_data_t, typename walker_data_t, typename query_data_t, typename response_data_t>
class SecondOrderTransitionConfig
{
public:
	// random walk setting
	std::function<real_t (Walker<walker_data_t>&, vertex_id_t)> extension_comp_func;
	std::function<real_t(vertex_id_t, AdjUnit<edge_data_t>*)> static_comp_func;
    std::function<void (Walker<walker_data_t>&, walker_id_t, vertex_id_t, AdjUnit<edge_data_t> *)> post_query_func;
    std::function<void (vertex_id_t, stateQuery<query_data_t> &, AdjList<edge_data_t>*)> respond_query_func;
    std::function<real_t (Walker<walker_data_t>&, stateResponse<response_data_t> &, vertex_id_t, AdjUnit<edge_data_t> *)> dynamic_comp_func;
	std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_upperbound_func;
	std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_lowerbound_func;
	std::function<void(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, real_t&, vertex_id_t&)> outlier_upperbound_func;
	std::function<AdjUnit<edge_data_t>*(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, vertex_id_t)> outlier_search_func;

    SecondOrderTransitionConfig()
    {
        this->extension_comp_func = nullptr;
        this->static_comp_func = nullptr;
        this->post_query_func = nullptr;
        this->respond_query_func = nullptr;
        this->dynamic_comp_func = nullptr;
        this->dcomp_upperbound_func = nullptr;
        this->dcomp_lowerbound_func = nullptr;
        this->outlier_upperbound_func = nullptr;
        this->outlier_search_func = nullptr;
    }

    SecondOrderTransitionConfig(
        std::function<real_t (Walker<walker_data_t>&, vertex_id_t)> extension_comp_func_param,
        std::function<real_t(vertex_id_t, AdjUnit<edge_data_t>*)> static_comp_func_param,
        std::function<void (Walker<walker_data_t>&, walker_id_t, vertex_id_t, AdjUnit<edge_data_t> *)> post_query_func_param,
        std::function<void (vertex_id_t, stateQuery<query_data_t> &, AdjList<edge_data_t>*)> respond_query_func_param,
        std::function<real_t (Walker<walker_data_t>&, stateResponse<response_data_t> &, vertex_id_t, AdjUnit<edge_data_t> *)> dynamic_comp_func_param,
        std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_upperbound_func_param,
        std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_lowerbound_func_param = nullptr,
        std::function<void(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, real_t&, vertex_id_t&)> outlier_upperbound_func_param = nullptr,
        std::function<AdjUnit<edge_data_t>*(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, vertex_id_t)> outlier_search_func_param = nullptr
    )
    {
        SecondOrderTransitionConfig();
        set_transition(
            extension_comp_func_param,
            static_comp_func_param,
            post_query_func_param,
            respond_query_func_param,
            dynamic_comp_func_param,
            dcomp_upperbound_func_param,
            dcomp_lowerbound_func_param,
            outlier_upperbound_func_param,
            outlier_search_func_param
        );
    }

    void set_transition (
        std::function<real_t (Walker<walker_data_t>&, vertex_id_t)> extension_comp_func_param,
        std::function<real_t(vertex_id_t, AdjUnit<edge_data_t>*)> static_comp_func_param,
        std::function<void (Walker<walker_data_t>&, walker_id_t, vertex_id_t, AdjUnit<edge_data_t> *)> post_query_func_param,
        std::function<void (vertex_id_t, stateQuery<query_data_t> &, AdjList<edge_data_t>*)> respond_query_func_param,
        std::function<real_t (Walker<walker_data_t>&, stateResponse<response_data_t> &, vertex_id_t, AdjUnit<edge_data_t> *)> dynamic_comp_func_param,
        std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_upperbound_func_param,
        std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_lowerbound_func_param = nullptr,
        std::function<void(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, real_t&, vertex_id_t&)> outlier_upperbound_func_param = nullptr,
        std::function<AdjUnit<edge_data_t>*(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, vertex_id_t)> outlier_search_func_param = nullptr
    )
    {
        if (dynamic_comp_func_param != nullptr)
        {
            assert(dcomp_upperbound_func_param != nullptr);
        } else
        {
            assert(dcomp_upperbound_func_param == nullptr);
            assert(dcomp_lowerbound_func_param == nullptr);
            assert(outlier_upperbound_func_param == nullptr);
            assert(outlier_search_func_param == nullptr);
        }
        assert(post_query_func_param != nullptr);
        assert(respond_query_func_param != nullptr);
        assert(outlier_upperbound_func_param == nullptr && outlier_search_func_param == nullptr || outlier_upperbound_func_param != nullptr && outlier_search_func_param != nullptr);

        this->extension_comp_func = extension_comp_func_param;
        this->static_comp_func = static_comp_func_param;
        this->post_query_func = post_query_func_param;
        this->respond_query_func = respond_query_func_param;
        this->dynamic_comp_func = dynamic_comp_func_param;
        this->dcomp_upperbound_func = dcomp_upperbound_func_param;
        this->dcomp_lowerbound_func = dcomp_lowerbound_func_param;
        this->outlier_upperbound_func = outlier_upperbound_func_param;
        this->outlier_search_func = outlier_search_func_param;
    }
};

template<typename edge_data_t, typename walker_data_t>
class WalkEngine : public GraphEngine<edge_data_t>
{
    StdRandNumGenerator* randgen;

#ifdef COLLECT_WALK_SEQUENCE
    std::vector<std::vector<Footprint> > footprints;
#endif
#ifdef COLLECT_WALKER_INIT_STATE
    std::vector<Walker<walker_data_t> > walker_init_state;
#endif

public:
    WalkEngine()
    {
        randgen = new StdRandNumGenerator[this->worker_num];
    }

    ~WalkEngine()
    {
        if (randgen != nullptr)
        {
            delete []randgen;
        }
    }

    void set_concurrency(int worker_num_param)
    {
        this->set_graph_engine_concurrency(worker_num_param);
        delete []randgen;
        randgen = new StdRandNumGenerator[worker_num_param];
    }

    StdRandNumGenerator* get_thread_local_rand_gen()
    {
        return &randgen[omp_get_thread_num()];
    }

    std::function<vertex_id_t (walker_id_t)> get_equal_dist_func()
    {
        auto equal_dist_func = [&] (walker_id_t w_id)
        {
            vertex_id_t start_v = w_id % this->v_num;
            return start_v;
        };
        return equal_dist_func;
    }

    std::function<vertex_id_t (walker_id_t)> get_uniform_dist_func()
    {
        auto uniform_dist_func = [&] (walker_id_t w_id)
        {
            vertex_id_t start_v = get_thread_local_rand_gen()->gen(this->v_num);
            return start_v;
        };
        return uniform_dist_func;
    }

    template<typename transition_config_t>
    void random_walk(WalkerConfig<edge_data_t, walker_data_t> *walker_config, transition_config_t *transition_config, WalkConfig *walk_config_param = nullptr)
    {
        WalkConfig* walk_config = walk_config_param;
        if (walk_config_param == nullptr)
        {
            walk_config = new WalkConfig();
        }
        internal_random_walk_wrap(walker_config, transition_config, walk_config);
        if (walk_config_param == nullptr)
        {
            delete walk_config;
        }
    }

private:

    walker_id_t init_walkers(
        Message<Walker<walker_data_t> >* &local_walkers,
        Message<Walker<walker_data_t> >* &local_walkers_bak,
        walker_id_t walker_begin,
        walker_id_t walker_end,
        std::function<vertex_id_t (walker_id_t)> walker_init_dist_func,
        std::function<void (Walker<walker_data_t>&, vertex_id_t)> walker_init_state_func
    )
    {
        typedef Walker<walker_data_t> walker_t;
        typedef Message<walker_t> walker_msg_t;

        walker_id_t local_walker_num = 0;

        this->template distributed_execute<walker_t>(
            [&] (void) {
                #pragma omp parallel for
                for (walker_id_t w_i = walker_begin / this->partition_num * this->partition_num + this->local_partition_id; w_i < walker_end; w_i += this->partition_num)
                {
                    if (w_i < walker_begin)
                    {
                        continue;
                    }
                    vertex_id_t start_v = walker_init_dist_func(w_i);
                    Walker<walker_data_t> walker;
                    walker.id = w_i;
                    walker.step = 0;
                    assert(start_v < this->v_num);
                    this->emit(start_v, walker, omp_get_thread_num());
                #ifdef COLLECT_WALK_SEQUENCE
                    footprints[omp_get_thread_num()].push_back(Footprint(walker.id, start_v, walker.step));
                #endif
                }
            },
            [&](walker_msg_t *begin, walker_msg_t *end)
            {
                local_walker_num = end - begin;
                std::swap(local_walkers, local_walkers_bak);
            },
            local_walkers_bak
        );
        if (walker_init_state_func != nullptr)
        {
            #pragma omp parallel for
            for (walker_id_t w_i = 0; w_i < local_walker_num; w_i ++)
            {
                walker_init_state_func(local_walkers[w_i].data, local_walkers[w_i].dst_vertex_id);
            }
        }
        #ifdef COLLECT_WALKER_INIT_STATE
        walker_init_state.clear();
        for (walker_id_t w_i = 0; w_i < local_walker_num; w_i ++)
        {
            walker_init_state.push_back(local_walkers[w_i].data);
        }
        #endif
        return local_walker_num;
    }

    void init_dcomp_upperbound(std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_upperbound_func, real_t* &dcomp_upperbound)
    {
        assert(dcomp_upperbound == nullptr);
        dcomp_upperbound = this->template alloc_vertex_array<real_t>();

        vertex_id_t local_v_begin = this->get_local_vertex_begin();
        vertex_id_t local_v_end = this->get_local_vertex_end();
#pragma omp parallel for
        for (vertex_id_t v_i = local_v_begin; v_i < local_v_end; v_i++)
        {
            dcomp_upperbound[v_i] = dcomp_upperbound_func(v_i, this->csr->adj_lists + v_i);
        }
    }

    void free_dcomp_upperbound(real_t *dcomp_upperbound)
    {
        assert(dcomp_upperbound != nullptr);
        this->dealloc_vertex_array(dcomp_upperbound);
    }

    void init_dcomp_lowerbound(std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_lowerbound_func, real_t* &dcomp_lowerbound)
    {
        assert(dcomp_lowerbound == nullptr);
        dcomp_lowerbound = this->template alloc_vertex_array<real_t>();
        vertex_id_t local_v_begin = this->get_local_vertex_begin();
        vertex_id_t local_v_end = this->get_local_vertex_end();
#pragma omp parallel for
        for (vertex_id_t v_i = local_v_begin; v_i < local_v_end; v_i++)
        {
            dcomp_lowerbound[v_i] = dcomp_lowerbound_func(v_i, this->csr->adj_lists + v_i);
        }
    }

    void free_dcomp_lowerbound(real_t *dcomp_lowerbound)
    {
        assert(dcomp_lowerbound != nullptr);
        this->dealloc_vertex_array(dcomp_lowerbound);
    }

    void init_alias_tables(std::function<real_t(vertex_id_t, AdjUnit<edge_data_t>*)> static_comp_func, AliasTableContainer<edge_data_t> *&alias_tables, real_t* regular_area = nullptr)
    {
        assert(static_comp_func != nullptr);
        assert(alias_tables == nullptr);

        Timer timer;
        vertex_id_t local_v_begin = this->get_local_vertex_begin();
        vertex_id_t local_v_end = this->get_local_vertex_end();
        if (alias_tables == nullptr)
        {
            alias_tables = new AliasTableContainer<edge_data_t>();
            alias_tables->index = this->template alloc_vertex_array<AliasBucket<edge_data_t>*>();
            alias_tables->buckets = this->template alloc_array<AliasBucket<edge_data_t> >(this->local_e_num);
        }
        edge_id_t local_edge_p = 0;
        vertex_id_t max_degree = 0;
        for (vertex_id_t v_i = local_v_begin; v_i != local_v_end; v_i++)
        {
            alias_tables->index[v_i] = alias_tables->buckets + local_edge_p;
            local_edge_p += this->vertex_out_degree[v_i];
            max_degree = std::max(max_degree, this->vertex_out_degree[v_i]);
        }
        vertex_id_t progress = local_v_begin;
        #pragma omp parallel
        {
            int worker_id = omp_get_thread_num();
            const vertex_id_t work_step_length = PARALLEL_CHUNK_SIZE;
            vertex_id_t next_workload;
            vertex_id_t *p_set = new vertex_id_t[max_degree];
            vertex_id_t *q_set = new vertex_id_t[max_degree];
            real_t *prob = new real_t[max_degree];
            while((next_workload =  __sync_fetch_and_add(&progress, work_step_length)) < local_v_end)
            {
                vertex_id_t v_begin = next_workload;
                vertex_id_t v_end = std::min(next_workload + work_step_length, local_v_end);
                for (vertex_id_t v_i = v_begin; v_i != v_end; v_i++)
                {
                    AdjUnit<edge_data_t> *begin = this->csr->adj_lists[v_i].begin;
                    AdjUnit<edge_data_t> *end = this->csr->adj_lists[v_i].end;
                    if (begin != end)
                    {
                        real_t sum = 0;
                        for (vertex_id_t e_i = 0; e_i < this->vertex_out_degree[v_i]; e_i++)
                        {
                            prob[e_i] = static_comp_func(v_i, begin + e_i); 
                            sum += prob[e_i];
                        }
                        if (regular_area != nullptr)
                        {
                            regular_area[v_i] = sum;
                        }
                        real_t ave_prob = sum / this->vertex_out_degree[v_i];
                        vertex_id_t p_set_head = 0;
                        vertex_id_t p_set_tail = 0;
                        vertex_id_t q_set_head = 0;
                        vertex_id_t q_set_tail = 0;
                        for (vertex_id_t e_i = 0; e_i < this->vertex_out_degree[v_i]; e_i++)
                        {
                            if (prob[e_i] < ave_prob)
                            {
                                p_set[p_set_tail++] = e_i;
                            } else
                            {
                                q_set[q_set_tail++] = e_i;
                            }
                        }
                        AliasBucket<edge_data_t>* alias_p = alias_tables->index[v_i];
                        while (p_set_tail != p_set_head)
                        {
                            vertex_id_t p_idx = p_set[p_set_head++];
                            if (q_set_tail != q_set_head)
                            {
                                vertex_id_t q_idx = q_set[q_set_head++];
                                alias_p->p_ptr = begin + p_idx;
                                alias_p->p = prob[p_idx] / ave_prob;
                                alias_p->q_ptr = begin + q_idx;
                                real_t diff = ave_prob - prob[p_idx];
                                prob[q_idx] -= diff;
                                if (prob[q_idx] < ave_prob)
                                {
                                    p_set[p_set_tail++] = q_idx;
                                } else
                                {
                                    q_set[q_set_tail++] = q_idx;
                                }
                            } else
                            {
                                alias_p->p_ptr = begin + p_idx;
                                alias_p->p = 1.0;
                                alias_p->q_ptr = nullptr;
                            }
                            alias_p ++;
                        }
                        while(q_set_tail != q_set_head)
                        {
                            vertex_id_t q_idx = q_set[q_set_head++];
                            alias_p->p_ptr = begin + q_idx;
                            alias_p->p = 1.0;
                            alias_p->q_ptr = nullptr;
                            alias_p++;
                        }
                        assert(q_set_tail <= this->vertex_out_degree[v_i]);
                        assert(p_set_tail <= this->vertex_out_degree[v_i]);
                        assert(alias_p == alias_tables->index[v_i] + this->vertex_out_degree[v_i]);
                    }
                }
            }
            delete []p_set;
            delete []q_set;
            delete []prob;
        }
#ifdef PERF_PROF
        printf("%d: init alias tables in %lf seconds\n", this->local_partition_id, timer.duration());
#endif
    }

    void free_alias_tables(AliasTableContainer<edge_data_t> *alias_tables)
    {
        assert(alias_tables != nullptr);
        if (alias_tables->index != nullptr)
        {
            this->dealloc_vertex_array(alias_tables->index);
        }
        if (alias_tables->buckets != nullptr)
        {
            this->dealloc_array(alias_tables->buckets, this->local_e_num);
        }
        delete alias_tables; 
    }

public:
    template<typename query_data_t, typename response_data_t>
    struct InternalWalkData
    {
        typedef Walker<walker_data_t> walker_t;
        typedef Message<walker_t> walker_msg_t;
        typedef stateQuery<query_data_t> query_t;
        typedef stateResponse<response_data_t> response_t;

        Timer timer;

        AliasTableContainer<edge_data_t>* alias_tables;
        real_t *dcomp_lowerbound;
        real_t *dcomp_upperbound;
        bool outlier_opt_flag;
        real_t *regular_area;

        response_t* remote_response_cache;
        SecondOrderCandidate<edge_data_t>* remote_fetch_candidate;
        Message<query_t>* cached_request;

        walker_msg_t *local_walkers;
        walker_msg_t *local_walkers_bak;
        walker_id_t local_walker_num;
        walker_id_t active_walker_num;

        bool collect_path_flag;
        PathCollector *pc;
    };

    template<typename query_data_t, typename response_data_t, typename transition_config_t>
    void internal_random_walk(WalkerConfig<edge_data_t, walker_data_t> *walker_config, transition_config_t *transition_config, WalkConfig* walk_config, int order)
    {
        typedef Walker<walker_data_t> walker_t;
        typedef Message<walker_t> walker_msg_t;
        typedef stateQuery<query_data_t> query_t;
        typedef stateResponse<response_data_t> response_t;

        walker_id_t walker_num = walker_config->walker_num;
        walker_id_t walker_per_iter = walker_num * walk_config->rate;
        if (walker_per_iter == 0) walker_per_iter = 1;
        if (walker_per_iter > walker_num) walker_per_iter = walker_num;
        size_t walker_array_size = walker_per_iter;
        walker_id_t remained_walker = walker_num;

        InternalWalkData<query_data_t, response_data_t> walk_data;

        walk_data.collect_path_flag = false;
        walk_data.pc = nullptr;
        if (walk_config->output_file_flag || walk_config->output_consumer_flag)
        {
            walk_data.collect_path_flag = true;
        }

        walk_data.outlier_opt_flag = false;
        if (transition_config->outlier_upperbound_func != nullptr)
        {
            walk_data.outlier_opt_flag = true;
        }

        walk_data.remote_response_cache = nullptr;
        walk_data.remote_fetch_candidate = nullptr;
        walk_data.cached_request = nullptr;
        if (order == 2)
        {
            walk_data.remote_response_cache = this->template alloc_array<response_t>(walker_array_size);
            walk_data.remote_fetch_candidate = this->template alloc_array<SecondOrderCandidate<edge_data_t> >(walker_array_size);
            walk_data.cached_request = this->template alloc_array<Message<query_t> > (walker_array_size);
        }

        walk_data.dcomp_upperbound = nullptr;
        if (transition_config->dcomp_upperbound_func != nullptr)
        {
            init_dcomp_upperbound(transition_config->dcomp_upperbound_func, walk_data.dcomp_upperbound);
        }

        walk_data.dcomp_lowerbound = nullptr;
        if (transition_config->dcomp_lowerbound_func != nullptr)
        {
            init_dcomp_lowerbound(transition_config->dcomp_lowerbound_func, walk_data.dcomp_lowerbound);
        }

        walk_data.regular_area = nullptr;
        if (walk_data.outlier_opt_flag)
        {
            walk_data.regular_area = this->template alloc_vertex_array<real_t>();
        }
        walk_data.alias_tables = nullptr;
        if (transition_config->static_comp_func != nullptr)
        {
            init_alias_tables(transition_config->static_comp_func, walk_data.alias_tables, walk_data.regular_area); 
            if (walk_data.outlier_opt_flag)
            {
                vertex_id_t local_v_begin = this->get_local_vertex_begin();
                vertex_id_t local_v_end = this->get_local_vertex_end();
#pragma omp parallel for
                for (vertex_id_t v_i = local_v_begin; v_i < local_v_end; v_i++)
                {
                    walk_data.regular_area[v_i] *= walk_data.dcomp_upperbound[v_i];
                }
            }
        } else
        {
            if (walk_data.outlier_opt_flag)
            {
                vertex_id_t local_v_begin = this->get_local_vertex_begin();
                vertex_id_t local_v_end = this->get_local_vertex_end();
#pragma omp parallel for
                for (vertex_id_t v_i = local_v_begin; v_i < local_v_end; v_i++)
                {
                    walk_data.regular_area[v_i] = this->vertex_out_degree[v_i] * walk_data.dcomp_upperbound[v_i];
                }
            }
        }

        if (walker_config->walker_init_dist_func == nullptr)
        {
            walker_config->walker_init_dist_func = this->get_equal_dist_func();
        }
        walk_data.local_walkers = this->template alloc_array<Message<Walker<walker_data_t> > >(walker_array_size);
        walk_data.local_walkers_bak = this->template alloc_array<Message<Walker<walker_data_t> > >(walker_array_size);

        size_t max_msg_size = 0;
        if (order == 1)
        {
            max_msg_size = sizeof(walker_msg_t);
        } else
        {
            max_msg_size = std::max(sizeof(walker_msg_t), std::max(sizeof(Message<query_t>), sizeof(Message<response_t>)));
        }
        this->set_msg_buffer(walker_array_size, max_msg_size);

#ifdef COLLECT_WALK_SEQUENCE
        footprints.resize(this->worker_num);
#endif
        int iter = 0;
        while (remained_walker != 0)
        {
            if (walk_data.collect_path_flag)
            {
                walk_data.pc = new PathCollector(this->worker_num);
            }
            walker_id_t walker_begin = walker_num - remained_walker;
            walk_data.active_walker_num = std::min(remained_walker, walker_per_iter);
            remained_walker -= walk_data.active_walker_num;
            walk_data.local_walker_num = init_walkers(walk_data.local_walkers, walk_data.local_walkers_bak, walker_begin, walker_begin + walk_data.active_walker_num, walker_config->walker_init_dist_func, walker_config->walker_init_state_func);

            if (walk_data.collect_path_flag)
            {
#pragma omp parallel for
                for (walker_id_t w_i = 0; w_i < walk_data.local_walker_num; w_i++)
                {
                    walk_data.pc->add_footprint(Footprint(walk_data.local_walkers[w_i].data.id, walk_data.local_walkers[w_i].dst_vertex_id, 0), omp_get_thread_num());
                }
            }

            internal_walk_epoch(&walk_data, walker_config, transition_config);

            if (walk_data.collect_path_flag)
            {
                auto* paths = walk_data.pc->assemble_path(walker_begin);
                if (walk_config->output_file_flag)
                {
                    std::string local_output_path = walk_config->output_path_prefix + "." + std::to_string(this->local_partition_id);
                    paths->dump(local_output_path.c_str(), iter == 0 ? "w": "a", walk_config->print_with_head_info);
                }
                if (walk_config->output_consumer_flag)
                {
                    walk_config->output_consumer_func(paths);
                } else
                {
                    delete paths;
                }
                delete walk_data.pc;
            }
            iter++;
        }

        if (walk_data.remote_response_cache != nullptr)
        {
            this->dealloc_array(walk_data.remote_response_cache, walker_array_size);
        }
        if (walk_data.remote_fetch_candidate != nullptr)
        {
            this->dealloc_array(walk_data.remote_fetch_candidate, walker_array_size);
        }
        if (walk_data.cached_request != nullptr)
        {
            this->dealloc_array(walk_data.cached_request, walker_array_size);
        }

        if (transition_config->static_comp_func != nullptr)
        {
            free_alias_tables(walk_data.alias_tables);
        }
        if (transition_config->dcomp_upperbound_func != nullptr)
        {
            free_dcomp_upperbound(walk_data.dcomp_upperbound);
        }
        if (transition_config->dcomp_lowerbound_func != nullptr)
        {
            free_dcomp_lowerbound(walk_data.dcomp_lowerbound);
        }
        if (walk_data.outlier_opt_flag)
        {
            this->dealloc_vertex_array(walk_data.regular_area);
        }

        this->dealloc_array(walk_data.local_walkers, walker_array_size);
        this->dealloc_array(walk_data.local_walkers_bak, walker_array_size);
    }

    void internal_random_walk_wrap (WalkerConfig<edge_data_t, walker_data_t> *walker_config, TransitionConfig<edge_data_t, walker_data_t> *transition_config, WalkConfig *walk_config)
    {
        internal_random_walk<EmptyData, EmptyData> (walker_config, transition_config, walk_config, 1);
    }

    template<typename query_data_t, typename response_data_t>
    void internal_random_walk_wrap (WalkerConfig<edge_data_t, walker_data_t> *walker_config, SecondOrderTransitionConfig<edge_data_t, walker_data_t, query_data_t, response_data_t> *transition_config, WalkConfig *walk_config)
    {
        internal_random_walk<query_data_t, response_data_t> (walker_config, transition_config, walk_config, 2);
    }

    template<typename query_data_t, typename response_data_t>
    void internal_walk_epoch(
        InternalWalkData<query_data_t, response_data_t> *walk_data,
        WalkerConfig<edge_data_t, walker_data_t> *walker_config,
        TransitionConfig<edge_data_t, walker_data_t> *transition_config
    )
    {
        typedef Walker<walker_data_t> walker_t;
        typedef Message<walker_t> walker_msg_t;

        auto extension_comp_func = transition_config->extension_comp_func;
        auto static_comp_func = transition_config->static_comp_func;
        auto dynamic_comp_func = transition_config->dynamic_comp_func;
        auto dcomp_upperbound_func = transition_config->dcomp_upperbound_func;
        auto dcomp_lowerbound_func = transition_config->dcomp_lowerbound_func;
        auto outlier_upperbound_func = transition_config->outlier_upperbound_func;
        auto outlier_search_func = transition_config->outlier_search_func;

        auto walker_update_state_func = walker_config->walker_update_state_func;

        auto* alias_tables = walk_data->alias_tables;
        auto* dcomp_lowerbound = walk_data->dcomp_lowerbound;
        auto* dcomp_upperbound = walk_data->dcomp_upperbound;
        auto outlier_opt_flag = walk_data->outlier_opt_flag;
        auto* regular_area = walk_data->regular_area;

        auto* local_walkers = walk_data->local_walkers;
        auto* local_walkers_bak = walk_data->local_walkers_bak;
        auto local_walker_num = walk_data->local_walker_num;

        auto output_flag = walk_data->collect_path_flag;
        auto* pc = walk_data->pc;

        auto active_walker_num = walk_data->active_walker_num;    
        int super_step = 0;
        while (active_walker_num != 0)
        {
            #ifndef UNIT_TEST
            if (this->local_partition_id == 0)
            {
                printf("step(%d), active(%u), time(%.3lf)\n", super_step++, active_walker_num, walk_data->timer.duration());
            }
            #endif
            bool use_parallel = (active_walker_num >= OMP_PARALLEL_THRESHOLD);
            active_walker_num = this->template distributed_execute<walker_t>(
                [&] (void) {
                    walker_id_t progress = 0;
                    walker_id_t data_amount = local_walker_num;
                    auto *data_begin = local_walkers;
                    const walker_id_t work_step_length = PARALLEL_CHUNK_SIZE;
                    #pragma omp parallel if (use_parallel)
                    {
                        int worker_id = omp_get_thread_num();
                        StdRandNumGenerator* gen = get_thread_local_rand_gen();
                        vertex_id_t next_workload;
                        while((next_workload =  __sync_fetch_and_add(&progress, work_step_length)) < data_amount)
                        {
                            walker_msg_t *begin = data_begin + next_workload;
                            walker_msg_t *end = data_begin + std::min(next_workload + work_step_length, data_amount);
                            for (walker_msg_t *p = begin; p != end; p++)
                            {
                                walker_t& walker = p->data;
                                vertex_id_t current_v = p->dst_vertex_id;
                                while (true)
                                {
                                    vertex_id_t degree = this->vertex_out_degree[current_v];
                                    bool terminate_flag = false;
                                    if (degree == 0)
                                    {
                                        terminate_flag = true;
                                    } else
                                    {
                                        real_t extension_comp = extension_comp_func(walker, current_v);
                                        // continue to walk at a probability of extension_comp
                                        if (extension_comp == 0.0 || (extension_comp < 1.0 && gen->gen_float(1.0) > extension_comp))
                                        {
                                            terminate_flag = true;
                                        }
                                    }
                                    if (terminate_flag)
                                    {
                                        break;
                                    }

                                    AdjList<edge_data_t> *adj = this->csr->adj_lists + current_v;
                                    AdjUnit<edge_data_t> *candidate = nullptr;
                                    AdjUnit<edge_data_t> *ac_edge = nullptr;
                                    
                                    while (ac_edge == nullptr)
                                    {
                                        if (outlier_opt_flag)
                                        {
                                            real_t outlier_overflow_upperbound;
                                            vertex_id_t outlier_num_upperbound;
                                            outlier_upperbound_func(p->data, current_v, adj, outlier_overflow_upperbound, outlier_num_upperbound);
                                            if (outlier_overflow_upperbound > 0 && outlier_num_upperbound > 0)
                                            {
                                                real_t randval = randgen[worker_id].gen_float(outlier_overflow_upperbound * outlier_num_upperbound + regular_area[current_v]) - regular_area[current_v];
                                                if (randval > 0)
                                                {
                                                    //fall in appendix area
                                                    vertex_id_t outlier_idx = randval / outlier_overflow_upperbound;
                                                    if (outlier_idx + 1 > outlier_num_upperbound)
                                                    {
                                                        //in case of round-off error
                                                        outlier_idx = outlier_num_upperbound - 1;
                                                    }
                                                    candidate = outlier_search_func(p->data, current_v, adj, outlier_idx);
                                                    if (candidate != nullptr)
                                                    {
                                                        real_t outlier_static_comp =  static_comp_func ? static_comp_func(current_v, candidate) : 1.0;
                                                        randval = (randval - outlier_idx * outlier_overflow_upperbound) / outlier_static_comp + dcomp_upperbound[current_v];
                                                        #ifdef UNIT_TEST
                                                        assert(outlier_overflow_upperbound / outlier_static_comp + dcomp_upperbound[current_v] + 1e-6 >= dynamic_comp_func(walker, current_v, candidate));
                                                        #endif
                                                        if (randval <= dynamic_comp_func(walker, current_v, candidate))
                                                        {
                                                            ac_edge = candidate;
                                                        }
                                                    }
                                                    continue;
                                                }
                                            }
                                        }

                                        if (alias_tables == nullptr)
                                        {
                                            candidate = adj->begin + gen->gen(degree);
                                        } else
                                        {
                                            AliasBucket<edge_data_t>* bucket = alias_tables->index[current_v] + gen->gen(degree);
                                            if (bucket->q_ptr == nullptr || gen->gen_float(1.0) < bucket->p)
                                            {
                                                candidate = bucket->p_ptr;
                                            } else
                                            {
                                                candidate = bucket->q_ptr;
                                            }
                                        }
                                        if (dynamic_comp_func)
                                        {
                                            real_t dart_height = gen->gen_float(dcomp_upperbound[current_v]);
                                            if ((dcomp_lowerbound != nullptr && dart_height <= dcomp_lowerbound[current_v]) || dart_height <= dynamic_comp_func(walker, current_v, candidate))
                                            {
                                                ac_edge = candidate;
                                            }
                                        } else
                                        {
                                            ac_edge = candidate;
                                        }
                                    }

                                    walker.step ++;
                                    if (walker_update_state_func != nullptr)
                                    {
                                        walker_update_state_func(walker, current_v, ac_edge);
                                    }

                                    if (output_flag)
                                    {
                                        pc->add_footprint(Footprint(walker.id, ac_edge->neighbour, walker.step), worker_id);
                                    }
                                    #ifdef COLLECT_WALK_SEQUENCE
                                    footprints[worker_id].push_back(Footprint(walker.id, ac_edge->neighbour, walker.step));
                                    #endif
                                    if (this->is_local_vertex(ac_edge->neighbour))
                                    {
                                        current_v = ac_edge->neighbour;
                                    } else
                                    {
                                        this->emit(ac_edge->neighbour, walker, worker_id);
                                        break;
                                    }
                                }
                            }
                            this->notify_progress(begin - data_begin, end - data_begin, data_amount, active_walker_num >= PHASED_EXECTION_THRESHOLD * this->partition_num);
                        }
                    }
                    local_walker_num = 0;
                },
                [&](Message<Walker<walker_data_t> > *begin, Message<Walker<walker_data_t> > *end)
                {
                    local_walker_num = end - begin;
                    std::swap(local_walkers, local_walkers_bak);
                },
                local_walkers_bak,
                active_walker_num >= PHASED_EXECTION_THRESHOLD * this->partition_num
            );
        }
    }

    template<typename query_data_t, typename response_data_t>
    void internal_walk_epoch(
        InternalWalkData<query_data_t, response_data_t> *walk_data,
        WalkerConfig<edge_data_t, walker_data_t> *walker_config,
        SecondOrderTransitionConfig<edge_data_t, walker_data_t, query_data_t, response_data_t> *transition_config
    )
    {
        typedef Walker<walker_data_t> walker_t;
        typedef Message<walker_t> walker_msg_t;
        typedef stateQuery<query_data_t> query_t;
        typedef stateResponse<response_data_t> response_t;

        auto extension_comp_func = transition_config->extension_comp_func;
        auto static_comp_func = transition_config->static_comp_func;
        auto post_query_func = transition_config->post_query_func;
        auto respond_query_func = transition_config->respond_query_func;
        auto dynamic_comp_func = transition_config->dynamic_comp_func;
        auto dcomp_upperbound_func = transition_config->dcomp_upperbound_func;
        auto dcomp_lowerbound_func = transition_config->dcomp_lowerbound_func;
        auto outlier_upperbound_func = transition_config->outlier_upperbound_func;
        auto outlier_search_func = transition_config->outlier_search_func;

        auto walker_update_state_func = walker_config->walker_update_state_func;

        auto* alias_tables = walk_data->alias_tables;
        auto* dcomp_lowerbound = walk_data->dcomp_lowerbound;
        auto* dcomp_upperbound = walk_data->dcomp_upperbound;
        auto outlier_opt_flag = walk_data->outlier_opt_flag;
        auto* regular_area = walk_data->regular_area;

        auto* local_walkers = walk_data->local_walkers;
        auto* local_walkers_bak = walk_data->local_walkers_bak;
        auto local_walker_num = walk_data->local_walker_num;

        auto output_flag = walk_data->collect_path_flag;
        auto* pc = walk_data->pc;

        auto* remote_response_cache = walk_data->remote_response_cache;
        auto* remote_fetch_candidate = walk_data->remote_fetch_candidate;
        auto* cached_request = walk_data->cached_request;
        auto* cached_request_end = walk_data->cached_request;

        walker_id_t active_walker_num = walk_data->active_walker_num;    
        int super_step = 0;
        while (active_walker_num != 0)
        {
            #ifndef UNIT_TEST
            if (this->local_partition_id == 0)
            {
                printf("step(%d), active(%u), time(%.3lf)\n", super_step++, active_walker_num, walk_data->timer.duration());
            }
            #endif
            bool use_parallel = (active_walker_num >= OMP_PARALLEL_THRESHOLD);
            this->template distributed_execute<query_t>(
                [&] (void) {
                    walker_id_t progress = 0;
                    walker_id_t data_amount = local_walker_num;
                    auto *data_begin = local_walkers;
                    const walker_id_t work_step_length = PARALLEL_CHUNK_SIZE;
                    #pragma omp parallel if (use_parallel)
                    {
                        int worker_id = omp_get_thread_num();
                        auto *gen = &randgen[worker_id];
                        vertex_id_t next_workload;
                        while((next_workload =  __sync_fetch_and_add(&progress, work_step_length)) < data_amount)
                        {
                            auto *begin = data_begin + next_workload;
                            auto *end = data_begin + std::min(next_workload + work_step_length, data_amount);
                            for (auto *p = begin; p != end; p++)
                            {
                                walker_id_t walker_idx = p - local_walkers;
                                auto *cd = &remote_fetch_candidate[walker_idx];
                                bool local_cal = true;
                                while (local_cal)
                                {
                                    cd->candidate = nullptr;
                                    cd->accepted = false;
                                    vertex_id_t current_v = p->dst_vertex_id;

                                    bool terminate_flag = false;
                                    if (this->vertex_out_degree[current_v] == 0)
                                    {
                                        terminate_flag = true;
                                    } else
                                    {
                                        real_t extension_comp = extension_comp_func(p->data, current_v);
                                        // continue to walk at a probability of extension_comp
                                        if (extension_comp == 0.0 || (extension_comp < 1.0 && gen->gen_float(1.0) > extension_comp))
                                        {
                                            terminate_flag = true;
                                        }
                                    }
                                    if (!terminate_flag)
                                    {
                                        while (local_cal)
                                        {
                                            AdjList<edge_data_t> *adj = this->csr->adj_lists + current_v;
                                            vertex_id_t degree = this->vertex_out_degree[current_v];
                                            bool in_appendix_area = false;
                                            if (outlier_opt_flag)
                                            {
                                                real_t outlier_overflow_upperbound;
                                                vertex_id_t outlier_num_upperbound;
                                                outlier_upperbound_func(p->data, current_v, adj, outlier_overflow_upperbound, outlier_num_upperbound);
                                                if (outlier_overflow_upperbound > 0 && outlier_num_upperbound > 0)
                                                {
                                                    real_t randval = randgen[worker_id].gen_float(outlier_overflow_upperbound * outlier_num_upperbound + regular_area[current_v]) - regular_area[current_v];
                                                    if (randval > 0)
                                                    {
                                                        in_appendix_area = true;
                                                        vertex_id_t outlier_idx = randval / outlier_overflow_upperbound;
                                                        if (outlier_idx + 1 > outlier_num_upperbound)
                                                        {
                                                            //in case of round-off error
                                                            outlier_idx = outlier_num_upperbound - 1;
                                                        }
                                                        cd->candidate = outlier_search_func(p->data, current_v, adj, outlier_idx);
                                                        if (cd->candidate != nullptr)
                                                        {
                                                            real_t outlier_static_comp =  static_comp_func ? static_comp_func(current_v, cd->candidate) : 1.0;
                                                            cd->randval = (randval - outlier_idx * outlier_overflow_upperbound) / outlier_static_comp + dcomp_upperbound[current_v];
                                                            post_query_func(p->data, walker_idx, current_v, cd->candidate);
                                                            local_cal = false;
                                                        }
                                                    }
                                                }
                                            }
                                            if (!in_appendix_area)
                                            {
                                                if (alias_tables == nullptr)
                                                {
                                                    cd->candidate = adj->begin + randgen[worker_id].gen(degree);
                                                } else
                                                {
                                                    AliasBucket<edge_data_t>* bucket = alias_tables->index[current_v] + randgen[worker_id].gen(degree);
                                                    if (bucket->q_ptr == nullptr || randgen[worker_id].gen_float(1.0) < bucket->p)
                                                    {
                                                        cd->candidate = bucket->p_ptr;
                                                    } else
                                                    {
                                                        cd->candidate = bucket->q_ptr;
                                                    }
                                                }
                                                cd->randval = randgen[worker_id].gen_float(dcomp_upperbound[current_v]);
                                                if (dcomp_lowerbound != nullptr && cd->randval <= dcomp_lowerbound[current_v])
                                                {
                                                    cd->accepted = true;
                                                } else
                                                {
                                                    post_query_func(p->data, walker_idx, current_v, cd->candidate);
                                                    local_cal = false;
                                                }
                                            }

                                            if (cd->accepted == true)
                                            {
                                                if (this->is_local_vertex(cd->candidate->neighbour))
                                                {
                                                    p->data.step ++;
                                                    walker_update_state_func(p->data, current_v, cd->candidate);
                                                    p->dst_vertex_id = cd->candidate->neighbour;

                                                    if (output_flag)
                                                    {
                                                        pc->add_footprint(Footprint(p->data.id, cd->candidate->neighbour, p->data.step), worker_id);
                                                    }
                                                    #ifdef COLLECT_WALK_SEQUENCE
                                                    footprints[worker_id].push_back(Footprint(p->data.id, cd->candidate->neighbour, p->data.step));
                                                    #endif
                                                } else
                                                {
                                                    local_cal = false;
                                                }
                                                //break inner loop to play next step
                                                break;
                                            }
                                        }
                                    } else
                                    {
                                        local_cal = false;
                                    }
                                }
                            }
                            this->notify_progress(begin - data_begin, end - data_begin, data_amount, active_walker_num >= PHASED_EXECTION_THRESHOLD * this->partition_num);
                        }
                    }
                },
                [&] (Message<query_t> *msg_begin, Message<query_t> *msg_end)
                {
                    cached_request_end = msg_end;
                },
                cached_request,
                active_walker_num >= PHASED_EXECTION_THRESHOLD * this->partition_num
            );

            this->template distributed_execute<response_t>(
                [&] (void) {
                    walker_id_t progress = 0;
                    walker_id_t data_amount = cached_request_end - cached_request;
                    auto *data_begin = cached_request;
                    const walker_id_t work_step_length = PARALLEL_CHUNK_SIZE;
                    #pragma omp parallel if (use_parallel)
                    {
                        int worker_id = omp_get_thread_num();
                        vertex_id_t next_workload;
                        while((next_workload =  __sync_fetch_and_add(&progress, work_step_length)) < data_amount)
                        {
                            auto *begin = data_begin + next_workload;
                            auto *end = data_begin + std::min(next_workload + work_step_length, data_amount);
                            for (auto *p = begin; p != end; p++)
                            {
                                vertex_id_t current_v = p->dst_vertex_id;
                                auto *adj_list = this->csr->adj_lists + current_v;
                                respond_query_func(current_v, p->data, adj_list);
                            }
                            this->notify_progress(begin - data_begin, end - data_begin, data_amount, active_walker_num >= PHASED_EXECTION_THRESHOLD * this->partition_num);
                        }
                    }
                },
                [&] (Message<response_t> *msg_begin, Message<response_t> *msg_end)
                {
                    walker_id_t progress = 0;
                    walker_id_t data_amount = msg_end - msg_begin;
                    auto *data_begin = msg_begin;
                    const walker_id_t work_step_length = PARALLEL_CHUNK_SIZE;
                    #pragma omp parallel if (use_parallel)
                    {
                        int worker_id = omp_get_thread_num();
                        vertex_id_t next_workload;
                        while((next_workload =  __sync_fetch_and_add(&progress, work_step_length)) < data_amount)
                        {
                            auto *begin = data_begin + next_workload;
                            auto *end = data_begin + std::min(next_workload + work_step_length, data_amount);
                            for (auto *p = begin; p != end; p++)
                            {
                                remote_response_cache[p->data.walker_idx] = p->data;
                            }
                        }
                    }
                },
                nullptr,
                active_walker_num >= PHASED_EXECTION_THRESHOLD * this->partition_num
            );

            active_walker_num = this->template distributed_execute<walker_t>(
                [&] (void) {
                    walker_id_t progress = 0;
                    walker_id_t data_amount = local_walker_num;
                    auto *data_begin = local_walkers;
                    const walker_id_t work_step_length = PARALLEL_CHUNK_SIZE;
                    #pragma omp parallel if (use_parallel)
                    {
                        int worker_id = omp_get_thread_num();
                        vertex_id_t next_workload;
                        while((next_workload =  __sync_fetch_and_add(&progress, work_step_length)) < data_amount)
                        {
                            walker_msg_t *begin = data_begin + next_workload;
                            walker_msg_t *end = data_begin + std::min(next_workload + work_step_length, data_amount);
                            for (walker_msg_t *p = begin; p != end; p++)
                            {
                                walker_t& walker = p->data;
                                vertex_id_t current_v = p->dst_vertex_id;
                                walker_id_t walker_idx = p - local_walkers;
                                auto *cd = &remote_fetch_candidate[walker_idx];
                                if (cd->candidate != nullptr)
                                {
                                    if (cd->accepted || cd->randval <= dynamic_comp_func(walker, remote_response_cache[walker_idx], current_v, cd->candidate))
                                    {
                                        walker.step ++;
                                        walker_update_state_func(walker, current_v, cd->candidate);
                                        this->emit(cd->candidate->neighbour, walker, worker_id);

                                        if (output_flag)
                                        {
                                            pc->add_footprint(Footprint(walker.id, cd->candidate->neighbour, walker.step), worker_id);
                                        }
                                        #ifdef COLLECT_WALK_SEQUENCE
                                        footprints[worker_id].push_back(Footprint(walker.id, cd->candidate->neighbour, walker.step));
                                        #endif
                                    } else
                                    {
                                        this->emit(current_v, walker, worker_id);
                                    }
                                }
                            }
                            this->notify_progress(begin - data_begin, end - data_begin, data_amount, active_walker_num >= PHASED_EXECTION_THRESHOLD * this->partition_num);
                        }
                    }
                    local_walker_num = 0;
                },
                [&](Message<Walker<walker_data_t> > *begin, Message<Walker<walker_data_t> > *end)
                {
                    local_walker_num = end - begin;
                    std::swap(local_walkers, local_walkers_bak);
                },
                local_walkers_bak,
                active_walker_num >= PHASED_EXECTION_THRESHOLD * this->partition_num
            );
        }
    }

#ifdef COLLECT_WALK_SEQUENCE
    void collect_walk_sequence(std::vector<std::vector<vertex_id_t> > &sequence, walker_id_t walker_num)
    {
        std::thread send_thread([&]() {
            for (int t_i = 0; t_i < this->worker_num; t_i++)
            {
                MPI_Send(footprints[t_i].data(), footprints[t_i].size() * sizeof(Footprint), get_mpi_data_type<char>(), 0, Tag_Msg, MPI_COMM_WORLD);
                footprints[t_i].clear();
            }
        });
        if (get_mpi_rank() == 0)
        {
            std::vector<step_t> walk_length(walker_num, 0);
            std::vector<Footprint> glb_footprints;
            for (partition_id_t p_i = 0; p_i < get_mpi_size(); p_i++)
            {
                for (int t_i = 0; t_i < this->worker_num; t_i++)
                {
                    int recv_size = 0;
                    MPI_Status recv_status;
                    MPI_Probe(p_i, Tag_Msg, MPI_COMM_WORLD, &recv_status);
                    MPI_Get_count(&recv_status, get_mpi_data_type<char>(), &recv_size);
                    int recv_n = recv_size / sizeof(Footprint);
                    std::vector<Footprint> recv_data(recv_n);
                    MPI_Recv(recv_data.data(), recv_size, get_mpi_data_type<char>(), p_i, Tag_Msg, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    for (int f_i = 0; f_i < recv_n; f_i++)
                    {
                        auto fp = recv_data[f_i];
                        glb_footprints.push_back(fp);
                        walk_length[fp.walker] = std::max(walk_length[fp.walker], fp.step);
                    }
                }
            }

            sequence.resize(walker_num);
            for (walker_id_t w_i = 0; w_i < walker_num; w_i++)
            {
                sequence[w_i].resize(walk_length[w_i] + 1);
            }
            for (auto fp : glb_footprints)
            {
                sequence[fp.walker][fp.step] = fp.vertex;
            }
        }
        send_thread.join();
    }
#endif
#ifdef COLLECT_WALKER_INIT_STATE
    void collect_walker_init_state(std::vector<Walker<walker_data_t> > &output)
    {
        if (get_mpi_rank() != 0)
        {
            MPI_Send(walker_init_state.data(), walker_init_state.size() * sizeof(Walker<walker_data_t>), get_mpi_data_type<char>(), 0, Tag_Msg, MPI_COMM_WORLD);
        } else
        {
            output = walker_init_state;
            for (partition_id_t p_i = 1; p_i < get_mpi_size(); p_i++)
            {
                int recv_size = 0;
                MPI_Status recv_status;
                MPI_Probe(p_i, Tag_Msg, MPI_COMM_WORLD, &recv_status);
                MPI_Get_count(&recv_status, get_mpi_data_type<char>(), &recv_size);
                int recv_n = recv_size / sizeof(Walker<walker_data_t>);
                std::vector<Walker<walker_data_t> > recv_data(recv_n);
                MPI_Recv(recv_data.data(), recv_size, get_mpi_data_type<char>(), p_i, Tag_Msg, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (auto x : recv_data)
                {
                    output.push_back(x);
                }
            }
            std::sort(output.begin(), output.end(), [](const Walker<walker_data_t> a, const Walker<walker_data_t> b) { return a.id < b.id; });
        }
    }
#endif
};
