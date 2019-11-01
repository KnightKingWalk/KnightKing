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

template<typename edge_data_t, typename walker_data_t>
class WalkEngine : public GraphEngine<edge_data_t>
{
    StdRandNumGenerator* randgen;

    //for setting walkers
    walker_id_t walker_num;
    std::function<vertex_id_t (walker_id_t)> walker_init_dist_func;
    std::function<void (Walker<walker_data_t>&, vertex_id_t)> walker_init_state_func;
    std::function<void (Walker<walker_data_t>&, vertex_id_t, AdjUnit<edge_data_t> *)> walker_update_state_func;

    //for outputs
    bool output_flag;
    PathSet path_data;
#ifdef COLLECT_WALK_SEQUENCE
    std::vector<std::vector<Footprint> > footprints;
#endif

    template<typename T>
    T* alloc_walker_array()
    {
        return this->template alloc_array<T>(walker_num);
    }

    template<typename T>
    void dealloc_walker_array(T* arr)
    {
        this->dealloc_array(arr, walker_num);
    }

public:
    WalkEngine()
    {
        walker_num = 0;
        walker_init_dist_func = nullptr;
        walker_init_state_func = nullptr;
        walker_update_state_func = nullptr;

        output_flag = false;

        randgen = new StdRandNumGenerator[this->worker_num];
    }
    ~WalkEngine()
    {
        if (randgen != nullptr)
        {
            delete []randgen;
        }
        _internal_free_path_data(path_data);
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

private:
    void init_walkers(
        Message<Walker<walker_data_t> >* &local_walkers,
        Message<Walker<walker_data_t> >* &local_walkers_bak,
        walker_id_t &local_walker_num
    )
    {
        typedef Walker<walker_data_t> walker_t;
        typedef Message<walker_t> walker_msg_t;

        assert(local_walkers == nullptr);
        assert(local_walkers_bak == nullptr);
        assert(local_walker_num == 0);

        //local_walkers = new Message<Walker<walker_data_t> > [walker_capacity];
        local_walkers = this->template alloc_walker_array<Message<Walker<walker_data_t> > >();
        //local_walkers_bak = new Message<Walker<walker_data_t> > [walker_capacity];
        local_walkers_bak = this->template alloc_walker_array<Message<Walker<walker_data_t> > >();
        local_walker_num = 0;
#ifdef COLLECT_WALK_SEQUENCE
        footprints.resize(this->worker_num);
#endif
        this->set_msg_buffer(walker_num, sizeof(walker_msg_t));
        this->template distributed_execute<walker_t>(
            [&] (void) {
                #pragma omp parallel for
                for (walker_id_t w_i = this->local_partition_id; w_i < walker_num; w_i += this->partition_num)
                {
                    vertex_id_t start_v = walker_init_dist_func(w_i);
                    Walker<walker_data_t> walker;
                    walker.id = w_i;
                    walker.step = 0;
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
    }

    void free_walkers(
        Message<Walker<walker_data_t> >* local_walkers,
        Message<Walker<walker_data_t> >* local_walkers_bak
    )
    {
        assert(local_walkers != nullptr);
        assert(local_walkers_bak != nullptr);
        this->dealloc_walker_array(local_walkers);
        this->dealloc_walker_array(local_walkers_bak);
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

    void init_alias_tables(std::function<real_t(vertex_id_t, AdjUnit<edge_data_t>*)> static_comp_func, AliasTableContainer<edge_data_t> *&alias_tables)
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

    void set_output()
    {
        output_flag = true;
    }

    PathSet get_path_data()
    {
        PathSet temp = path_data;
        path_data = PathSet();
        return temp;
    }

    void dump_path_data(PathSet ps, std::string output_path_root)
    {
        std::string local_output_path = output_path_root + std::string("_") + std::to_string(this->local_partition_id);
        _internal_write_path_data(ps, local_output_path.c_str());
    }

    void free_path_data(PathSet &ps)
    {
        _internal_free_path_data(ps);
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

    void set_walkers(
        walker_id_t walker_num_param,
        std::function<void (Walker<walker_data_t>&, vertex_id_t)> walker_init_state_func_param = nullptr,
        std::function<void (Walker<walker_data_t>&, vertex_id_t, AdjUnit<edge_data_t> *)> walker_update_state_func_param = nullptr,
        std::function<vertex_id_t (walker_id_t)> walker_init_dist_func_param = nullptr
    )
    {
        typedef Walker<walker_data_t> walker_t;
        typedef Message<walker_t> walker_msg_t;

        this->walker_num = walker_num_param;
        this->walker_init_state_func = walker_init_state_func_param;
        this->walker_update_state_func = walker_update_state_func_param;
        if (walker_init_dist_func_param == nullptr)
        {
            walker_init_dist_func_param = get_equal_dist_func();
        }
        this->walker_init_dist_func = walker_init_dist_func_param;
    }

    void random_walk(
        std::function<real_t (Walker<walker_data_t>&, vertex_id_t)> extension_comp_func,
        std::function<real_t(vertex_id_t, AdjUnit<edge_data_t>*)> static_comp_func = nullptr,
        std::function<real_t (Walker<walker_data_t>&, vertex_id_t, AdjUnit<edge_data_t> *)> dynamic_comp_func = nullptr,
        std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_upperbound_func = nullptr,
        std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_lowerbound_func = nullptr
    )
    {
        typedef Walker<walker_data_t> walker_t;
        typedef Message<walker_t> walker_msg_t;

        Timer timer;

        PathCollector *pc = nullptr;
        if (output_flag)
        {
            pc = new PathCollector(this->worker_num);
        }

        walker_msg_t *local_walkers = nullptr;
        walker_msg_t *local_walkers_bak = nullptr;
        walker_id_t local_walker_num = 0;
        init_walkers(local_walkers, local_walkers_bak, local_walker_num);
        if (output_flag)
        {
#pragma omp parallel for
            for (walker_id_t w_i = 0; w_i < local_walker_num; w_i++)
            {
                pc->add_footprint(Footprint(local_walkers[w_i].data.id, local_walkers[w_i].dst_vertex_id, 0), omp_get_thread_num());
            }
        }

        AliasTableContainer<edge_data_t> *alias_tables = nullptr;
        if (static_comp_func != nullptr)
        {
            init_alias_tables(static_comp_func, alias_tables); 
        }

        real_t* dcomp_upperbound = nullptr;
        if (dcomp_lowerbound_func != nullptr)
        {
            init_dcomp_upperbound(dcomp_upperbound_func, dcomp_upperbound);
        }

        real_t* dcomp_lowerbound = nullptr;
        if (dcomp_lowerbound_func != nullptr)
        {
            init_dcomp_lowerbound(dcomp_lowerbound_func, dcomp_lowerbound);
        }

        this->set_msg_buffer(walker_num, sizeof(walker_msg_t));

        walker_id_t active_walker_num = walker_num;    
        int super_step = 0;
        while (active_walker_num != 0)
        {
            #ifndef UNIT_TEST
            if (this->local_partition_id == 0)
            {
                printf("step(%d), active(%u), time(%.3lf)\n", super_step++, active_walker_num, timer.duration());
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

                                    if (walker_update_state_func != nullptr)
                                    {
                                        walker_update_state_func(walker, current_v, ac_edge);
                                    }
                                    walker.step ++;

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

        if (output_flag)
        {
            path_data = pc->assemble_path();
            delete pc;
        }
        free_walkers(local_walkers, local_walkers_bak);
        if (static_comp_func != nullptr)
        {
            free_alias_tables(alias_tables);
        }
        if (dcomp_upperbound_func != nullptr)
        {
            free_dcomp_upperbound(dcomp_upperbound);
        }
        if (dcomp_lowerbound_func != nullptr)
        {
            free_dcomp_lowerbound(dcomp_lowerbound);
        }
    }

    template<typename query_data_t, typename response_data_t>
    void second_order_random_walk (
        std::function<real_t (Walker<walker_data_t>&, vertex_id_t)> extension_comp_func,
        std::function<real_t(vertex_id_t, AdjUnit<edge_data_t>*)> static_comp_func,
        std::function<void (Walker<walker_data_t>&, walker_id_t, vertex_id_t, AdjUnit<edge_data_t> *)> post_query_func,
        std::function<void (vertex_id_t, stateQuery<query_data_t> &, AdjList<edge_data_t>*)> respond_query_func,
        std::function<real_t (Walker<walker_data_t>&, stateResponse<response_data_t> &, vertex_id_t, AdjUnit<edge_data_t> *)> dynamic_comp_func,
        std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_upperbound_func,
        std::function<real_t(vertex_id_t, AdjList<edge_data_t>*)> dcomp_lowerbound_func = nullptr,
        std::function<bool (Walker<walker_data_t>&, vertex_id_t, AdjUnit<edge_data_t> *, real_t&)> try_local_compute_func = nullptr
    )
    {
        typedef Walker<walker_data_t> walker_t;
        typedef Message<walker_t> walker_msg_t;
        typedef stateQuery<query_data_t> query_t;
        typedef stateResponse<response_data_t> response_t;

        Timer timer;

        PathCollector *pc = nullptr;
        if (output_flag)
        {
            pc = new PathCollector(this->worker_num);
        }

        assert(post_query_func != nullptr);
        assert(respond_query_func != nullptr);
        response_t* remote_response_cache = this->template alloc_walker_array<response_t>();
        SecondOrderCandidate<edge_data_t>* remote_fetch_candidate = this->template alloc_walker_array<SecondOrderCandidate<edge_data_t> >();
        Message<query_t>* cached_request = this->template alloc_walker_array<Message<query_t> > ();
        Message<query_t> *cached_request_end = cached_request;

        walker_msg_t *local_walkers = nullptr;
        walker_msg_t *local_walkers_bak = nullptr;
        walker_id_t local_walker_num = 0;
        init_walkers(local_walkers, local_walkers_bak, local_walker_num);

        if (output_flag)
        {
#pragma omp parallel for
            for (walker_id_t w_i = 0; w_i < local_walker_num; w_i++)
            {
                pc->add_footprint(Footprint(local_walkers[w_i].data.id, local_walkers[w_i].dst_vertex_id, 0), omp_get_thread_num());
            }
        }

        AliasTableContainer<edge_data_t> *alias_tables = nullptr;
        if (static_comp_func != nullptr)
        {
            init_alias_tables(static_comp_func, alias_tables); 
        }

        real_t* dcomp_upperbound = nullptr;
        if (dcomp_lowerbound_func != nullptr)
        {
            init_dcomp_upperbound(dcomp_upperbound_func, dcomp_upperbound);
        }

        real_t* dcomp_lowerbound = nullptr;
        if (dcomp_lowerbound_func != nullptr)
        {
            init_dcomp_lowerbound(dcomp_lowerbound_func, dcomp_lowerbound);
        }

        size_t max_msg_size = std::max(sizeof(walker_msg_t), std::max(sizeof(Message<query_t>), sizeof(Message<response_t>)));
        this->set_msg_buffer(walker_num, max_msg_size);

        walker_id_t active_walker_num = walker_num;    
        int super_step = 0;
        while (active_walker_num != 0)
        {
            #ifndef UNIT_TEST
            if (this->local_partition_id == 0)
            {
                printf("step(%d), active(%u), time(%.3lf)\n", super_step++, active_walker_num, timer.duration());
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
                                            AdjUnit<edge_data_t> *candidate;
                                            if (alias_tables == nullptr)
                                            {
                                                candidate = adj->begin + randgen[worker_id].gen(degree);
                                            } else
                                            {
                                                AliasBucket<edge_data_t>* bucket = alias_tables->index[current_v] + randgen[worker_id].gen(degree);
                                                if (bucket->q_ptr == nullptr || randgen[worker_id].gen_float(1.0) < bucket->p)
                                                {
                                                    candidate = bucket->p_ptr;
                                                } else
                                                {
                                                    candidate = bucket->q_ptr;
                                                }
                                            }
                                            cd->randval = randgen[worker_id].gen_float(dcomp_upperbound[current_v]);
                                            cd->candidate = candidate;
                                            if (dcomp_lowerbound != nullptr && cd->randval <= dcomp_lowerbound[current_v])
                                            {
                                                cd->accepted = true;
                                            } else
                                            {
                                                real_t try_prob;
                                                if (try_local_compute_func != nullptr && try_local_compute_func(p->data, current_v, candidate, try_prob))
                                                {
                                                    if (cd->randval <= try_prob)
                                                    {
                                                        cd->accepted = true;
                                                    }
                                                } else
                                                {
                                                    post_query_func(p->data, walker_idx, current_v, candidate);
                                                    local_cal = false;
                                                }
                                            }
                                            if (cd->accepted == true)
                                            {
                                                if (this->is_local_vertex(cd->candidate->neighbour))
                                                {
                                                    walker_update_state_func(p->data, current_v, cd->candidate);
                                                    p->data.step ++;
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
                                        walker_update_state_func(walker, current_v, cd->candidate);
                                        walker.step ++;
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

        if (output_flag)
        {
            path_data = pc->assemble_path();
            delete pc;
        }

        if (remote_response_cache != nullptr)
        {
            this->dealloc_walker_array(remote_response_cache);
        }
        if (remote_fetch_candidate != nullptr)
        {
            this->dealloc_walker_array(remote_fetch_candidate);
        }
        if (cached_request != nullptr)
        {
            this->dealloc_walker_array(cached_request);
        }

        free_walkers(local_walkers, local_walkers_bak);
        if (static_comp_func != nullptr)
        {
            free_alias_tables(alias_tables);
        }
        if (dcomp_upperbound_func != nullptr)
        {
            free_dcomp_upperbound(dcomp_upperbound);
        }
        if (dcomp_lowerbound_func != nullptr)
        {
            free_dcomp_lowerbound(dcomp_lowerbound);
        }
    }

#ifdef COLLECT_WALK_SEQUENCE
    void collect_walk_sequence(std::vector<std::vector<vertex_id_t> > &sequence)
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
};
