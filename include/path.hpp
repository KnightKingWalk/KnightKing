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

#include <string.h>
#include <vector>
#include <mutex>
#include <thread>

#include "type.hpp"
#include "util.hpp"
#include "constants.hpp"
#include "mpi_helper.hpp"

struct Footprint
{
    walker_id_t walker;
    vertex_id_t vertex;
    step_t step;
    Footprint () {}
    Footprint (walker_id_t _walker, vertex_id_t _vertex, step_t _step) : walker(_walker), vertex(_vertex), step(_step) {}
};

struct PathSet
{
    int seg_num;
    vertex_id_t **path_set;
    walker_id_t **walker_id;
    vertex_id_t ***path_begin;
    vertex_id_t ***path_end;
    step_t **path_length;
    walker_id_t *path_num;
    PathSet ()
    {
        seg_num = 0;
        path_set = nullptr;
        walker_id = nullptr;
        path_begin = nullptr;
        path_end = nullptr;
        path_length = nullptr;
        path_num = nullptr;
    }
    ~PathSet ()
    {
        if (seg_num != 0)
        {
            for (int s_i = 0; s_i < seg_num; s_i++)
            {
                delete []path_set[s_i];
                delete []walker_id[s_i];
                delete []path_begin[s_i];
                delete []path_end[s_i];
                delete []path_length[s_i];
            }
            delete []path_set;
            delete []walker_id;
            delete []path_begin;
            delete []path_end;
            delete []path_length;
            delete []path_num;
        }
    }
    void dump(const char* output_path, const char* fopen_mode, bool with_head_info)
    {
        Timer timer;
        FILE* f = fopen(output_path, fopen_mode);
        assert(f != NULL);
        for (int wo_i = 0; wo_i < seg_num; wo_i++)
        {
            for (walker_id_t wa_i = 0; wa_i < path_num[wo_i]; wa_i++)
            {
                if (with_head_info)
                {
                    fprintf(f, "%u %u", walker_id[wo_i][wa_i], path_length[wo_i][wa_i]);
                }
                for (step_t p_i = 0; p_i < path_length[wo_i][wa_i]; p_i++)
                {
                    fprintf(f, " %u", *(path_begin[wo_i][wa_i] + p_i));
                }
                fprintf(f, "\n");
            }
        }
        fclose(f);
#ifndef UNIT_TEST
        printf("finish write path data in %lf seconds \n", timer.duration());
#endif
    }
};

class PathCollector
{
    int worker_num;
    std::vector<MessageBuffer*> node_local_fp;
    std::mutex node_local_fp_lock;
    MessageBuffer** thread_local_fp;
public:
    PathCollector(int worker_num_param)
    {
        this->worker_num = worker_num_param;
        thread_local_fp = new MessageBuffer*[worker_num];
        for (int w_i = 0; w_i < worker_num; w_i++)
        {
            thread_local_fp[w_i] = nullptr;
        }
    }
    ~ PathCollector()
    {
        for (int w_i = 0; w_i < worker_num; w_i++)
        {
            if (thread_local_fp[w_i] != nullptr)
            {
                delete thread_local_fp[w_i];
            }
        }
        delete []thread_local_fp;
        for (auto buf : node_local_fp)
        {
            if (buf != nullptr)
            {
                delete buf;
            }
        }
    }

    void add_footprint(Footprint ft, int worker)
    {
        MessageBuffer* &fp_buffer = thread_local_fp[worker];
        if (fp_buffer == nullptr)
        {
            fp_buffer = new MessageBuffer(sizeof(Footprint) * FOOT_PRINT_CHUNK_SIZE);
        }
        fp_buffer->write(&ft);
        if (fp_buffer->count == FOOT_PRINT_CHUNK_SIZE)
        {
            node_local_fp_lock.lock();
            node_local_fp.push_back(fp_buffer);
            fp_buffer = nullptr;
            node_local_fp_lock.unlock();
        }
    }

    PathSet* assemble_path(walker_id_t walker_begin)
    {
        Timer timer;
        //flush thread local fp
        for (int w_i = 0; w_i < worker_num; w_i++)
        {
            if (thread_local_fp[w_i] != nullptr)
            {
                node_local_fp.push_back(thread_local_fp[w_i]);
                thread_local_fp[w_i] = nullptr;
            }
        }

        //shuffle the footprints
        partition_id_t partition_num = get_mpi_size();
        partition_id_t partition_id = get_mpi_rank();
        auto get_walker_partition_id = [&] (walker_id_t walker)
        {
            return walker % partition_num;
        };
        size_t send_fp_num[partition_num];
        std::fill(send_fp_num, send_fp_num + partition_num, 0);
        size_t progress = 0;
#pragma omp parallel
        {
            int worker_id = omp_get_thread_num();
            size_t next_workload;
            size_t local_counter[partition_num];
            std::fill(local_counter, local_counter + partition_num, 0);
            while ((next_workload =  __sync_fetch_and_add(&progress, 1)) < node_local_fp.size())
            {
                MessageBuffer* buf = node_local_fp[next_workload];
                Footprint* begin = (Footprint*)buf->data; 
                Footprint* end = begin + buf->count;
                for (Footprint *fp = begin; fp < end; fp++)
                {
                    local_counter[get_walker_partition_id(fp->walker)]++;
                }
            }
            for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
            {
                __sync_fetch_and_add(&send_fp_num[p_i], local_counter[p_i]);
            }
        }
        size_t tot_send_fp_num = 0;
        for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
        {
            tot_send_fp_num += send_fp_num[p_i];
        }
        Footprint* send_buf = new Footprint[tot_send_fp_num];
        size_t partition_progress[partition_num];
        partition_progress[0] = 0;
        for (partition_id_t p_i = 1; p_i < partition_num; p_i++)
        {
            partition_progress[p_i] = partition_progress[p_i - 1] + send_fp_num[p_i - 1];
        }
        size_t send_fp_pos[partition_num];
        for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
        {
            send_fp_pos[p_i] = partition_progress[p_i];
        }
        progress = 0;
#pragma omp parallel
        {
            int worker_id = omp_get_thread_num();
            size_t next_workload;
            const size_t LOCAL_BUF_SIZE = THREAD_LOCAL_BUF_CAPACITY;
            MessageBuffer* local_buf[partition_num];
            for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
            {
                local_buf[p_i] = new MessageBuffer(LOCAL_BUF_SIZE * sizeof(Footprint));
            }
            auto flush_local_buf = [&] (partition_id_t p_i)
            {
                if (local_buf[p_i]->count != 0)
                {
                    size_t start_pos = __sync_fetch_and_add(&partition_progress[p_i], (size_t)local_buf[p_i]->count);
                    memcpy(send_buf + start_pos, local_buf[p_i]->data, local_buf[p_i]->count * sizeof(Footprint));
                    local_buf[p_i]->count = 0;
                }
            };
            while ((next_workload =  __sync_fetch_and_add(&progress, 1)) < node_local_fp.size())
            {
                MessageBuffer* &buf = node_local_fp[next_workload];
                Footprint* begin = (Footprint*)buf->data; 
                Footprint* end = begin + buf->count;
                for (Footprint* fp = begin; fp < end; fp++)
                {
                    partition_id_t p_i = get_walker_partition_id(fp->walker);
                    local_buf[p_i]->write(fp);
#ifdef UNIT_TEST
                    local_buf[p_i]->self_check<Footprint>();
#endif
                    if (local_buf[p_i]->count == LOCAL_BUF_SIZE)
                    {
                        flush_local_buf(p_i);
                    }
                }
                delete buf;
                buf = nullptr;
            }
            for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
            {
                flush_local_buf(p_i);
                delete local_buf[p_i];
            }
        }
        node_local_fp.clear();
        size_t glb_send_fp_num[partition_num];
        MPI_Allreduce(send_fp_num, glb_send_fp_num, partition_num, get_mpi_data_type<size_t>(), MPI_SUM, MPI_COMM_WORLD);
        size_t total_recv_fp_num = glb_send_fp_num[partition_id];
		size_t recv_fp_pos = 0;
        Footprint* recv_buf = new Footprint[total_recv_fp_num];
        std::thread send_thread([&](){
            for (partition_id_t step = 0; step < partition_num; step++)
            {
                partition_id_t dst = (partition_id + step) % partition_num;
                size_t tot_send_sz = send_fp_num[dst] * sizeof(Footprint);
#ifdef UNIT_TEST
                const int max_single_send_sz = (1 << 12) / sizeof(Footprint) * sizeof(Footprint);
#else
                const int max_single_send_sz = (1 << 28) / sizeof(Footprint) * sizeof(Footprint);
#endif
                void* send_data = send_buf + send_fp_pos[dst];
                while (true)
                {
                    MPI_Send(&tot_send_sz, 1, get_mpi_data_type<size_t>(), dst, 0, MPI_COMM_WORLD);
                    if (tot_send_sz == 0)
                    {
                        break;
                    }
                    size_t send_sz = std::min((size_t)max_single_send_sz, tot_send_sz);
                    tot_send_sz -= send_sz;
                    MPI_Send(send_data, send_sz, get_mpi_data_type<char>(), dst, 0, MPI_COMM_WORLD);
                    send_data = (char*)send_data + send_sz;
                }
            }
        });
        std::thread recv_thread([&](){
            for (partition_id_t step = 0; step < partition_num; step++)
            {
                partition_id_t src = (partition_id + partition_num - step) % partition_num;
                while (true)
                {
                    size_t remained_sz;
                    MPI_Recv(&remained_sz, 1, get_mpi_data_type<size_t>(), src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    if (remained_sz == 0)
                    {
                        break;
                    }
                    MPI_Status recv_status;
                    MPI_Probe(src, 0, MPI_COMM_WORLD, &recv_status);
                    int sz;
                    MPI_Get_count(&recv_status, get_mpi_data_type<char>(), &sz);

                    MPI_Recv(recv_buf + recv_fp_pos, sz, get_mpi_data_type<char>(), src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    recv_fp_pos += sz / sizeof(Footprint);
                }
            }
        });
        send_thread.join();
        recv_thread.join();
        delete []send_buf;

        //shuffle task locally
        const size_t large_step_size = PARALLEL_CHUNK_SIZE * worker_num;
        const size_t task_num = (total_recv_fp_num + large_step_size - 1) / large_step_size;
        Footprint** task_begin[worker_num];
        Footprint** task_end[worker_num];
        for (int w_i = 0; w_i < worker_num; w_i++)
        {
            task_begin[w_i] = new Footprint*[task_num];
            task_end[w_i] = new Footprint*[task_num];
        }
        auto get_task_partition = [&] (walker_id_t walker)
        {
            return (walker - walker_begin) / partition_num % worker_num;
        };
        progress = 0;
#pragma omp parallel
        {
            int worker_id = omp_get_thread_num();
            size_t next_workload;
            size_t task_counter[worker_num];
            Footprint* task_pos[worker_num];
            Footprint* local_fp = new Footprint[large_step_size];
            while ((next_workload =  __sync_fetch_and_add(&progress, 1)) < task_num)
            {
                Footprint* begin = recv_buf + next_workload * large_step_size; 
                Footprint* end = recv_buf + std::min((next_workload + 1) * large_step_size, total_recv_fp_num);
                size_t work_load_size = end - begin;
                Footprint* local_begin = local_fp;
                Footprint* local_end = local_fp + work_load_size;
                memcpy(local_begin, begin, work_load_size * sizeof(Footprint));
                std::fill(task_counter, task_counter + worker_num, 0);
                for (Footprint *fp = local_begin; fp < local_end; fp++)
                {
                    task_counter[get_task_partition(fp->walker)]++;
                }
                Footprint* current_task_pos = begin;
                for (partition_id_t p_i = 0; p_i < worker_num; p_i++)
                {
                    task_begin[p_i][next_workload] = current_task_pos;
                    task_pos[p_i] = current_task_pos;
                    current_task_pos += task_counter[p_i];
                    task_end[p_i][next_workload] = current_task_pos;
                }
                for (Footprint *fp = local_begin; fp < local_end; fp++)
                {
                    partition_id_t task_p = get_task_partition(fp->walker);
                    *task_pos[task_p] = *fp;
                    task_pos[task_p]++;
                }
            }
            delete []local_fp;
        }

        //do tasks
        vertex_id_t* path_set[worker_num];
        vertex_id_t** path_begin[worker_num];
        vertex_id_t** path_end[worker_num];
        walker_id_t* walker_id[worker_num];
        step_t* path_length[worker_num];
        walker_id_t path_num[worker_num];
        auto get_walker_local_idx = [&] (walker_id_t walker)
        {
            return (walker - walker_begin) / partition_num / worker_num;
        };
#pragma omp parallel
        {
            int worker_id = omp_get_thread_num();
            walker_id_t max_walker_id = 0;
            bool has_any_walker = false;
            size_t step_num = 0;
            for (size_t t_i = 0; t_i < task_num; t_i++)
            {
                Footprint* begin = task_begin[worker_id][t_i];
                Footprint* end = task_end[worker_id][t_i];
                for (Footprint* fp = begin; fp < end; fp++)
                {
                    has_any_walker = true;
                    max_walker_id = std::max(max_walker_id, fp->walker);
                    step_num++;
                }
            }
            walker_id_t max_walker_idx = get_walker_local_idx(max_walker_id);
            walker_id_t thread_walker_num = has_any_walker ? max_walker_idx + 1 : 0;
            path_set[worker_id] = new vertex_id_t[step_num];
            path_begin[worker_id] = new vertex_id_t*[thread_walker_num];
            path_end[worker_id] = new vertex_id_t*[thread_walker_num];
            walker_id[worker_id] = new walker_id_t[thread_walker_num];
            path_length[worker_id] = new step_t[thread_walker_num];
            path_num[worker_id] = thread_walker_num;
            for (walker_id_t w_i = 0; w_i < thread_walker_num; w_i++)
            {
                path_length[worker_id][w_i] = 0;
                walker_id[worker_id][w_i] = w_i * partition_num * worker_num + partition_num * worker_id + partition_id + walker_begin;
            }
            for (size_t t_i = 0; t_i < task_num; t_i++)
            {
                Footprint* begin = task_begin[worker_id][t_i];
                Footprint* end = task_end[worker_id][t_i];
                for (Footprint* fp = begin; fp < end; fp++)
                {
                    walker_id_t idx = get_walker_local_idx(fp->walker);
                    path_length[worker_id][idx]++;
                }
            }
            size_t step_counter = 0;
            for (walker_id_t w_i = 0; w_i < thread_walker_num; w_i++)
            {
                path_begin[worker_id][w_i] = path_set[worker_id] + step_counter;
                step_counter += path_length[worker_id][w_i];
                path_end[worker_id][w_i] = path_set[worker_id] + step_counter;
            }
            for (size_t t_i = 0; t_i < task_num; t_i++)
            {
                Footprint* begin = task_begin[worker_id][t_i];
                Footprint* end = task_end[worker_id][t_i];
                for (Footprint* fp = begin; fp < end; fp++)
                {
                    walker_id_t idx = get_walker_local_idx(fp->walker);
                    *(path_begin[worker_id][idx] + fp->step) = fp->vertex;
                }
            }
        }
        for (int w_i = 0; w_i < worker_num; w_i++)
        {
            delete []task_begin[w_i];
            delete []task_end[w_i];
        }
        delete []recv_buf;
        PathSet* ps = new PathSet();
        ps->seg_num = worker_num;
        ps->path_set = new vertex_id_t*[worker_num];
        ps->walker_id = new walker_id_t*[worker_num];
        ps->path_begin = new vertex_id_t**[worker_num];
        ps->path_end = new vertex_id_t**[worker_num];
        ps->path_length = new step_t*[worker_num];
        ps->path_num = new walker_id_t[worker_num];
        for (int w_i = 0; w_i < worker_num; w_i++)
        {
            ps->path_set[w_i] = path_set[w_i];
            ps->walker_id[w_i] = walker_id[w_i];
            ps->path_begin[w_i] = path_begin[w_i];
            ps->path_end[w_i] = path_end[w_i];
            ps->path_length[w_i] = path_length[w_i];
            ps->path_num[w_i] = path_num[w_i];
        }
#ifndef UNIT_TEST
        printf("finish assembling in %lf seconds\n", timer.duration());
#endif
        return ps;
    }
};
