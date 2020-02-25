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

#include <stdio.h>
#include <assert.h>
#include <unistd.h>

#include <vector>

#include "type.hpp"

template<typename T>
void read_graph(const char* fname, int partition_id, int partition_num, T* &edge, size_t &e_num)
{
    FILE *f = fopen(fname, "r");
    assert(f != NULL);
    fseek(f, 0, SEEK_END);
    size_t total_size = ftell(f);
    size_t total_e_num = total_size / sizeof(T);
    e_num = total_e_num / partition_num;
    if (partition_id == partition_num -1)
    {
        e_num += total_e_num % partition_num;
    }
    size_t f_offset = sizeof(T) * (total_e_num / partition_num) * partition_id; 

    edge = new T[e_num];
    fseek(f, f_offset, SEEK_SET);
    auto ret = fread(edge, sizeof(T), e_num, f);
    assert(ret == e_num);
    fclose(f);
}

template<typename T>
void write_graph(const char* fname, const T* es, const size_t e_num)
{
    FILE *out_f = fopen(fname, "w");
    assert(out_f != NULL);

    auto ret = fwrite(es, sizeof(T), e_num, out_f);
    assert(ret == e_num);
    fclose(out_f);
}

size_t next_endline_pos(FILE *f)
{
    size_t current_pos = ftell(f);
    while (true)
    {
        char ch;
        auto ret = fread(&ch, 1, 1, f);
        if (ret != 1 || ch == '\n')
        {
            break;
        }
        current_pos++;
    }
    return current_pos;
}

std::vector<size_t> partition_text_file(const char* fname, int partition_num)
{
    std::vector<size_t> partition_end;
    FILE *f = fopen(fname, "r");
    assert(f != NULL);
    fseek(f, 0, SEEK_END);
    size_t total_size = ftell(f);
    for (int p_i = 0; p_i < partition_num; p_i++)
    {
        size_t f_offset = total_size / partition_num * (p_i + 1);
        if (f_offset >= total_size)
        {
            partition_end.push_back(f_offset);
        } else
        {
            fseek(f, f_offset, SEEK_SET);
            partition_end.push_back(next_endline_pos(f));
        }
    }
    fclose(f);
    return partition_end;
}

bool read_edge_txt(FILE* f, Edge<EmptyData>* edge)
{
    return (2 == fscanf(f, "%u %u", &edge->src, &edge->dst));
}

bool read_edge_txt(FILE* f, Edge<real_t>* edge)
{
    return (3 == fscanf(f, "%u %u %f", &edge->src, &edge->dst, &edge->data));
}

template<typename T>
bool read_edge_txt(FILE* f, Edge<T>* edge)
{
    fprintf(stderr, "Edge type doesn't support reading from text\n");
    exit(1);
}

template<typename T>
void read_edgelist(const char* fname, int partition_id, int partition_num, Edge<T>* &edge, size_t &e_num)
{
    std::vector<size_t> partition_end = partition_text_file(fname, partition_num);
    size_t begin = (partition_id == 0 ? 0 : partition_end[partition_id - 1]);
    size_t end = partition_end[partition_id];
    FILE *f = fopen(fname, "r");
    assert(f != NULL);

    Edge<T> temp;
    e_num = 0;
    fseek(f, begin, SEEK_SET);
    while (ftell(f) < end)
    {
        if (read_edge_txt(f, &temp))
        {
            e_num++;
        }
    }
    edge = new Edge<T>[e_num];

    size_t e_i = 0;
    fseek(f, begin, SEEK_SET);
    while (ftell(f) < end)
    {
        if (read_edge_txt(f, &temp))
        {
            edge[e_i] = temp;
            e_i++;
        }
    }

    fclose(f);
}

void print_edge(FILE* f, const Edge<EmptyData>* edge)
{
   fprintf(f, "%u %u\n", edge->src, edge->dst);
}

void print_edge(FILE* f, const Edge<real_t>* edge)
{
   fprintf(f, "%u %u %f\n", edge->src, edge->dst, edge->data);
}

template<typename T>
void print_edge(FILE* f, const Edge<T>* edge)
{
    fprintf(stderr, "Edge type doesn't support writing to text\n");
    exit(1);
}

template<typename T>
void write_edgelist(const char* fname, const Edge<T>* es, const size_t e_num)
{
    FILE *out_f = fopen(fname, "w");
    assert(out_f != NULL);
    for (size_t e_i = 0; e_i < e_num; e_i++)
    {
        print_edge(out_f, es + e_i);
    }
    fclose(out_f);
}
