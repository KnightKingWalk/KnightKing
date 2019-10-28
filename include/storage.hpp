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
