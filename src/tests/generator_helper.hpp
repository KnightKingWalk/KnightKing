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

#include <type_traits>
#include <stdlib.h>
#include <vector>
#include <random>

#include "type.hpp"
#include "util.hpp"

enum GraphType
{
    DirectedGraph,
    UndirectedGraph
};

enum EdgeDataType
{
    TypeEmptyData,
    TypeFloat
};

template <typename T>
void gen_rand_edge_data(typename std::enable_if<std::is_same<EmptyData, T>::value, T>::type &t, uint32_t range = 0, RandNumGenerator *gen = nullptr)
{
}

template <typename T>
void gen_rand_edge_data(typename std::enable_if<std::is_integral<T>::value, T>::type &t, uint32_t range = 5, RandNumGenerator *gen = nullptr)
{
    assert(range > 1);
    if (gen == nullptr)
    {
        t = rand() % range + 1;
    } else
    {
        t = gen->gen(range - 1) + 1;
    }
}

template <typename T>
void gen_rand_edge_data(typename std::enable_if<std::is_same<float, T>::value, T>::type &t, uint32_t range = 5, RandNumGenerator *gen = nullptr)
{
    assert(range > 1);
    if (gen == nullptr)
    {
        t = (T) (rand() % range + 1);
    } else
    {
        t = gen->gen_float(range - 1.0) + 1.0;
    }
}

template <typename T>
void gen_rand_edge_data(typename std::enable_if<std::is_class<T>::value && !std::is_same<EmptyData, T>::value, T>::type &t)
{
    printf("[error] Please define edge data function generator\n");
    exit(1);
}
