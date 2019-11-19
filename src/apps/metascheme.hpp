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
#include <vector>

#include "type.hpp"

typedef int scheme_id_t;
typedef int meta_state_t;

struct MetapathState
{
    scheme_id_t scheme_id;
    meta_state_t state;
};

struct WeightedMetaData
{
    real_t weight;
    meta_state_t meta_info;
    meta_state_t get_meta()
    {
        return meta_info;
    }
    real_t get_weight()
    {
        return weight;
    }
};

struct UnweightedMetaData
{
    meta_state_t meta_info;
    meta_state_t get_meta()
    {
        return meta_info;
    }
    real_t get_weight()
    {
        return 1.0;
    }
};

std::vector<std::vector<std::vector<bool> > > read_metapath_schemes(const char* path)
{
    FILE *f = fopen(path, "r");
    assert(f != NULL);
    int scheme_num, state_num;
    assert(2 == fscanf(f, "%d %d", &scheme_num, &state_num));
    std::vector<std::vector<std::vector<bool> > > schemes(scheme_num);
    for (int sc_i = 0; sc_i < scheme_num; sc_i++)
    {
        int scheme_length;
        assert(1 == fscanf(f, "%d", &scheme_length));
        schemes[sc_i].resize(scheme_length); 
        for (int l_i = 0; l_i < scheme_length; l_i++)
        {
            schemes[sc_i][l_i].resize(state_num);
            for (int st_i = 0; st_i < state_num; st_i++)
            {
                int st;
                assert(1 == fscanf(f, "%d", &st));
                schemes[sc_i][l_i][st_i] = (st != 0);
            }
        }
    }
    return schemes;
    fclose(f);
};

void write_metapath_schemes(std::vector<std::vector<std::vector<bool> > > schemes, const char* path)
{
    FILE *f = fopen(path, "w");
    assert(f != NULL);
    int scheme_num = schemes.size();
    int state_num = schemes[0][0].size();
    fprintf(f, "%d %d\n", scheme_num, state_num);
    for (int sc_i = 0; sc_i < scheme_num; sc_i++)
    {
        fprintf(f, "%d\n", (int)schemes[sc_i].size());
        for (auto trans : schemes[sc_i])
        {
            for (auto s : trans)
            {
                fprintf(f, "%d ", (int)s);
            }
            fprintf(f, "\n");
        }
    }
    fclose(f);
}
