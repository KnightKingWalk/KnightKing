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

#define L1_CACHE_LINE_SIZE 64

#ifdef UNIT_TEST

#define THREAD_LOCAL_BUF_CAPACITY 16
#define OMP_PARALLEL_THRESHOLD 10
#define PARALLEL_CHUNK_SIZE 4
#define PHASED_EXECTION_THRESHOLD 100
#define DISTRIBUTEDEXECUTIONCTX_PHASENUM 5
#define FOOT_PRINT_CHUNK_SIZE 16

#else

#define THREAD_LOCAL_BUF_CAPACITY 1024
#define OMP_PARALLEL_THRESHOLD 4000
#define PARALLEL_CHUNK_SIZE 128
#define PHASED_EXECTION_THRESHOLD 500000
#define DISTRIBUTEDEXECUTIONCTX_PHASENUM 16
#define FOOT_PRINT_CHUNK_SIZE 65536

#endif
