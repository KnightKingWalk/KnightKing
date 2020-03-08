# KnightKing

**KnightKing** is a general-purpose, distributed graph random walk engine. It provides:

- Unified edge transition probability definition
- Rejection-based, fast and exact edge sampling
- Walker-centric design and programming model
- Common optimizations for different random walk algorithms.

This repository gives the reference implementation of the KnightKing graph random walk engine (SOSP '19 paper [[PDF]](resources/sosp2019_paper.pdf), 16-minute conference talk slides [[PDF]](resources/sosp2019_slides.pdf)). It provides the same APIs as introduced in the paper, with a semi-asynchronous system design.

Contributors: Ke Yang<sup>1 </sup>, Mingxing Zhang<sup>1, 2</sup>, Kang Chen<sup>1</sup>, Xiaosong Ma<sup>3</sup>, Yang Bai<sup>4</sup>, Yong Jiang<sup>1</sup>, Yongwei Wu<sup>1</sup>, Shuke Wang<sup>1 </sup>, Junqi Huang<sup>1 </sup>

<sup>1 </sup> Tsinghua University, <sup>2 </sup>Sangfor, <sup>3 </sup>QCRI, <sup>4 </sup>4 Paradigm.

## Update

In this update, KnightKing has new functionalities:

- KnightKing now is stateless. Users can set the walks via **WalkConfig**, **WalkerConfig**, and **TransitionConfig**. Each time **random_walk** is invoked, a new random walk execution is launched, which is configured by the above 3 configurations. As a result, users now can launch multiple different random walks on the same engine.
- KnightKing supports breaking random walk into multiple epochs. In each epoch only a fraction of walkers are processed. This can save a lot of memory.

But its APIs also change greatly. For its old version, please refer to the v0_1_0 tag. Sorry for the inconvenience.

## Content

- [Quick Start](#Quick-Start)
- [Create Your Own Applications](#Create-Your-Own-Applications)
- [APIs](#APIs)
- [Publication](#Publication)

## Quick Start

This section gives a guide about how to compile KnightKing and how to use the built-in applications.

### Tested Environment

- Ubuntu 14.04
- gcc 4.8.4
- cmake 2.8.12

### Download and Compile

KnightKing uses OpenMPI for inter-node communications, OpenMP for multi-threading, Google Test for unit test, and CMake for compilation. 

To use KnightKing, first download this repository from GitHub:

```
git clone https://github.com/KnightKingWalk/KnightKing.git --recurse-submodules
cd KnightKing
```

Then compile it with CMake:

```
mkdir build && cd build

cmake ..

make
```

The compiled executable files will be installed at the "bin" directory:

```
ls ./bin
```

[Optional] Run unit tests:

```
cmake -DWITH_TESTS=on ..
ctest
```
note that the unit test may take about one hour.

### Run Built-in Applications

Here we take node2vec as an example to show how to run the built-in applications. The usage of other three applications is quite similar.

Use "-h" option to print out node2vec's options

```
~/graph/knightking/build$./bin/node2vec -h
 ./bin/node2vec {OPTIONS}


  OPTIONS:
      -h, --help                        Display this help menu
      -v[vertex]                        vertex number
      -g[graph]                         graph data path
      --make-undirected                 load graph and treat each edge as
                                        undirected edge
      -w[walker]                        walker number
      -o[output]                        [optional] the output path. Omit this
                                        option for pure random walk performance
                                        testing without output.
      -r[rate]                          Set this option will break random walk
                                        into multiple iterations to save memory.
                                        Each iteration has rate% walkers.
      -l[length]                        walk length
      -s[static_comp]                   [weighted | unweighted] a weighted graph
                                        usually indicates a non-trivial static
                                        component.
      -p[p]                             hyperparameter p
      -q[q]                             hyperparameter q
```

A graph random walk application usually takes a graph as input, then setups a group of walkers to wander around the graph. So we need to specify the path of the graph file, the number of vertices, and the number of walkers.

"-v" option specifies how many vertices the graph have. The range of vertex ID is \[0, vertex_num).

"-g" option specifies the path of input graph file.

"--make-undirected" option will treat each edge as undirected edge. For each edge (u, v) in the file, a new edge of (v, u) will be created.

"-w" option specifies how many walkers there are.

"-o" option specifies the prefix of the name of output files. This is an optional paramemter. If this parameter is set, then after random walk, the walking paths will be dumped to the disk.

"-r" option will break random walk into multiple epochs. Each iteration has rate% walkers. Use this option if memory is not enough to hold all walkers and their paths.

"-l" option specifies the walk length. Node2vec is a truncated random walk, which means each walker walks a pre-defined length.

"-s" option specifies whether the graph is a weighted graph.

"-p" and "-q" option specify the hyper-parameters p and q, which controls the walking strategy. Please read the [node2vec paper](https://dl.acm.org/citation.cfm?id=2939754) for more information.

There is a text file containing a sample graph, but node2vec takes a binary file as input. So first we convert the text file to binary format:

```
./bin/gconverter -i ../dataset/karate.txt -o ./karate.data -s weighted
```

Then we can invoke node2vec:

```
mkdir out
./bin/node2vec -g ./karate.data -v 34 -w 5 -s weighted -l 10 -p 2 -q 0.5 -o ./out/walks.txt
```

See the random walk output:

```
~/graph/KnightKing/build$ cat ./out/walks.txt.0
0 5 16 6 4 6 16 6 4 6 16
1 21 0 2 8 30 33 15 32 8 32
2 13 0 10 0 31 0 6 16 5 0
3 12 0 11 0 17 1 2 9 33 30
4 0 1 21 1 17 1 2 8 32 2
```

There are 5 lines, each representing a path for one walker.

Note that you may have an output different from above example. Since this is a random walk, so the output is also random.

### Run in Distributed Environment

First, copy the graph file to the same path of each node, or simply place it to a shared file system. Second, write each node's IP address to a text file (e.g. ./hosts). Then use MPI to run the application. Suppose the graph file is placed at ./karate.data. For OpenMPI:

```
mpiexec -npernode 1 -hostfile ./hosts ./bin/node2vec -g ./karate.data -v 34 -w 34 -s weighted -l 80 -p 2 -q 0.5 -o ./out/walks.txt
```

For MPICH:

```
mpiexec -np 1 -hostfile ./hosts ./bin/node2vec -g ./karate.data -v 34 -w 34 -s weighted -l 80 -p 2 -q 0.5 -o ./out/walks.txt
```

The "-npernode 1" or -"np 1" setting is recommended, which tells MPI to instantiate one instance per node. KnightKing will automatically handle the concurrency within each node. Instantiating more than one instances per node may make the graph more fragmented and thus hinge the performance.

See the random walk output:

```
cat ./out/walks.txt.*
```

## Create Your Own Applications

This section introduces how to use KnightKing's API to write your own apps.

### A Simple Case

This sub-section demonstrates the simplest case: the walkers randomly walk around the graph, without any preference on which edges to follow.

**Step 1:** Create a random walk engine instance and load the graph into memory:

```c++
WalkEngine<real_t, EmptyData> graph;
graph.load_graph(34 /*vertex number*/, "./karate.data" /*graph file path*/);
```

**Step 2:** Set walk configuration:

```c++
WalkConfig walk_conf;
walk_conf.set_output_file("./out/walks.txt");
```

**Step 3:** Set walker configuration:

```c++
WalkerConfig<real_t, EmptyData> walker_conf(34 /*walker number*/);
```

**Step 4**: Define the extension component, which controls the termination condition and sets the probability to continue to walk:

```c++
auto extension_comp = [&] (Walker<EmptyData>& walker, vertex_id_t current_v)
{
    return 0.875; /*the probability to continue the walk*/
};
TransitionConfig<real_t, EmptyData> tr_conf(extension_comp);
```

In this case, at each step, each active walker will terminate at a probability of 0.125.

After the configuration, invoke the function *random_walk*, then KnightKing will begin the random walk process until all walkers terminate:

```c++
graph.random_walk(&walker_conf, &tr_conf, &walk_conf);
```

The complete sample code can be found in [src/examples/simple_walk.cpp](src/examples/simple_walk.cpp).

### Biased Graph Random Walk

Now let's consider a more complicated case: a truncated biased graph random walk.

First, unlike the above case where each walker terminates at a given probability at each step, in truncated walk, all walkers terminates at a given step.

Second, the edges are no longer identical. In real world, it is common that some edges are more important than others, and are more likely to be visited by walkers. In this case, each edge is associated to a real number called edge weight, denoting its un-normalized transition probability. For example, suppose a vertex **v** has four outgoing edges with edge weight of 1, 2, 1, 4, respectively. When a walker walks to **v**, the probability of being selected by the walker to follow is 0.125, 0.25, 0.125, 0.5, respectively.

In this case, we re-define the extension component, which sets the termination probability:

```c++
auto extension_comp = [&] (Walker<EmptyData>& walker, vertex_id_t current_v)
{
   return walker.step >= 10 ? 0.0 : 1.0; /*walk 10 steps then terminate*/
};
```

Then we define the static component, which tells how static properties such as edge weight affect the preference of walkers:

```c++
auto static_comp = [&] (vertex_id_t v, AdjUnit<real_t> *edge)
{
    return edge->data; /*edge->data is a real number denoting edge weight*/
};
TransitionConfig<real_t, EmptyData> tr_conf(extension_comp, static_comp);
```

Then invoke the random walk:

```
graph.random_walk(&walker_conf, &tr_conf, &walk_conf);
```

The complete sample code can be found in [src/examples/biased_walk.cpp](src/examples/biased_walk.cpp).

### Dynamic Random Walk

So far until now, the examples we meet all have fixed edge transition probability. But there are many cases where we need dynamic edge transition probability to create more flexible walking strategy. For example, if we expect the walkers to prefer exploring locally, a simple strategy is to increase the transition probability of the *return edge*, i.e., the edge just traveled.

To do so, first add an additional property to walker definition, so that the walkers can remember where they come from:

```c++
struct WalkState
{
    vertex_id_t last_vertex;
};
```

Then specify the walker properties when create the graph engine:

```c++
WalkEngine<real_t, WalkState> graph;
```

Next, define how to initiate and how to update the state:

```c++
auto init_walker_func = [&] (Walker<WalkState> &walker, vertex_id_t start_vertex)
{
    /*At first, the last vertex is not defined*/
    walker.data.last_vertex = UINT_MAX;
};
auto update_walker_func = [&] (Walker<WalkState> &walker, vertex_id_t current_v, AdjUnit<real_t> *edge)
{
    walker.data.last_vertex = current_v;
};
WalkerConfig<real_t, WalkState> walker_conf(34, init_walker_func, update_walker_func);
```

Then define dynamic component and its upper bound:

```c++
auto dynamic_comp = [&] (Walker<WalkState> &walker, vertex_id_t current_v, AdjUnit<real_t> *edge)
{
    if (walker.step == 0)
    {
        /*No return edge for the first step*/
        return 1.0;
    } else if (edge->neighbour == walker.data.last_vertex)
    {
        /*if return edge, double the un-normalized transition probability*/
        return 2.0;
    } else
    {
        /*if not return edge*/
        return 1.0;
    }
};
auto dynamic_comp_upperbound = [&] (vertex_id_t v_id, AdjList<real_t> *adj_lists)
{
    return 2.0;
};
TransitionConfig<real_t, WalkState> tr_conf(extension_comp, static_comp, dynamic_comp, dynamic_comp_upperbound);
```

Then invoke the random walk:

```c++
graph.random_walk(&walker_conf, &tr_conf, &walk_conf);
```

The complete sample code can be found in [src/examples/dynamic_walk.cpp](src/examples/dynamic_walk.cpp).

## APIs

This section lists KnightKing's APIs and their usage.

### Walkers

#### Walker Definition

In KnightKing, walkers are defined as a template:

```c++
template<typename walker_data_t>
struct Walker
{
public:
    walker_id_t id;
    step_t step;
    walker_data_t data;
};
```

 In default, it has two properties: id and step.

**id**: each walker is assigned to a unique id.

**step**: this property records how many steps the walker has walked.

These two properties are automatically maintained by KnightKing.

Users can define additional properties using the **data** field.

#### Walker Configuration 

The class *WalkerConfig* defines walkers' properties. Its constructor is:

```c++
WalkerConfig(
    walker_id_t walker_num,
    std::function<void (Walker<walker_data_t>&, vertex_id_t)> walker_init_state_func = nullptr,
    std::function<void (Walker<walker_data_t>&, vertex_id_t, AdjUnit<edge_data_t> *)> walker_update_state_func = nullptr,
    std::function<vertex_id_t (walker_id_t)> walker_init_dist_func = nullptr
)
```

**walker_num**: This parameter is an integer which tells how many walkers there are.

**walker_init_state_func**: This parameter is a function that takes a walker and its initial vertex as input, and initialize that walker. Note that users are only responsible to initialize the *data* field. If the *data* field needs not to be initialized, just give a *nullptr* to this parameter.

**walker_update_state_func**: This parameter is a function that takes a walker, its current residing vertex, and an edge as input.  At each step, once an edge is selected for the walker's next step, this function is invoked to update the walker's states. Note that users are only responsible to update the *data* field. If the *data* field needs not to be updated, just give a *nullptr* to this parameter.

**walker_init_dist_func**: This parameter is a function that takers a walker's id as input, and return the id of the vertex it starts to walk from. KnightKing provides two implementations for this function:

- **std::function<vertex_id_t (walker_id_t)> get_equal_dist_func()**: This API returns a function that assigns i_th walker to (i % vertex_num)_th vertex.
- **std::function<vertex_id_t (walker_id_t)> get_uniform_dist_func()**: This API returns a function that for each walker it randomly assign the walker to all vertices with equal probability.

 Besides the above two functions, users can define custom functions as well.

### Random Walk

#### First Order Random Walk

To define the edge transition probability for random walk, users need to define the extension component, static component, and dynamic component, respectively. These components are wrapped in class **TransitionConfig**.

```c++
TransitionConfig(
    std::function<real_t (Walker<walker_data_t>&, vertex_id_t)> extension_comp_func,
    std::function<real_t (vertex_id_t, AdjUnit<edge_data_t>*)> static_comp_func = nullptr,
    std::function<real_t (Walker<walker_data_t>&, vertex_id_t, AdjUnit<edge_data_t> *)> dynamic_comp_func = nullptr,
    std::function<real_t (vertex_id_t, AdjList<edge_data_t>*)> dcomp_upperbound_func = nullptr,
    std::function<real_t (vertex_id_t, AdjList<edge_data_t>*)> dcomp_lowerbound_func = nullptr,
    std::function<void(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, real_t&, vertex_id_t&)> outlier_upperbound_func = nullptr,
    std::function<AdjUnit<edge_data_t>*(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, vertex_id_t)> outlier_search_func = nullptr
)
```

**extension_comp_func**: This function returns a real number representing the extension component.

**static_comp_func**: This function returns a real number representing the static component.

**dynamic_comp_func**: This function returns a real number representing the dynamic component.

If *static_comp_func* is not defined (left as *nullptr*),  then KnightKing assumes a trivial static component. Otherwise, KnightKing uses alias method to deal with static component.

If *dynamic_comp_func* is not defined (left as *nullptr*), then KnightKing assumes a trivial dynamic component, and the rejection-sampling-based solution degenerates to a static sampling algorithm (alias method in our implementation).

If *dynamic_comp_func* is defined, then *dcomp_upperbound_func* must also be defined. While *dcomp_lowerbound_func* is just an optional optimization.

There is no default terminate condition, so *extension_comp_func* cannot be left as undefined.

**dcomp_upperbound_func**: This function returns a real number representing the upper bound of the dynamic component. Each vertex has its own upper bound.

**dcomp_lowerbound_func**: This function returns a real number representing the lower bound of the dynamic component. Each vertex has its own lower bound.

**outlier_upperbound_func**: This function sets the upperbound for how much area an outlier can occupy at the appendix area, and the upperbound of outlier number. These two values are set to the last two parameters of this function.

**outlier_search_func**: This function returns the i-th outlier. The number *i* is the last parameter of this function. If there are less than i outliers, just return nullptr.

For example, in node2vec with p = 0.5 and q = 2, all edges have dynamic component no more than 1 except the return edge, whose dynamic component is 2. In this case, setting the upper bound of dynamic component as 2 will lead to high rejection probability for rejection sampling. Instead, we can set the dynamic component upper bound as 1. Then before each sampling, give an estimated upperbound for the overflow area of the return edge. Suppose we find that the return edge has static component set as *S*, then the overflow area (area beyong dynamic component upper bound) is *S × (2 - 1)*, i.e. *S*. The example usage can be found at [src/apps/node2vec.hpp](src/apps/node2vec.hpp).

At the beginning, *static_comp_func* are invoked to initiate the alias table, and *dcomp_upperbound_func*, and *dcomp_lowerbound_func* are invoked, too. Their returned values will be stored for future use.

The *extension_comp_func*, *dynamic_comp_func*, *outlier_upperbound_func* and *outlier_search_func* are invoked ad-hocly during the walk.

For more information, please refer to the [KnightKing slides and paper](#Additional-Documentation).

### Second Order Random Walk

With first-order walk algorithms, walkers are oblivious to the vertices visited before the current one, while with second-order algorithms, in selecting the next stop, a walker considers the previous vertex visited, from which it transited to the current one.

In distributed environment, the previously visited vertex and currently residing vertex may not locate at the same node. So access previously visited vertex information from currently residing vertex requires communications. Thus to run second order random walk, two additional functions are required for communication.

```c++
SecondOrderTransitionConfig(
    std::function<real_t (Walker<walker_data_t>&, vertex_id_t)> extension_comp_func,
    std::function<real_t (vertex_id_t, AdjUnit<edge_data_t>*)> static_comp_func,
    std::function<void (Walker<walker_data_t>&, walker_id_t, vertex_id_t, AdjUnit<edge_data_t> *)> post_query_func,
    std::function<void (vertex_id_t, stateQuery<query_data_t> &, AdjList<edge_data_t>*)> respond_query_func,
    std::function<real_t (Walker<walker_data_t>&, vertex_id_t, AdjUnit<edge_data_t> *)> dynamic_comp_func,
    std::function<real_t (vertex_id_t, AdjList<edge_data_t>*)> dcomp_upperbound_func,
    std::function<real_t (vertex_id_t, AdjList<edge_data_t>*)> dcomp_lowerbound_func = nullptr,
    std::function<void(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, real_t&, vertex_id_t&)> outlier_upperbound_func = nullptr,
    std::function<AdjUnit<edge_data_t>*(Walker<walker_data_t>&, vertex_id_t, AdjList<edge_data_t>*, vertex_id_t)> outlier_search_func = nullptr
)
```

**post_query_func**: This function sends a query message from the walker's currently residing vertex to another vertex. The query message should be attached with the source vertex id and the walker's walker index. i.e. its second parameter. These two items jointly tell the source location of the message, so that latter the destination vertex knows where to send back the response.

**respond_query_func**: This function receives a query message, generates corresponding response and sends it back.

## Walk Configuration

Some miscellaneous parameters are put into the **WalkConfig**, which supports following settings:

```c++
void set_output_file(const char* output_path_prefix_param, bool with_head_info = false);
void set_output_consumer(std::function<void (PathSet*)> output_consumer_func_param);
void set_walk_rate(double rate_param);
```

**set_output_file**: This function tells KnightKing to record the walking paths for the walkers, and dump the paths into files.

**set_output_consumer**: This function tells KnightKing to record the walking paths for the walkers, and return them to users via a callback function provided by users. The paths are placed at a class *PathSet*, the detail of which will be given later.

Note that output file and consumer can be set at the same time. In this case, paths will first be flushed to files, then be given to users.

**set_walk_rate**: As introduced before, to save memory, this function sets how many walkers can simutaneously walk. For example, set walk rate to 0.5, then there will be 2 epochs. Each time 50% walkers will walk. After they stop, if the output file is set, then the walks will be flushed to files before the next epoch. If the output consumer is set, then the walks will be given to users before the next epoch.

Instead of dumping the paths, users can directly pass the in-memory output data to other applications. The paths are stored in *PathSet* object:

```c++
struct PathSet
{
    int seg_num;
    vertex_id_t **path_set;
    walker_id_t **walker_id;
    vertex_id_t ***path_begin;
    vertex_id_t ***path_end;
    step_t **path_length;
    walker_id_t *path_num;
}
```

The paths are stored distributedly, i.e., each KnightKing instance holds a different part of the paths. Further more, in each instance, the paths are stored in several segments. But it is guaranteed that the path of the same walker is located at the continuous space.

Each *PathSet* stores the paths for one KnightKing instance.

The *seg_num* specifies how many segments there are;

The *path_num\[i\]* specifies how many paths the i-th segment have;

The *walker_id\[i\]\[j\]* specifies the walker id of the j-th path of the i-th segment;

The *path_length\[i\]\[j\]* specifies the path length of the j-th path of the i-th segment;

The *path_begin\[i\]\[j\]* and *path_end\[i\]\[j\]* specifies the begin and end of the location of the j-th path of the i-th segment.

A sample traversing code:

```c++
PathSet *ps = nullptr;
WalkConfig walk_conf;
walk_conf.set_output_consumer(
    [&] (PathSet* ps_param) {
        // Assume only has one iteration
        ps = ps_param;
    }
);
graph.random_walk(&walker_conf, &tr_conf, &walk_conf);
for (int i = 0; i < ps->seg_num; i++)
{
    for (walker_id_t j = 0; j < ps->path_num[wo_i]; j++)
    {
        printf("%u %u", ps->walker_id[i][j], ps->path_length[i][j]);
        for (step_t k = 0; k < ps->path_length[i][j]; k++)
        {
            printf(" %u", *(ps->path_begin[i][j] + k));
        }
        printf("\n");
    }
}
```

### Miscellaneous

#### Load Graph

**load_graph**: This function takes the vertex number and graph data file path as input, and loads the graph into memory.

```c++
void load_graph(vertex_id_t vertex_num, const char* graph_file_path, bool load_as_undirected = false);
```

#### Set Concurrency
**set_concurrency**: This function sets the number of threads for concurrent task execution. The default setting is std::thread::hardware_concurrency() - 1, i.e. the number of cores minus one. This number is set slightly lower than core number in default, since there are two additional threads that are responsible for message sending and receiving.

```c++
void set_concurrency(int worker_num);
```

### Constants

There are several numerical constants that can be adjusted for performance tuning. They are defined in *include/constants.hpp*.

**L1_CACHE_LINE_SIZE**: This constant tells KnightKing the size of L1 cache line. KnightKing uses this information to decide the size of paddings that are used to avoid performance degenerating caused by simultaneously writing to close positions by different threads.

**THREAD_LOCAL_BUF_CAPACITY**: Each thread has its own local message buffer. Messages are first written to that buffer, then flushed to global buffer when the local buffer gets full. This constant defines how many messages the local buffer can hold.

**OMP_PARALLEL_THRESHOLD**: This constant defines the threshold of the number of active walkers. Beyond the threshold KnightKing executes tasks concurrently using multi-threading, while below the threshold KnightKing executes the tasks without multi-threading

**PARALLEL_CHUNK_SIZE**: The tasks, oftentimes the walkers or queries, are grouped as chunks, then put into a task pool. This constant defines the multi-thread scheduling granularity.

## Publication

Ke Yang, Mingxing Zhang, Kang Chen, Xiaosong Ma, Yang Bai, Yong Jiang. KnightKing: A Fast Distributed Graph Random Walk Engine. In ACM SIGOPS 27th Symposium on Operating Systems Principles (SOSP ’19).
