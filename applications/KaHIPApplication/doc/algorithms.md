# KaHIP Algorithms

## Overview

KaHIP implements a multilevel framework for graph partitioning. The core idea is to solve the NP-hard partitioning problem on a hierarchy of progressively smaller (coarser) graphs, then project and refine the solution back to the original. This "compress → solve → expand" approach achieves near-optimal quality in practice.

```
Original graph G₀   (n nodes, m edges)
      │ coarsening
      ▼
Coarsened graph G₁  (n₁ << n nodes)
      │ coarsening
      ▼
      ⋮
      │ coarsening
      ▼
Coarsest graph Gₗ   (~2k nodes)
      │ initial partitioning
      ▼
Partition of Gₗ
      │ projection + refinement
      ▼
Partition of Gₗ₋₁
      │ projection + refinement
      ▼
      ⋮
      │ projection + refinement
      ▼
Partition of G₀    (final result)
```

---

## Phase 1: Coarsening

### Edge Rating

Before matching, each edge is assigned an importance score. Highly-rated edges are more likely to be contracted (their endpoints merged). Available ratings:

| Rating           | Description                                                          |
|------------------|----------------------------------------------------------------------|
| `EXPANSIONSTAR`  | `c(u,v) / (deg(u) * deg(v))` — favors edges between low-degree nodes |
| `EXPANSIONSTAR2` | Variation of EXPANSIONSTAR                                           |
| `WEIGHT`         | Raw edge weight `c(u,v)`                                             |
| `PSEUDOGEOM`     | Geometry-inspired heuristic for mesh graphs                          |
| `SEPARATOR_*`    | Variants tuned for node separator computation                        |

### Matching Algorithms

After rating edges, a matching is computed: a set of disjoint edges whose endpoints will be merged.

| Algorithm            | Description                                                                                                           |
|----------------------|-----------------------------------------------------------------------------------------------------------------------|
| `RANDOM_MATCHING`    | Random matching — each node randomly picks a neighbor                                                                 |
| `GPA_MATCHING`       | Global Path Algorithm — computes a near-optimal maximum weighted matching by finding augmenting paths in a path graph |
| `RANDOM_GPA`         | Randomized GPA — combines randomness with GPA quality                                                                 |
| `CLUSTER_COARSENING` | Label propagation with size constraints — assigns nodes to clusters rather than pairs; better for social networks     |

**GPA (Global Path Algorithm)**: The key insight is that paths and cycles in the matching graph can be handled optimally in linear time. GPA first computes a maximal path partition of the graph and then extracts a maximum weight matching from each path. This yields near-optimal matchings without the cubic cost of general maximum matching.

### Contraction

After matching, matched node pairs (u, v) are merged into a single supernode with:
- Weight = w(u) + w(v)
- Adjacency = union of neighbor lists (with merged parallel edges)

The coarsening stops when the graph is small enough for initial partitioning (typically when it has ~2k nodes, configurable via `maximum_rounds` and `contraction_factor`).

---

## Phase 2: Initial Partitioning

At the coarsest level (~2k nodes), an exact or high-quality partition is computed. Multiple attempts are made; the best is kept.

### Bipartition (2-way)

**BFS Bipartition**: Grows two regions from random seeds using BFS until balance is achieved. Fast but low quality.

**FM Bipartition**: Applies Fiduccia-Mattheyses local search to a BFS seed. Iteratively moves nodes across the boundary to reduce edge cut while maintaining balance.

### Recursive Partitioning

For k > 2: recursively bisect the coarsest graph into k parts. Each bisection uses the bipartition methods above.

### Multiple Repetitions

The `initial_partitioning_repetitions` parameter controls how many independent initial partitions are computed. The best (minimum edge cut) is kept. STRONG mode uses 64 repetitions; ECO uses 16; FAST uses 1.

---

## Phase 3: Uncoarsening and Refinement

### Projection

The partition of the coarser graph Gᵢ₊₁ is projected to the finer graph Gᵢ: each node in Gᵢ inherits the block assignment of the supernode it was merged into.

### 2-Way FM Refinement (Fiduccia-Mattheyses)

The FM algorithm is the classical local search for bisection. At each step:
1. Compute gain of moving each boundary node
2. Move the highest-gain node (even if gain is negative)
3. Lock moved nodes
4. After all nodes moved, revert to the best prefix of moves

KaHIP uses a **bucket priority queue** for O(1) gain updates, giving O(n) time per pass.

**Active Block Scheduling**: Rather than refining all pairs of adjacent blocks, KaHIP maintains a priority queue of block pairs ordered by their potential gain. Only the most promising pairs are refined, reducing total refinement time.

### Flow-Based Refinement

For each pair of adjacent blocks, the algorithm:
1. Extracts a small subgraph around the shared boundary
2. Solves a max-flow problem on this subgraph (using push-relabel)
3. The min-cut of this flow problem gives the locally optimal cut

This is more expensive than FM but finds improvements that FM misses. Flow refinement is used selectively in ECO and STRONG modes.

**Most Balanced Minimum Cut**: Among all minimum cuts of the flow problem, the algorithm selects the one that minimizes imbalance. This leads to better balance with no edge cut penalty.

### K-Way FM Refinement

Extends the 2-way FM idea to k-way partitioning. Each boundary node can be moved to any adjacent block. Gain = reduction in edge cut from moving node v from block A to block B.

KaHIP implements several k-way FM variants:
- **Simple stop rule**: Stop after `fm_search_limit` consecutive non-improving moves
- **Adaptive stop rule**: Dynamically adjust the search limit based on graph size

### Label Propagation Refinement

Each boundary node is re-assigned to the block that maximizes the node's connectivity (weighted sum of edges to neighboring blocks). Multiple rounds of label propagation are applied. Faster than FM but may miss complex improvements.

### Cycle Refinement (KaBaPE)

KaBaPE (KaFFPa Beyond a Partition Equilibrium) searches for improving cycles in the **quotient graph** (the graph of blocks). Moving nodes along a negative-weight cycle in the quotient graph reduces edge cut while maintaining balance.

This operator can find improvements that no local 2-way or k-way FM move finds, because it considers simultaneous coordinated moves across multiple blocks.

### Tabu Search

An optional tabu search layer (enabled by `--mh_enable_tabu_search` in kaffpaE) performs a sequence of moves with a tabu list to escape local optima.

---

## Preconfigurations

### FAST

Optimized for speed with acceptable quality loss:
- Matching: RANDOM_GPA (fast approximation)
- Edge rating: EXPANSIONSTAR
- 1 coarsening level per multigrid cycle
- 1 initial partitioning attempt
- FM refinement with simple stop rule and small search limit
- No flow refinement

### ECO

Balanced quality/speed (recommended default):
- Matching: RANDOM_GPA
- Edge rating: EXPANSIONSTAR
- V-cycle multilevel scheme
- 16 initial partitioning repetitions
- FM + flow refinement with active block scheduling
- Most balanced minimum cut enabled

### STRONG

Maximum quality:
- Matching: GPA (near-optimal matching)
- Edge rating: EXPANSIONSTAR
- Full multigrid (W-cycles): 2 iterations
- 64 initial partitioning repetitions
- Full FM + flow refinement
- Most balanced minimum cut
- All refinement operators enabled

### FASTSOCIAL / ECOSOCIAL / STRONGSOCIAL

Same as FAST/ECO/STRONG but with settings tuned for social/scale-free graphs:
- Cluster coarsening instead of matching (handles high-degree hubs better)
- Different edge rating suited for power-law degree distributions
- Adjusted coarsening stop criteria

---

## ParHIP: Distributed Parallel Partitioning

ParHIP extends the multilevel framework to distributed memory using MPI.

### Distributed Coarsening

Replaces matching-based coarsening with **parallel label propagation with size constraints**:

1. Each node is initially its own cluster
2. Iteratively, each node adopts the label of the most connected neighbor (by edge weight)
3. Size constraints prevent any cluster from growing too large (≤ C times average node weight)
4. After convergence, nodes in the same cluster are contracted

This works well for social networks' hierarchical cluster structure. For mesh graphs, the label propagation is adapted with different size constraints.

### Parallel Refinement

After projecting the coarsest partition back, distributed label propagation with balancing is used for refinement. Each MPI rank refines its local boundary.

### Coarsest Level

When the distributed graph is small enough, it is gathered on rank 0 and partitioned using the full serial KaFFPa (or a parallel memetic algorithm in high-quality modes). The partition is then scattered back.

---

## Node Ordering

Fill-in minimizing orderings for sparse direct solvers (Cholesky, LU factorization).

### Reduction Rules

Before nested dissection, degree-1 and simplicial vertex elimination rules are applied exhaustively:
- **Degree-1 reduction**: Remove pendant nodes (they can be placed at the end without fill-in)
- **Simplicial reduction**: Remove vertices whose neighborhood forms a clique (no fill-in introduced)

These reductions significantly shrink the graph before applying the more expensive nested dissection.

### Nested Dissection (ND)

1. Find a separator S in the graph
2. Recursively order both sides A and B
3. Place S at the end: `order(A) ++ order(B) ++ order(S)`

KaHIP uses its multilevel partitioner to find the separator, which produces better separators than traditional BFS-based methods, leading to smaller fill-in.

With METIS installed: `reduced_nd_fast` applies reductions then calls METIS ND, which is typically the best combination.

---

## Edge Partitioning (SPAC)

The SPAC (Split-and-Connect) algorithm partitions edges rather than nodes:

1. **Split phase**: Each edge (u, v) is split into two directed "edge nodes"
2. **Connect phase**: Edge nodes are connected based on shared endpoints (edges sharing a node compete for the same processing unit)
3. **Partition**: The resulting bipartite graph is partitioned using standard node partitioning
4. **Map back**: Each edge is assigned to the block of its edge node

Edge partitioning is particularly useful for graph algorithms where edge-centric computation is preferable (e.g., link analysis, triangle counting on large-scale social networks).

---

## Process Mapping

Models the problem as a **Quadratic Assignment Problem (QAP)**:

Given:
- Communication graph C = (V, E_c, w_c): task graph with communication volumes
- Processor network P = (U, E_p, d_p): processor topology with distances

Find a bijection f: V → U that minimizes:
```
∑_{(i,j) ∈ E_c} w_c(i,j) · d_p(f(i), f(j))
```

KaHIP solves this by:
1. **Multisection**: Recursively partition C along the processor hierarchy
2. **Local search**: FM-like moves that swap node assignments

The processor hierarchy is described by `hierarchy_parameter` (e.g., `[4, 8, 8]` = 4 nodes × 8 sockets × 8 cores) and `distance_parameter` (e.g., `[1, 10, 100]` = intra-core, intra-socket, intra-node communication costs).

---

## ILP Improvements

For the highest possible quality on small graphs:

### ILP Improver (`ilp_improve`)

Given an existing partition, defines a reduced model:
1. Extract the **gain graph**: nodes are boundary vertices, edges model possible moves
2. Solve the partitioning problem on this small model exactly using Gurobi ILP
3. Apply the optimal solution as a post-processing step

Symmetry breaking constraints avoid the exponential symmetry inherent in the graph partitioning problem, making the ILP tractable for moderate-sized instances.

### ILP Exact Solver (`ilp_exact`)

Solves the partitioning problem exactly from scratch using ILP. Only feasible for small graphs (~100 nodes).

---

## Connected Blocks (Experimental)

When `--connected_blocks` is used with the `strong` preconfiguration:

1. A standard multilevel partition is computed
2. **Connectivity enforcement**: After uncoarsening each level, any disconnected block is detected using BFS/DFS
3. Small components are greedily reassigned to connected neighboring blocks, maintaining balance
4. Connectivity-aware FM refinement only makes moves that preserve block connectivity

**Limitation**: The input graph must be connected. Only the `strong` preconfiguration is supported. This is an experimental feature introduced in v3.25.
