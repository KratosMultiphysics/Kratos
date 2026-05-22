# KaHIP Configuration Reference

## Preconfigurations

Preconfigurations are the recommended entry point. Each sets a bundle of `PartitionConfig` parameters tuned for a specific trade-off.

| Mode | Flag | Quality | Speed | Graph type |
|---|---|---|---|---|
| Fast | `--preconfiguration=fast` | ★★☆☆☆ | ★★★★★ | Meshes |
| Eco | `--preconfiguration=eco` | ★★★☆☆ | ★★★☆☆ | Meshes |
| Strong | `--preconfiguration=strong` | ★★★★★ | ★★☆☆☆ | Meshes |
| Fast Social | `--preconfiguration=fsocial` | ★★☆☆☆ | ★★★★★ | Social/scale-free |
| Eco Social | `--preconfiguration=esocial` | ★★★☆☆ | ★★★☆☆ | Social/scale-free |
| Strong Social | `--preconfiguration=ssocial` | ★★★★★ | ★★☆☆☆ | Social/scale-free |

---

## `kaffpa` Command-Line Options

### Required

| Option | Description |
|---|---|
| `graph_filename` | Path to the graph file (first positional argument) |
| `--k N` | Number of blocks |

### Core Options

| Option | Default | Description |
|---|---|---|
| `--preconfiguration=MODE` | `strong` | Preconfiguration mode (see above) |
| `--imbalance=FLOAT` | 3 | Allowed imbalance in percent (e.g. 3 = 3%) |
| `--seed=INT` | 0 | Random seed |
| `--time_limit=FLOAT` | 0 | Time limit in seconds (0 = no limit) |
| `--output_filename=FILE` | — | Write partition to file |
| `--input_partition=FILE` | — | Start from existing partition |
| `--suppress_output` | off | Suppress verbose output |
| `--version` | — | Print version and exit |
| `--mmap_io` | off | Use memory-mapped I/O (faster for large files) |

### Feature Flags

| Option | Default | Description |
|---|---|---|
| `--connected_blocks` | off | Enforce connected partition blocks (strong only, experimental) |
| `--enable_mapping` | off | Enable process mapping (use with `--hierarchy_parameter_string`) |
| `--perfectly_balance` | off | Enforce exact balance (no imbalance tolerance) |
| `--balance_edges` | off | Also balance edge count across blocks |

### Process Mapping Options

| Option | Description |
|---|---|
| `--enable_mapping` | Enable process mapping |
| `--hierarchy_parameter_string=S` | Colon-separated hierarchy sizes, e.g. `4:8:8` |
| `--distance_parameter_string=S` | Colon-separated distances, e.g. `1:10:100` |
| `--mapping_mode=MODE` | `multisection` or `bisection` |

---

## `kaffpaE` Command-Line Options

kaffpaE accepts all kaffpa options plus:

| Option | Default | Description |
|---|---|---|
| `--time_limit=FLOAT` | required | Wall-clock time limit (seconds) |
| `--mh_enable_tabu_search` | off | Enable tabu search improvement |
| `--mh_enable_kabapE` | off | Enable KaBaPE combine operator |
| `--mh_enable_gal_combine` | off | Enable Galinier combine operator |
| `--mh_initial_population_fraction=F` | 1.0 | Fraction of ranks used for initial population |
| `--mh_flip_coin=INT` | 1 | Probability control for variation |

---

## `parhip` Command-Line Options

| Option | Default | Description |
|---|---|---|
| `graph_filename` | required | Graph file (METIS or BGF) |
| `--k N` | required | Number of blocks |
| `--preconfiguration=MODE` | `fast` | `ultrafastmesh`, `fastmesh`, `ecomesh`, `ultrafastsocial`, `fastsocial`, `ecosocial` |
| `--imbalance=FLOAT` | 3 | Allowed imbalance percent |
| `--seed=INT` | 0 | Random seed |
| `--output_filename=FILE` | — | Write partition to file |
| `--suppress_output` | off | Suppress output |

---

## Internal Configuration Parameters

The `PartitionConfig` struct (in `lib/partition/partition_config.h`) exposes all internal parameters. These are not directly CLI-accessible but are set by preconfigurations.

### Coarsening Parameters

| Parameter | Type | FAST | ECO | STRONG | Description |
|---|---|---|---|---|---|
| `matching_type` | enum | RANDOM_GPA | RANDOM_GPA | GPA | Matching algorithm |
| `edge_rating` | enum | EXPANSIONSTAR | EXPANSIONSTAR | EXPANSIONSTAR | Edge importance metric |
| `contraction_factor` | double | 0.8 | 0.75 | 0.75 | Stop when graph shrinks by this factor |
| `maximum_rounds` | int | unlimited | unlimited | unlimited | Max coarsening rounds |
| `aggressive_random_levels` | int | 3 | 3 | 3 | Random levels before GPA activates |
| `permutation_quality` | enum | FAST | FAST | GOOD | Permutation quality for matching |
| `edge_rating_tiebreaking` | bool | false | false | true | Break ties in edge rating |

### Initial Partitioning Parameters

| Parameter | Type | FAST | ECO | STRONG | Description |
|---|---|---|---|---|---|
| `initial_partitioning_repetitions` | int | 1 | 16 | 64 | Number of IP attempts |
| `bipartition_tries` | int | 10 | 20 | 64 | Bipartition attempts per repeat |
| `minipreps` | int | 0 | 6 | 10 | Mini-preparation steps |
| `bipartition_algorithm` | enum | BIPARTITION | BIPARTITION | BIPARTITION | Algorithm type |

### Refinement Parameters

| Parameter | Type | FAST | ECO | STRONG | Description |
|---|---|---|---|---|---|
| `refinement_type` | enum | FM | FM_FLOW | FM_FLOW | Refinement algorithm |
| `refinement_scheduling_algorithm` | enum | FAST | ACTIVE_BLOCKS | ACTIVE_BLOCKS_REF_KWAY | Scheduling strategy |
| `fm_search_limit` | int | 5 | 15 | 25 | FM search depth limit |
| `kway_rounds` | int | 1 | 3 | 6 | K-way refinement rounds |
| `kway_stop_rule` | enum | SIMPLE | ADAPTIVE | ADAPTIVE | K-way stop criterion |
| `corner_refinement_enabled` | bool | false | false | true | Enable corner refinement |
| `flow_region_factor` | double | 4 | 8 | 8 | Size of flow refinement region |
| `most_balanced_minimum_cuts` | bool | false | false | true | Find most balanced min-cut |
| `toposort_iterations` | int | 0 | 0 | 2 | Topological sort iterations |

### W-Cycle / Multigrid Parameters

| Parameter | Type | FAST | ECO | STRONG | Description |
|---|---|---|---|---|---|
| `global_cycle_iterations` | int | 1 | 1 | 2 | Number of multilevel cycles |
| `use_wcycles` | bool | false | false | true | Use W-cycle vs V-cycle |
| `use_fullmultigrid` | bool | false | false | true | Full multigrid scheme |
| `level_split` | int | 0 | 0 | 0 | Level at which to split into subproblems |

### Node Separator Parameters

| Parameter | Type | Description |
|---|---|---|
| `mode_node_separators` | bool | Enable node separator computation |
| `sep_flows_disabled` | bool | Disable flow-based separator improvement |
| `sep_fm_disabled` | bool | Disable FM separator improvement |
| `sep_loc_fm_disabled` | bool | Disable localized FM |
| `sep_fm_unsucc_steps` | int | FM unsuccessful steps limit |

### Parallel Memetic (kaffpaE) Parameters

| Parameter | Type | Description |
|---|---|---|
| `time_limit` | double | Wall-clock time limit |
| `mh_enable_gal_combine` | bool | Enable Galinier combine |
| `mh_enable_tabu_search` | bool | Enable tabu search |
| `mh_enable_kabapE` | bool | Enable KaBaPE |
| `mh_initial_population_fraction` | double | Population size fraction |
| `mh_flip_coin` | int | Mutation probability |
| `mh_num_ncs_to_compute` | int | Number of non-connected samples |

### Special Features

| Parameter | Type | Description |
|---|---|---|
| `connected_blocks` | bool | Enforce connectivity in each block |
| `balance_edges` | bool | Balance edge count in addition to node weight |
| `kaffpa_perfectly_balance` | bool | Enforce exact node weight balance |
| `disable_reductions` | bool | Skip nested dissection reduction rules |

---

## Quality vs. Speed Trade-offs

### Single run vs. multiple seeds

Running kaffpa multiple times with different seeds and taking the best result is a simple way to improve quality:

```bash
best_cut=999999
for seed in 1 2 3 4 5; do
    cut=$(./kaffpa graph.graph --k 8 --preconfiguration=eco --seed=$seed \
          --suppress_output --output_filename=tmp_part 2>&1 | grep "cut" | awk '{print $NF}')
    if [ "$cut" -lt "$best_cut" ]; then
        best_cut=$cut
        cp tmp_part best_part
    fi
done
```

### Using kaffpaE for maximum quality

```bash
mpirun -n 8 ./kaffpaE graph.graph \
    --k 8 \
    --time_limit=300 \
    --preconfiguration=strong \
    --mh_enable_tabu_search \
    --mh_enable_kabapE \
    --output_filename=best_part
```

### Decision chart

```
Need a partition? Start here:
│
├─ graph > available RAM?  →  Use ParHIP (parhip)
│
├─ Need highest quality?   →  Use kaffpaE with time budget
│
├─ Balanced quality/speed? →  Use kaffpa --preconfiguration=eco
│
├─ Need fast result?       →  Use kaffpa --preconfiguration=fast
│
└─ Social/scale-free?      →  Use *social preconfigurations
```
