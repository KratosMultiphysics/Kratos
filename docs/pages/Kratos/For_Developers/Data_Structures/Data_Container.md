---
title: DataContainer
keywords:
tags: [DataContainer DataChunk DataAccessor DataValuePolicy DataHistoryPolicy]
sidebar: kratos_for_developers
summary:
---

# DataContainer

A chunked, type-erased, variable-keyed storage system for entity data, intended as the future backend for the data held by `Node`, geometric entities and constraints (Phase II). In this first phase it is a standalone module: nothing in the core uses it yet.

## Storage model

A `DataContainer` maps each added `Variable` to one **`DataChunk`**: a single contiguous block of memory holding the values of that variable for *all* entities (struct-of-arrays layout), instead of scattering per-entity values across entity objects. A chunk additionally holds one such block *per stored step slot*, laid out step-by-step.

```
DataContainer
 ├── DataChunk<double>   (PRESSURE)    [step0: e0 e1 e2 ...][step1: ...][step2: ...]
 ├── DataChunk<int>      (STEP)        [e0 e1 e2 ...]
 └── DataChunk<array_1d> (VELOCITY)    [e0 e1 e2 ...]
```

How a chunk behaves is decided by two orthogonal policies passed to `DataContainer::Add`:

- **Value policy** (`DataValuePolicy<T>` hierarchy) — *what* is stored and how. `DataValuePolicy<T>` stores one dense value of type `T` per entity and also owns the **zero value** used to initialize entries (note: the zero belongs to the policy, not to the `Variable`; e.g. `DataValuePolicy<int>(-1)` for sparse index variables). `LayeredDataValuePolicy<T>` stores a fixed number of layers (`std::vector<T>`) per entity. `SparseDataValuePolicy<T>` stores values for only a subset of entities (see below). The base class also carries the type-erased operations (`Clone`/`Copy`/`Assign`/`AssignZero`/`Delete`/`Destruct`/`Print`/`Allocate`) that let `DataChunkBase` be manipulated without knowing `T`.
- **History policy** (`DataHistoryPolicyBase` hierarchy) — *how many steps* are kept. `NonHistoricalDataPolicy` keeps a single slot, no history. `HistoricalDataPolicy(category, n)` keeps a ring buffer of `n` step slots bound to one `StepCategory` (time step, iteration step, sub-step); `DataContainer::CloneStepData(category)` advances the ring only for chunks of that category and copies the current data onto the new current slot, so time-step history can be cloned without touching iteration-level data.

Policy comparison semantics matter for lookups: `IsSameType` compares only the dynamic type (used by `Has`/`GetAccessor` — the *configuration*, e.g. number of layers or buffer size, does not need to match), while `IsSame` also compares the configuration (used by `Add` to detect conflicting re-additions and by the variable+policy `GetDataSpan` overloads). Note that `SparseDataValuePolicy` does not override `IsSameType`: a sparse and a dense policy of the same value type compare as the same *type* (sparseness itself is queried through `IsSparse`).

## Accessors

`Add`/`GetAccessor` return a **`DataAccessor<T>`**: a small, trivially copyable handle caching the *index* of the chunk inside the container, so `GetDataSpan(accessor)` addresses the chunk in O(1) — no search by variable — resolving only the step slot through the history policy. An accessor also records which step it refers to (`StepCategory` + steps before current); `GetStepAccessor` derives an accessor to another step of the same data. Accessors do not own anything and are only valid for the container that created them (the variable identity is verified on access).

## Variable canonicalization

`Add` looks the variable up in the Kratos `Registry` (`variables.all.<name>`) and binds the chunk and accessors to the *registry-owned* variable object, so the stored data is independent of the lifetime of the `Variable` object passed by the caller. Consequently only registered variables can be added (all standard Kratos variables are; locally created ones must call `Variable::Register()` first) — `Add` does not register variables itself, it only canonicalizes already-registered ones.

## Sparse storage

Sparse data (values existing for only a few entities) is driven by a dense *index variable* (`Variable<int>`) stored in the same container: entry `i` holds the sparse position of entity `i`, or a negative value if the entity has no data. A chunk created with `SparseDataValuePolicy<T>(index_accessor)` starts **empty** and grows on demand: `UpdateSparseStorage(index_accessor)` rebuilds every sparse chunk keyed on that index to the number of active entries (existing values are discarded, all entries re-zero-initialized), while `AddToSparseStorage(index_accessor, entity_indices)` assigns the next free sparse indices to currently inactive entities and grows the chunks, **preserving** existing values and zero-initializing only the appended entries. Both operations reallocate chunk storage, so previously obtained spans (and NumPy views) are invalidated and must be re-fetched.

## Header/source split

Only the non-template classes (`DataValuePolicyBase`, `DataHistoryPolicyBase` and its two concrete policies, `DataChunkBase`, and the non-template methods of `DataContainer`) have their logic-bearing bodies in a matching `.cpp` file, mirroring how `VariableData` (non-template) has `variable_data.cpp` while `Variable<T>` (template) stays fully inline in `variable.h`. `DataAccessor<T>`, `DataValuePolicy<T>`, `LayeredDataValuePolicy<T>`, `SparseDataValuePolicy<T>` and `DataChunk<T>` are templates and stay header-only, since the module must remain generic over arbitrary user value types (as exercised by the tests); the same applies to `DataContainer::Add`/`GetAccessor`/`GetDataSpan`/`Has`, which are templates for the same reason.

## Python interface

`DataContainer`, the policies and accessors are exposed per value type (`double`, `int`, `array_1d<double,3>`, `std::string`) following the `DoubleVariable` naming precedent: `DoubleDataValuePolicy`, `IntegerSparseDataValuePolicy`, `Array1DDataValuePolicy3`, `StringDataAccessor`, ... `GetDataSpan` returns a zero-copy NumPy view (shape `(n,)`, or `(n, 3)` for `array_1d<double,3>`) whose base object keeps the container alive; string data is exposed through a `StringDataSpan` proxy instead. The same invalidation caveat as in C++ applies: re-fetch views after sparse storage updates. `LayeredDataValuePolicy` is not bound in this phase.

## Semantics worth knowing

- Copying a `DataContainer` is **shallow**: copies share the chunks (`std::shared_ptr`).
- Thread-safety: only `Add` is internally locked; everything else needs external synchronization.
- `DataChunk` owns its storage manually (`new[]`/`delete[]`) and is non-copyable; chunks with zero entities are valid (that is how sparse chunks start).
- `HistoricalDataPolicy::Clone()` resets the running step index — freshly created chunks (`CreateNew`, `Initialize`) start their history at slot 0.
- `Initialize(other, chunk_size)` mirrors another container's chunk *structure* with fresh zero-initialized chunks (dense at `chunk_size`, sparse empty); it does not copy data.
- `Resize(n)` grows/shrinks all **dense** chunks to `n` entities per step, preserving the first `min(old, new)` values of every step slot and zero-filling the tails (sparse chunks are skipped; previously obtained spans are invalidated).