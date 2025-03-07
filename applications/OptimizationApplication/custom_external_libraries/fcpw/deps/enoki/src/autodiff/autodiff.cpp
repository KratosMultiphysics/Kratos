/*
    src/autodiff/autodiff.cpp -- Reverse mode automatic differentiation

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyrighe (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#include <enoki/dynamic.h>
#if defined(ENOKI_CUDA)
#include <enoki/cuda.h>
#endif
#include <enoki/autodiff.h>

#include <unordered_map>
#include <set>
#include <sstream>
#include <iomanip>
#include <algorithm>

#if defined(NDEBUG)
#  define ENOKI_AUTODIFF_DEFAULT_LOG_LEVEL 0
#else
#  define ENOKI_AUTODIFF_DEFAULT_LOG_LEVEL 1
#endif

/// Max. allowed cost in number of arithmetic operations that a simplification can do
#define ENOKI_AUTODIFF_MAX_SIMPLIFICATION_COST 10

NAMESPACE_BEGIN(enoki)

using Index = uint32_t;

template <typename Value>
Value safe_mul(const Value &value1, const Value &value2);
template <typename Value>
Value safe_fmadd(const Value &value1, const Value &value2, const Value &value3);

template <typename Value> struct Tape<Value>::Node {
    /// Descriptive label
    std::string label;

    /// Gradient value
    Value grad;

    /// Pointer to incident edge linked list
    std::vector<Edge> edges;

    /// Reverse edge list
    std::vector<uint32_t> edges_rev;

    /// External (i.e. by Enoki) reference count
    uint32_t ref_count_ext = 0;

    /// Internal (i.e. within the computation graph) reference count
    uint32_t ref_count_int = 0;

    /// Size of the variable
    uint32_t size = 0;

    Node(size_t size, const char *label)
        : label(label ? label : ""), size((uint32_t) size) { }

    bool is_scalar() const {
        return size == 1;
    }

    bool collapse_allowed() const {
        return !edges.empty() && !edges_rev.empty();
    }

    Index score() const {
        return (uint32_t) (edges.size() * edges_rev.size());
    }

    Edge *edge(Index source) {
        for (auto &e: edges) {
            if (e.source == source)
                return &e;
        }
        return nullptr;
    }

    Edge remove_edge(Index source) {
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            if (it->source == source) {
                Edge temp(std::move(*it));
                edges.erase(it);
                return temp;
            }
        }
        throw std::runtime_error("Node::remove_edge(): not found!");
    }

    Node() = default;
    Node(const Node &) = delete;
    Node(Node&&) = default;
    Node& operator=(const Node &) = delete;
    Node& operator=(Node&&) = default;
};

template <typename Value> struct Tape<Value>::Edge {
    /// Source node ID associated with this edge
    Index source;

    /// Edge weight
    Value weight;

    /// Optional: special operation (scatter/gather/reduction)
    std::unique_ptr<Special> special;

    /// Pointer to next edge
    std::unique_ptr<Edge> next;

    Edge(Index source, const Value &weight)
        : source(source), weight(weight) { }

    Edge(Index source, Special *special)
        : source(source), special(special) { }

    bool is_special() const { return special != nullptr; }

    Edge() = default;
    Edge(const Edge &) = delete;
    Edge(Edge&&) = default;
    Edge& operator=(const Edge &) = delete;
    Edge& operator=(Edge&&) = default;
};

template <typename Value> struct Tape<Value>::Special {
    virtual void backward(Detail *detail, Index target_idx, const Edge &edge) const {
        throw std::runtime_error("Special::backward(): not implemented!");
    }

    virtual void forward(Detail *detail, Index target_idx, const Edge &edge) const {
        throw std::runtime_error("Special::forward(): not implemented!");
    }

    virtual ~Special() = default;
};

template <typename Value> struct Tape<Value>::Detail {
    Index node_counter = 1,
          node_counter_last = 1;

    std::unordered_map<Index, Node> nodes;
    std::vector<std::string> prefix;
    Index *scatter_gather_index = nullptr;
    size_t scatter_gather_size = 0;
    bool scatter_gather_permute = false;
    uint32_t log_level = ENOKI_AUTODIFF_DEFAULT_LOG_LEVEL;
    bool graph_simplification = true,
         is_simplified = true;

    /// Set of indices selected for next backward pass
    std::set<uint32_t> scheduled;

    Node &node(Index index) {
        auto it = nodes.find(index);
        if (it == nodes.end())
            throw std::runtime_error("autodiff: Detail::node(): Unknown index " +
                                     std::to_string(index));
        return it->second;
    }

    void dfs(Index k, bool backward, bool clear_grad) {
        if (scheduled.find(k) != scheduled.end())
            return;
        scheduled.insert(k);

        Node &n = node(k);
        if (clear_grad) {
            if (is_dynamic_v<Value>)
                n.grad = Value();
            else
                n.grad = zero<Value>();
        }

        if (backward) {
            for (const Edge &edge : n.edges)
                dfs(edge.source, backward, clear_grad);
        } else {
            for (Index k2: n.edges_rev)
                dfs(k2, backward, clear_grad);
        }
    }
};

template <typename Value> struct Tape<Value>::SimplificationLock {
    SimplificationLock(Tape &tape) : tape(tape) {
        std::swap(state, tape.d->graph_simplification);
    }

    ~SimplificationLock() {
        std::swap(state, tape.d->graph_simplification);
    }

    Tape &tape;
    bool state = false;
};

template <typename Value> std::unique_ptr<Tape<Value>> Tape<Value>::s_tape;
template <typename Value> ENOKI_PURE Tape<Value> *Tape<Value>::get() {
    if (ENOKI_UNLIKELY(!s_tape))
        s_tape = std::unique_ptr<Tape>(new Tape());
    return s_tape.get();
}

template <typename Value> Tape<Value>::Tape() {
    d = new Detail();

#if defined(ENOKI_CUDA)
    if constexpr (is_cuda_array_v<Value>)
        cuda_register_callback((void (*)(void *)) & Tape::cuda_callback, this);
#endif
}

template <typename Value> Tape<Value>::~Tape() {
#if !defined(NDEBUG)
    if (d->log_level >= 1) {
        if (d->node_counter != 1)
            std::cerr << "autodiff: shutdown." << std::endl;
        size_t n_live = 0;
        for (const auto &it : d->nodes) {
            if (n_live < 10)
                std::cerr << "autodiff: variable " << it.first
                          << " still live at shutdown. (ref_count_int="
                          << it.second.ref_count_int
                          << ", ref_count_ext=" << it.second.ref_count_ext << ")"
                          << std::endl;
            if (n_live == 9)
                std::cerr << "(skipping remainder)" << std::endl;
            n_live++;
        }
        if (n_live > 0)
            std::cerr << "autodiff: " << n_live
                      << " variables were still live at shutdown." << std::endl;
    }
#endif
    delete d;
}

template <typename Value> void Tape<Value>::cuda_callback(void *ptr) {
    Tape *tape = (Tape *) ptr;
    if (tape->d->graph_simplification)
        tape->simplify_graph();
}

template <typename Value> void Tape<Value>::set_log_level(uint32_t level) {
    d->log_level = level;
}

template <typename Value> uint32_t Tape<Value>::log_level() const {
    return d->log_level;
}

template <typename Value> void Tape<Value>::set_graph_simplification(bool value) {
    d->graph_simplification = value;
}

template <typename Value>
Index Tape<Value>::append(const char *label, size_t size, Index i1, const Value &w1) {
    if (i1 == 0)
        return 0;
    Index idx = append_node(size, label);
#if !defined(NDEBUG)
    if (d->log_level >= 3)
        std::cerr << "autodiff: append(\"" << (label ? label : "") << "\", " << idx
                  << " <- " << i1 << ")" << std::endl;
#endif
    append_edge(i1, idx, w1);
    return idx;
}

template <typename Value>
Index Tape<Value>::append(const char *label, size_t size, Index i1, Index i2,
                          const Value &w1, const Value &w2) {
    if (i1 == 0 && i2 == 0)
        return 0;
    Index idx = append_node(size, label);
#if !defined(NDEBUG)
    if (d->log_level >= 3)
        std::cerr << "autodiff: append(\"" << (label ? label : "") << "\", " << idx
                  << " <- [" << i1 << ", " << i2 << "])" << std::endl;
#endif
    append_edge(i1, idx, w1);
    append_edge(i2, idx, w2);
    return idx;
}

template <typename Value>
Index Tape<Value>::append(const char *label, size_t size, Index i1, Index i2, Index i3,
                          const Value &w1, const Value &w2, const Value &w3) {
    if (i1 == 0 && i2 == 0 && i3 == 0)
        return 0;
    Index idx = append_node(size, label);
#if !defined(NDEBUG)
    if (d->log_level >= 3)
        std::cerr << "autodiff: append(\"" << (label ? label : "") << "\", " << idx
                  << " <- [" << i1 << ", " << i2 << ", " << i3 << "])" << std::endl;
#endif
    append_edge(i1, idx, w1);
    append_edge(i2, idx, w2);
    append_edge(i3, idx, w3);
    return idx;
}

template <typename Value>
Index Tape<Value>::append_node(size_t size, const char *label) {
    Index idx = d->node_counter++;
    auto result = d->nodes.emplace(std::make_pair(idx, Node(size, label)));

    Node &node = result.first->second;
    for (auto it = d->prefix.rbegin(); it != d->prefix.rend(); ++it)
        node.label = *it + '/' + node.label;

#if !defined(NDEBUG)
    if (d->log_level >= 3)
        std::cerr << "autodiff: append_node(\"" << (label ? label : "")
                  << "\", size=" << size << ") -> " << idx << std::endl;
#endif
    inc_ref_ext(idx);
    d->is_simplified = false;
    return idx;
}

template <typename Value>
Index Tape<Value>::append_leaf(size_t size) {
    Index idx = append_node(size, "'unnamed'");
    Node &n = d->node(idx);
    n.grad = zero<Value>(n.size);
    return idx;
}

template <typename Value>
void Tape<Value>::set_label(Index idx, const char *label) {
    if (idx == 0)
        return;
#if !defined(NDEBUG)
    if (d->log_level >= 3)
        std::cerr << "autodiff: set_label(" << idx << ") -> " << label << std::endl;
#endif
    std::string name = "'" + std::string(label) + "'";
    Node &n = d->node(idx);
    n.label = name;
    enoki::set_label(n.grad, (label + std::string(".grad")).c_str());
}

template <typename Value>
Index Tape<Value>::append_gather(const Int64 &offset, const Mask &mask) {
    if constexpr (is_dynamic_v<Value>) {
        if (d->scatter_gather_index == nullptr ||
           *d->scatter_gather_index == 0)
            return 0;
        Index source = *d->scatter_gather_index;

        struct Gather : Special {
            Int64 offset;
            Mask mask;
            size_t size;
            bool permute;

            void forward(Detail *detail, Index target_idx,
                         const Edge &edge) const override {
                const Value &grad_source = detail->node(edge.source).grad;
                Value &grad_target = detail->node(target_idx).grad;

                if (grad_source.size() != size)
                    throw std::runtime_error("Internal error in Gather::forward()!");

                Value value = gather<Value>(grad_source, offset, mask);

                if (grad_target.empty())
                    grad_target = value;
                else
                    grad_target += value;
            }

            void backward(Detail *detail, Index target_idx,
                          const Edge &edge) const override {
                const Value &grad_target = detail->node(target_idx).grad;
                Value &grad_source = detail->node(edge.source).grad;

                if (grad_source.empty())
                    grad_source = zero<Value>(size);
                else if (grad_source.size() != size)
                    throw std::runtime_error("Internal error in Gather::backward()!");

                if (permute)
                    scatter(grad_source, grad_target, offset, mask);
                else
                    scatter_add(grad_source, grad_target, offset, mask);
            }
        };

        Gather *gather = new Gather();
        gather->offset = offset;
        gather->mask = mask;
        gather->size = d->scatter_gather_size;
        gather->permute = d->scatter_gather_permute;

        Index target = append_node(slices(offset), "gather");
        d->node(target).edges.emplace_back(source, gather);
        inc_ref_int(source, target);

#if !defined(NDEBUG)
        if (d->log_level >= 3)
            std::cerr << "autodiff: append_gather(" << target << " <- " << source << ")"
                      << std::endl;
#endif

        return target;
    } else {
        return 0;
    }
}

template <typename Value>
Index Tape<Value>::append_reverse(Index source) {
    if (source == 0)
        return 0;

    if constexpr (is_dynamic_v<Value>) {
        struct Reverse : Special {
            void forward(Detail *detail, Index target_idx, const Edge &edge) const override {
                const Value &grad_source = detail->node(edge.source).grad;
                Value &grad_target = detail->node(target_idx).grad;

                Value result = reverse(grad_source);

                if (grad_target.empty())
                    grad_target = result;
                else
                    grad_target += result;
            }

            void backward(Detail *detail, Index target_idx, const Edge &edge) const override {
                const Value &grad_target = detail->node(target_idx).grad;
                Value &grad_source = detail->node(edge.source).grad;

                Value result = reverse(grad_target);

                if (grad_source.empty())
                    grad_source = result;
                else
                    grad_source += result;
            }
        };

        Reverse *s = new Reverse();

        Index target = append_node(d->node(source).size, "reverse");
        d->node(target).edges.emplace_back(source, s);
        inc_ref_int(source, target);

#if !defined(NDEBUG)
        if (d->log_level >= 3)
            std::cerr << "autodiff: append_reverse(" << target << " <- "
                      << source << ")" << std::endl;
#endif

        return target;
    } else {
        throw std::runtime_error("append_reverse(): internal error!");
    }
}

template <typename Value>
Index Tape<Value>::append_psum(Index source) {
    if (source == 0)
        return 0;

    if constexpr (is_dynamic_v<Value>) {
        struct PrefixSum : Special {
            void forward(Detail *detail, Index target_idx, const Edge &edge) const override {
                const Value &grad_source = detail->node(edge.source).grad;
                Value &grad_target = detail->node(target_idx).grad;

                Value result = psum(grad_source);

                if (grad_target.empty())
                    grad_target = result;
                else
                    grad_target += result;
            }

            void backward(Detail *detail, Index target_idx, const Edge &edge) const override {
                const Value &grad_target = detail->node(target_idx).grad;
                Value &grad_source = detail->node(edge.source).grad;

                Value result = reverse(psum(reverse(grad_target)));

                if (grad_source.empty())
                    grad_source = result;
                else
                    grad_source += result;
            }
        };

        PrefixSum *s = new PrefixSum();

        Index target = append_node(d->node(source).size, "psum");
        d->node(target).edges.emplace_back(source, s);
        inc_ref_int(source, target);

#if !defined(NDEBUG)
        if (d->log_level >= 3)
            std::cerr << "autodiff: append_psum(" << target << " <- "
                      << source << ")" << std::endl;
#endif

        return target;
    } else {
        throw std::runtime_error("append_psum(): internal error!");
    }
}

template <typename Value>
void Tape<Value>::append_scatter(Index source, const Int64 &offset, const Mask &mask, bool scatter_add) {
    if constexpr (is_dynamic_v<Value>) {
        SimplificationLock lock(*this);

        if (d->scatter_gather_index == nullptr || source == 0)
            return;
        Index target_orig = *d->scatter_gather_index;

        struct Scatter : Special {
            Int64 offset;
            Mask mask;
            size_t size;
            bool scatter_add;

            void forward(Detail *detail, Index target_idx, const Edge &edge) const override {
                const Value &grad_source = detail->node(edge.source).grad;
                Value &grad_target = detail->node(target_idx).grad;

                if (grad_target.empty())
                    grad_target = zero<Value>(size);
                else if (grad_target.size() != size)
                    throw std::runtime_error("Internal error in Scatter::forward()!");

                if (scatter_add)
                    enoki::scatter_add(grad_target, grad_source, offset, mask);
                else
                    enoki::scatter(grad_target, grad_source, offset, mask);
            }

            void backward(Detail *detail, Index target_idx, const Edge &edge) const override {
                Node &source = detail->node(edge.source);
                const Value &grad_target = detail->node(target_idx).grad;
                Value &grad_source = source.grad;

                if (grad_target.size() != size)
                    throw std::runtime_error("Internal error in Scatter::backward()!");

                Value result = gather<Value>(grad_target, offset, mask);
                if (source.size == 1)
                    result = hsum(result);
                else if (result.size() == 1 && source.size != 1)
                    set_slices(result, source.size);

                if (grad_source.empty())
                    grad_source = result;
                else
                    grad_source += result;
            }
        };

        Scatter *s = new Scatter();
        s->offset = offset;
        s->mask = mask;
        s->size = d->scatter_gather_size;
        s->scatter_add = scatter_add;

        Index target_new = append_node(d->scatter_gather_size,
                                       scatter_add ? "scatter_add" : "scatter");
        d->node(target_new).edges.emplace_back(source, s);
        inc_ref_int(source, target_new);

        if (target_orig != 0) {
            Index sa_node = target_new;

            Value weight = 1.f;
            if (!scatter_add && !d->scatter_gather_permute) {
                weight = full<Value>(1.f, d->scatter_gather_size);
                scatter(weight, Value(0), offset, mask);
            }
            target_new = append("scatter_combine", d->scatter_gather_size,
                                target_new, target_orig, 1, weight);
            dec_ref_ext(sa_node);
            dec_ref_ext(target_orig);
        }

        *d->scatter_gather_index = target_new;

#if !defined(NDEBUG)
        if (d->log_level >= 3)
            std::cerr << "autodiff: append_scatter(" << target_orig << " <- "
                      << source << ", scatter_add=" << scatter_add << ") -> "
                      << target_new << std::endl;
#endif
    }
}

template <typename Value>
void Tape<Value>::append_edge(Index source_idx, Index target_idx,
                              const Value &weight) {
    if (source_idx == 0)
        return;
    assert(target_idx != 0);

#if !defined(NDEBUG)
    if (d->log_level >= 4)
        enoki::set_label(weight, ("edge[" + std::to_string(source_idx) + " -> " +
                                  std::to_string(target_idx) + "]").c_str());
#endif

    Node &target = d->node(target_idx);
    if (Edge *edge = target.edge(source_idx); edge != nullptr) {
#if !defined(NDEBUG)
        if (d->log_level >= 4)
            std::cerr << "autodiff: append_edge(" << target_idx << " <- "
                      << source_idx << "): merging."
                      << std::endl;
#endif
        SimplificationLock lock(*this);
        edge->weight += weight;
    } else {
#if !defined(NDEBUG)
        if (d->log_level >= 4)
            std::cerr << "autodiff: append_edge(" << target_idx << " <- "
                      << source_idx << "): creating."
                      << std::endl;
#endif
        target.edges.emplace_back(source_idx, weight);
        inc_ref_int(source_idx, target_idx);
    }
}

template <typename Value>
void Tape<Value>::append_edge_prod(Index source_idx, Index target_idx,
                                   const Value &weight1, const Value &weight2) {
    if (source_idx == 0)
        return;
    assert(target_idx != 0);

    Node &target = d->node(target_idx);
    if (Edge *edge = target.edge(source_idx); edge != nullptr) {
        Value weight = safe_fmadd(weight1, weight2, edge->weight);
#if !defined(NDEBUG)
        if (d->log_level >= 4) {
            std::cerr << "autodiff: append_edge_prod(" << target_idx << " <- "
                      << source_idx << "): merging."
                      << std::endl;
            enoki::set_label(weight, ("edge[" + std::to_string(source_idx) + " -> " +
                                      std::to_string(target_idx) + "]").c_str());
        }
#endif
        edge->weight = weight;
    } else {
        Value weight = safe_mul(weight1, weight2);
#if !defined(NDEBUG)
        if (d->log_level >= 4) {
            std::cerr << "autodiff: append_edge_prod(" << target_idx << " <- "
                      << source_idx << "): creating."
                      << std::endl;
            enoki::set_label(weight, ("edge[" + std::to_string(source_idx) + " -> " +
                                      std::to_string(target_idx) + "]").c_str());
        }
#endif
        target.edges.emplace_back(source_idx, weight);
        inc_ref_int(source_idx, target_idx);
    }
}

template <typename Value> void Tape<Value>::inc_ref_int(Index index, Index from) {
    Node &node = d->node(index);

#if !defined(NDEBUG)
    if (d->log_level >= 4)
        std::cerr << "autodiff: inc_ref_int(" << index << ", " << from
                  << ") -> " << (node.ref_count_int + 1) << std::endl;
#endif

    auto it = std::find(node.edges_rev.begin(), node.edges_rev.end(), from);
    if (it != node.edges_rev.end())
        throw std::runtime_error("inc_ref_int(): internal error -- edge already exists!");

    node.edges_rev.push_back(from);
    node.ref_count_int++;
}

template <typename Value> void Tape<Value>::dec_ref_int(Index index, Index from) {
    if (index == 0)
        return;
    Node &node = d->node(index);

#if !defined(NDEBUG)
    if (d->log_level >= 4)
        std::cerr << "autodiff: dec_ref_int(" << index << ", " << from
                  << ") -> " << (node.ref_count_int - 1) << std::endl;

    if (node.ref_count_int == 0)
        throw std::runtime_error("autodiff: dec_ref_int(): Node " +
                                 std::to_string(index) +
                                 " has no internal references!");
#endif
    --node.ref_count_int;

    auto it = std::find(node.edges_rev.begin(), node.edges_rev.end(), from);
    if (ENOKI_UNLIKELY(it == node.edges_rev.end()))
        throw std::runtime_error("dec_ref_int(): internal error -- edge not found!");

    node.edges_rev.erase(it);

    if (node.ref_count_int == 0 && node.ref_count_ext == 0)
        free_node(index);
}

template <typename Value> void Tape<Value>::inc_ref_ext(Index index) {
    if (index == 0)
        return;
    Node &node = d->node(index);
    node.ref_count_ext++;

#if !defined(NDEBUG)
    if (d->log_level >= 4)
        std::cerr << "autodiff: inc_ref_ext(" << index << ") -> " << node.ref_count_ext << std::endl;
#endif
}

template <typename Value> void Tape<Value>::dec_ref_ext(Index index) {
    if (index == 0)
        return;
    Node &node = d->node(index);

#if !defined(NDEBUG)
    if (d->log_level >= 4)
        std::cerr << "autodiff: dec_ref_ext(" << index << ") -> " << (node.ref_count_ext - 1) << std::endl;

    if (node.ref_count_ext == 0)
        throw std::runtime_error("autodiff: dec_ref_ext(): Node " +
                                 std::to_string(index) +
                                 " has no external references!");
#endif

    --node.ref_count_ext;

    if (node.ref_count_int == 0 && node.ref_count_ext == 0)
        free_node(index);
}

template <typename Value> void Tape<Value>::free_node(Index index) {
#if !defined(NDEBUG)
    if (d->log_level >= 4)
        std::cerr << "autodiff: free_node(" << index << ")" << std::endl;
#endif

    auto it = d->nodes.find(index);
    if (it == d->nodes.end())
        throw std::runtime_error("autodiff: free_node(): Unknown index " +
                                 std::to_string(index));

    Node &node = it->second;
    for (const Edge &edge : node.edges)
        dec_ref_int(edge.source, index);

    d->nodes.erase(it);
}

template <typename Value> void Tape<Value>::push_prefix(const char *value) {
    d->prefix.push_back(value);
}

template <typename Value> void Tape<Value>::pop_prefix() {
    if (d->prefix.empty())
        throw std::runtime_error("pop_prefix(): prefix list is already empty!");
    d->prefix.pop_back();
}

template <typename Value>
void Tape<Value>::set_scatter_gather_operand(Index *index, size_t size,
                                             bool permute) {
    if (ENOKI_UNLIKELY(index != nullptr && d->scatter_gather_index != 0))
        throw std::runtime_error("set_scatter_gather_operand(): attempted to override an existing operand!");
    d->scatter_gather_index = index;
    d->scatter_gather_size = size;
    d->scatter_gather_permute = permute;
}

template <typename Value> const Value &Tape<Value>::gradient(Index index) {
    if (index == 0)
        throw std::runtime_error(
            "No gradient was computed for this variable! (a call to "
            "requires_gradient() is necessary.)");
    return d->node(index).grad;
}

template <typename Value>
void Tape<Value>::backward(Index index, bool free_graph) {
    using Scalar = scalar_t<Value>;

    SimplificationLock lock(*this);
    set_gradient(index, Scalar(1), true);
    backward(free_graph);
}

template <typename Value>
void Tape<Value>::forward(Index index, bool free_graph) {
    using Scalar = scalar_t<Value>;

    SimplificationLock lock(*this);
    set_gradient(index, Scalar(1), false);
    forward(free_graph);
}

template <typename Value>
void Tape<Value>::set_gradient(Index index, const Value &value, bool backward) {
    if (index == 0)
        throw std::runtime_error(
            "set_gradient(): no gradients are associated with this variable (a "
            "prior call to requires_gradient() is required.) ");

    d->dfs(index, backward, true);
    Node &node = d->node(index);
    node.grad = value;
    if constexpr (is_dynamic_v<Value>) {
        if (node.size > 1 && node.grad.size() == 1)
            set_slices(node.grad, node.size);
    }
}

template <typename Value>
void Tape<Value>::backward(bool free_graph) {
    auto &scheduled = d->scheduled;

    if (free_graph) {
        for (auto it = scheduled.begin(); it != scheduled.end(); ++it)
            inc_ref_ext(*it);
    }

    for (auto it = scheduled.rbegin(); it != scheduled.rend(); ++it) {
        Index target_idx = *it;
        Node &target = d->node(target_idx);

        if constexpr (is_dynamic_v<Value>) {
            if (ENOKI_UNLIKELY(target.size != target.grad.size())) {
                if (target.grad.size() == 1)
                    set_slices(target.grad, target.size);
                else
                    throw std::runtime_error(
                        "backward(): gradient sizes don't match: expected " +
                        std::to_string(target.size) + ", got " +
                        std::to_string(target.grad.size()));
            }
        }

        for (Edge &edge : target.edges) {
            Node &source = d->node(edge.source);
            if (ENOKI_LIKELY(!edge.is_special())) {
                if constexpr (is_dynamic_v<Value>) {
                    if (source.size == 1 && (edge.weight.size() != 1 || target.grad.size() != 1)) {
                        if (source.grad.empty())
                            source.grad = hsum(safe_mul(edge.weight, target.grad));
                        else
                            source.grad += hsum(safe_mul(edge.weight, target.grad));
                    } else {
                        if (source.grad.empty())
                            source.grad = safe_mul(edge.weight, target.grad);
                        else
                            source.grad = safe_fmadd(edge.weight, target.grad, source.grad);
                    }
                } else {
                    source.grad = safe_fmadd(edge.weight, target.grad, source.grad);
                }
            } else {
                edge.special->backward(d, target_idx, edge);
            }
            if (free_graph) {
                dec_ref_int(edge.source, target_idx);
                edge.source = 0;
            }
        }
        if (free_graph) {
            if (target.edges.size() > 0) {
                target.edges.clear();
                target.grad = Value();
            }
            dec_ref_ext(target_idx);
        } else {
            if (target.ref_count_int > 0)
                target.grad = Value();
        }
    }

    if (d->log_level >= 1)
        std::cerr << "autodiff: backward(): processed " << scheduled.size() << "/"
                  << (d->node_counter - d->node_counter_last) << " nodes."
                  << std::endl;

    if (free_graph)
        d->node_counter_last = d->node_counter;

    scheduled.clear();
}

template <typename Value>
void Tape<Value>::forward(bool free_graph) {
    auto &scheduled = d->scheduled;

    if (free_graph) {
        for (auto it = scheduled.begin(); it != scheduled.end(); ++it)
            inc_ref_ext(*it);
    }

    for (auto it = scheduled.begin(); it != scheduled.end(); ++it) {
        Index source_idx = *it;
        Node &source = d->node(source_idx);

        if constexpr (is_dynamic_v<Value>) {
            if (source.size == 1 && source.grad.size() > 1)
                source.grad = hsum(source.grad);
        }

        for (Index target_idx : source.edges_rev) {
            Node &target = d->node(target_idx);
            Edge *edge = target.edge(source_idx);
            if (edge == nullptr)
                throw std::runtime_error("forward(): invalid graph structure!");

            if (ENOKI_LIKELY(!edge->is_special())) {
                if constexpr (is_dynamic_v<Value>) {
                    if (target.size == 1 && (edge->weight.size() != 1 || source.grad.size() != 1)) {
                        if (target.grad.empty())
                            target.grad = hsum(safe_mul(edge->weight, source.grad));
                        else
                            target.grad += hsum(safe_mul(edge->weight, source.grad));
                    } else {
                        if (target.grad.empty())
                            target.grad = safe_mul(edge->weight, source.grad);
                        else
                            target.grad = safe_fmadd(edge->weight, source.grad, target.grad);
                    }
                } else {
                    target.grad = safe_fmadd(edge->weight, source.grad, target.grad);
                }
            } else {
                edge->special->forward(d, target_idx, *edge);
            }
            if constexpr (is_dynamic_v<Value>) {
                if (ENOKI_UNLIKELY(target.size != target.grad.size())) {
                    if (target.grad.size() == 1)
                        set_slices(target.grad, target.size);
                    else
                        throw std::runtime_error(
                            "forward(): gradient sizes don't match: expected " +
                            std::to_string(target.size) + ", got " +
                            std::to_string(target.grad.size()));
                }
            }
        }
        if (source.ref_count_int > 0)
            source.grad = Value();
        if (free_graph) {
            auto edges_rev = source.edges_rev;
            for (Index target_idx : edges_rev) {
                dec_ref_int(source_idx, target_idx);
                d->node(target_idx).remove_edge(source_idx);
            }
            dec_ref_ext(source_idx);
        }
    }

    if (d->log_level >= 1)
        std::cerr << "autodiff: forward(): processed " << scheduled.size() << "/"
                  << (d->node_counter - d->node_counter_last) << " nodes."
                  << std::endl;

    if (free_graph)
        d->node_counter_last = d->node_counter;

    scheduled.clear();
}

template <typename Value> void Tape<Value>::simplify_graph() {
    if (d->is_simplified)
        return;

    SimplificationLock lock(*this);
    if (d->log_level >= 2)
        std::cerr << "autodiff: simplify_graph(): starting.." << std::endl;

    std::set<std::pair<Index, Index>> todo;
    std::vector<std::pair<Index, Index>> update;
    std::vector<Index> edges_rev;
    for (const auto &it : d->nodes)
        todo.emplace(it.second.score(), it.first);
    size_t cost = 0;

    while (!todo.empty()) {
        auto it = todo.begin();
        uint32_t score = it->first, index = it->second;
        todo.erase(it);
        Node &node = d->node(index);
        if (!node.collapse_allowed())
            continue;
        if (score > ENOKI_AUTODIFF_MAX_SIMPLIFICATION_COST) {
            if (d->log_level >= 2)
                std::cerr << "autodiff: simplify_graph(): cost of next simplification = " << cost << ", giving up." << std::endl;
            break;
        }

        update.clear();
        bool skip = false;

        /* Collect predecessors and successors */ {
            for (Index k : node.edges_rev) {
                Edge *e = d->node(k).edge(index);
                assert(e != nullptr);
                if (e->is_special())
                    skip = true;
                update.emplace_back(d->node(k).score(), k);
            }
            for (const Edge &edge : node.edges) {
                const Node &node2 = d->node(edge.source);
                update.emplace_back(node2.score(), edge.source);
                if ((node.size == 1 && (node2.size != node.size)) || edge.is_special())
                    skip = true;
            }
        }

        if (skip)
            continue;

        if (d->log_level >= 3)
            std::cerr << "autodiff: simplify_graph(): collapsing node " << index << ", cost = " << score << std::endl;

        /* Remove node and create edges */ {
            edges_rev = node.edges_rev;
            for (Index other : edges_rev) {
                Edge edge1 = d->node(other).remove_edge(index);

                for (auto const &edge2 : node.edges) {
                    append_edge_prod(edge2.source, other, edge1.weight, edge2.weight);
                    cost++;
                }

                dec_ref_int(index, other);
            }
        }

        /* Update costs (if changed) */ {
            for (auto [old_score, id] : update) {
                auto it = todo.find({ old_score, id });
                if (it == todo.end())
                    continue;
                uint32_t new_score = d->node(id).score();
                if (old_score != new_score) {
                    todo.erase(it);
                    todo.emplace(new_score, id);
                }
            }
        }
    }

    if (d->log_level >= 2)
        std::cerr << "autodiff: simplify_graph(): done. (cost = " << cost << ")" << std::endl;
    d->is_simplified = true;
}

template <typename Value>
std::string Tape<Value>::graphviz(const std::vector<Index> &indices_) {
    std::ostringstream oss;
    oss << "digraph {" << std::endl
        << "  rankdir=BT;" << std::endl // BT or RL
        << "  fontname=Consolas;" << std::endl
        << "  node [shape=record fontname=Consolas];" << std::endl;

    for (Index index : indices_)
        d->dfs(index, true, false);

    auto &indices = d->scheduled;

    int current_depth = 0;
    auto hasher = std::hash<std::string>();
    std::string current_path = "";

    for (Index index : indices) {
        const Node &node = d->node(index);
        if (node.label.empty())
            continue;
        std::string label = node.label;

        auto sepidx = label.rfind("/");
        std::string path, suffix;
        if (sepidx != std::string::npos) {
            path = label.substr(0, sepidx);
            label = label.substr(sepidx + 1);
        }

        if (current_path != path) {
            for (int i = 0; i < current_depth; ++i)
                oss << "  }" << std::endl;
            current_depth = 0;
            current_path = path;

            do {
                sepidx = path.find('/');
                std::string graph_label = path.substr(0, sepidx);
                if (graph_label.empty())
                    break;

                oss << "  subgraph cluster"
                    << std::to_string(hasher(graph_label)) << " {" << std::endl
                    << "  label=\"" << graph_label << "\";" << std::endl;
                ++current_depth;

                if (sepidx == std::string::npos)
                    break;
                path = path.substr(sepidx + 1, std::string::npos);
            } while (true);
        }

        oss <<  "  " << std::to_string(index) << " [label=\"" + label;
        if (node.is_scalar())
            oss << " [s]";

        oss << "\\n#" << std::to_string(index) << " [E/I: "
            << std::to_string(node.ref_count_ext) << "/"
            << std::to_string(node.ref_count_int) << "]"
            << "\"";
        if (node.label[0] == '\'')
            oss << " fillcolor=salmon style=filled";
        oss << "];" << std::endl;
    }
    for (int i = 0; i < current_depth; ++i)
        oss << "  }\n";

    for (Index index : indices) {
        const Node &node = d->node(index);
        for (const Edge &edge : node.edges) {
            oss << "  " << std::to_string(index) << " -> "
                << std::to_string(edge.source) << ";" << std::endl;

            if (edge.is_special())
                oss << "  " << std::to_string(index)
                    << " [shape=doubleoctagon];" << std::endl;
        }
    }

    for (Index idx : indices_)
        oss << "  " << std::to_string(idx)
            << " [fillcolor=cornflowerblue style=filled];" << std::endl;

    oss << "}";
    indices.clear();
    return oss.str();
}

template <typename Value> std::string Tape<Value>::whos() const {
    std::ostringstream oss;
    oss << std::endl
        << "  ID      E/I Refs   Size        Label" << std::endl
        << "  ====================================" << std::endl;

    std::vector<uint32_t> indices;
    indices.reserve(d->nodes.size());
    for (const auto& it : d->nodes)
        indices.push_back(it.first);
    std::sort(indices.begin(), indices.end());

    for (uint32_t id : indices) {
        const Node &n = d->node(id);
        oss << "  " << std::left << std::setw(7) << id << " ";
        oss << std::left << std::setw(10) << (std::to_string(n.ref_count_ext) + " / " + std::to_string(n.ref_count_int)) << " ";
        oss << std::left << std::setw(12) << n.size;
        oss << n.label;
        oss << std::endl;
    }

    oss << "  ====================================" << std::endl << std::endl;

    return oss.str();
}

template <typename Value> Value safe_mul(const Value &value1, const Value &value2) {
    Value tentative = value1 * value2;
    if constexpr (!is_cuda_array_v<Value>) {
        Value zero = scalar_t<Value>(0);
        mask_t<Value> is_zero = eq(value1, zero) || eq(value2, zero);
        return select(is_zero, zero, tentative);
    } else {
#if defined(ENOKI_CUDA)
        using Mask = mask_t<Value>;
        Mask m1 = Mask::from_index_(cuda_trace_append(EnokiType::Bool, "setp.eq.f32 $r1, $r2, 0.0", value1.index_())),
             m2 = Mask::from_index_(cuda_trace_append(EnokiType::Bool, "setp.eq.or.f32 $r1, $r2, 0.0, $r3", value2.index_(), m1.index_()));
        return Value::from_index_(cuda_trace_append(Value::Type, "selp.$t1 $r1, 0.0, $r2, $r3", tentative.index_(), m2.index_()));
#endif
    }
}

template <typename Value> Value safe_fmadd(const Value &value1, const Value &value2, const Value &value3) {
    Value tentative = fmadd(value1, value2, value3);
    if constexpr (!is_cuda_array_v<Value>) {
        Value zero = scalar_t<Value>(0);
        mask_t<Value> is_zero = eq(value1, zero) || eq(value2, zero);
        return select(is_zero, value3, tentative);
    } else {
#if defined(ENOKI_CUDA)
        using Mask = mask_t<Value>;
        Mask m1 = Mask::from_index_(cuda_trace_append(EnokiType::Bool, "setp.eq.f32 $r1, $r2, 0.0", value1.index_())),
             m2 = Mask::from_index_(cuda_trace_append(EnokiType::Bool, "setp.eq.or.f32 $r1, $r2, 0.0, $r3", value2.index_(), m1.index_()));
        return Value::from_index_(cuda_trace_append(Value::Type, "selp.$t1 $r1, $r2, $r3, $r4", value3.index_(), tentative.index_(), m2.index_()));
#endif
    }
}

template struct ENOKI_EXPORT Tape<float>;
template struct ENOKI_EXPORT DiffArray<float>;

template struct ENOKI_EXPORT Tape<double>;
template struct ENOKI_EXPORT DiffArray<double>;

template struct ENOKI_EXPORT Tape<DynamicArray<Packet<float>>>;
template struct ENOKI_EXPORT DiffArray<DynamicArray<Packet<float>>>;

template struct ENOKI_EXPORT Tape<DynamicArray<Packet<double>>>;
template struct ENOKI_EXPORT DiffArray<DynamicArray<Packet<double>>>;

#if defined(ENOKI_CUDA)
    template struct ENOKI_EXPORT Tape<CUDAArray<float>>;
    template struct ENOKI_EXPORT DiffArray<CUDAArray<float>>;

    template struct ENOKI_EXPORT Tape<CUDAArray<double>>;
    template struct ENOKI_EXPORT DiffArray<CUDAArray<double>>;
#endif

NAMESPACE_END(enoki)
