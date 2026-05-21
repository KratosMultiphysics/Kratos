#ifndef CONTEXT_SETUP_HPP
#define CONTEXT_SETUP_HPP

#include <memory>
#include <iostream>
#include <vexcl/devlist.hpp>
#include "random_vector.hpp"

struct ContextSetup {
    ContextSetup() {
        try {
            context = std::make_shared<vex::Context>(vex::Filter::DoublePrecision && vex::Filter::Env);
        } catch(const vex::backend::error &e) {
            std::cerr << "Failed to initialize compute context" << std::endl;
            std::cerr << "Error: " << e << std::endl;
            throw;
        }

        unsigned seed = static_cast<uint>(time(0));
        std::cout << "seed: " << seed << std::endl;

        srand(seed);

#ifndef VEXCL_BACKEND_JIT
        // If there is only one device in context, duplicate the command queues
        // in order to properly test multi-device capabilities.
        if (context->queue().size() == 1) {
            vex::Context second(vex::Filter::DoublePrecision && vex::Filter::Env);

            std::vector<vex::backend::context>       c = context->context();
            std::vector<vex::backend::command_queue> q = context->queue();

            c.push_back(second.context(0));
            q.push_back(second.queue(0));

            context = std::make_shared<vex::Context>(c, q);
            vex::StaticContext<>::set(*context);
        }
#endif

        std::cout << *context << std::endl;
    }

    std::shared_ptr<vex::Context> context;
};

struct ContextReference {
    ContextReference() : ctx( vex::current_context() ) { }

    const vex::Context &ctx;
};

#define SAMPLE_SIZE 32

template<class V, class F>
void check_sample(const V &v, F f) {
    for(size_t i = 0; i < SAMPLE_SIZE; ++i) {
        size_t idx = rand() % v.size();
        f(idx, v[idx]);
    }
}

template<class V1, class V2, class F>
void check_sample(const V1 &v1, const V2 &v2, F f) {
    BOOST_REQUIRE(v1.size() == v2.size());
    for(size_t i = 0; i < SAMPLE_SIZE; ++i) {
        size_t idx = rand() % v1.size();
        f(idx, v1[idx], v2[idx]);
    }
}

template<class V1, class V2, class V3, class F>
void check_sample(const V1 &v1, const V2 &v2, const V3 &v3, F f) {
    BOOST_REQUIRE(v1.size() == v2.size());
    BOOST_REQUIRE(v1.size() == v3.size());

    for(size_t i = 0; i < SAMPLE_SIZE; ++i) {
        size_t idx = rand() % v1.size();
        f(idx, v1[idx], v2[idx], v3[idx]);
    }
}

BOOST_GLOBAL_FIXTURE( ContextSetup );
BOOST_FIXTURE_TEST_SUITE(cr, ContextReference)

BOOST_AUTO_TEST_CASE(context_ready)
{
    BOOST_REQUIRE(ctx);
}

#endif
