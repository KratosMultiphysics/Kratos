#define BOOST_TEST_MODULE Image
#include <boost/test/unit_test.hpp>
#include <vexcl/vector.hpp>
#include <vexcl/function.hpp>
#include <vexcl/element_index.hpp>
#include <vexcl/image.hpp>
#include "context_setup.hpp"

#if defined(VEXCL_BACKEND_CUDA)
BOOST_AUTO_TEST_CASE(image1d)
{
    if (!vex::Filter::CC(3, 0)(ctx.device(0))) return;

    VEX_FUNCTION(float, imread, (vex::cuda_tex_object, tex)(int, i),
            return tex1Dfetch<float>(tex, i);
            );

    std::vector<vex::backend::command_queue> q1(1, ctx.queue(0));

    const int n = 1024;
    const int m = 16;

    vex::vector<float> imdata(q1, m);
    imdata = 42;

    CUDA_RESOURCE_DESC res_desc; memset(&res_desc, 0, sizeof(res_desc));

    res_desc.resType = CU_RESOURCE_TYPE_LINEAR;
    res_desc.res.linear.devPtr = imdata(0).raw();
    res_desc.res.linear.format = CU_AD_FORMAT_FLOAT;
    res_desc.res.linear.numChannels = 1;
    res_desc.res.linear.sizeInBytes = m * sizeof(float);

    CUDA_TEXTURE_DESC tex_desc; memset(&tex_desc, 0, sizeof(tex_desc));

    CUtexObject tex;
    cuTexObjectCreate(&tex, &res_desc, &tex_desc, NULL);

    vex::vector<int> p(q1, n);
    p = vex::element_index() % m;

    vex::vector<float> x(q1, n);

    x = imread(vex::wrap(tex), p);

    check_sample(x, [](size_t, float a) {
            BOOST_CHECK_EQUAL(a, 42);
            });

    cuTexObjectDestroy(tex);
}
#elif defined(CL_VERSION_1_2)
#if defined(VEXCL_BACKEND_OPENCL)
BOOST_AUTO_TEST_CASE(image1d)
{
    struct imread_t : vex::UserFunction<imread_t, cl_float4(cl::Image1D, int)> {
        imread_t() {}

        static std::string name() { return "imread"; }

        static void define(vex::backend::source_generator &src) {
            src << VEX_STRINGIZE_SOURCE(
                    __constant sampler_t sampler =
                          CLK_NORMALIZED_COORDS_FALSE
                        | CLK_ADDRESS_CLAMP_TO_EDGE
                        | CLK_FILTER_NEAREST;

                    float4 imread(__read_only image1d_t im, int i) {
                        return read_imagef(im, sampler, i);
                    }
                    );
        }
    } const imread;

    std::vector<vex::backend::command_queue> q1(1, ctx.queue(0));

    const int n = 1024;
    const int m = 16;
    std::vector<float> imdata;
    for(int i = 0; i < m; ++i)
        for(int j = 0; j < 4; ++j)
            imdata.push_back(static_cast<float>(j));

    cl::Image1D image(ctx.context(0), CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            cl::ImageFormat(CL_RGBA, CL_FLOAT), m, imdata.data());

    vex::vector<int> p(q1, n);
    p = vex::element_index() % m;

    vex::vector<cl_float4> x(q1, n);

    x = imread(image, p);

    check_sample(x, [](size_t, cl_float4 a) {
            BOOST_CHECK_EQUAL(a.s[0], 0);
            BOOST_CHECK_EQUAL(a.s[1], 1);
            BOOST_CHECK_EQUAL(a.s[2], 2);
            BOOST_CHECK_EQUAL(a.s[3], 3);
            });
}
#elif defined(VEXCL_BACKEND_COMPUTE)
BOOST_AUTO_TEST_CASE(image1d)
{
    struct imread_t : vex::UserFunction<imread_t, cl_float4(boost::compute::image1d, int)> {
        imread_t() {}

        static std::string name() { return "imread"; }

        static void define(vex::backend::source_generator &src) {
            src << VEX_STRINGIZE_SOURCE(
                    __constant sampler_t sampler =
                          CLK_NORMALIZED_COORDS_FALSE
                        | CLK_ADDRESS_CLAMP_TO_EDGE
                        | CLK_FILTER_NEAREST;

                    float4 imread(__read_only image1d_t im, int i) {
                        return read_imagef(im, sampler, i);
                    }
                    );
        }
    } const imread;

    std::vector<vex::backend::command_queue> q1(1, ctx.queue(0));

    const int n = 1024;
    const int m = 16;
    std::vector<float> imdata;
    for(int i = 0; i < m; ++i)
        for(int j = 0; j < 4; ++j)
            imdata.push_back(j);

    boost::compute::image1d image(
            ctx.context(0), m,
            boost::compute::image_format(CL_RGBA, CL_FLOAT),
            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            imdata.data());

    vex::vector<int> p(q1, n);
    p = vex::element_index() % m;

    vex::vector<cl_float4> x(q1, n);

    x = imread(image, p);

    check_sample(x, [](size_t, cl_float4 a) {
            BOOST_CHECK_EQUAL(a.s[0], 0);
            BOOST_CHECK_EQUAL(a.s[1], 1);
            BOOST_CHECK_EQUAL(a.s[2], 2);
            BOOST_CHECK_EQUAL(a.s[3], 3);
            });
}
#endif
#endif

BOOST_AUTO_TEST_SUITE_END()
