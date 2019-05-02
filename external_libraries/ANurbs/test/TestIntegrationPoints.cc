#include "catch.hpp"

#include <ANurbs/src/IntegrationPoints.h>

using namespace ANurbs;

TEST_CASE( "Throws exception for unsupported degree", "[IntegrationPoints]" ) {
    Interval<double> domain(0, 1);

    CHECK_THROWS( IntegrationPoints<double>::Points1(  0, domain) );
    CHECK_THROWS( IntegrationPoints<double>::Points1(100, domain) );
}

TEST_CASE( "Integration points of domain [0, 1]", "[IntegrationPoints]" ) {
    int degree = 2;
    Interval<double> domain(0, 1);

    auto integrationPoints = IntegrationPoints<double>::Points1(degree, domain);

    REQUIRE( integrationPoints[0].t == 0.21132486540518713 );
    REQUIRE( integrationPoints[1].t == 0.7886751345948129  );

    REQUIRE( integrationPoints[0].weight == 0.5 );
    REQUIRE( integrationPoints[1].weight == 0.5 );
}

TEST_CASE( "Integration points of domain [-1, 1]", "[IntegrationPoints]" ) {
    int degree = 2;
    Interval<double> domain(-1, 1);

    auto integrationPoints = IntegrationPoints<double>::Points1(degree, domain);

    REQUIRE( integrationPoints[0].t == -0.5773502691896257 );
    REQUIRE( integrationPoints[1].t ==  0.5773502691896257 );

    REQUIRE( integrationPoints[0].weight == 1.0 );
    REQUIRE( integrationPoints[1].weight == 1.0 );
}

TEST_CASE( "Integration points of domain [[0, 1], [0, 1]]",
    "[IntegrationPoints]" ) {
    int degreeU = 2;
    int degreeV = 3;
    Interval<double> domainU(0, 1);
    Interval<double> domainV(0, 1);

    auto integrationPoints = IntegrationPoints<double>::Points2(degreeU,
        degreeV, domainU, domainV);

    REQUIRE( integrationPoints[0].u == 0.21132486540518713 );
    REQUIRE( integrationPoints[1].u == 0.21132486540518713 );
    REQUIRE( integrationPoints[2].u == 0.21132486540518713 );
    REQUIRE( integrationPoints[3].u == 0.7886751345948129  );
    REQUIRE( integrationPoints[4].u == 0.7886751345948129  );
    REQUIRE( integrationPoints[5].u == 0.7886751345948129  );

    REQUIRE( integrationPoints[0].v == 0.1127016653792583  );
    REQUIRE( integrationPoints[1].v == 0.5                 );
    REQUIRE( integrationPoints[2].v == 0.8872983346207417  );
    REQUIRE( integrationPoints[3].v == 0.1127016653792583  );
    REQUIRE( integrationPoints[4].v == 0.5                 );
    REQUIRE( integrationPoints[5].v == 0.8872983346207417  );

    REQUIRE( integrationPoints[0].weight == 0.13888888888888892 );
    REQUIRE( integrationPoints[1].weight == 0.2222222222222222  );
    REQUIRE( integrationPoints[2].weight == 0.13888888888888892 );
    REQUIRE( integrationPoints[3].weight == 0.13888888888888892 );
    REQUIRE( integrationPoints[4].weight == 0.2222222222222222  );
    REQUIRE( integrationPoints[5].weight == 0.13888888888888892 );
}

TEST_CASE( "Integration points of domain [[-1, 1], [-1, 1]]",
    "[IntegrationPoints]" ) {
    int degreeU = 2;
    int degreeV = 3;
    Interval<double> domainU(-1, 1);
    Interval<double> domainV(-1, 1);

    auto integrationPoints = IntegrationPoints<double>::Points2(degreeU,
        degreeV, domainU, domainV);

    REQUIRE( integrationPoints[0].u == -0.5773502691896257 );
    REQUIRE( integrationPoints[1].u == -0.5773502691896257 );
    REQUIRE( integrationPoints[2].u == -0.5773502691896257 );
    REQUIRE( integrationPoints[3].u ==  0.5773502691896257 );
    REQUIRE( integrationPoints[4].u ==  0.5773502691896257 );
    REQUIRE( integrationPoints[5].u ==  0.5773502691896257 );

    REQUIRE( integrationPoints[0].v == -0.7745966692414834 );
    REQUIRE( integrationPoints[1].v ==  0.0                );
    REQUIRE( integrationPoints[2].v ==  0.7745966692414834 );
    REQUIRE( integrationPoints[3].v == -0.7745966692414834 );
    REQUIRE( integrationPoints[4].v ==  0.0                );
    REQUIRE( integrationPoints[5].v ==  0.7745966692414834 );

    REQUIRE( integrationPoints[0].weight == 0.5555555555555557 );
    REQUIRE( integrationPoints[1].weight == 0.8888888888888888 );
    REQUIRE( integrationPoints[2].weight == 0.5555555555555557 );
    REQUIRE( integrationPoints[3].weight == 0.5555555555555557 );
    REQUIRE( integrationPoints[4].weight == 0.8888888888888888 );
    REQUIRE( integrationPoints[5].weight == 0.5555555555555557 );
}
