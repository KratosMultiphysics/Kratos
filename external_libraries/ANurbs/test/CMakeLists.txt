add_executable(run_tests
    main.cc
    TestKnots.cc
    TestPoint.cc
    TestCurveShapeEvaluator.cc
    TestCurveGeometry.cc
    TestCurve.cc
    TestCurveTessellation.cc
    TestSurfaceShapeEvaluator.cc
    TestSurfaceGeometry.cc
    TestSurface.cc
    TestCurveOnSurface.cc
    TestIntegrationPoints.cc
    TestCurveSpanIntersection.cc
    TestPointOnCurveProjection.cc

    "${PROJECT_SOURCE_DIR}/external_libraries/ClipperLib/clipper.cpp"
)

install(TARGETS run_tests DESTINATION bin)
