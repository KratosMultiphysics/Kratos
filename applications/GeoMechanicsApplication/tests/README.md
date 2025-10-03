This subdirectory contains geotechnical tests.


# Result plots

Some tests can plot their results.  For instance, some settlement test cases (see `test_dsettlement_validation.py`) can plot the evolution of settlement over time.  And some of them can plot the effective vertical stress distribution as well as the water pressure distribution over the depth.  To let the test scripts generate such plots, run them in a shell where the environment variable `KRATOS_GEO_MAKE_TEST_PLOTS` is set to `ON`.
