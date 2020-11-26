# These tests should *only* be run in parallel. E.g.
# runcompss parallel_tests.py

import unittest

if __name__ == "__main__":
    parallel_loader = unittest.TestLoader()
    parallel_loader.testMethodPrefix = "parallel_test_"
    parallel_tests = parallel_loader.discover(".", pattern="*Tests.py")
    parallel_tests.addTest(parallel_loader.discover(".", pattern="test*.py"))
    runner = unittest.runner.TextTestRunner(verbosity=1)
    runner.run(parallel_tests)
