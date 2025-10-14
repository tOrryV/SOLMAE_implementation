import unittest
import time
from io import StringIO

from tests import test_modular_big, test_poly, test_ntt, test_cfft


def run_tests_with_timing(test):
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromModule(test)

    result_stream = StringIO()
    runner = unittest.TextTestRunner(stream=result_stream, verbosity=2)

    start_time = time.perf_counter()
    result = runner.run(suite)
    total_time = (time.perf_counter() - start_time) * 1000

    print(result_stream.getvalue())
    print("=" * 70)
    print(f"Complete {result.testsRun} tests")
    print(f"Failed: {len(result.failures) + len(result.errors)}")
    print(f"Total execution time: {total_time:.3f} ms")
    print("=" * 70)


def main():
    # print(f'============================== TEST BIG MODULAR OPERATIONS ==========================')
    # run_tests_with_timing(test_modular_big)

    # print(f'============================== TEST POLYNOM OPERATIONS ==========================')
    # run_tests_with_timing(test_poly)

    # print(f'============================== TEST NTT ==========================')
    # run_tests_with_timing(test_ntt)

    # print(f'============================== TEST NTT ==========================')
    # run_tests_with_timing(test_cfft)

if __name__ == "__main__":
    main()