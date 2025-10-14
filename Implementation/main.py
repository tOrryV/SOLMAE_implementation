import time
import unittest
import sys

from tests import test_modular_big, test_poly


def run_tests_with_timing(test_module) -> bool:
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromModule(test_module)

    print(f"\n===== Running {test_module.__name__} =====")
    start = time.perf_counter()

    runner = unittest.TextTestRunner(stream=sys.stdout, verbosity=2, buffer=True)
    result = runner.run(suite)

    elapsed_ms = (time.perf_counter() - start) * 1000.0
    print("=" * 70)
    print(f"Complete {result.testsRun} tests")
    print(f"Failures: {len(result.failures)}  Errors: {len(result.errors)}  Skipped: {len(result.skipped)}")
    print(f"Total execution time: {elapsed_ms:.3f} ms")
    print("=" * 70)
    return result.wasSuccessful()


def main():
    ok = True

    # ok &= run_tests_with_timing(test_modular_big)

    ok &= run_tests_with_timing(test_poly)

    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
