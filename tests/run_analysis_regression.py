#!/usr/bin/env python3
"""End-to-end numerical regression checks for representative AMDAT analyses."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import subprocess
import tempfile


REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_FILE = REPO_ROOT / "tests" / "analysis_regression.in"
EXPECTED_DIR = REPO_ROOT / "tests" / "expected" / "analysis_regression"
REL_TOL = 5e-6
ABS_TOL = 1e-10

ANALYSES = {
    "msd": ("msd.dat", "Mean square displacement data", 2),
    "md": ("md.dat", "Mean displacement data", 4),
    "msd2d": ("msd2d.dat", "2-D mean square displacement data", 2),
    "baf": ("baf.dat", "Bond autocorrelation function data", 2),
}


def read_table(path: Path, columns: int, *, has_header: bool) -> list[list[float]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if has_header:
        if not lines:
            raise AssertionError(f"{path}: output is empty")
        lines = lines[1:]

    rows: list[list[float]] = []
    for line_number, line in enumerate(lines, start=2 if has_header else 1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != columns:
            raise AssertionError(
                f"{path}:{line_number}: expected {columns} columns, got {len(fields)}"
            )
        row = [float(field) for field in fields]
        if not all(math.isfinite(value) for value in row):
            raise AssertionError(f"{path}:{line_number}: non-finite value")
        rows.append(row)
    return rows


def compare_tables(name: str, actual: list[list[float]], expected: list[list[float]]) -> tuple[float, float]:
    if len(actual) != len(expected):
        raise AssertionError(
            f"{name}: expected {len(expected)} rows, got {len(actual)}"
        )

    max_absolute = 0.0
    max_relative = 0.0
    for row_number, (actual_row, expected_row) in enumerate(
        zip(actual, expected), start=1
    ):
        for column, (actual_value, expected_value) in enumerate(
            zip(actual_row, expected_row), start=1
        ):
            absolute = abs(actual_value - expected_value)
            relative = absolute / max(abs(expected_value), 1e-30)
            max_absolute = max(max_absolute, absolute)
            max_relative = max(max_relative, relative)
            if not math.isclose(
                actual_value,
                expected_value,
                rel_tol=REL_TOL,
                abs_tol=ABS_TOL,
            ):
                raise AssertionError(
                    f"{name}: row {row_number}, column {column}: "
                    f"expected {expected_value:.17g}, got {actual_value:.17g}"
                )
    return max_absolute, max_relative


def run_regression(binary: Path, threads: int) -> None:
    with tempfile.TemporaryDirectory(prefix=f"amdat-regression-{threads}-") as output:
        output_dir = Path(output)
        result = subprocess.run(
            [
                str(binary),
                "-q",
                "-n",
                str(threads),
                "-c",
                "OUTPUT_DIR",
                str(output_dir),
                "-i",
                str(INPUT_FILE),
            ],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        if result.returncode != 0:
            raise AssertionError(
                f"AMDAT failed with {threads} thread(s), exit {result.returncode}\n"
                f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
            )

        for name, (filename, header_prefix, columns) in ANALYSES.items():
            actual_path = output_dir / filename
            if not actual_path.is_file():
                raise AssertionError(f"{name}: AMDAT did not create {actual_path}")

            actual_lines = actual_path.read_text(encoding="utf-8").splitlines()
            if not actual_lines or not actual_lines[0].startswith(header_prefix):
                raise AssertionError(f"{name}: missing or unexpected output header")

            actual = read_table(actual_path, columns, has_header=True)
            expected = read_table(
                EXPECTED_DIR / f"{name}.expected", columns, has_header=False
            )
            max_absolute, max_relative = compare_tables(name, actual, expected)
            print(
                f"PASS {name}: {threads} thread(s), {len(actual)} rows, "
                f"max abs={max_absolute:.3g}, max rel={max_relative:.3g}"
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--binary",
        type=Path,
        default=REPO_ROOT / "AMDAT",
        help="AMDAT executable to test (default: repository AMDAT binary)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        nargs="+",
        default=[1, 4],
        help="thread counts to test (default: 1 4)",
    )
    args = parser.parse_args()

    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error(f"binary does not exist: {binary}; run make first")
    if any(threads < 1 for threads in args.threads):
        parser.error("thread counts must be positive")

    for threads in args.threads:
        run_regression(binary, threads)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
