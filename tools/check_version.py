#!/usr/bin/env python3
"""Validate AMDAT's VERSION against reachable and release tags."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


SEMVER_RE = re.compile(r"^(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)$")


class VersionPolicyError(ValueError):
    """Raised when VERSION or a tag violates the release policy."""


def normalize_version(value: str, *, label: str, allow_v_prefix: bool = False) -> tuple[str, tuple[int, int, int]]:
    normalized = value.strip()
    if allow_v_prefix and normalized.startswith("v"):
        normalized = normalized[1:]

    match = SEMVER_RE.fullmatch(normalized)
    if match is None:
        raise VersionPolicyError(
            f"{label} must use MAJOR.MINOR.PATCH with no leading zeroes; found {value!r}"
        )

    return normalized, tuple(int(part) for part in match.groups())


def validate_version_policy(
    version: str,
    *,
    latest_tag: str | None = None,
    release_tag: str | None = None,
) -> str:
    normalized_version, version_key = normalize_version(version, label="VERSION")

    if release_tag is not None:
        normalized_release, _ = normalize_version(
            release_tag, label="release tag", allow_v_prefix=True
        )
        if normalized_version != normalized_release:
            raise VersionPolicyError(
                f"VERSION {normalized_version} does not match release tag {release_tag}"
            )

    if latest_tag is not None:
        normalized_latest, latest_key = normalize_version(
            latest_tag, label="latest tag", allow_v_prefix=True
        )
        if version_key < latest_key:
            raise VersionPolicyError(
                f"VERSION {normalized_version} is older than latest reachable tag {normalized_latest}"
            )

    return normalized_version


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version-file", type=Path, default=Path("VERSION"))
    parser.add_argument("--latest-tag")
    parser.add_argument("--release-tag")
    args = parser.parse_args()

    try:
        version = args.version_file.read_text(encoding="utf-8")
    except FileNotFoundError as error:
        parser.error(str(error))

    try:
        normalized = validate_version_policy(
            version,
            latest_tag=args.latest_tag,
            release_tag=args.release_tag,
        )
    except VersionPolicyError as error:
        parser.error(str(error))

    if args.release_tag is not None:
        print(f"VERSION {normalized} matches release tag {args.release_tag}.")
    elif args.latest_tag is not None:
        print(f"VERSION {normalized} is not older than {args.latest_tag}.")
    else:
        print(f"VERSION {normalized} is valid.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
