import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from check_version import VersionPolicyError, validate_version_policy


class VersionPolicyTests(unittest.TestCase):
    def test_branch_may_match_latest_tag(self):
        self.assertEqual(
            validate_version_policy("1.1.0", latest_tag="v1.1.0"),
            "1.1.0",
        )

    def test_branch_may_advance_beyond_latest_tag(self):
        self.assertEqual(
            validate_version_policy("1.2.0", latest_tag="1.1.0"),
            "1.2.0",
        )

    def test_branch_may_not_move_behind_latest_tag(self):
        with self.assertRaisesRegex(VersionPolicyError, "older than"):
            validate_version_policy("1.0.9", latest_tag="1.1.0")

    def test_release_tag_must_match_exactly(self):
        self.assertEqual(
            validate_version_policy(
                "1.1.0", latest_tag="1.1.0", release_tag="v1.1.0"
            ),
            "1.1.0",
        )

    def test_release_tag_mismatch_is_rejected(self):
        with self.assertRaisesRegex(VersionPolicyError, "does not match"):
            validate_version_policy("1.2.0", release_tag="v1.1.0")

    def test_non_semantic_version_is_rejected(self):
        with self.assertRaisesRegex(VersionPolicyError, "MAJOR.MINOR.PATCH"):
            validate_version_policy("release-1.1")

    def test_leading_zero_is_rejected(self):
        with self.assertRaisesRegex(VersionPolicyError, "leading zeroes"):
            validate_version_policy("01.1.0")


if __name__ == "__main__":
    unittest.main()
