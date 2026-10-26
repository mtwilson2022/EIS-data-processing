"""Run with: python -m unittest discover -s tests -p test_smart_converter.py"""

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from eis_processing import smart_dat_to_csv as converter


def trial_data(real, imag=None):
    imag = [-value for value in real] if imag is None else imag
    return pd.DataFrame({
        converter.SWEEP: range(1, len(real) + 1),
        converter.FREQUENCY: [1000.0] * len(real),
        converter.REAL: real,
        converter.IMAG: imag,
    })


def write_dat(path, n_trials=1):
    """Small native-format fixture with known impedance and two frequencies."""
    header = b"SMaRT_Data".ljust(32, b" ") + b"2".ljust(16, b" ")
    header = header.ljust(converter.HEADER_SIZE, b"\x00")
    records = []
    for sweep in range(n_trials):
        for point, frequency in enumerate((1000.0, 100.0)):
            records.append(converter.RECORD.pack(
                sweep, point, 0,
                0.0, frequency, 0.1, 0.0, 0.0, 0.0, 0.0,
                10.0 + sweep + point, -5.0 - sweep - point, 2, 0,
            ))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header + b"".join(records) + b"SDEM")


class ConverterTests(unittest.TestCase):
    def test_single_trial_has_only_raw_columns_and_never_screens(self):
        for enabled in (False, True):
            with self.subTest(enabled=enabled), patch.object(
                converter, "paired_outlier_mask", side_effect=AssertionError("screened")
            ):
                result = converter.combine_trials(trial_data([10]), exclude_outliers=enabled)
                self.assertEqual(list(result), [
                    converter.FREQUENCY, "Z_real_trial1", "Z_imag_trial1",
                ])
                self.assertEqual(result.loc[0, "Z_imag_trial1"], -10)

    def test_default_uses_all_trials_without_audit_columns(self):
        with patch.object(converter, "paired_outlier_mask", side_effect=AssertionError("screened")):
            result = converter.combine_trials(trial_data([1, 1, 1, 1, 100]))
        self.assertEqual(len(result.columns), 15)
        self.assertAlmostEqual(result.loc[0, "Z_real_avg"], 20.8)
        self.assertAlmostEqual(result.loc[0, "Z_real_std"], np.std([1, 1, 1, 1, 100], ddof=1))
        self.assertNotIn("has_outliers", result)

    def test_enabled_screening_excludes_one_pair_and_keeps_original(self):
        result = converter.combine_trials(trial_data([1, 1, 1, 1, 100]), exclude_outliers=True)
        self.assertEqual(result.loc[0, "excluded_trials"], "5")
        self.assertEqual(result.loc[0, "n_trials_used"], 4)
        self.assertEqual(result.loc[0, "Z_real_avg"], 1)
        self.assertEqual(result.loc[0, "Z_imag_avg"], -1)
        self.assertEqual(result.loc[0, "Z_real_std"], 0)
        self.assertEqual(result.loc[0, "Z_real_trial5"], 100)

    def test_conflicting_components_require_review_without_exclusion(self):
        result = converter.combine_trials(
            trial_data([1, 1, 1, 1, 100], [-100, -1, -1, -1, -1]),
            exclude_outliers=True,
        )
        self.assertTrue(result.loc[0, "outlier_review_required"])
        self.assertTrue(result.loc[0, "has_outliers"])
        self.assertEqual(result.loc[0, "n_trials_used"], 5)
        self.assertEqual(result.loc[0, "excluded_trials"], "")

    def test_two_to_four_trials_are_averaged_without_detection(self):
        for n in (2, 3, 4):
            with self.subTest(n=n):
                result = converter.combine_trials(trial_data(list(range(n))), exclude_outliers=True)
                self.assertFalse(result.loc[0, "has_outliers"])
                self.assertEqual(result.loc[0, "n_trials_used"], n)
                self.assertAlmostEqual(result.loc[0, "Z_real_std"], np.std(range(n), ddof=1))

    def test_holm_correction_is_across_both_components(self):
        np.testing.assert_array_equal(
            converter.holm_rejections(np.array([[0.01, 0.04], [0.03, 0.5]]), 0.05),
            [[True, False], [False, False]],
        )

    def test_mismatched_frequency_grids_fail(self):
        data = trial_data([1, 2])
        data.loc[1, converter.FREQUENCY] = 10
        with self.assertRaisesRegex(ValueError, "mismatched frequency"):
            converter.combine_trials(data)

    def test_conversion_paths_schema_and_overwrite_protection(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "raw" / "500C" / "test.dat"
            write_dat(source)
            original = source.read_bytes()
            kwargs = dict(smart_root=root / "raw", output_root=root / "new-output")
            destination = converter.convert_dat_to_csv("500C/test.dat", **kwargs)
            self.assertEqual(destination, root / "new-output/500C/combined_impedance_data_test.csv")
            result = pd.read_csv(destination)
            self.assertEqual(result.shape, (2, 3))
            self.assertEqual(result[converter.FREQUENCY].tolist(), [100, 1000])
            self.assertEqual(result["Z_real_trial1"].tolist(), [11, 10])
            with self.assertRaises(FileExistsError):
                converter.convert_dat_to_csv("500C/test.dat", **kwargs)
            converter.convert_dat_to_csv("500C/test.dat", overwrite=True, **kwargs)
            self.assertEqual(source.read_bytes(), original)

    def test_cli_defaults_and_overrides_ignore_toml(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_dat(root / "raw/test.DAT", 5)
            (root / "config").mkdir()
            (root / "config/paths.local.toml").write_text("invalid TOML!!!")
            with patch("os.getcwd", return_value=str(root)):
                status = converter.main([
                    "--smart-root", str(root / "raw"),
                    "--output-root", str(root / "override"),
                    "--no-handle-outliers",
                ], output_root=root / "unused", handle_outliers=True)
            self.assertEqual(status, 0)
            csv = pd.read_csv(root / "override/combined_impedance_data_test.csv")
            self.assertNotIn("has_outliers", csv)
            status = converter.main([
                "--handle-outliers", "--overwrite",
            ], smart_root=root / "raw", output_root=root / "override")
            self.assertEqual(status, 0)
            self.assertIn("has_outliers", pd.read_csv(root / "override/combined_impedance_data_test.csv"))


if __name__ == "__main__":
    unittest.main()
