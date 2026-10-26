r"""Reusable Solartron SMaRT DAT -> combined impedance CSV conversion library.

This module exposes the functions used by the CLI wrapper and is importable from
packages installed under the src/ layout.

Requires Python 3.11+ and NumPy, pandas, SciPy:
    python -m pip install numpy pandas scipy

Version: 2026-10-26
"""

import argparse
import math
from pathlib import Path
import struct
import sys

import numpy as np
import pandas as pd


# Native file reader

HEADER_SIZE = 66

RECORD = struct.Struct("<HHI9dHH")

FREQUENCY = "Frequency (Hz)"
REAL = "Z_real (Ω)"
IMAG = "Z_imag (Ω)"
SWEEP = "Sweep Number"
POINT = "Point Number"


class SmartFormatError(ValueError):
    """The input is damaged or does not match the supported native layout."""


def read_zplot_ascii(path: str | Path) -> pd.DataFrame:
    """Read one stored spectrum from a ZPLOT2 ASCII export."""
    path = Path(path)
    lines = path.read_text(encoding="utf-8-sig").splitlines()

    if not lines or lines[0].strip() != "ZPLOT2 ASCII":
        raise SmartFormatError(f"{path.name}: not a ZPLOT2 ASCII file.")

    try:
        end = next(
            i for i, line in enumerate(lines)
            if line.strip() == "End Comments"
        )
        columns = next(
            line.split() for line in lines[:end]
            if line.strip().startswith("Freq(Hz)")
        )
        indices = [
            columns.index(name)
            for name in ("Freq(Hz)", "Z'(a)", "Z''(b)")
        ]
        expected = int(next(
            line.split(":", 1)[1] for line in lines[:end]
            if line.strip().startswith("Data Points:")
        ))
    except (StopIteration, ValueError) as error:
        raise SmartFormatError(
            f"{path.name}: unsupported or incomplete ASCII header."
        ) from error

    rows = []
    for number, line in enumerate(lines[end + 1:], start=end + 2):
        if not line.strip():
            continue

        fields = line.split()
        try:
            if len(fields) != len(columns):
                raise ValueError("incorrect number of columns")

            values = [float(field) for field in fields]
            frequency, real, imag = [values[i] for i in indices]

            if (
                not all(math.isfinite(x) for x in (frequency, real, imag))
                or frequency <= 0
            ):
                raise ValueError("invalid frequency or impedance")

        except ValueError as error:
            raise SmartFormatError(
                f"{path.name}: invalid ASCII data on line {number}: {error}"
            ) from error

        rows.append((1, len(rows) + 1, frequency, real, imag))

    if not rows or len(rows) != expected:
        raise SmartFormatError(
            f"{path.name}: expected {expected} ASCII points; "
            f"found {len(rows)}."
        )

    frame = pd.DataFrame(
        rows, columns=[SWEEP, POINT, FREQUENCY, REAL, IMAG]
    )
    if frame[FREQUENCY].duplicated().any():
        raise SmartFormatError(
            f"{path.name}: repeated frequencies; "
            "multiple ASCII spectra are not supported."
        )

    frame.attrs["source_file"] = str(path.resolve())
    return frame


def read_smart_dat(path: str | Path) -> pd.DataFrame:
    """Read stored components in ohms, with one-based sweep/point identifiers.

    The third prefix field, ancillary floating-point fields, and range flags
    are intentionally not interpreted. Negative imaginary impedance is kept.
    """
    path = Path(path)
    data = path.read_bytes()

    if data.startswith((b"ZPLOT2 ASCII", b"\xef\xbb\xbfZPLOT2 ASCII")):
        return read_zplot_ascii(path)

    if data[:32].rstrip(b" ") != b"SMaRT_Data":
        raise SmartFormatError(f"{path.name}: not a native SMaRT_Data file.")
    if data[32:48].strip() != b"2":
        raise SmartFormatError(f"{path.name}: unsupported SMaRT version marker.")
    if len(data) < HEADER_SIZE + RECORD.size + 4 or not data.endswith(b"SDEM"):
        raise SmartFormatError(f"{path.name}: empty/truncated file or missing SDEM marker.")
    payload = memoryview(data)[HEADER_SIZE:-4]
    if len(payload) % RECORD.size:
        raise SmartFormatError(f"{path.name}: unsupported layout or truncated 84-byte record.")

    rows = []
    for number, record in enumerate(RECORD.iter_unpack(payload), start=1):
        # Prefix: sweep, point, unresolved uint32. Then nine doubles:
        # time, frequency, AC, DC, three ancillary values, Z_real, Z_imag.
        sweep, point = record[:2]
        frequency, real, imag = record[4], record[10], record[11]
        if record[12] != 2:
            raise SmartFormatError(f"{path.name}: unsupported measurement tag in record {number}.")
        if not all(math.isfinite(x) for x in (frequency, real, imag)) or frequency <= 0:
            raise SmartFormatError(f"{path.name}: invalid frequency/impedance in record {number}.")
        rows.append((sweep + 1, point + 1, frequency, real, imag))

    frame = pd.DataFrame(rows, columns=[SWEEP, POINT, FREQUENCY, REAL, IMAG])
    sweeps = sorted(frame[SWEEP].unique())
    if sweeps != list(range(1, len(sweeps) + 1)):
        raise SmartFormatError(f"{path.name}: unsupported or incomplete native sweep indices.")
    for sweep, trial in frame.groupby(SWEEP):
        if sorted(trial[POINT].tolist()) != list(range(1, len(trial) + 1)):
            raise SmartFormatError(f"{path.name}: missing/duplicate point indices in sweep {sweep}.")
    frame.attrs["source_file"] = str(path.resolve())
    return frame


# Outlier screening

MIN_TRIALS = 5


def grubbs_candidate(values: np.ndarray) -> tuple[int, float]:
    """Return the most extreme index and a two-sided Grubbs p-value.

    The 2*n factor accounts for selecting the most extreme of n measurements.
    Inputs are finite scalar repetitions, assumed approximately normal under
    the null. Constant data has no candidate. The formula inverts NIST's
    two-sided Grubbs critical value, with n-2 Student-t degrees of freedom.
    """
    values = np.asarray(values, dtype=float)
    n = len(values)
    if n < MIN_TRIALS:
        return -1, 1.0
    sd = values.std(ddof=1)
    if sd == 0:
        return -1, 1.0
    deviations = np.abs(values - values.mean())
    index = int(np.argmax(deviations))
    g = deviations[index] / sd
    denominator = (n - 1) ** 2 - n * g ** 2
    if denominator <= 0:
        # The limiting case: n-1 identical values and one different value.
        return index, 0.0
    from scipy.stats import t

    statistic = np.sqrt(n * (n - 2) * g ** 2 / denominator)
    return index, float(min(1.0, 2 * n * t.sf(statistic, df=n - 2)))


def holm_rejections(pvalues: np.ndarray, alpha: float) -> np.ndarray:
    """Holm step-down family-wise rejection mask; valid with dependent tests."""
    pvalues = np.asarray(pvalues, dtype=float)
    rejected = np.zeros(pvalues.shape, dtype=bool)
    flat = pvalues.ravel()
    for rank, index in enumerate(np.argsort(flat, kind="stable")):
        if flat[index] > alpha / (len(flat) - rank):
            break
        rejected.flat[index] = True
    return rejected


def paired_outlier_mask(
    real: np.ndarray, imag: np.ndarray, *, enabled: bool = True, alpha: float = 0.05
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return retained-pair mask, detected-frequency flags, and review flags.

    Arrays are frequency x trial, already validated and aligned. At most one
    complex point is excluded per row. If significant real/imaginary tests
    identify different trials, retain all and mark review_required instead.
    Fewer than five trials disables automatic detection, not averaging.
    """
    if not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("outlier_alpha must be between 0 and 1.")
    retained = np.ones(real.shape, dtype=bool)
    detected = np.zeros(real.shape[0], dtype=bool)
    review = np.zeros(real.shape[0], dtype=bool)
    if not enabled or real.shape[1] < MIN_TRIALS:
        return retained, detected, review
    candidates = np.full((real.shape[0], 2), -1, dtype=int)
    pvalues = np.ones(candidates.shape, dtype=float)
    for row in range(real.shape[0]):
        for component, values in enumerate((real[row], imag[row])):
            candidates[row, component], pvalues[row, component] = grubbs_candidate(values)
    significant = holm_rejections(pvalues, alpha)
    for row in range(real.shape[0]):
        flagged = set(candidates[row, significant[row]]) - {-1}
        detected[row] = bool(flagged)
        if len(flagged) == 1:
            retained[row, flagged.pop()] = False
        elif len(flagged) > 1:
            review[row] = True
    return retained, detected, review


# Trial statistics and CSV writer

DEFAULT_SMART_ROOT = Path(".")
DEFAULT_OUTPUT_ROOT = Path(".")


def resolve_input_path(path: str | Path, *, smart_root: str | Path | None = None) -> Path:
    """Resolve a relative input beneath the supplied raw SMaRT directory."""
    root = Path(DEFAULT_SMART_ROOT if smart_root is None else smart_root).expanduser().resolve()
    path = Path(path).expanduser()
    return (path if path.is_absolute() else root / path).resolve()


def combine_trials(
    data: pd.DataFrame, *, exclude_outliers: bool = False, outlier_alpha: float = 0.05
) -> pd.DataFrame:
    """Preserve trial columns; summarize retained pairs with sample SD (ddof=1).

    Frequency grids must match exactly. No interpolation, frequency exclusion,
    or partial-trial averaging is performed. Native sweep IDs name the columns.
    Optional screening: one Grubbs candidate per component/frequency, Holm
    correction across both components and all frequencies, alpha=0.05 per file.
    Detection needs >=5 trials; at most one pair is excluded at any frequency.
    Conflicting component candidates are retained and marked for review.
    A single trial exports only frequency and raw trial columns. With multiple
    trials, summary statistics are added; audit columns require exclude_outliers.
    """
    required = [SWEEP, FREQUENCY, REAL, IMAG]
    if data.empty or any(column not in data for column in required):
        raise ValueError(f"Nonempty data with columns {required} is required.")
    values = data[required].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Sweep, frequency, and impedance values must be finite.")
    ids = data[SWEEP].to_numpy(dtype=float)
    if (ids < 1).any() or not np.equal(ids, np.floor(ids)).all():
        raise ValueError("Sweep identifiers must be positive integers.")
    if (data[FREQUENCY] <= 0).any():
        raise ValueError("Frequencies must be positive.")
    if data.duplicated([SWEEP, FREQUENCY]).any():
        raise ValueError("Duplicate frequency within a trial; cannot combine unambiguously.")

    real = data.pivot(index=FREQUENCY, columns=SWEEP, values=REAL).sort_index().sort_index(axis=1)
    imag = data.pivot(index=FREQUENCY, columns=SWEEP, values=IMAG).reindex_like(real)
    if real.isna().any().any() or imag.isna().any().any():
        raise ValueError("Trials have mismatched frequency grids or missing points.")

    result = pd.DataFrame({FREQUENCY: real.index.to_numpy()})
    for sweep in real.columns:
        result[f"Z_real_trial{int(sweep)}"] = real[sweep].to_numpy()
    for sweep in imag.columns:
        result[f"Z_imag_trial{int(sweep)}"] = imag[sweep].to_numpy()
    if real.shape[1] == 1:
        return result

    re, im = real.to_numpy(), imag.to_numpy()
    if exclude_outliers:
        retained, detected, review = paired_outlier_mask(re, im, alpha=outlier_alpha)
    else:
        retained = np.ones(re.shape, dtype=bool)
    for label, values in (("real", re), ("imag", im)):
        result[f"Z_{label}_avg"] = [row[mask].mean() for row, mask in zip(values, retained)]
    for label, values in (("real", re), ("imag", im)):
        result[f"Z_{label}_std"] = [
            row[mask].std(ddof=1) if mask.sum() > 1 else np.nan
            for row, mask in zip(values, retained)
        ]
    if exclude_outliers:
        result["has_outliers"] = detected
        result["excluded_trials"] = [
            ";".join(str(int(trial)) for trial in real.columns[~mask]) for mask in retained
        ]
        result["n_trials_used"] = retained.sum(axis=1)
        result["outlier_review_required"] = review
    return result


def convert_dat_to_csv(
    source: str | Path,
    output: str | Path | None = None,
    *,
    overwrite: bool = False,
    smart_root: str | Path | None = None,
    output_root: str | Path | None = None,
    exclude_outliers: bool = False,
    outlier_alpha: float = 0.05,
) -> Path:
    """Convert one input; relative source/output paths use the supplied roots.

    Roots default to the working directory. Without an explicit output filename,
    preserve source subfolders beneath output_root. Inputs outside smart_root
    write directly beneath output_root. Output folders are created as needed.
    Existing CSVs require overwrite=True. No configuration files are read.
    """
    raw_root = Path(DEFAULT_SMART_ROOT if smart_root is None else smart_root).expanduser().resolve()
    output_root = Path(DEFAULT_OUTPUT_ROOT if output_root is None else output_root).expanduser().resolve()
    source = resolve_input_path(source, smart_root=raw_root)
    if source.suffix.lower() != ".dat":
        raise ValueError("The source must be a .dat file.")
    if output is not None:
        destination = Path(output).expanduser()
        if not destination.is_absolute():
            destination = output_root / destination
    else:
        try:
            relative_folder = source.parent.relative_to(raw_root)
        except ValueError:
            relative_folder = Path()
        destination = output_root / relative_folder / f"combined_impedance_data_{source.stem}.csv"
    destination = destination.resolve()
    if destination.suffix.lower() != ".csv":
        raise ValueError("The output must have a .csv extension.")
    if destination.resolve() == source.resolve():
        raise ValueError("The output cannot replace the source file.")
    combined = combine_trials(
        read_smart_dat(source), exclude_outliers=exclude_outliers, outlier_alpha=outlier_alpha
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation protects an existing CSV, including concurrent writes.
    with destination.open("w" if overwrite else "x", encoding="utf-8", newline="") as stream:
        combined.to_csv(stream, index=False)
    return destination


# Command-line interface


def main(
    argv: list[str] | None = None,
    *,
    smart_root: str | Path = DEFAULT_SMART_ROOT,
    output_root: str | Path = DEFAULT_OUTPUT_ROOT,
    handle_outliers: bool = False,
    outlier_alpha: float = 0.05,
    recursive: bool = False,
    overwrite: bool = False,
) -> int:
    """Run the CLI; keyword defaults let the wrapper keep settings in one place."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs", nargs="*", type=Path, help="DAT files/folders relative to --smart-root. Absolute paths are accepted. With no inputs, scan the raw SMaRT root.")
    parser.add_argument("--smart-root", type=Path, default=smart_root, help="Raw input root (default: wrapper setting or working directory).")
    parser.add_argument("--output-root", type=Path, default=output_root, help="Output root (default: wrapper setting or working directory).")
    parser.add_argument("--recursive", action=argparse.BooleanOptionalAction, default=recursive, help="Search input folders and their subfolders.")
    parser.add_argument("-o", "--output", type=Path, help="Output CSV relative to --output-root, or an absolute path (single DAT input only).")
    parser.add_argument("--overwrite", action=argparse.BooleanOptionalAction, default=overwrite, help="Replace existing combined CSVs.")
    parser.add_argument("--handle-outliers", action=argparse.BooleanOptionalAction, default=handle_outliers, help="Enable Grubbs/Holm outlier handling and audit columns (default: wrapper setting; off unless changed).")
    parser.add_argument("--outlier-alpha", type=float, default=outlier_alpha, help="Grubbs/Holm significance level per input file (default: 0.05).")
    args = parser.parse_args(argv)
    if not 0 < args.outlier_alpha < 1:
        parser.error("--outlier-alpha must be between 0 and 1.")
    try:
        inputs = [resolve_input_path(path, smart_root=args.smart_root) for path in (args.inputs or [Path('.')])]
    except (OSError, ValueError) as error:
        parser.error(str(error))
    if args.output is not None:
        if len(inputs) != 1:
            parser.error(f"--output requires one input DAT file; received {len(inputs)} inputs.")
        if inputs[0].suffix.lower() != ".dat":
            parser.error(f"--output requires a DAT file, not a folder or another file type: {inputs[0]}")
        if not inputs[0].is_file():
            parser.error(
                f"Input DAT file was not found or is not a regular file: {inputs[0]}\n"
                "Relative input arguments resolve beneath the SMaRT root. "
                "Check SMART_ROOT at the top of the wrapper, "
                "pass --smart-root, or use an absolute input path."
            )

    sources = []
    seen = set()
    for path in inputs:
        if path.is_dir():
            candidates = sorted(
                (p for p in (path.rglob("*") if args.recursive else path.iterdir())
                 if p.is_file() and p.suffix.lower() == ".dat"),
                key=lambda p: str(p).casefold(),
            )
        elif path.is_file() and path.suffix.lower() == ".dat":
            candidates = [path]
        else:
            parser.error(f"Input is not a DAT file or directory: {path}")
        for candidate in candidates:
            resolved = candidate.resolve()
            if resolved not in seen:
                sources.append(candidate)
                seen.add(resolved)
    if not sources:
        parser.error("No DAT files found. Use --recursive to search subfolders.")

    failures = 0
    for source in sources:
        try:
            output = convert_dat_to_csv(
                source, args.output, overwrite=args.overwrite,
                smart_root=args.smart_root, output_root=args.output_root,
                exclude_outliers=args.handle_outliers, outlier_alpha=args.outlier_alpha,
            )
        except (OSError, ValueError) as error:
            failures += 1
            print(f"FAILED: {source}: {error}", file=sys.stderr)
        else:
            print(f"Saved: {output}")
    print(f"Converted {len(sources) - failures} of {len(sources)} file(s).")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
