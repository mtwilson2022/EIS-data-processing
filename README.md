# EIS Data Processing

Python tools for converting Solartron SMaRT impedance measurements into CSV
tables for plotting and analysis. Raw DAT files remain unchanged. Each DAT file
produces its own CSV; trials from different files are never averaged together.

## Project layout

| Location | Purpose |
| --- | --- |
| `src/eis_processing/` | Reusable readers, processing, and analysis code. |
| `src/eis_processing/smart_dat_to_csv.py` | Native DAT reader, trial statistics, outlier handling, and converter CLI. |
| `scripts/smart_dat_to_csv_standalone.py` | Converter entry point; edit the settings near the top here. |
| `tests/` | Automated checks, if included in your checkout. |
| `pyproject.toml` | Package metadata and dependencies. |

The wrapper depends on the installed package despite its historical
`standalone` filename. Keep both Python files in their indicated folders.

## Installation

Use Python 3.11 or newer. Open a terminal in the project root, the folder
containing `pyproject.toml`, `src/`, and `scripts/`.

On Windows (PowerShell):

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -e .
```

On macOS/Linux:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -e .
```

In VS Code, select this environment as the Python interpreter. The examples
below use `python`; use `.\.venv\Scripts\python.exe` on Windows or
`.venv/bin/python` on macOS/Linux if the environment is not activated.
The converter uses NumPy and pandas, plus SciPy for outlier screening.
Matplotlib is used by the project's plotting tools. Project installation uses
the dependencies declared in your existing `pyproject.toml`.

## Running the file converter

### 1. Set your folders and options

Open `scripts/smart_dat_to_csv_standalone.py` and find
`EDIT YOUR SETTINGS HERE` near the top. For example:

```python
SMART_ROOT = Path("D:/SMaRT measurements")
OUTPUT_ROOT = Path("C:/Users/yourname/results")
HANDLE_OUTLIERS = False
OUTLIER_ALPHA = 0.05
RECURSIVE = False
OVERWRITE = False
```

`SMART_ROOT` is the folder containing raw data. `OUTPUT_ROOT` is any destination
folder, including a synced folder if desired. It is created if needed. Use
forward slashes for Windows paths, or raw strings such as `Path(r"D:\SMaRT")`.
Relative settings resolve from the terminal's working directory. The supplied
defaults read the current folder and write into its `csv/` subfolder.

No configuration files are read. You do not need a TOML settings file.

### 2. Run a conversion

Run these commands from the project root. Quote paths containing spaces.

Convert DAT files directly in `SMART_ROOT`:

```powershell
python scripts/smart_dat_to_csv_standalone.py
```

Convert one file relative to `SMART_ROOT`:

```powershell
python scripts/smart_dat_to_csv_standalone.py "500C/measurement.dat"
```

Convert a folder and all its subfolders:

```powershell
python scripts/smart_dat_to_csv_standalone.py . --recursive
```

Override the roots for one run without editing the wrapper:

```powershell
python scripts/smart_dat_to_csv_standalone.py "500C/measurement.dat" --smart-root "D:/SMaRT" --output-root "C:/results"
```

Enable outlier handling for one run:

```powershell
python scripts/smart_dat_to_csv_standalone.py . --recursive --handle-outliers
```

Choose an exact output filename for one DAT input:

```powershell
python scripts/smart_dat_to_csv_standalone.py "D:/SMaRT/test.dat" -o "C:/results/test.csv"
```

To replace existing CSVs, add `--overwrite`. Otherwise existing files are
protected and reported as failures. Multiple files or folders can be supplied
in one command. Batch conversion avoids repeating Python's package startup for
each file. `--output` accepts exactly one input DAT file.

### 3. Find the results

The default filename is `combined_impedance_data_<input-stem>.csv`. For example:

```text
Input:  D:/SMaRT/500C/measurement.dat
Output: C:/results/500C/combined_impedance_data_measurement.csv
```

Relative inputs resolve under `SMART_ROOT`. Relative `--output` paths resolve
under `OUTPUT_ROOT`. Absolute input/output paths are used directly. Without
`--output`, subfolders beneath `SMART_ROOT` are preserved beneath `OUTPUT_ROOT`;
an absolute input outside `SMART_ROOT` writes directly into `OUTPUT_ROOT`.
Relative command-line root overrides resolve from the working directory.

### Available options

Command-line values override wrapper settings for that run.

| Option | Effect |
| --- | --- |
| `--smart-root PATH` | Raw input base folder. |
| `--output-root PATH` | CSV output base folder. |
| `-o PATH`, `--output PATH` | Exact CSV destination for one input DAT file. |
| `--recursive` / `--no-recursive` | Include / skip input subfolders. |
| `--handle-outliers` / `--no-handle-outliers` | Enable / disable outlier handling. |
| `--outlier-alpha VALUE` | Significance level between 0 and 1; default 0.05. |
| `--overwrite` / `--no-overwrite` | Allow / prevent replacing existing CSVs. |
| `-h`, `--help` | Show command-line help. |

### CSV columns and statistics

Every CSV contains one row per frequency in ascending order. Impedance values
are in ohms; the original imaginary sign is preserved. Original trial values
are retained even if excluded from summary statistics.

| Input and selected mode | Output columns |
| --- | --- |
| One trial, either mode | `Frequency (Hz)`, `Z_real_trial1`, `Z_imag_trial1` only. |
| Multiple trials, outlier handling off | Frequency; all `Z_real_trialN` and `Z_imag_trialN` columns; `Z_real_avg`, `Z_imag_avg`, `Z_real_std`, `Z_imag_std`. |
| Multiple trials, outlier handling on | All the above plus `has_outliers`, `excluded_trials`, `n_trials_used`, `outlier_review_required`. |

With one trial, averaging, standard deviations, and outlier screening are
skipped entirely. With multiple trials and outlier handling off (the default),
all trials contribute to ordinary arithmetic means and sample standard
deviations (`ddof=1`). No audit columns are added.

With outlier handling on, the existing screening method is used:

- At least five trials are required for detection. With two to four trials,
  all trials are retained and audit columns show no detections.
- Two-sided Grubbs tests select one candidate per component and frequency.
  Holm correction covers both components and all frequencies within that file.
- At most one complex trial point is excluded at any frequency. The same
  retained-trial mask is used for real and imaginary statistics.
- If significant real and imaginary tests identify different trials, all
  trials at that frequency are retained and `outlier_review_required` is true.
- Exclusions are not repeated iteratively. `excluded_trials` lists the native
  trial ID excluded at that frequency; it is blank when none is excluded.

Screening assumes approximately normal repeat measurements. A detection alone
does not establish instrument error. Filtered standard deviations describe the
retained measurements. Save the options and code version with your analysis;
the CSV does not record the selected alpha.

Frequency grids must match exactly across trials. Missing points, duplicate
frequencies, nonfinite values, and unsupported native layouts raise errors;
the converter does not interpolate measurements.

## Calling the converter from Python

```python
from eis_processing.smart_dat_to_csv import convert_dat_to_csv

csv_path = convert_dat_to_csv(
    "500C/measurement.dat",
    smart_root="D:/SMaRT",
    output_root="C:/results",
    exclude_outliers=False,  # True enables screening and audit columns.
    overwrite=False,
)
```

Python calls use the roots passed to them, defaulting to the working directory.
Settings in the wrapper affect wrapper runs only.

## Plotting and further analysis

Converted CSVs can be used for impedance plots and subsequent conductivity or
Arrhenius analysis. Multi-trial average columns keep their existing names.
Code that reads averages unconditionally must handle single-trial files by
using `Z_real_trial1` and `Z_imag_trial1` instead. Code that reads audit columns
must account for their absence when outlier handling is off.

## Troubleshooting and migration

- **`No module named eis_processing`:** run `python -m pip install -e .` from
  the project root using the same environment that runs the wrapper.
- **Installation says `src` or `README.md` is missing:** check the folder
  layout and that you are installing the project root rather than `scripts/`.
- **Input not found:** check `SMART_ROOT` and the relative input path, or use
  an absolute input path. Use `--recursive` to include nested folders.
- **Existing output:** add `--overwrite` if replacement is intended.
- **Unsupported DAT layout:** this reader handles the verified native
  `SMaRT_Data` version-2 binary layout, not every file with a `.dat` extension.
- **Old commands:** replace `--box-root` with `--output-root`, and replace
  `--no-outlier-filter` with `--no-handle-outliers`. Screening now defaults off;
  use `--handle-outliers` to enable it.
- **Old configuration files:** this converter ignores `config/paths.local.toml`.
  Set paths in the wrapper or command line. Other project scripts may still
  depend on older configuration code; check their imports before deleting it.

To apply this update, replace the converter module and wrapper in their listed
folders and put this README at the project root. Keep your existing
`pyproject.toml`, package initializer, and other analysis modules.

## Converter checks

Run the included focused checks after installing the project:

```powershell
python -m unittest discover -s tests -p test_smart_converter.py
```
