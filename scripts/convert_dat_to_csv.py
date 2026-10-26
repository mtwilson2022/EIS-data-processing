r"""Convert Solartron SMaRT measurement files (.dat) into combined CSV tables.

WHAT THIS SCRIPT DOES
    Reads the repeated measurement sweeps (trials) in each native SMaRT file
    and creates one CSV per input file. Each row is one frequency. Columns
    contain every trial's real and imaginary impedance, their averages,
    sample standard deviations, and any outlier flags. Impedance is in ohms.

    This file is the command-line wrapper. The conversion code lives in
    src/eis_processing/smart_dat_to_csv.py. Keep both files in the project.

FIRST-TIME SETUP
    Use Python 3.11 or newer. Open a terminal in the project root: the folder
    containing pyproject.toml, src/, and scripts/. Select/activate your Python
    environment, then install the project and its dependencies:

        python -m pip install -e .

    This makes the eis_processing module available to the wrapper. Copying
    only this wrapper to another computer is not sufficient.

QUICK START
    The examples below assume this wrapper is named
    scripts/smart_dat_to_csv_standalone.py and the terminal is at project root.
    Put paths containing spaces in quotes. Replace the example paths with
    actual files/folders on your computer.

    1. Convert one file, specifying both paths explicitly:

        python scripts/smart_dat_to_csv_standalone.py "C:/measurements/test.dat" -o "C:/results/test.csv"

    2. Use separate base folders for input and output:

        python scripts/smart_dat_to_csv_standalone.py "500C/test.dat" --smart-root "C:/measurements" --box-root "C:/results"

       Reads C:/measurements/500C/test.dat. Writes
       C:/results/500C/combined_impedance_data_test.csv.

    3. Convert all DAT files in an input folder, including its subfolders:

        python scripts/smart_dat_to_csv_standalone.py . --smart-root "C:/measurements" --box-root "C:/results" --recursive

       Each input file gets its own CSV; different files are not averaged
       together. Processing a batch also avoids repeated package startup.

    4. Replace a CSV that already exists:

        python scripts/smart_dat_to_csv_standalone.py "C:/measurements/test.dat" -o "C:/results/test.csv" --overwrite

HOW PATHS WORK
    --smart-root is the base folder for relative input paths.
    --box-root is the base folder for relative output paths. It can be any
    output folder; using Box is optional.

    Omit either option to use the source module's configured default for that
    root. Explicit root arguments take precedence over optional local config,
    which takes precedence over built-in defaults. The TOML settings, when
    used, are [paths].smart_root and [paths].box_root.

    Absolute input/output paths are used directly. Relative root overrides
    resolve from the terminal's current folder. With no input argument, the
    script scans the selected raw-data root. --recursive includes subfolders.

    Without -o, the filename is combined_impedance_data_<input-stem>.csv,
    preserving the input's subfolders relative to the raw-data root under the
    output root. An absolute input outside the raw-data root writes directly
    into the output root. The stem is the input filename without .dat.

    The output root must exist when used. Output subfolders are created as
    needed. -o accepts one input DAT file, not a folder or multiple files.

OUTLIERS AND STATISTICS
    Every original trial value is preserved in the CSV and the input DAT
    file is never changed. Outlier exclusions affect only means and standard
    deviations at the corresponding frequency.

    With at least five trials, the default uses two-sided Grubbs tests with
    Holm correction across real/imaginary components and all frequencies in
    that input file (significance level 0.05). At most one complex trial point
    is excluded per frequency, using the same mask for real and imaginary
    values. There is no repeated removal. If components flag different trials,
    all trials are retained and the frequency is marked for review.

    Sample standard deviation (ddof=1) is always used, even with filtering
    disabled. With only one trial, standard deviation is undefined and blank.
    Filtered variability describes the retained data, not the full experiment.
    The statistical test assumes approximately normal repeat measurements;
    a flag alone does not prove an instrument error.

    Audit columns:
        has_outliers             A statistically significant point was found.
        excluded_trials          Trial ID actually excluded; blank if none.
        n_trials_used            Trial count used for both component statistics.
        outlier_review_required  Components disagree; no automatic exclusion.

OPTIONS
    -h, --help             Show the source module's command-line help.
    -o, --output PATH      Choose the output CSV for one input file.
    --smart-root PATH      Override the raw SMaRT input base folder.
    --box-root PATH        Override the output base folder.
    --recursive            Include subfolders when searching for DAT files.
    --overwrite            Allow existing output CSVs to be replaced.
    --no-outlier-filter    Include all trial points in the summary statistics.
    --outlier-alpha VALUE  Change the significance level (default 0.05).

    Save the command/settings and code version with an analysis, especially
    when changing outlier options. The CSV records exclusions, not the alpha
    or whether filtering was disabled.

IF SOMETHING GOES WRONG
    "No module named eis_processing": install the project using the setup
    command above, with the same Python environment used to run this wrapper.

    Input not found: check the resolved input path, raw-data root, and filename.
    Existing output: use --overwrite if you intend to replace it.
    Unsupported DAT layout or mismatched trial frequencies: inspect the input;
    the reader supports the verified native SMaRT format, not arbitrary .dat
    files. It does not silently interpolate or discard inconsistent points.
"""

from eis_processing.smart_dat_to_csv import main


if __name__ == "__main__":
    raise SystemExit(main())
