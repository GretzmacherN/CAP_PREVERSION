# Environment and setup (Windows)

This guide describes the current `CAP_Analysis` layout and the Windows setup used for the reported run. Keep the repository, external data root and individual run outputs separate.

## Verification status

- **Confirmed:** notebook 005 (`05_three_proteome`) completed successfully in the Windows `cap-rerun` environment, using the same `CAP_RUN_ID` as the other stages.
- **Windows R setup recorded:** R 4.6.2 was made available to the Python process through `R_HOME` and `PATH`; `SingleCellExperiment` was installed in the R user library visible to `rpy2`.
- **Not yet confirmed:** a fresh Windows environment installation and a complete run of all R-dependent single-cell notebooks using the new CAP paths.
- The checked-in `environment.yml` and `environment_R.json` are records from Apple Silicon/macOS. They are not a validated Windows lockfile or Windows R package inventory.

## Repository, data and local configuration

The Python project root is `CAP_Analysis/`, identified by its `pyproject.toml`. Large inputs and run outputs are kept outside the repository in `CAP_DATA_ROOT`. The path helper reads `CAP_DATA_ROOT` and `CAP_RUN_ID` from the process environment or the local `.env` file beside `pyproject.toml`.

Keep `.env` local and Git-ignored. Commit `.env.example` with placeholders only. For example, a local `.env` can contain:

```text
CAP_DATA_ROOT=C:/path/to/CAP_Data
CAP_RUN_ID=your_run_id
```

Use the same run ID for all stages in one run. Choose a new ID for a new run; a non-empty stage is not an overwrite target. Do not commit machine-specific paths or data under `CAP_DATA_ROOT`.

## Use the Windows Conda environment

The already-used environment is named `cap-rerun`. From PowerShell, activate it, install the current repository in editable mode, and register the kernel:

```powershell
conda activate cap-rerun
Set-Location C:\path\to\CAP_Analysis
python -m pip install -e .
python -m ipykernel install --user --name cap-rerun --display-name "Python (cap-rerun)"
```

Start Jupyter or the IDE from that environment. Restart any existing notebook kernel after changing the package or environment. Because the current Python YAML was captured on macOS, do not treat a fresh `conda env create -f environment.yml` on Windows as validated until it has been tested separately.

Check that the kernel imports the intended checkout and resolves paths correctly:

```powershell
python -c "import cap.paths as p; x=p.get_paths(); print(p.__file__); print(x.root); print(x.data_root); print(x.run_id); print(x.stage_dir('05_three_proteome'))"
python -m pip check
```

The repository root printed by the first command should be `CAP_Analysis`, while `data_root` should be the external `CAP_Data` directory.

## Make R available to `rpy2` on Windows

Set `R_HOME` and add the matching R executable directory to `PATH` in the same PowerShell session that starts Jupyter or the IDE. Adjust the versioned directory if R is installed elsewhere:

```powershell
$env:R_HOME = 'C:\Program Files\R\R-4.6.2'
$env:Path = "$env:R_HOME\bin\x64;$env:Path"
Test-Path "$env:R_HOME\bin\x64\R.exe"
```

`Test-Path` should return `True`. If R packages are installed in the user library, also set `R_LIBS_USER` so the R session started by `rpy2` can find them:

```powershell
$env:R_LIBS_USER = "$env:USERPROFILE\Documents\R\win-library\4.6"
```

Start Jupyter from this same PowerShell session, or configure these variables in the IDE's kernel settings. A server already running will not inherit variables added later in another terminal. For troubleshooting and verification inside a notebook, see `R on Windows - Notebook Startup.md`.

Verify R from the active notebook kernel:

```powershell
python -c "import rpy2.robjects as ro; print(ro.r('R.version.string')[0])"
```

Install Bioconductor packages into a library visible to this R installation. The Windows troubleshooting record confirms `SingleCellExperiment` plus the `R_HOME`/`PATH` setup; record the actual Windows versions of `miQC`, `scDblFinder`, `scry` and ALRA separately rather than treating the macOS `environment_R.json` as a Windows record.

## Compute accelerator

The recorded Apple Silicon notebooks used `accelerator="mps"`. MPS is not available on Windows; use a supported Windows accelerator such as CUDA when available, or CPU, and verify that the relevant training notebook completes. Successful execution of the Python-only three-proteome notebook does not validate scVI training or the full single-cell workflow on Windows.

## Environment files

| File | Purpose |
|---|---|
| `environment.yml` | Direct Python dependency record captured on Apple Silicon/macOS; not a Windows lockfile |
| `environment_R.json` | R package inventory from the recorded macOS environment |
| `.env.example` | Shareable template with placeholder data root and run ID |
| `.env` | Local ignored configuration; never commit it |
