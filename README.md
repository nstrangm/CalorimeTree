# CalorimeTree

Analysis framework for isolated photons, neutral pions, jets, and correlation observables from ALICE tree-based inputs (AliPhysics/O2Physics style).

## Quick start

1. Edit one of the run configs (for example `RunConfig.yaml`) so dataset paths match your environment.
2. Run the main driver:
   - `python3 start.py --config RunConfig.yaml`
3. Enable workflow steps in the config (`doMakeHistos`, `doPurity`, `doPlotting`, ...) as needed.

## Available flags at a glance

| Flag | Macro/script run by `start.py` | Minimal description |
|---|---|---|
| `doMakeHistos` | `Analysis/makeHistosFromTree.C` | Produces `HistosFromTree*.root` from input tree files. |
| `doPlotting` | `Analysis/plotHistosFromTree.C` | Creates standard QA/physics plots from histogram outputs. |
| `doAnalysisExclGammaJet` | `Analysis/analyseExclGammaJet.C` | Runs exclusive gamma-jet analysis from produced histograms. |
| `doPlottingExclGammaJet` | `Analysis/plotExclGammaJet.C` | Plots and summarizes exclusive gamma-jet analysis results. |
| `doCombineExclGammaJet` | `Analysis/combineExclGammaJet.C` | Merges exclusive gamma-jet outputs across datasets/cuts. |
| `doPurity` | `Analysis/calculatePurity.C` | Computes purity-related quantities and derived histograms. |
| `doPlottingPurity` | `Analysis/plotPurity.C` | Produces plots for purity study outputs. |
| `doIsoGamma` | *(no standalone macro; used inside `makeHistosFromTree.C`)* | Enables isolated-photon observables during histogram production. |
| `doJets` | *(no standalone macro; used inside `makeHistosFromTree.C`)* | Enables jet-related observables during histogram production. |
| `doGGPi0` | *(no standalone macro; used inside `makeHistosFromTree.C`)* | Enables gamma-gamma neutral-pion channel filling. |
| `domPi0` | *(no standalone macro; used inside `makeHistosFromTree.C`)* | Enables invariant-mass neutral-pion channel filling. |
| `doSubstructure` | *(no standalone macro; used inside `makeHistosFromTree.C`)* | Enables jet-substructure observable filling. |

## Main entry points

- `start.py`: orchestrates ROOT macro execution and optional SLURM submission.
- `Analysis/`: core C++ analysis and plotting macros plus shared headers.
- `Cuts.yaml`: central cut definitions used by analysis macros.
- `RunConfig*.yaml`: predefined campaign-specific run configurations.
- `ML/`: preparation, training, and evaluation scripts/macros for ML-based studies.
- `tools/`: conversion, download, job recovery, and full-chain helper scripts.
- `DummyDataSet/`: minimal example input structure for local testing.

## Typical workflow

1. Prepare or convert input trees (helper scripts in `tools/`).
2. Configure datasets/train configs in a `RunConfig*.yaml` file.
3. Run `start.py` to produce histograms and optional derived outputs.
4. Post-process with plotting/purity macros in `Analysis/` or ML scripts in `ML/`.