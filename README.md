# sni-simulation

Reproducible simulation studies in R, managed with [renv](https://rstudio.github.io/renv/).

## Getting Started

### RStudio

1. Open `sni-simulation.Rproj` — renv activates automatically via `.Rprofile`.
2. Restore packages:

   ```r
   renv::restore()
   ```

3. Run a simulation:

   ```bash
   Rscript run_simulation.R simulations/example_sim
   ```

### VS Code

1. Install the [R extension](https://marketplace.visualstudio.com/items?itemName=REditorSupport.r) (`REditorSupport.r`).

2. Make sure `R` is on your `PATH`, or point VS Code to your R installation.
   Open **Settings** (`Ctrl+,`) and set:

   ```json
   "r.rterm.windows": "C:\\Program Files\\R\\R-4.5.2\\bin\\R.exe"
   ```

3. Open this project folder in VS Code. The `.Rprofile` will bootstrap renv
   when you launch an R terminal (`Ctrl+Shift+`` → R Terminal).

4. In the R terminal, restore the environment:

   ```r
   renv::restore()
   ```

5. Run a simulation from the VS Code integrated terminal:

   ```bash
   Rscript run_simulation.R simulations/example_sim
   ```

> **Tip:** If renv doesn't activate automatically, run `source("renv/activate.R")`
> in your R terminal, or ensure `.Rprofile` is not being overridden by a
> user-level profile (`~/.Rprofile`).

## Project Structure

```
R/                        Shared helper functions
simulations/              One subfolder per simulation study
  example_sim/
    config.R              Parameters (sample size, replications, DGP settings)
    dgp.R                 Data-generating process
    estimator.R           Estimator(s) under evaluation
    run.R                 Orchestration script
    results/              Output (gitignored)
run_simulation.R          CLI entry point
```

## Adding a New Simulation

1. Create a new folder under `simulations/`, e.g. `simulations/my_sim/`.
2. Add `config.R`, `dgp.R`, `estimator.R`, and `run.R` following the pattern in `example_sim`.
3. Run it:

   ```bash
   Rscript run_simulation.R simulations/my_sim
   ```

## Reproducibility

- `renv.lock` pins all package versions. Run `renv::restore()` on a fresh clone.
- Each simulation's `config.R` contains the random seed.
- Results are saved as both `.rds` and `.csv`.
