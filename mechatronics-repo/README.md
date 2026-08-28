# Mechatronics Class Repo

MATLAB/Simulink coursework and projects.

## Structure

```
├── models/     # .slx / .mdl Simulink models
├── src/        # .m scripts and functions
├── data/       # input datasets, sensor logs, etc.
├── results/    # exported plots, output data, generated reports
├── docs/       # lab writeups, notes, assignment PDFs
└── assets/     # images, diagrams used in docs/README
```

## Setup (MATLAB)

1. Clone the repo, then open MATLAB and `cd` into the project root.
2. Add the project to your MATLAB path:
   ```matlab
   addpath(genpath(pwd));
   savepath; % optional, persists across sessions
   ```
3. Open models from `models/`, scripts from `src/`.

## Conventions

- Name files by lab/assignment, e.g. `lab3_pid_control.slx`, `lab3_analysis.m`.
- Keep raw data in `data/` and generated plots/output in `results/` — don't commit large generated files unless needed for grading.
- One commit per meaningful change; write commit messages that describe *what* changed (e.g. `Add PID tuning script for Lab 4` not `update`).
