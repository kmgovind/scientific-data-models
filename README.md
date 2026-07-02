# Scientific Data Models

## Overview

Many sustainably powered robots are designed to collect scientific information about the world. Examples include environmental monitoring, ocean sampling, weather observation, and ecological surveys.

This repository contains examples of different types of scientific data, along with:

- Data sources and download scripts
- Mathematical models used to describe the data
- Visualization notebooks
- Supporting Julia code

The goal is to provide a collection of scientific data models that can be explored, visualized, and compared.

---

# Getting Started

## Prerequisites

Before using this repository, install the following software:

1. **Julia** – the programming language used throughout this project.
   - https://julialang.org/downloads/

2. **Git** – used to download the repository.
   - https://git-scm.com/downloads

3. **Jupyter Notebook** – used to view and run the visualization notebooks.
   - https://jupyter.org/install

---

# Downloading the Repository

The repository is hosted on GitHub.

Clone the repository using:

```bash
git clone <repository-url>
```

Move into the repository directory:

```bash
cd <repository-name>
```

---

# Setting Up the Julia Environment

This repository includes a `Project.toml` file that specifies all required Julia packages.

Start Julia from the repository root:

```bash
julia --project
```

You should see a prompt similar to:

```julia
julia>
```

Press `]` to enter the Julia package manager. The prompt will change to something like:

```julia
(@v1.x) pkg>
```

Run:

```julia
instantiate
```

The first installation may take several minutes while Julia downloads and installs all required packages.

When installation is complete, press **Backspace** (or **Ctrl+C**) to return to the standard Julia prompt.

> **Note:** You only need to run `instantiate` once after cloning the repository (or whenever the project dependencies change).

---

# Opening the Visualization Notebooks
## VS Code
Ensure that the **Jupyter** extension is installed in **VS Code**. This allows `.ipynb` notebooks to open and run directly within VS Code. When opening a notebook, make sure the **Julia kernel** is selected from the kernel picker in the upper-right corner.

## Using Jupyter in the browser
Each data type includes one or more Jupyter notebooks (`.ipynb`) containing:

- Background information
- Executable code
- Figures and plots
- Interactive visualizations

Launch Jupyter with:

```bash
jupyter notebook
```

A browser window should open automatically.

Navigate to the notebook you wish to explore and open the corresponding `.ipynb` file.

To execute a notebook:

- Click inside a cell.
- Press **Shift + Enter**.

This executes the current cell and advances to the next one.


# Running Julia Scripts
There are two common ways to run a Julia script.

## Option 1: From the Julia REPL

Navigate to the appropriate directory and start Julia:

```bash
julia --project
```

Then execute the script:

```julia
include("filename.jl")
```

## Option 2: From the Command Line

Run the script directly:

```bash
julia --project filename.jl
```


---

# Repository Structure

The repository is organized by scientific data type.

```text
repository/
│
├── src/
│   └── Reusable Julia helper functions
│
├── DataTypeA/
│   ├── download_data.jl
│   ├── visualization.ipynb
│   ├── README.md
│   └── data/
│
├── DataTypeB/
│   ├── ...
│
├── Project.toml
└── Manifest.toml
```

---

# Data Type Directories

Each data type has its own directory containing the files needed to download, understand, and visualize that dataset.

## Download Script

The download script retrieves any required datasets that are not already present.

Example:

```text
download_data.jl
```

---

## Visualization Notebook

A Jupyter notebook used to explore and visualize the dataset.

Example:

```text
visualization.ipynb
```

---

## Documentation

Each directory includes a `README.md` describing:

- What the dataset represents
- Where the data comes from
- Relevant scientific background
- Mathematical models used to describe the data

---

## Data Files

The `data/` directory contains the raw data used by the notebooks and scripts.

Some datasets may not be stored directly in GitHub because of file size limitations. In those cases, use the provided download script to retrieve the required files.

---

# Helper Functions

Reusable Julia utilities are stored in the `src/` directory.

Load them into a script or notebook using `include`.

For example:

```julia
include("../src/plotting.jl")
```

This makes the functions defined in that file available for use.

---

# Suggested Workflow

If this is your first time using the repository:

1. Clone the repository.
2. Start Julia:

   ```bash
   julia --project
   ```

3. Install the project dependencies:

   ```julia
   ]
   instantiate
   ```

4. Navigate to the directory for the dataset you want to explore.
5. Read the `README.md`.
6. Run the data download script if necessary.
7. Open the visualization notebook.
8. Run the notebook from top to bottom.
9. Experiment with the code and visualizations.

---

# Troubleshooting

## Julia Cannot Find a Package

Ensure that Julia was started with the project environment:

```bash
julia --project
```

Then install the required packages:

```julia
]
instantiate
```

---

## Notebook Fails to Run

Try restarting the notebook kernel and then run all cells from the beginning.

---

## Data Files Are Missing

Consult the `README.md` in the corresponding data directory.

Some datasets must be downloaded separately using the provided `download_data.jl` script.