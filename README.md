# Scientific Data Models

## What is this Repository?

Many sustainably powered robots are designed to collect scientific information about the world. Examples include environmental monitoring, ocean sampling, weather observation, and ecological surveys.

This repository contains examples of different types of scientific data, along with:

* Data sources and download scripts
* Mathematical models used to describe the data
* Visualization notebooks
* Supporting Julia code

The goal is to provide a collection of scientific data models that can be explored, visualized, and compared.

---

# Getting Started

## What You Need

Before using this repository, install:

1. **Julia** (the programming language used in this project)

   * Download: https://julialang.org/downloads/

2. **Git** (used to download the repository)

   * Download: https://git-scm.com/downloads

3. **Jupyter Notebook** (used to view and run visualizations)

   * Installation instructions: https://jupyter.org/install


---

# Downloading the Repository

The repository is hosted on GitHub.

To download it, open a terminal and run:

```bash
git clone <repository-url>
```

This will create a local copy of the repository on your computer.

Move into the repository directory:

```bash
cd <repository-name>
```

---

# Setting Up the Julia Environment

This repository includes a file called `Project.toml`.

This file tells Julia which packages are required for the project.

From the repository directory, start Julia:

```bash
julia --project
```

You should see a prompt that looks something like:

```julia
julia>
```

Press the `]` key to enter Julia's package manager mode. The prompt should change to:

```julia
(@v1.x) pkg>
```

Now run:

```julia
instantiate
```

This may take several minutes the first time. Julia will automatically download and install all required packages.

When installation finishes, press Backspace to return to the normal Julia prompt.

You only need to run `instantiate` once after cloning the repository.

---

# Opening the Visualization Notebooks

Each data type includes one or more Jupyter notebooks (`.ipynb` files).

A notebook contains:

* Explanations
* Code
* Figures and plots
* Interactive exploration tools

To start Jupyter:

```bash
jupyter notebook
```

A browser window should open automatically.

Navigate to the notebook you want to explore and open the `.ipynb` file.

To run a notebook:

* Click inside a cell.
* Press **Shift + Enter**.

This executes the code and moves to the next cell.

---

# Repository Structure

The repository is organized by data type.

```text
repository/
│
├── src/
│   └── Helper Julia functions
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

Each data type has its own directory.

Inside you will typically find:

### Download Script

A Julia script that downloads the data.

Example:

```text
download_data.jl
```

Run this script if the required data files are not already present.

### Visualization Notebook

A Jupyter notebook used to explore and visualize the data.

Example:

```text
visualization.ipynb
```

### Documentation

A Markdown file explaining:

* What the data represents
* Where it comes from
* Relevant scientific background
* Mathematical models used

Example:

```text
README.md
```

### Data Files

The raw data used by the notebooks and scripts.

Large datasets may not be stored directly in GitHub due to size limits.

---

# Helper Functions

Reusable Julia functions are stored in the `src/` directory.

These functions can be loaded into a script or notebook using:

```julia
include("../src/filename.jl")
```

For example:

```julia
include("../src/plotting.jl")
```

This makes the functions defined in that file available for use.

---

# Suggested Workflow

If this is your first time using the repository:

1. Clone the repository.

2. Start Julia with:

   ```bash
   julia --project
   ```

3. Run:

   ```julia
   ]
   instantiate
   ```

4. Open a data type directory.

5. Read the documentation file.

6. Open the visualization notebook.

7. Run the notebook from top to bottom.

8. Experiment with the code and plots.

---

# Troubleshooting

## Julia Cannot Find a Package

Make sure you started Julia with:

```bash
julia --project
```

and have run:

```julia
]
instantiate
```

---

## Notebook Fails to Run

Try restarting the notebook kernel and running all cells again from the beginning.

---

## Data Files Are Missing

Check the documentation in the corresponding data directory. Some datasets must be downloaded separately using the provided download script.

---
