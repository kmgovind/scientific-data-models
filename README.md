# Scientific Data Models
In my research, I look at the persistent optimal control of sustaainbly powered robot, and in general, these robots are operating with the goal of collecting some form of scientific information.

In order to better understand and characterize the various types of scientific information that exist, this repository hosts a variety of information types, with associated models and data visualizations.

# To Use the Repository
Ensure that you have the Julia programming language installed ([installation guide here](https://julialang.org/downloads/)).

Clone the repository. The `Project.toml` file contains the list of all necessary packages for the project and their associated versions.

Once the repository has been cloned, navigate to the repository's directory and run
```bash
julia --project
```

This will start the Julia REPL with the associated project. The `]` key can be used to put the REPL into package manager mode. Once here, run `instantiate`. This will install all the necessary packages.

# Repository Organization
## Data Type
Each type of data will have it's own directory. Inside that directory there will be a script for data download, a notebook (.ipynb) for data visualization, and a markdown (.md) file explaining the data and associated models. The datafiles will also be contained within that directory (though for space reasons, may not be uploaded to github).

## Helper Function
Any helper functions will go into julia (.jl) files in the `src/` directory. The functions can be called by importing the files by:
```julia
include("../src/[filename].jl)
```