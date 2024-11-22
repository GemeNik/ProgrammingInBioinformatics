# Programming in Bioinformatics

## Overview
This project maps query sequences to reference sequences in a genome.

## Requirements
- Python 3.10+
- `flake8` for linting
- `pytest` for testing

## Setup Instructions
1. Clone the repository:
    ```bash
    git clone <repository-url>
    ```
2. Install dependencies using Conda:
    ```bash
    conda env create -f environment.yml
    conda activate myenv
    ```

OR build the Docker image:
```bash
docker build -t ProgrammingInBioinformatics
 .