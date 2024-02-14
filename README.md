# CHIPS_wrappers
Scripts to run CHIPS power spectrum estimator, and plot the outputs. This repo is under development, and does not have unit/integration tests as of yet.

## Installation
Real basic at the mo, you need to do
```bash
git clone https://github.com/JLBLine/CHIPS_wrappers.git
cd CHIPS_wrappers
pip install .
```
maybe we can make it fancier later

## Usage
You'll need a working `CHIPS` installation, and to set a few environment variables to use `run_CHIPS.py`. The plotting code should work without env variables. *Fill in later*

### Docker

```bash
docker run --rm -it -w /app --entrypoint python d3vnull0/chips_wrappers:latest scripts/plotchips_all.py --help
```

### singularity

```bash
singularity exec --cleanenv -B $PWD --home $PWD docker://d3vnull0/chips_wrappers:latest python /app/scripts/plotchips_all.py --help
```