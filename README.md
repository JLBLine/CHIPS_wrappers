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

### chips1D_tsv example

rather than plotting the 1D power, you can export to a tsv file for use in another program.

this example assumes:
- eor high-band data averaged to 80kHZ
- lssa_fg_simple bias_mode=10 (upper channels)
- dat files are in `./crosspower_${pol}_${bias_mode}.iter.${tag}.dat`
```bash
# export tag=...
# export pol=... xx or yy
export bias_mode=10
export n_chan=192
export lowerfreq=182515000

chips1D_tsv.py \
    --basedir "$PWD/"  \
    --chips_tag ${group}_${name}  \
    --polarisation $pol  \
    --lowerfreq $lowerfreq \
    --umax 300 \
    --N_chan $n_chan \
    --ktot_bin_edges ktot_bin_edges_cmt.txt \
    --plot_delta \
    --kperp_max 0.03 \
    --kperp_min 0.02 \
    --kparra_min 0.11 \
    --wedge_factor 2. \
    --bias_mode ${bias_mode}
```