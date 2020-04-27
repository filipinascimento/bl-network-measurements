[![Abcdspec-compliant](https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg)](https://github.com/brain-life/abcd-spec)
[![Run on Brainlife.io](https://img.shields.io/badge/Brainlife-bl.app.1-blue.svg)](https://doi.org/10.25663/brainlife.app.321)

# Network Measurements
App to calculate several basic statistics for networks and their respective null model distributions.

## List of measurements:

### Node-based measurements:
 - *Degree*
 - *In-Degree* (`directed`)
 - *Out-Degree* (`directed`)
 - *Strength* (`weighted`)
 - *In-Strength* (`directed`,`weighted`)
 - *Out-Strength* (`directed`,`weighted`)
 - *Clustering* Coefficient
 - *Coreness*
 - *Match Index*
 - *Betweeness Centrality*
 - *Betweeness CentralityWeighted* (`weighted`)

### Network-based measurements:
 - *Avg. Degree*
 - *Avg. In-Degree* (`directed`)
 - *Avg. Out-Degree* (`directed`)
 - *Avg. Strength* (`weighted`)
 - *Avg. In-Strength* (`directed`,`weighted`)
 - *Avg. Out-Strength* (`directed`,`weighted`)
 - *Avg. Clustering* Coefficient
 - *Betweenness Centralization*
 - *Rich-Club Coefficient*
 - *Degree Assortativity*
 - *Diameter*
 - *Module Degree Z-Score* (`communities`)
 - *Participation Coefficient* (`communities`)
 - *Modularity* (`communities`)

> TODO: Add references for all statistics

Some statistics are calculated only for `directed` or `weighted` networks. Those requiring `communities` are only calculated if the Communities App was executed on this data before.


### Authors
- [Filipi N. Silva](filsilva@iu.edu)

<!-- ### Contributors
- Franco Pestilli (franpest@indiana.edu) -->

<!-- ### Funding  -->
<!-- [![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1734853) -->


### Citations

1. Bassett, Danielle S., and Olaf Sporns. "Network neuroscience." Nature neuroscience 20, no. 3 (2017): 353. [https://doi.org/10.1038/nn.4502](https://doi.org/10.1038/nn.4502)

2. Costa, L. da F., Francisco A. Rodrigues, Gonzalo Travieso, and Paulino Ribeiro Villas Boas. "Characterization of complex networks: A survey of measurements." Advances in physics 56, no. 1 (2007): 167-242.[https://doi.org/10.1080/00018730601170527](https://doi.org/10.1080/00018730601170527)

## Running the App 

### On Brainlife.io

You can submit this App online at [https://doi.org/10.25663/brainlife.app.321](https://doi.org/10.25663/brainlife.app.321) via the "Execute" tab.

### Running Locally (on your machine)
Singularity is required to run the package locally.

1. git clone this repo.

```bash
git clone <repository URL>
cd <repository PATH>
```

2. Inside the cloned directory, edit `config-sample.json` with your data or use the provided data.

3. Rename `config-sample.json` to `config.json` .

```bash
mv config-sample.json config.json
```

4. Launch the App by executing `main`

```bash
./main
```

### Sample Datasets

A sample dataset is provided in folder `data` and `config-sample.json`

## Output

The output is a `conmat` enriched with the measurements. Network measurements and null model statistics can be accessed directly from the files ending with `_measurements.txt` in the output `csv` directory. Node measurements are named as `<Name of Network>_prop_<Name of Property>.txt`. If `generatePlots` option is enabled, the distributions of statistics are plotted together with null model if present in a new directory `figures`.

<!-- #### Product.json

The secondary output of this app is `product.json`. This file allows web interfaces, DB and API calls on the results of the processing.  -->

### Dependencies

This App only requires [singularity](https://www.sylabs.io/singularity/) to run. If you don't have singularity, you will need to install the python packages defined in `environment.yml`, then you can run the code directly from python using:  

```bash
./main.py config.json
```

