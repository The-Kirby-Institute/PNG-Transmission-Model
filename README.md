# PNG-Transmission-Model
This repository contains a compartmental transmission model based on ordinary differential equations, developed to evaluate the impact of scaling up viral load and drug resistance testing in Papua New Guinea. The model also assesses the cost-effectiveness of these interventions.

## Model description 
This compartmental model simulates the HIV epidemics in PNG since around 1993/94 - the first year when seroprevalence studies were conducted in the country.
The model reflects the time-varying effect of diagnosis and treatment scale-up on the HIV epidemics. It also reflects the rapid scale-up of dolutegravir in PNG since 2020. The model was calibrated using Bayesian methods informed by multiple epidemiological and modelling data sources.

While primarily designed to evaluate the scale-up of viral load and drug resistance testing, the model can also be adapted to explore broader HIV epidemic trends and the impact of other interventions.
## Project organisation
PNG-TRANSMISSION-MODEL/
├── data/
├── model/
├── output/
│   ├── figures
│   └── script2.py
├── projects/
│   ├── script1.py
│   └── script2.py
└── Calibration.py
└── README.md (optional)

## Data sources
Several open source data were used to calibrate this model. 
1. [HIV estimates with uncertainty bounds 1990-Present](https://www.unaids.org/en/resources/documents/2024/HIV_estimates_with_uncertainty_bounds_1990-present). Treatment cascade estimates from UNAIDS were used. 
2. [World Bank Open Data](https://data.worldbank.org). PNG population estimates and all-cause mortality were based on the Open data. 
We thank UNAIDS and World Bank for their generosity in making important modelling datasets available for research and public use. 
In addition, epidemiological estimates from PNG's literature and reports were also used. 
## Calibration
The model is calibrated using Bayesian methods, using the Python package [PyMC](https://www.pymc.io/welcome.html). The calibration process takes up about 20 hours on a Macbook Pro M1 16-inch with 32Gb of integrated memory. 
To save your time in reproducing the results, an Excel sheet with 1,000 posterior samples were provided together with this model. `model priors.xlsx` in the `data/priors` folder 
## Using the model. 
After cloning this repository, several software packages are required to run this model:
1. Python. It is preferable to install Anaconda Distribution which contains most of the scientific packages used in this model 
2. The Python packages used in the model include: 
`pymc >=20.0`, `numpy`, `scipy`, `matplotlib`
## Citation
Please cite Quang's PhD thesis for this repository. 

Quang Nguyen. Impact and Cost-Effectiveness of Simplified HIV Monitoring and Vertical Transmission Interventions in Papua New Guinea, *UNSW's PhD thesis repository* __2025__

## Contact 
This model in Python is designed, developed and maintained by [Quang Nguyen](https://github.com/DucQuang1); ORCiD ID: [0000-0001-6919-5211](https://orcid.org/0000-0001-6919-5211)

Please contact Quang at qnguyen@kirby.unsw.edu.au for collaboration, queries or methodology questions/ concerns related to this model. We will respond promptly to your questions and concerns in future model updates. 
## License 
This model is distributed under the GNU Affero General Public License v3.0 (AGPL-3.0). You are free to use, modify, and distribute the code, including for commercial purposes, provided that any derivative work is also made publicly available under the same license. For more details, see the LICENSE file.
