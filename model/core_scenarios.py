### different types of monitoring CD4, Routine and POC, NoDolutegravir
# Monitoring_Type = "Routine"
# scenario = "2"
mode = "actual"

# DR_testing_scenario = "Baseline"

import json
import os

def get_config():
    # Get the directory where this file is located
    base_dir = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(base_dir, 'config.json')
    with open(config_path, 'r') as file:
        return json.load(file)

config = get_config()
Monitoring_Type = config['Monitoring_Type']
scenario = config['scenario']
DR_testing_scenario = config['DR_testing_scenario']
treatment = config['treatment']
running_mode = config['running_mode']
analysis_mode = config['analysis_mode']
year_OffInfs = config["year_OffInfs"]
year_OffInterv = config["year_OffInterv"]

drugresistance_testingrate = 0.00
if DR_testing_scenario == "Baseline":
    drugresistance_testingrate = 0.00
elif DR_testing_scenario == "Medium":
    drugresistance_testingrate = 0.30
elif DR_testing_scenario == "High":
    drugresistance_testingrate = 0.90