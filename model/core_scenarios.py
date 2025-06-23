### different types of monitoring CD4, Routine and POC, NoDolutegravir
# Monitoring_Type = "Routine"
# scenario = "2"
mode = "actual"

# DR_testing_scenario = "Baseline"

import json

def get_config():
    with open('model/config.json', 'r') as file:
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