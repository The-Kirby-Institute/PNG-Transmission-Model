### line of code to add the parent directory to the system path
import sys, os; sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# Define base path as the parent directory of the current script
base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


from matplotlib.ticker import MultipleLocator

small_df_size = 37
starting_epidemic_year = 1994
reference_year = 2000
adj =  reference_year - starting_epidemic_year
tspan = np.arange(0, small_df_size, 1)


# Function to extract F_W1_median from an Excel file
def extract_excel_files(excel_file, small_df_size):
    df = pd.read_excel(excel_file)
    number_of_small_dfs = len(df) // small_df_size
    writing_sols = [df.iloc[i*small_df_size:(i+1)*small_df_size] for i in range(number_of_small_dfs)]
    return writing_sols

def setCorrectxAxis(frequency_ticks =5,starting_position=0):
    plt.figure(figsize=(7, 5))
    default_x_ticks = list(range(small_df_size))
    new_x_ticks = list(range(starting_epidemic_year,starting_epidemic_year+small_df_size ))
    plt.xticks(default_x_ticks,new_x_ticks)
    plt.xticks(np.append(np.arange(starting_position, small_df_size+1, frequency_ticks),default_x_ticks[-1]))
    plt.xticks(fontsize=8, rotation=45)

folder_path = "output/Chp2_scenarios"
img_path = "output/Chp2-figures"



def process_excel_files(excel_files, folder_path, small_df_size, tspan):
    # Initialize DataFrames
    dataframe_viralsuppression = pd.DataFrame()
    dataframe_DR = pd.DataFrame()
    dataframe_VL_ADR = pd.DataFrame()
    dataframe_VL_TDR = pd.DataFrame()
    dataframe_suppressionperc = pd.DataFrame()
    dataframe_newDFfailure = pd.DataFrame()
    dataframe_incidences = pd.DataFrame()

    # Generate scenario names
    scenario_names = ["Scenario_" + str(i) for i in range(1, len(excel_files) + 1)]

    # Loop through each Excel file and scenario
    for file, scenario in zip(excel_files, scenario_names):
        full_file_path = os.path.join(base_path,folder_path, file)
        writing_sols = extract_excel_files(full_file_path, small_df_size)

        # Calculate medians and populate DataFrames
        dataframe_viralsuppression[scenario] = np.median(
            np.array([(sol["F_W1"] + sol["F_W2"] + sol["F_TDR1"] + sol["F_ADR1"]) for sol in writing_sols]), axis=0
        )

        dataframe_DR[scenario] = np.median(
            np.array([(sol.D_DR + sol.I_DR) / (sol["D"] + sol["I"]) * 100 for sol in writing_sols]), axis=0
        )

        dataframe_VL_ADR[scenario] = np.median(
            np.array([sol["F_ADR1"] for sol in writing_sols]), axis=0
        )

        dataframe_VL_TDR[scenario] = np.median(
            np.array([sol["F_TDR1"] for sol in writing_sols]), axis=0
        )

        dataframe_suppressionperc[scenario] = np.median(
            np.array([sol["T"] / (sol["T"] + sol["F"]) * 100 for sol in writing_sols]), axis=0
        )

        dataframe_newDFfailure[scenario] = np.median(
            np.array([np.gradient(sol['newAcquiredDR'], tspan) + np.gradient(sol['newTransmittedDR'], tspan) for sol in writing_sols]), axis=0
        )

        dataframe_incidences[scenario] = np.median(
            np.array([np.gradient(sol['incidence'], tspan)  for sol in writing_sols]), axis=0
        )


    # Return all DataFrames as a dictionary
    results = {
        "dataframe_viralsuppression": dataframe_viralsuppression,
        "dataframe_DR": dataframe_DR,
        "dataframe_VL_ADR": dataframe_VL_ADR,
        "dataframe_VL_TDR": dataframe_VL_TDR,
        "dataframe_suppressionperc": dataframe_suppressionperc,
        "dataframe_newDFfailure": dataframe_newDFfailure,
        "dataframe_incidences":dataframe_incidences
    }

    return results






#### viral suppression levels in PNG - baseline

excel_files = ["Baseline_Routine_2_WithDolutegravir.xlsx","Baseline_POC_2_WithDolutegravir.xlsx","Baseline_POC_3_WithDolutegravir.xlsx","Baseline_POC_4_WithDolutegravir.xlsx","Baseline_Routine_2_WithDolutegravir.xlsx","Baseline_Routine_2_NoDolutegravir.xlsx"]  

results = process_excel_files(excel_files, folder_path, small_df_size, tspan)

# Access individual DataFrames
dataframe_viralsuppression = results["dataframe_viralsuppression"]
dataframe_DR = results["dataframe_DR"]
dataframe_VL_ADR = results["dataframe_VL_ADR"]
dataframe_VL_TDR = results["dataframe_VL_TDR"]
dataframe_suppressionperc = results["dataframe_suppressionperc"]
dataframe_newDFfailure = results["dataframe_newDFfailure"]
dataframe_incidences = results["dataframe_incidences"]

setCorrectxAxis(5,1)
plt.plot (dataframe_viralsuppression['Scenario_1'], label="- Central lab VL testing",c = "black")
plt.plot (dataframe_viralsuppression['Scenario_2'], label="-Current(20%) POC levels in ACT-UP study",c = "blue")
plt.plot (dataframe_viralsuppression['Scenario_3'], label="-40% POC levels ",c = "green")
plt.plot (dataframe_viralsuppression['Scenario_4'], label="-60% POC levels ",c = "red")
plt.plot (dataframe_viralsuppression['Scenario_6'], c = "purple",alpha=0.0)
plt.plot (dataframe_viralsuppression['Scenario_1'], c = "black")

plt.ylabel('Unsuppressed VL cases (Total)',fontsize=13)
plt.ylim(ymin=0)  # this line
plt.grid(True)
plt.xlim(xmin = adj+10)

plt.xticks(fontsize=10)
plt.yticks(fontsize=12)
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.savefig(os.path.join(base_path, 'output','Chp2-figures', 'Viral Failures cases_withDolutegravir_NoDR.png'), dpi=500)
plt.show()











setCorrectxAxis(5,1)
plt.plot (dataframe_VL_ADR['Scenario_1'], label="-Central lab VL testing",c = "black")
plt.plot (dataframe_VL_ADR['Scenario_2'], label="-Current(20%) POC levels in ACT-UP study",c = "blue")
plt.plot (dataframe_VL_ADR['Scenario_3'], label="-40% POC levels ",c = "green")
plt.plot (dataframe_VL_ADR['Scenario_4'], label="-60% POC levels ",c = "red")
plt.plot (dataframe_VL_ADR['Scenario_6'], c = "purple",alpha=0.0)
plt.plot (dataframe_VL_ADR['Scenario_1'], c = "black")


plt.ylabel('Treatment Failure (ADR)',fontsize=13)
plt.ylim(ymin=0, ymax = 6000)  # this line
plt.grid(True)
plt.xlim(xmin = adj+10)

plt.xticks(fontsize=10)
plt.yticks(fontsize=12)
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.savefig(os.path.join(base_path, 'output','Chp2-figures', 'VF of ADR_withDolutegravir_NoDR.png'), dpi=500)
plt.show()






setCorrectxAxis(5,1)
plt.plot (dataframe_VL_TDR['Scenario_1'], label="-Central lab VL testing",c = "black")
plt.plot (dataframe_VL_TDR['Scenario_2'], label="-Current(20%) POC levels in ACT-UP study",c = "blue")
plt.plot (dataframe_VL_TDR['Scenario_3'], label="-40% POC levels ",c = "green")
plt.plot (dataframe_VL_TDR['Scenario_4'], label="-60% POC levels ",c = "red")
plt.plot (dataframe_VL_TDR['Scenario_6'], c = "purple",alpha=0.0)
plt.plot (dataframe_VL_TDR['Scenario_1'], c = "black")

plt.ylabel('Treatment Failure (TDR)',fontsize=13)
plt.ylim(ymin=0,ymax = 3000)  
plt.grid(True)
plt.xlim(xmin = adj+10)

plt.xticks(fontsize=10)
plt.yticks(fontsize=12)
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.savefig(os.path.join(base_path, 'output','Chp2-figures', 'VF of TDR_withDolutegravir_NoDR.png'), dpi=500)
plt.show()









setCorrectxAxis(5,1)
plt.plot (dataframe_newDFfailure['Scenario_1'], label="-Central lab VL testing",c = "black")
plt.plot (dataframe_newDFfailure['Scenario_2'], label="-Current(20%) POC levels in ACT-UP study",c = "blue")
plt.plot (dataframe_newDFfailure['Scenario_3'], label="-40% POC levels ",c = "green")
plt.plot (dataframe_newDFfailure['Scenario_4'], label="-60% POC levels ",c = "red")
plt.plot (dataframe_newDFfailure['Scenario_6'], c = "purple",alpha=0.0)
plt.plot (dataframe_newDFfailure['Scenario_1'], c = "black")


plt.xlim(xmin=adj+10)
plt.ylabel('New Treatment Failure (Annual)')
plt.grid(True)

plt.xticks(fontsize=10)
plt.yticks(fontsize=12)
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.savefig(os.path.join(base_path, 'output','Chp2-figures', 'NewViralLoad_DFfailure_withDolutegravir_NoDR.png'), dpi=500)
plt.show()











############################################################################
############################################################################
############################################################################
############################################################################

#### New set of results 
excel_files = ["Baseline_POC_2_WithDolutegravir.xlsx","Timevarying_POC_AcquiredDR_2_WithDolutegravir.xlsx","Timevarying_POC_AcquiredDR_3_WithDolutegravir.xlsx","Timevarying_POC_AcquiredDR_4_WithDolutegravir.xlsx","Baseline_Routine_2_WithDolutegravir.xlsx","Baseline_Routine_2_NoDolutegravir.xlsx"]  

results = process_excel_files(excel_files, folder_path, small_df_size, tspan)

# Access individual DataFrames
dataframe_viralsuppression = results["dataframe_viralsuppression"]
dataframe_DR = results["dataframe_DR"]
dataframe_VL_ADR = results["dataframe_VL_ADR"]
dataframe_VL_TDR = results["dataframe_VL_TDR"]
dataframe_suppressionperc = results["dataframe_suppressionperc"]
dataframe_newDFfailure = results["dataframe_newDFfailure"]
dataframe_incidences = results["dataframe_incidences"]



setCorrectxAxis(5,1)
plt.plot (dataframe_viralsuppression['Scenario_5'], label="- Central lab VL testing ",c = "black")
plt.plot (dataframe_viralsuppression['Scenario_1'], label="-Current(20%) POC levels in ACT-UP study/ no DR testing",c = "blue")
plt.plot (dataframe_viralsuppression['Scenario_2'], label="-Current(20%) POC levels in ACT-UP study/ with DR testing",c = "blue",linestyle=":")
plt.plot (dataframe_viralsuppression['Scenario_3'], label="-40% POC levels / with DR testing",c = "green",linestyle=":")
plt.plot (dataframe_viralsuppression['Scenario_4'], label="-60% POC levels / with DR testing",c = "red",linestyle=":")
plt.plot (dataframe_viralsuppression['Scenario_6'], c = "purple",alpha=0.0)
plt.plot (dataframe_viralsuppression['Scenario_5'], c = "black")

plt.ylim(ymin=0)  # this line
plt.grid(True)
plt.xlim(xmin = adj+10)

plt.xticks(fontsize=10)
plt.yticks(fontsize=12)
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.savefig(os.path.join(base_path, 'output','Chp2-figures', 'Viral Failures cases_withDolutegravir_with DR testing.png'), dpi=500)
plt.show()















setCorrectxAxis(5,1)
plt.plot (dataframe_VL_ADR['Scenario_1'], label="-Current(20%) POC levels in ACT-UP study/ no DR testing",c = "blue")
plt.plot (dataframe_VL_ADR['Scenario_2'], label="-Current(20%) POC levels in ACT-UP study/ with DR testing",c = "blue",linestyle=":")
plt.plot (dataframe_VL_ADR['Scenario_3'], label="-40% POC levels / with DR testing",c = "green",linestyle=":")
plt.plot (dataframe_VL_ADR['Scenario_4'], label="-60% POC levels / with DR testing",c = "red",linestyle=":")
plt.plot (dataframe_VL_ADR['Scenario_5'], label="- Central lab VL testing ",c = "black")
plt.plot (dataframe_VL_ADR['Scenario_6'], c = "purple",alpha=0.0)
# plt.plot (dataframe_VL_ADR['Scenario_1'], c = "blue")





# plt.legend(loc='lower center', prop={'size': 8})
# plt.ylabel('Virological failure cases due to ADR\nwith DR testing')
plt.ylim(ymin=0 , ymax =6000)  # this line
plt.grid(True)
plt.xlim(xmin = adj+10)

plt.xticks(fontsize=10)
plt.yticks(fontsize=12)
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.savefig(os.path.join(base_path, 'output','Chp2-figures', 'VF of ADR_withDolutegravir_with DR.png'), dpi=500)
plt.show()







setCorrectxAxis(5,1)
plt.plot (dataframe_VL_TDR['Scenario_1'], label="-Current(20%) POC levels in ACT-UP study/ no DR testing",c = "blue")
plt.plot (dataframe_VL_TDR['Scenario_2'], label="-Current(20%) POC levels in ACT-UP study/ with DR testing",c = "blue",linestyle=":")
plt.plot (dataframe_VL_TDR['Scenario_3'], label="-40% POC levels / with DR testing",c = "green",linestyle=":")
plt.plot (dataframe_VL_TDR['Scenario_4'], label="-60% POC levels / with DR testing",c = "red",linestyle=":")
plt.plot (dataframe_VL_TDR['Scenario_5'], label="- Central lab VL testing ",c = "black")
plt.plot (dataframe_VL_TDR['Scenario_6'], c = "purple",alpha=0.0)

plt.ylim(ymin=0, ymax=3000)  
plt.grid(True)
plt.xlim(xmin = adj+10)

plt.xticks(fontsize=10)
plt.yticks(fontsize=12)
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.savefig(os.path.join(base_path, 'output','Chp2-figures', 'VF of TDR_withDolutegravir_with DR.png'), dpi=500)

plt.show()






setCorrectxAxis(5,1)
plt.plot (dataframe_newDFfailure['Scenario_1'], label="-Current(20%) POC levels in ACT-UP study/ no DR testing",c = "blue")
plt.plot (dataframe_newDFfailure['Scenario_2'], label="-Current(20%) POC levels in ACT-UP study/ with DR testing",c = "blue",linestyle=":")
plt.plot (dataframe_newDFfailure['Scenario_3'], label="-40% POC levels / with DR testing",c = "green",linestyle=":")
plt.plot (dataframe_newDFfailure['Scenario_4'], label="-60% POC levels / with DR testing",c = "red",linestyle=":")
plt.plot (dataframe_newDFfailure['Scenario_6'], c = "purple",alpha=0.0)
plt.plot (dataframe_newDFfailure['Scenario_5'], label="-Central lab VL testing",c = "black")


# plt.minorticks_on()

# plt.grid(which='major', linestyle='-', linewidth='0.5')
# plt.grid(which='minor', linestyle=':', linewidth='0.2')

# plt.legend(loc="lower right",prop={'size': 8})
# plt.ylabel('New Treatment Failure cases')
plt.xlim(xmin=adj+10)
plt.grid(True)

plt.xticks(fontsize=10)
plt.yticks(fontsize=12)
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.savefig(os.path.join(base_path, 'output','Chp2-figures', 'NewViralLoad_DFfailure_with DR and access levels.png'), dpi=500)

plt.show()











############################################################################
############################################################################
############################################################################
############################################################################



#### viral suppression levels in PNG - Drug Resistance testing 

excel_files = ["Baseline_POC_2_WithDolutegravir.xlsx","Timevarying_POC_2_WithDolutegravir.xlsx","Baseline_POC_AcquiredDR_2_WithDolutegravir.xlsx","Timevarying_POC_AcquiredDR_2_WithDolutegravir.xlsx","Baseline_Routine_2_WithDolutegravir.xlsx","Baseline_Routine_2_NoDolutegravir.xlsx"]  

results = process_excel_files(excel_files, folder_path, small_df_size, tspan)

# Access individual DataFrames
dataframe_viralsuppression = results["dataframe_viralsuppression"]
dataframe_DR = results["dataframe_DR"]
dataframe_VL_ADR = results["dataframe_VL_ADR"]
dataframe_VL_TDR = results["dataframe_VL_TDR"]
dataframe_suppressionperc = results["dataframe_suppressionperc"]
dataframe_newDFfailure = results["dataframe_newDFfailure"]
dataframe_incidences = results["dataframe_incidences"]



setCorrectxAxis(5,1)
plt.plot (dataframe_viralsuppression['Scenario_1'], label="- Current(20%) POC levels in ACT-UP study/ no DR testing",c = "blue")
plt.plot (dataframe_viralsuppression['Scenario_2'], label="-Current(20%) POC levels in ACT-UP study/ with only pre-treatment DR testing",c = "blue",linestyle="--")
plt.plot (dataframe_viralsuppression['Scenario_3'], label="-Current(20%) POC levels in ACT-UP study/ with only post-treatment DR testing",c = "blue",linestyle="-.")
plt.plot (dataframe_viralsuppression['Scenario_4'], label="-Current(20%) POC levels in ACT-UP study/ with both DR testing",c = "blue",linestyle=":")
plt.plot (dataframe_viralsuppression['Scenario_5'], label="- Central lab VL testing ",c = "black")
plt.plot (dataframe_viralsuppression['Scenario_6'], c = "purple",alpha=0.0)
# plt.plot (dataframe_viralsuppression['Scenario_1'], c = "blue")




# plt.legend(loc="upper left",prop={'size': 8})
# plt.ylabel('Virological failure cases \nwith DR testing')
plt.ylim(ymin=0)  # this line
plt.grid(True)
plt.xlim(xmin = adj+10)

plt.xticks(fontsize=10)
plt.yticks(fontsize=12)
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.savefig(os.path.join(base_path, 'output','Chp2-figures', 'Viral Failures cases_withDolutegravir_with ACT-UP levels+DR test.png'), dpi=500)
plt.show()









setCorrectxAxis(5,1)
plt.plot (dataframe_VL_ADR['Scenario_1'], label="- Current(20%) POC levels in ACT-UP study/ no DR testing",c = "blue")
plt.plot (dataframe_VL_ADR['Scenario_2'], label="-Current(20%) POC levels in ACT-UP study/ with only pre-treatment DR testing",c = "blue",linestyle="--")
plt.plot (dataframe_VL_ADR['Scenario_3'], label="-Current(20%) POC levels in ACT-UP study/ with only post-treatment DR testing",c = "blue",linestyle="-.")
plt.plot (dataframe_VL_ADR['Scenario_4'], label="-Current(20%) POC levels in ACT-UP study/ with both DR testing",c = "blue",linestyle=":")
plt.plot (dataframe_VL_ADR['Scenario_5'], label="- Central lab VL testing ",c = "black")
plt.plot (dataframe_VL_ADR['Scenario_6'], c = "purple",alpha=0.0)
# plt.plot (dataframe_VL_ADR['Scenario_1'], c = "blue")




# plt.legend(loc="lower right",prop={'size': 8})
# plt.ylabel('Virological failure cases due to ADR \nwith DR testing')
plt.ylim(ymin=0, ymax=6000)  # this line
plt.grid(True)
plt.xlim(xmin = adj+10)

plt.xticks(fontsize=10)
plt.yticks(fontsize=12)
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.savefig(os.path.join(base_path, 'output','Chp2-figures', 'VF of ADR_withDolutegravir_with  ACT-UP levels+DR test.png'), dpi=500)
plt.show()







setCorrectxAxis(5,1)
plt.plot (dataframe_VL_TDR['Scenario_1'], label="- Current(20%) POC levels in ACT-UP study/ no DR testing",c = "blue")
plt.plot (dataframe_VL_TDR['Scenario_2'], label="-Current(20%) POC levels in ACT-UP study/ with only pre-treatment DR testing",c = "blue",linestyle="--")
plt.plot (dataframe_VL_TDR['Scenario_3'], label="-Current(20%) POC levels in ACT-UP study/ with only post-treatment DR testing",c = "blue",linestyle="-.")
plt.plot (dataframe_VL_TDR['Scenario_4'], label="-Current(20%) POC levels in ACT-UP study/ with both DR testing",c = "blue",linestyle=":")
plt.plot (dataframe_VL_TDR['Scenario_5'], label="- Central lab VL testing ",c = "black")
plt.plot (dataframe_VL_TDR['Scenario_6'], c = "purple",alpha=0.0)
# plt.plot (dataframe_VL_TDR['Scenario_1'], c = "blue")




# plt.legend(loc="upper left",prop={'size': 8})
# plt.ylabel('Virological failure cases due to TDR \nwith DR testing')
plt.ylim(ymin=0,ymax=3000)  
plt.grid(True)
plt.xlim(xmin = adj+10)

plt.xticks(fontsize=10)
plt.yticks(fontsize=12)
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.savefig(os.path.join(base_path, 'output','Chp2-figures', 'VF of TDR_withDolutegravir_with ACT-UP levels+DR test.png'), dpi=500)
plt.show()






setCorrectxAxis(5,1)

plt.plot (dataframe_newDFfailure['Scenario_1'], label="- Current(20%) POC levels in ACT-UP study/ no DR testing",c = "blue")
plt.plot (dataframe_newDFfailure['Scenario_2'], label="-Current(20%) POC levels in ACT-UP study/ with only pre-treatment DR testing",c = "blue",linestyle="--")
plt.plot (dataframe_newDFfailure['Scenario_3'], label="-Current(20%) POC levels in ACT-UP study/ with only post-treatment DR testing",c = "blue",linestyle="-.")
plt.plot (dataframe_newDFfailure['Scenario_4'], label="-Current(20%) POC levels in ACT-UP study/ with both DR testing",c = "blue",linestyle=":")
plt.plot (dataframe_newDFfailure['Scenario_6'], c = "purple",alpha=0.0)
plt.plot (dataframe_newDFfailure['Scenario_5'], label="-Central lab VL testing",c = "black")


plt.xlim(xmin=adj+10)
plt.grid(True)

plt.xticks(fontsize=10)
plt.yticks(fontsize=12)
plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}'))
plt.savefig(os.path.join(base_path, 'output','Chp2-figures', 'NewViralLoad_DFfailure_ACTUP level +DR tests.png'), dpi=500)
plt.show()





listB = [
        "Viral Failures cases_withDolutegravir_NoDR.png", 
        "Viral Failures cases_withDolutegravir_with DR testing.png",
        "Viral Failures cases_withDolutegravir_with ACT-UP levels+DR test.png",

        "VF of ADR_withDolutegravir_NoDR.png",
        "VF of ADR_withDolutegravir_with DR.png",
        "VF of ADR_withDolutegravir_with  ACT-UP levels+DR test.png",
        
        "VF of TDR_withDolutegravir_NoDR.png",
        "VF of TDR_withDolutegravir_with DR.png",
        "VF of TDR_withDolutegravir_with ACT-UP levels+DR test.png",

        "NewViralLoad_DFfailure_withDolutegravir_NoDR.png",
        "NewViralLoad_DFfailure_with DR and access levels.png",
        "NewViralLoad_DFfailure_ACTUP level +DR tests.png"
        ]

img_list = [mpimg.imread(os.path.join(base_path,img_path,filename)) for filename in listB]

# Create a 4x3 grid for displaying the images:
fig, axes = plt.subplots(4, 3, figsize=(20,20))

for ax, img, title in zip(axes.ravel(), img_list, listB):
    ax.imshow(img)
    ax.axis('off')

# Add vertical lines to separate the 3 columns
for ax in axes[:, 0]:
    ax.axvline(x=ax.get_xlim()[1], color='black', linewidth=3, linestyle='--')  # Line between first and second column
for ax in axes[:, 1]:
    ax.axvline(x=ax.get_xlim()[1], color='black', linewidth=3, linestyle='--')  # Line between second and third column

# Add text annotations for column labels
fig.text(0.16, 0.01, '(a)', fontsize=18, ha='center')  
fig.text(0.5, 0.01, '(b)', fontsize=18, ha='center')   
fig.text(0.83, 0.01, '(c)', fontsize=18, ha='center')  


plt.tight_layout(h_pad = -8)  # Added vertical padding between subplots
plt.savefig(os.path.join(base_path, 'output','Chp2-figures', 'All scenarios figure 3_production_revised.png'), dpi=500)

plt.show()

