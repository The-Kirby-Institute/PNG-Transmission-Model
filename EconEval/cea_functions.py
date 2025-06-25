import pandas as pd
import numpy as np
years = 2030-1994+1
tspan = np.arange(0, years, 1)

def calculate_discounted_value(compartment_data, values, discount_rate, start_index, end_index, tspan=tspan):    # Initialize total discounted value
    total_discounted_value = 0
    
    # Initialize a DataFrame to store discounted values
    period_index = compartment_data.index[start_index:end_index + 1]
    discounted_values = pd.DataFrame(index=period_index)
    
    # Iterate over each time period
    for t in range(start_index, end_index + 1):
        # Calculate the discount factor for the current time period
        discount_factor = 1 / ((1 + discount_rate) ** t)
        
        # Apply values and discount to each compartment
        for compartment, (value, use_gradient) in values.items():
            # Determine whether to use the gradient or the compartment itself
            if use_gradient:
                compartment_values = np.gradient(compartment_data[compartment].values,tspan)
            else:
                compartment_values = compartment_data[compartment].values
            
            # Calculate discounted value for the current time period
            discounted_compartment_value = compartment_values[t] * value * discount_factor
            
            # Prepare the column name
            col_name = compartment + ('_gradient_discounted_value' if use_gradient else '_discounted_value')
            
            # Use .at[] to avoid SettingWithCopyWarning
            discounted_values.at[compartment_data.index[t], col_name] = discounted_compartment_value
            
            # Add to the total discounted value
            total_discounted_value += discounted_compartment_value
    
    return total_discounted_value, discounted_values

def meta_calculation_all_populations(sols, values, discount_rate,start_index=0, end_index=2030-1994):
    final_total_discounted_value = []
    
    for df in sols:
        # Calculate the discounted value for this particular DataFrame
        total_discounted_value, _ = calculate_discounted_value(df, values, discount_rate,start_index,end_index)
        
        # Sum the discounted value to the final total
        final_total_discounted_value.append(total_discounted_value)
    
    return final_total_discounted_value



def calculate_discounted_value_sp(compartment_data, values, discount_rate, start_index, end_index, tspan=tspan):    # Initialize total discounted value
    total_discounted_value = 0
    
    
    # Iterate over each time period
    for t in range(start_index, end_index + 1):
        # Calculate the discount factor for the current time period
        discount_factor = 1 / ((1 + discount_rate) ** t)
        
        # Apply values and discount to each compartment
        for compartment, (value, use_gradient) in values.items():
            # Determine whether to use the gradient or the compartment itself
            if use_gradient:
                compartment_values = np.gradient(compartment_data[compartment].values,tspan)
            else:
                compartment_values = compartment_data[compartment].values
            
            # Calculate discounted value for the current time period
            discounted_compartment_value = compartment_values[t] * value * discount_factor
                        
            
            # Add to the total discounted value
            total_discounted_value += discounted_compartment_value
    
    return total_discounted_value


def probabilistic_all_populations(sols, values, discount_rate,start_index=0, end_index=2030-1994,tspan=tspan):
    final_total_discounted_value = []
    
    for df,val in zip(sols,values):
        # Calculate the discounted value for this particular DataFrame
        total_discounted_value = calculate_discounted_value_sp(df, val, discount_rate,start_index,end_index,tspan)
        
        # Sum the discounted value to the final total
        final_total_discounted_value.append(total_discounted_value)
    
    return final_total_discounted_value



def calculate_discounted_value_aslist(compartment_data, values, discount_rate, start_index, end_index, tspan=tspan):    # Initialize total discounted value
    
    list_discounted_values = []
    # Iterate over each time period
    for t in range(start_index, end_index + 1):
        total_discounted_value = 0
        # Calculate the discount factor for the current time period
        discount_factor = 1 / ((1 + discount_rate) ** t)
        
        # Apply values and discount to each compartment
        for compartment, (value, use_gradient) in values.items():
            # Determine whether to use the gradient or the compartment itself
            if use_gradient:
                compartment_values = np.gradient(compartment_data[compartment].values,tspan)
            else:
                compartment_values = compartment_data[compartment].values
            
            # Calculate discounted value for the current time period
            discounted_compartment_value = compartment_values[t] * value * discount_factor
                        
            
            # Add to the total discounted value
            total_discounted_value += discounted_compartment_value
        

        list_discounted_values.append(total_discounted_value)
    
    return list_discounted_values


def probabilistic_all_populations_aslist(sols, values, discount_rate,start_index=0, end_index=2030-1994,tspan=tspan):
    final_total_discounted_value = []
    
    for df,val in zip(sols,values):
        # Calculate the discounted value for this particular DataFrame
        total_discounted_value = calculate_discounted_value_aslist(df, val, discount_rate,start_index,end_index,tspan)
        
        # Sum the discounted value to the final total
        final_total_discounted_value.append(total_discounted_value)
    
    return np.array(final_total_discounted_value)
