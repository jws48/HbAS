# Python functions for the Haase Lab


## Files
### heatmap.py
File containing heatmapping functions.

Author: Sophie Campione






## Functions
#### heatmap.py:
**heatmap_max(data, gene_list, first_period)**

input: 
1. data frame with time series data (index = gene names)
2. gene list with the genes to be plotted
3. first period

Z-score normalizes the dataframe and plots a heatmap ordered based off of the location of the maximum in the first period. 
Returns the order of the genes as plotted.


**heatmap_order(data, order)**

input:
1. data frame with the time series data (index = gene names)
2. gene list in the desired order

Z-score normalizes the data frame and plots a heat map ordered based off of the ordered gene list.
Returns the order as plotted. 


**heatmap_LE(data, gene_list, first_period)**

input:
1. data frame with the time series data (index = gene names)
2. gene list with the genes to be plotted
3. first period

Z-score normalizes the data frame and plots a heat map ordered based off of the left edge in the first period.
Converted partially from R function heatOrderLE_firstper.  
Returns the order of the genes as plotted.


**Optional inputs for functions:**
1. yticks: gives option to include or remove y-ticks in the figure
default = False
2. cbar_bool: gives option to include or remove the color bar in the figure
default = True
3. axis: gives the option to include a subplot axis to the heatmap 
default = None
