function output_standardized = standardize_data(input_data)

% Code for standardizing a data set
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

output_standardized = (input_data - mean(input_data))/std(input_data);