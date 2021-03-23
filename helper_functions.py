from scipy import stats
import pandas as pd


def clean_column_names(dd_list):
	if(type(dd_list) == str):
		dd_list = [dd_list]

	return_list = []
	for dd in dd_list:
		for value in dd.split():
			if(value):
				return_list.append(value)

	return return_list

def find_group_names(dd_list, group_dict):
	if(type(dd_list) == str):
		dd_list = [dd_list]

	return_list = []
	for key in group_dict.keys():
		if(" ".join(group_dict[key]) == " ".join(dd_list)):
			return key

	return dd_list

def row_t_test(row, x_column_list, y_column_list):
    '''
    Performs a t-test for a row of data of group 1 columns vs group 2 columns
    '''
    x_distribution = list(row[x_column_list])
    y_distribution = list(row[y_column_list])

    t_statistic, p_value = stats.ttest_ind(x_distribution, y_distribution)

    return abs(t_statistic)

def clean_gene_list(dd_list):
	if(type(dd_list) == str):
		dd_list = [dd_list]

	return dd_list