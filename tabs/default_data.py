import pandas as pd
import pathlib
from models.expression import expression

default_groups = {
    'CONTROL':['CTRL_0', 'CTRL_1', 'CTRL_2'],
    'P4':['P4_0', 'P4_1', 'P4_2'],
    'IL1B':['IL1B_0', 'IL1B_1', 'IL1B_2'],
    'P4_IL1B':['P4_IL1B_0', 'P4_IL1B_1', 'P4_IL1B_2'],
    'FSK':['FSK_0', 'FSK_1', 'FSK_2'],
    'FSK_P4':['FSK_P4_0', 'FSK_P4_1', 'FSK_P4_2'],
    'FSK_IL1B': ['FSK_IL1B_0','FSK_IL1B_1', 'FSK_IL1B_2',],
    'FSK_P4_IL1B': ['FSK_P4_IL1B_0', 'FSK_P4_IL1B_1','FSK_P4_IL1B_2']
}

df = pd.read_csv(str(pathlib.Path("data/GSE134896_FSKrepExpTable.csv")))

default_exp = expression(df)
default_exp.organize_from_initial_upload(default_groups)
