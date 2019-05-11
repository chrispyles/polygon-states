import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

state_polygons = pd.read_csv("../data/state-polygons-cleaned.csv")

def plot_state(state):
    """
    Plots state as a polygon
    
    Arguments:
        state - name of state to plot
    """
    polygon = state_polygons[state_polygons["state"] == state]
    plt.plot("lon", "lat", data=polygon)
    plt.axis(option = "equal");

def write_to_utils(function):
    """
    Writes function to utils.py file
    
    Arguments:
        function - name of function to be written
    """
    func = inspect.getsource(function)
    with open("utils.py", "r+") as f:
        if func.split(" ")[1] not in f.read():
            with open("utils.py", "a") as f:
                f.write(func + "\n")

def get_state(state):
    """
    Returns portion of `state_polygons` df for specific state
    
    Arguments:
        state - name of state to select
    """
    return state_polygons[state_polygons["state"] == state]
