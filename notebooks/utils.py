import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

state_polygons = pd.read_csv("../data/state-polygons-cleaned.csv")

def get_state(state):
    return state_polygons[state_polygons["state"] == state]

def plot_state(state):
    polygon = state_polygons[state_polygons["state"] == state]
    plt.plot("lon", "lat", data=polygon)
    plt.axis(option = "equal");
