from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

example_directory = Path(__file__).resolve().parent
csv_path = example_directory / "artifacts" / "trajectory.csv"

df = pd.read_csv(csv_path)

frame_times = np.linspace(df["time"].iloc[0], df["time"].iloc[-1], 5000)

x = np.interp(frame_times, df["time"], df["x"])
y = np.interp(frame_times, df["time"], df["y"])
z = np.interp(frame_times, df["time"], df["z"])

trajectory = go.Scatter3d(
    x=x,
    y=y,
    z=z,
    mode="lines",
    line={"color": "lightgray", "width": 2},
    name="trajectory",
)

moving_point = go.Scatter3d(
    x=[x[0]],
    y=[y[0]],
    z=[z[0]],
    mode="markers",
    marker={"color": "red", "size": 5},
    name="current state",
)

frames = [
    go.Frame(
        name=f"{time:.2f}",
        data=[
            go.Scatter3d(
                x=[x_value],
                y=[y_value],
                z=[z_value],
                mode="markers",
                marker={"color": "red", "size": 5},
            )
        ],
        traces=[1],
    )
    for time, x_value, y_value, z_value in zip(frame_times, x, y, z)
]

figure = go.Figure(
    data=[trajectory, moving_point],
    frames=frames,
)

figure.update_layout(
    scene={
        "xaxis_title": "x",
        "yaxis_title": "y",
        "zaxis_title": "z",
    },
    updatemenus=[
        {
            "type": "buttons",
            "buttons": [
                {
                    "label": "Play",
                    "method": "animate",
                    "args": [
                        None,
                        {
                            "frame": {"duration": 30, "redraw": True},
                            "fromcurrent": True,
                        },
                    ],
                }
            ],
        }
    ],
)

figure.show()
