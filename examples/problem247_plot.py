from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

root = Path(__file__).resolve().parent
csv_path = root / "problem247_trajectory_to_l4_artifacts" / "trajectory.csv"

df = pd.read_csv(csv_path)

earth_moon_distance = 3.844e5
earth_mass = 5.974e24
moon_mass = 7.348e22

pi1 = earth_mass / (earth_mass + moon_mass)
pi2 = moon_mass / (earth_mass + moon_mass)

earth_x = -pi2 * earth_moon_distance
moon_x = pi1 * earth_moon_distance

l4_x = (0.5 - pi2) * earth_moon_distance
l4_y = 3.0**0.5 * earth_moon_distance / 2.0

figure = go.Figure()

figure.add_trace(
    go.Scatter(
        x=df["x"],
        y=df["y"],
        mode="lines",
        name="Spacecraft trajectory",
        line={"color": "royalblue", "width": 2},
    )
)

figure.add_trace(
    go.Scatter(
        x=[earth_x],
        y=[0.0],
        mode="markers+text",
        text=["Earth"],
        textposition="bottom center",
        name="Earth",
        marker={"size": 12, "color": "deepskyblue"},
    )
)

figure.add_trace(
    go.Scatter(
        x=[moon_x],
        y=[0.0],
        mode="markers+text",
        text=["Moon"],
        textposition="bottom center",
        name="Moon",
        marker={"size": 8, "color": "lightgray"},
    )
)

figure.add_trace(
    go.Scatter(
        x=[0.0],
        y=[0.0],
        mode="markers+text",
        text=["Barycenter"],
        textposition="top center",
        name="Barycenter",
        marker={"size": 7, "color": "black", "symbol": "x"},
    )
)

figure.add_trace(
    go.Scatter(
        x=[l4_x],
        y=[l4_y],
        mode="markers+text",
        text=["L4"],
        textposition="top center",
        name="L4",
        marker={"size": 10, "color": "red", "symbol": "star"},
    )
)

limit = 1.1 * earth_moon_distance

figure.update_layout(
    title="Earth–Moon L4 transfer in the rotating barycentric frame",
    xaxis={
        "title": "x (km)",
        "range": [-limit, limit],
    },
    yaxis={
        "title": "y (km)",
        "range": [-limit, limit],
        "scaleanchor": "x",
        "scaleratio": 1,
    },
    width=900,
    height=900,
)

figure.show()
