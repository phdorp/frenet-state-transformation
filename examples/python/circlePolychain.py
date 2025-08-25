from typing import List

import matplotlib.pyplot as plt
import numpy as np
from frenetTransform import _core as ft
from matplotlib.axes import Axes


def main():
    radius = 10.0
    lengthsCircle = np.linspace(0, 2 * np.pi, 101)
    circlePointsX = radius * np.cos(lengthsCircle)
    circlePointsY = radius * np.sin(lengthsCircle)

    bound = 15.0
    grid = np.linspace(-bound, bound, 61)
    posMeshX, posMeshY = np.meshgrid(grid, grid)

    cartesPoints = ft.Points(posMeshX.ravel(), posMeshY.ravel())
    cartesVels = ft.Points(
        np.ones(cartesPoints.numPoints()) / 2, np.ones(cartesPoints.numPoints()) / 2
    )
    cartesAccs = ft.Points(
        3 * np.ones(cartesPoints.numPoints()) / 4,
        -3 * np.ones(cartesPoints.numPoints()) / 4,
    )

    # Generate polyline along circle
    circlePoly = ft.Polychain(circlePointsX, circlePointsY)
    transform = ft.Transform(circlePoly)

    # Transform query points to Frenet frame
    frenetPointsTf = transform.posFrenet(cartesPoints)
    frenetVelsTf = transform.velFrenet(cartesVels, frenetPointsTf)
    frenetAccsTf = transform.accFrenet(cartesAccs, frenetVelsTf, frenetPointsTf)

    # Transform query points back to Cartesian frame
    cartesPointsTf = transform.posCartes(frenetPointsTf)
    cartesVelsTf = transform.velCartes(frenetVelsTf, frenetPointsTf)
    cartesAccsTf = transform.accCartes(frenetAccsTf, frenetVelsTf, frenetPointsTf)

    # Plot velocity vector field in Cartesian coordinates
    _, axs = plt.subplots(3, 2)
    axs: List[Axes] = axs.ravel().tolist()
    axs[0].plot(circlePointsX, circlePointsY, label="Circle")
    axs[0].quiver(cartesPoints.x(), cartesPoints.y(), cartesVels.x(), cartesVels.y())
    axs[0].set_xlabel("Cartesian x-axis/m")
    axs[0].set_ylabel("Cartesian y-axis/m")
    axs[0].set_title("Velocity transformation")

    # # Plot acceleration vector field in Cartesian coordinates
    axs[1].plot(circlePointsX, circlePointsY, label="Circle")
    axs[1].quiver(cartesPoints.x(), cartesPoints.y(), cartesAccs.x(), cartesAccs.y())
    axs[1].set_xlabel("Cartesian x-axis/m")
    axs[1].set_ylabel("Cartesian y-axis/m")
    axs[1].set_title("Acceleration transformation")

    # # Draw velocity vector field Frenet coordinates
    axs[2].quiver(
        frenetPointsTf.x(), frenetPointsTf.y(), frenetVelsTf.x(), frenetVelsTf.y()
    )
    axs[2].set_xlabel("Frenet x-axis/m")
    axs[2].set_ylabel("Frenet y-axis/m")

    # # Draw acceleration field Frenet coordinates
    axs[3].quiver(
        frenetPointsTf.x(), frenetPointsTf.y(), frenetAccsTf.x(), frenetAccsTf.y()
    )
    axs[3].set_xlabel("Frenet x-axis/m")
    axs[3].set_ylabel("Frenet y-axis/m")

    # # Draw velocity vector field Cartesian coordinates (transformed)
    axs[4].quiver(
        cartesPointsTf.x(), cartesPointsTf.y(), cartesVelsTf.x(), cartesVelsTf.y()
    )
    axs[4].set_xlabel("Cartesian x-axis/m")
    axs[4].set_ylabel("Cartesian y-axis/m")

    # # Draw acceleration field Cartesian coordinates (transformed)
    axs[5].quiver(
        cartesPointsTf.x(), cartesPointsTf.y(), cartesAccsTf.x(), cartesAccsTf.y()
    )
    axs[5].set_xlabel("Cartesian x-axis/m")
    axs[5].set_ylabel("Cartesian y-axis/m")

    plt.show()


if __name__ == "__main__":
    main()
