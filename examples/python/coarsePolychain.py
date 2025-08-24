import numpy as np
import matplotlib.pyplot as plt
from frenetTransform import _core as ft


def main():
    # create circle with radius 10 m
    radius = 10.0
    lengths_circle = np.linspace(0.0, 2 * np.pi, 101)
    circle_points_x = radius * np.cos(lengths_circle)
    circle_points_y = radius * np.sin(lengths_circle)

    # generate polyline from 4 points at 90° angles
    idcs_sub = np.arange(0, 101, 25)
    circle_poly = ft.Polychain(circle_points_x[idcs_sub], circle_points_y[idcs_sub])

    # get points along the polychain
    lengths_poly = np.linspace(0.0, 2 * np.pi * radius, 500)
    poly_points = circle_poly(lengths_poly)

    # instantiate transform
    transform = ft.Transform(circle_poly)

    # setup query points in Cartesian frame
    cartes_points = ft.Points([5.0, 12.0, -2.5, 0.3], [0.0, 12.0, 3.0, 11.5])

    # transform query points to Frenet frame
    frenet_points = transform.posFrenet(cartes_points)

    # get next points along polyline to query points
    proj_points = circle_poly(frenet_points.x())

    # plot next points
    plt.scatter(proj_points.x(), proj_points.y(), label="Projected Points")

    # get vectors from next points to query points
    normals = circle_poly.normal(frenet_points.x())
    normals = normals * frenet_points.y()

    # plot circle
    plt.plot(circle_points_x, circle_points_y, label="Circle")
    plt.gca().set_aspect("equal")

    # plot polychain
    plt.plot(poly_points.x(), poly_points.y(), label="Polychain")

    # plot query points
    plt.scatter(cartes_points.x(), cartes_points.y(), label="Query Points")

    # draw vectors
    plt.quiver(
        proj_points.x(),
        proj_points.y(),
        normals.x(),
        normals.y(),
        angles="xy",
        scale_units="xy",
        scale=1,
        width=0.005,
        color="r",
    )

    plt.legend(loc="lower right")
    plt.show()


if __name__ == "__main__":
    main()
