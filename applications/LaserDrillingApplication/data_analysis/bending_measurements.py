import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splprep, splev, griddata
from skimage.measure import EllipseModel
import csv
from math import ceil, sqrt
import argparse
from pathlib import Path



# ==== Ellipse Fitting Function ====
def fit_ellipse_least_squares(points_2d):
    """
    Fit an ellipse to all given 2D points using least-squares method.

    Parameters:
    - points_2d: Nx2 array of (x, y) points.

    Returns:
    - Tuple: (cx, cy, a, b, eccentricity, theta) or None if fitting fails.
    """
    if len(points_2d) < 5:
        raise ValueError("Not enough points to fit the ellipse")

    model = EllipseModel()
    if not model.estimate(points_2d):
        raise ValueError("Couldn't fit the ellipse to the data")

    cx, cy, a, b, theta = model.params  # (x, y, a, b, theta)
    if b > a:
        a, b = b, a
        theta += np.pi / 2

    eccentricity = np.sqrt(1 - (b**2 / a**2))
    return (cx, cy, a, b, eccentricity, theta)


# ==== Plotting Functions ====
def plot_3d_geometry(x, y, z, spline, cx, cy, cz, ellipse_spline, ecx, ecy, ecz, slice_bounds, plot_planes, sample_limits):
    """
    Plot 3D point cloud, medial spline, ellipse center spline, centroids, ellipse centers, and optional slicing planes.

    Parameters:
    - x, y, z: Arrays of point coordinates.
    - spline: Tuple of (x, y, z) medial spline coordinates (through centroids).
    - cx, cy, cz: Arrays of centroid coordinates.
    - ellipse_spline: Tuple of (x, y, z) ellipse center spline coordinates.
    - ecx, ecy, ecz: Arrays of ellipse center coordinates.
    - slice_bounds: Array of z-boundaries for slices.
    - plot_planes: Boolean to toggle slicing plane plotting.
    - sample_limits: Coordinates of the corners of the sample
    """
    sample_x_min, sample_x_max, sample_y_min, sample_y_max = sample_limits

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z, s=1, alpha=0.3, label="Point Cloud")
    ax.plot(*spline, c="r", lw=2, label="Medial Spline")
    ax.scatter(cx, cy, cz, c="k", label="Centroids")
    if ellipse_spline is not None:
        ax.plot(*ellipse_spline, c="g", lw=2, label="Ellipse Center Spline")
        ax.scatter(ecx, ecy, ecz, c="g", marker="*", s=100, label="Ellipse Centers")

    if plot_planes:
        for bound in slice_bounds:
            xx, yy = np.meshgrid(
                np.linspace(sample_x_min, sample_x_max, 10), np.linspace(sample_y_min, sample_y_max, 10)
            )
            zz = np.full_like(xx, bound)
            ax.plot_surface(xx, yy, zz, alpha=0.1, color="gray")

    ax.set_title("3D Hole Geometry with Splines\n" + filename)
    ax.set_xlabel("X (um)")
    ax.set_ylabel("Y (um)")
    ax.set_zlabel("Z (um)")
    ax.legend()
    plt.show()


def plot_outliers(x, y, z, surface_outliers, deep_outliers, plot_outlier_types):
    """
    Plot used points and optional surface/deep outliers in 3D.

    Parameters:
    - x, y, z: Arrays of used point coordinates.
    - surface_outliers: Array of surface outlier points.
    - deep_outliers: Array of deep outlier points.
    - plot_outlier_types: String ("all", "surface", "deep") to select outlier types.
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z, s=1, label="Used Points", alpha=0.4)
    if plot_outlier_types in ["all", "surface"] and len(surface_outliers) > 0:
        ax.scatter(
            surface_outliers[:, 0],
            surface_outliers[:, 1],
            surface_outliers[:, 2],
            c="orange",
            label="Surface Outliers",
            s=5,
        )
    if plot_outlier_types in ["all", "deep"] and len(deep_outliers) > 0:
        ax.scatter(
            deep_outliers[:, 0],
            deep_outliers[:, 1],
            deep_outliers[:, 2],
            c="purple",
            label="Deep Outliers",
            s=5,
        )
    ax.set_title("Filtered vs. Kept Points\n" + filename)
    ax.set_xlabel("X (um)")
    ax.set_ylabel("Y (um)")
    ax.set_zlabel("Z (um)")
    ax.legend()
    plt.show()


def plot_xy_centers_path(
    plot_xy_centroids, plot_xy_ellipse_centers, cx, cy, cz, centroid_ids, ellipse_centers, sample_limits):
    """
    Plot XY projection of centroid and ellipse center paths with annotations.

    Parameters:
    - plot_xy_centroids, plot_xy_ellipse_centers: toggle the polotting of the centroids and of the ellipse centers resp.
    - cx, cy, cz: Arrays of centroid coordinates.
    - centroid_ids: List of slice indices (1-based).
    - ellipse_centers: Array of ellipse center coordinates.
    """
    sample_x_min, sample_x_max, sample_y_min, sample_y_max = sample_limits

    if plot_xy_centroids:
        if plot_xy_ellipse_centers:
            title = "XY Projection of centroid and ellipse center paths\n"
        else:
            title = "XY Projection of centroid path\n"
    else:
        if plot_xy_ellipse_centers:
            title = "XY Projection of ellipse center path\n"
        else:
            print(
                "Are you sure you wanted an empty plot? (plot_xy_centers_path is being passed a False to both plotting centroids and ellipses' centers)"
            )
            return  # Nothing to plot

    title += filename

    fig, ax = plt.subplots(figsize=(8, 8))

    if plot_xy_centroids:
        ax.plot(cx, cy, "-o", color="blue", label="Centroid Path")
        for i, (x_, y_, z_, idx) in enumerate(zip(cx, cy, cz, centroid_ids)):
            label = f"{idx}"
            ax.annotate(label, (x_, y_), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=8)
    

    if plot_xy_ellipse_centers:
        if len(ellipse_centers) > 0:
            ax.plot(ellipse_centers[:, 0], ellipse_centers[:, 1], "-*", color="green", label="Ellipse Centers")
            for i, (x_, y_, z_, idx) in enumerate(
                zip(ellipse_centers[:, 0], ellipse_centers[:, 1], ellipse_centers[:, 2], centroid_ids)
            ):
                label = f"{idx}"
                ax.annotate(label, (x_, y_), textcoords="offset points", xytext=(0, -12), ha="center", fontsize=8)
        else:
            raise ValueError("There are no ellipse centers")

    ax.set_title(title)
    ax.set_xlabel("X (um)")
    ax.set_ylabel("Y (um)")
    ax.axis("square")
    ax.set_xlim(sample_x_min, sample_x_max)
    ax.set_ylim(sample_y_min, sample_y_max)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.show()


def plot_superposed_slices_xy(
    slices,
    centroids,
    centroid_ids,
    ellipses,
    sample_limits,
    plot_centroids_in_all_slices_xy,
    plot_ellipses_in_all_slices_xy,
    plot_ellipse_centers
):
    """
    Plot XY projection of all slices with optional centroids and ellipses.

    Parameters:
    - slices: List of point arrays for each slice (slices[i] for slice i+1).
    - centroids: List of centroid coordinates for each slice.
    - centroid_ids: List of slice indices (1-based).
    - ellipses: List of ellipses.
    - num_slices: Total number of slices.
    - plot_centroids_in_all_slices_xy: Boolean to toggle centroid plotting.
    - plot_ellipses_in_all_slices_xy: Boolean to toggle ellipse plotting.
    - sample_limits: (sample_x_min, sample_x_max, sample_y_min, sample_y_max)
    """

    num_slices = len(slices)
    sample_x_min, sample_x_max, sample_y_min, sample_y_max = sample_limits

    fig, ax = plt.subplots(figsize=(8, 8))
    colors = plt.cm.viridis(np.linspace(0, 1, num_slices))  # Distinct colors for slices
    for i in range(num_slices):
        slice_points = slices[i]
        centroid = centroids[i]
        idx = centroid_ids[i]
        if len(slice_points) > 0 and centroid is not None:
            ax.scatter(slice_points[:, 0], slice_points[:, 1], s=3, c=[colors[i]], label=f"Slice {idx}")
            if plot_centroids_in_all_slices_xy:
                ax.plot(centroid[0], centroid[1], "o", c=colors[i], markersize=8)
                label = f"Slice {idx}\nZ={centroid[2]:.2f}"
                ax.annotate(
                    label,
                    (centroid[0], centroid[1]),
                    textcoords="offset points",
                    xytext=(0, 6),
                    ha="center",
                    fontsize=8,
                )

            ellipse = ellipses[i]
            if plot_ellipses_in_all_slices_xy and ellipse is not None:
                center = ellipse["center"]
                cx, cy = center[0], center[1]
                a = ellipse["a"]
                b = ellipse["b"]
                theta = ellipse["angle"]

                t = np.linspace(0, 2 * np.pi, 100)

                x_ellipse = cx + a * np.cos(t) * np.cos(theta) - b * np.sin(t) * np.sin(theta)
                y_ellipse = cy + a * np.cos(t) * np.sin(theta) + b * np.sin(t) * np.cos(theta)
                ax.plot(x_ellipse, y_ellipse, "-", color="black", linewidth=3.5)
                ax.plot(x_ellipse, y_ellipse, "-", color="white", linewidth=2)
                # Annotate ellipse at a point (e.g., t=0) with slice index
                ax.annotate(
                    f"{idx}",
                    (x_ellipse[0], y_ellipse[0]),
                    textcoords="offset points",
                    xytext=(5, 5),
                    ha="center",
                    fontsize=8,
                    color="black",
                    bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
                )
                if plot_ellipse_centers:
                    ax.plot(cx, cy, "g*", markersize=10)
                    ax.annotate(
                        f"{idx}",
                        (cx, cy),
                        textcoords="offset points",
                        xytext=(5, 5),
                        ha="center",
                        fontsize=8,
                        color="black",
                        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
                    )
    ax.set_title("XY Projection of All Slices\n" + filename)
    ax.set_xlabel("X (um)")
    ax.set_ylabel("Y (um)")
    ax.axis("square")
    ax.set_xlim(sample_x_min, sample_x_max)
    ax.set_ylim(sample_y_min, sample_y_max)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.show()


def plot_ellipses_and_axes(
    ellipses,
    sample_limits,
    shift_to_common_center=False
):
    """
    Plot all ellipses, their centers, and major axes in the same XY plot, with optional centering.

    Parameters:
    - ellipses: List of ellipses.
    - num_slices: Total number of slices.
    - shift_to_common_center: Boolean to shift all ellipses to a common center (default: False).
    - sample_limits = (sample_x_min, sample_x_max, sample_y_min, sample_y_max)
    """
    num_ellipses = len(ellipses)
    sample_x_min, sample_x_max, sample_y_min, sample_y_max = sample_limits

    fig, ax = plt.subplots(figsize=(8, 8))
    colors = plt.cm.viridis(np.linspace(0, 1, num_ellipses))  # Distinct colors for slices
    common_center_x, common_center_y = ((sample_x_max + sample_x_min) / 2, (sample_y_max + sample_y_min) / 2)

    for i in range(num_ellipses):
        ellipse = ellipses[i]
        if ellipse is None:
            raise ValueError("An ellipse is None")

        center = ellipse["center"]
        cx, cy = center[0], center[1]
        a = ellipse["a"]
        b = ellipse["b"]
        theta = ellipse["angle"]

        idx = i + 1

        # Adjust coordinates if shifting to common center
        plot_cx, plot_cy = (common_center_x, common_center_y) if shift_to_common_center else (cx, cy)

        # Plot ellipse
        t = np.linspace(0, 2 * np.pi, 100)
        x_ellipse = plot_cx + a * np.cos(t) * np.cos(theta) - b * np.sin(t) * np.sin(theta)
        y_ellipse = plot_cy + a * np.cos(t) * np.sin(theta) + b * np.sin(t) * np.cos(theta)
        ax.plot(x_ellipse, y_ellipse, "-", color=colors[i], linewidth=2, label=f"Slice {idx}")

        # Plot ellipse center
        ax.plot(plot_cx, plot_cy, "*", color=colors[i], markersize=10)

        # Plot major axis (line through center along theta with length 2*a)
        major_axis_x = [plot_cx - a * np.cos(theta), plot_cx + a * np.cos(theta)]
        major_axis_y = [plot_cy - a * np.sin(theta), plot_cy + a * np.sin(theta)]
        ax.plot(major_axis_x, major_axis_y, "--", color=colors[i], linewidth=1.5)

        # Annotate ellipse with slice index
        ax.annotate(
            f"{idx}",
            (x_ellipse[0], y_ellipse[0]),
            textcoords="offset points",
            xytext=(5, 5),
            ha="center",
            fontsize=8,
            color="black",
            bbox=dict(facecolor="white", alpha=0.9, edgecolor="none"),
        )

    # Set title and optional subtitle
    title = "XY projection of all ellipses with centers and major axes\n" + filename
    if shift_to_common_center:
        title += "\nEllipses shifted to common center"
    ax.set_title(title)
    ax.set_xlabel("X (um)")
    ax.set_ylabel("Y (um)")
    ax.axis("square")
    ax.set_xlim(sample_x_min, sample_x_max)
    ax.set_ylim(sample_y_min, sample_y_max)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.show()


def plot_contour(x, y, z, contour_levels, sample_limits):
    """
    Plot filled contour plot of hole depth in XY projection.

    Parameters:
    - x, y, z: Arrays of point coordinates.
    - contour_levels: Number of contour levels.
    - sample_limits = (sample_x_min, sample_x_max, sample_y_min, sample_y_max)
    """
    
    sample_x_min, sample_x_max, sample_y_min, sample_y_max = sample_limits

    x_grid = np.linspace(sample_x_min, sample_x_max, 100)
    y_grid = np.linspace(sample_y_min, sample_y_max, 100)
    X, Y = np.meshgrid(x_grid, y_grid)
    points = np.vstack((x, y)).T
    Z = griddata(points, z, (X, Y), method="linear", fill_value=np.nan)

    fig, ax = plt.subplots(figsize=(8, 8))
    contour = ax.contourf(X, Y, Z, levels=contour_levels, cmap="viridis_r")
    ax.set_title("Filled Contour Plot of Hole Depth (XY Projection)\n" + filename)
    ax.set_xlabel("X (um)")
    ax.set_ylabel("Y (um)")
    ax.axis("square")
    ax.set_xlim(sample_x_min, sample_x_max)
    ax.set_ylim(sample_y_min, sample_y_max)
    fig.colorbar(contour, ax=ax, label="Depth (um)")
    ax.grid(True, alpha=0.3)
    plt.show()


def plot_ellipse_metrics(ellipses):
    """
    Plot eccentricity and major axis angle vs. depth.

    Parameters:
    - ellipse_params: List of ellipse parameters or None.
    """
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 6), sharex=True)

    eccentricities = [ellipse["eccentricity"] for ellipse in ellipses]
    thetas = [np.degrees(ellipse["angle"]) % 360 for ellipse in ellipses]  # Convert to degrees
    depths = [ellipse["center"][2] for ellipse in ellipses]

    ax1.plot(depths, eccentricities, ".-", color="blue", markersize=12, label="Eccentricity")
    ax1.set_ylabel("Eccentricity")
    ax1.set_title("Ellipse eccentricity vs. depth\n" + filename)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 1)
    ax1.legend(loc="upper right")

    ax2.plot(depths, thetas, ".-", color="red", markersize=12, label="Major axis angle (0 to 180º)")
    ax2.set_xlabel("Depth (Z, um)")
    ax2.set_ylabel("Major axis angle (degrees)")
    ax2.set_title("Ellipse major axis angle vs. depth\n" + filename)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 180)
    ax2.set_yticks([30 * i for i in range(7)])
    ax2.legend(loc="upper right")

    shifted_thetas = [thetas[i] if thetas[i] <= 90 else thetas[i] - 180 for i in range(len(thetas))]
    ax3.plot(depths, shifted_thetas, ".-", color="green", markersize=12, label="Major axis angle (-90º to 90º)")
    ax3.set_xlabel("Depth (Z, um)")
    ax3.set_ylabel("Major axis angle (degrees)")
    ax3.set_title("Ellipse major axis angle vs. depth\n" + filename)
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(-90, 90)
    ax3.set_yticks([-90 + 30 * i for i in range(7)])
    ax3.legend(loc="upper right")

    plt.tight_layout()
    plt.show()


def plot_individual_slices_grid(slices, slice_bounds, centroids, ellipses, sample_limits, plot_ellipses_in_slices):
    """
    Plot individual slices in a subplot grid with points, centroids, and optional ellipses.

    Parameters:
    - slices: List of point arrays for each slice (slices[i] for slice i+1).
    - slice_bounds: Array of z-boundaries for slices.
    - centroids: List of centroid coordinates for each slice.
    - ellipses: List of ellipses.
    - plot_ellipses_in_slices: Boolean to toggle ellipse plotting.
    """
    num_slices = len(slices)
    valid_slices = [i for i in range(num_slices) if len(slices[i]) > 0]
    n_plots = len(valid_slices)
    n_cols = ceil(sqrt(n_plots))
    n_rows = ceil(n_plots / n_cols)

    
    x_min, x_max, y_min, y_max = sample_limits

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
    # plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=10)
    axes = axes.flatten() if n_plots > 1 else [axes]
    handles, labels = [], []  # For collecting legend handles and labels

    plot_idx = 0
    for i in valid_slices:
        slice_points = slices[i]
        if len(slice_points) == 0:
            continue

        z_start, z_end = slice_bounds[i], slice_bounds[i + 1]
        centroid = centroids[i]
        ax = axes[plot_idx]
        # Plot slice points
        sc = ax.scatter(slice_points[:, 0], slice_points[:, 1], s=3, c="gray")
        if plot_idx == 0:
            handles.append(sc)
            labels.append("Slice Points")
        # Plot centroid
        pt = ax.plot(centroid[0], centroid[1], "ro")[0]
        if plot_idx == 0:
            handles.append(pt)
            labels.append("Centroid")
        # Plot ellipse and center if enabled
        ellipse = ellipses[i]
        if plot_ellipses_in_slices and ellipse is not None:
            center = ellipse["center"]
            cx, cy = center[0], center[1]
            a = ellipse["a"]
            b = ellipse["b"]
            theta = ellipse["angle"]

            t = np.linspace(0, 2 * np.pi, 100)
            x_ellipse = cx + a * np.cos(t) * np.cos(theta) - b * np.sin(t) * np.sin(theta)
            y_ellipse = cy + a * np.cos(t) * np.sin(theta) + b * np.sin(t) * np.cos(theta)
            el = ax.plot(x_ellipse, y_ellipse, "g-")[0]
            ec = ax.plot(cx, cy, "g*", markersize=10)[0]
            if plot_idx == 0:
                handles.append(el)
                labels.append("Fitted Ellipse")
                handles.append(ec)
                labels.append("Ellipse Center")
        ax.set_title(f"Slice {valid_slices[plot_idx] + 1} (z ∈ [{z_start:.2f}, {z_end:.2f}] um)\n" + filename)
        ax.set_xlabel("X (um)")
        ax.set_ylabel("Y (um)")
        ax.axis("square")
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.grid(True, alpha=0.3)
        plot_idx += 1

    for j in range(plot_idx, len(axes)):
        axes[j].axis("off")

    if handles:  # Only add legend if there are valid handles
        fig.legend(handles, labels, loc="lower center", ncol=len(handles), bbox_to_anchor=(0.5, 0.01))

    plt.suptitle("XY Projection of Slices with Centroids and Ellipses\n" + filename, fontsize=14)
    plt.subplots_adjust(top=0.87, bottom=0.1)  # Adjust bottom to accommodate legend
    plt.tight_layout()
    plt.show()


def plot_individual_slices_separately(
    slices,
    slice_bounds,
    centroids,
    centroid_ids,
    ellipses,
    sample_limits,
    plot_ellipses_in_slices
):
    """
    Plot each valid slice in a separate figure with points, centroid, and optional ellipse.

    Parameters:
    - slices: List of point arrays for each slice (slices[i] for slice i+1).
    - slice_bounds: Array of z-boundaries for slices.
    - centroids: List of centroid coordinates for each slice.
    - centroid_ids: List of slice indices (1-based).
    - ellipses: List of ellipse parameters or None.
    - sample_limits = (sample_x_min, sample_x_max, sample_y_min, sample_y_max).
    - plot_ellipses_in_slices: Boolean to toggle ellipse plotting.
    """
    
    sample_x_min, sample_x_max, sample_y_min, sample_y_max = sample_limits

    num_slices = len(slices)
    for i in range(num_slices):
        slice_points = slices[i]
        if len(slice_points) == 0:
            continue

        z_start, z_end = slice_bounds[i], slice_bounds[i + 1]
        centroid = centroids[i]
        idx = centroid_ids[i]

        # Create new figure for each slice
        fig, ax = plt.subplots(figsize=(6, 6))
        handles, labels = [], []  # For legend

        # Plot slice points
        sc = ax.scatter(slice_points[:, 0], slice_points[:, 1], s=3, c="gray")
        handles.append(sc)
        labels.append("Slice Points")

        # Plot centroid
        pt = ax.plot(centroid[0], centroid[1], "ro")[0]
        handles.append(pt)
        labels.append("Centroid")

        # Plot ellipse and center if enabled
        ellipse = ellipses[i]
        if plot_ellipses_in_slices and ellipse is not None:
            center = ellipse["center"]
            cx, cy = center[0], center[1]
            a = ellipse["a"]
            b = ellipse["b"]
            theta = ellipse["angle"]


            t = np.linspace(0, 2 * np.pi, 100)
            x_ellipse = cx + a * np.cos(t) * np.cos(theta) - b * np.sin(t) * np.sin(theta)
            y_ellipse = cy + a * np.cos(t) * np.sin(theta) + b * np.sin(t) * np.cos(theta)
            el = ax.plot(x_ellipse, y_ellipse, "g-")[0]
            ec = ax.plot(cx, cy, "g*", markersize=10)[0]
            handles.append(el)
            labels.append("Fitted Ellipse")
            handles.append(ec)
            labels.append("Ellipse Center")

        # Set plot properties
        ax.set_title(f"Slice {idx} (z ∈ [{z_start:.2f}, {z_end:.2f}] um)\n" + filename)
        ax.set_xlabel("X (um)")
        ax.set_ylabel("Y (um)")
        ax.axis("square")
        ax.set_xlim(sample_x_min, sample_x_max)
        ax.set_ylim(sample_y_min, sample_y_max)
        ax.grid(True, alpha=0.3)

        # Add legend
        if handles:
            ax.legend(handles, labels, loc="lower center", ncol=len(handles), bbox_to_anchor=(0.5, -0.1))

        plt.tight_layout()
        plt.show()


def calculate_chord_length(a):
    """
    Calculate the chord length of a curve, that is, the distance from its endpoints,
    given by points in array a of shape (N,3) where the columns represent
    the (x,y,z) coordinates of the points
    """
    ax, ay, az = a[:, 0], a[:, 1], a[:, 2]
    return np.linalg.norm([ax[-1] - ax[0], ay[-1] - ay[0], az[-1] - az[0]])


def calculate_arc_length(a):
    """
    Calculate the arc-length of a curve given by points in array a of shape (N,3)
    where the columns represent the (x,y,z) coordinates of the points
    """
    return np.sum(np.linalg.norm(np.diff(np.array(a), axis=0), axis=1))


def calculate_tortuosity(a):
    """
    Calculate the tortuosity of a curve given by points in array a of shape (N,3)
    where the columns represent the (x,y,z) coordinates of the points
    """
    arc_length = calculate_arc_length(a)
    chord_length = calculate_chord_length(a)
    return arc_length / chord_length


def fit_spline_3d(x, y, z, t_arr, spline_order=3):
    """
    Fits a spline through 3D points defined by x, y, z coordinates.

    Parameters:
        x, y, z (array-like): Coordinates of the points.
        t_arr (array of float): Values of the spline parameter where to evaluate the spline.
        spline_order (int): The order of the spline (default: 3).

    Returns:
        spline_points (np.ndarray): Array of shape (num_points, 3) with the fitted spline.
        tck (tuple): The B-spline representation returned by splprep.

    Raises:
        ValueError: If fewer than 2 points are provided or fitting fails.
    """
    if len(x) < 2:
        raise ValueError(f"Need at least 2 points to fit a spline, got {len(x)}.")

    try:
        k = min(spline_order, len(x) - 1)
        tck, _ = splprep([x, y, z], s=0, k=k)
        spline = splev(t_arr, tck)
        spline_points = np.array(spline).T
        return spline, spline_points, tck
    except Exception as e:
        raise ValueError(f"Failed to fit spline: {str(e)}")


def compute_derivatives_spline(t_arr, tck):
    """
    Given a list of points and the knots of a spline, it computes the first,
    second and third derivatives of the spline at said points.
    """
    dx, dy, dz = splev(t_arr, tck, der=1)
    ddx, ddy, ddz = splev(t_arr, tck, der=2)
    dddx, dddy, dddz = splev(t_arr, tck, der=3)

    # Each row represents the (n-th order) derivative at a point in the spline (column 1 = x, column 2 = y, column 3 = z)
    D1 = np.vstack((dx, dy, dz)).T
    D2 = np.vstack((ddx, ddy, ddz)).T
    D3 = np.vstack((dddx, dddy, dddz)).T

    return D1, D2, D3


def compute_curvature(D1, D2):
    """
    Computes the local curvature of a curve given its first (D1) and second (D2) derivatives
    """

    cross = np.cross(D1, D2)
    numerator = np.linalg.norm(cross, axis=1)
    denominator = np.linalg.norm(D1, axis=1) ** 3
    curvature = numerator / denominator
    curvature = np.where(np.isfinite(curvature), curvature, 0)  # Handle division by zero

    return curvature


def compute_torsion(D1, D2, D3):
    """
    Computes the local torsion of a curve given its first (D1), second (D2) and third (D3) derivatives
    """
    cross = np.cross(D1, D2)
    torsion_numerator = np.einsum("ij,ij->i", cross, D3)
    torsion_denominator = np.linalg.norm(cross, axis=1) ** 2
    torsion = torsion_numerator / torsion_denominator
    torsion = np.where(np.isfinite(torsion), torsion, 0)  # Handle division by zero

    return torsion


def read_single_bore(filepath):
    """
    Reads the data of a single bore (from files like those named 10W5P02.dat)
    """
    data = []
    try:
        with open(filepath, "r") as file:
            reader = csv.reader(file, delimiter=" ")
            for row in reader:
                try:
                    point = [float(coord) for coord in row if coord.strip()]
                    if len(point) == 3 and all(np.isfinite(point)):
                        data.append(point)
                except ValueError:
                    continue
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: File {filepath} not found.")

    data = np.array(data)

    return data


def clean_data(data):
    """
    Removes invalid data
    """
    total_points = len(data)

    # Remove non-numerical or NaN entries
    valid_mask = np.isfinite(data).all(axis=1)
    data_clean = data[valid_mask]
    num_nan = total_points - len(data)
    print(f"Discarded {num_nan} non-numerical/NaN points ({100 * num_nan / total_points:.2f}%)")

    return data_clean


def remove_surface_and_outliers(data, shallow_z_threshold, deep_outlier_quantile_threshold):
    """
    Removes the points corresponding to the surface by removing the points shallower than
    shallow_z_threshold (points closer to the surface than shallow_z_threshold), and removes
    outliers that lie too deep in the sample that correspond to wrongly measured data by
    removing the deep_outlier_quantile_threshold deepest points
    """

    z_vals = data[:, 2]

    # If shallow_z_threshold is None, do not exclude any points from the surface
    surface_mask = np.zeros_like(z_vals, dtype=bool) if shallow_z_threshold is None else z_vals > shallow_z_threshold

    # Similarly, to not exclude deep outliers
    if deep_outlier_quantile_threshold is None:
        deep_mask = np.zeros_like(z_vals, dtype=bool)
    else:
        z_deep_cutoff = np.quantile(z_vals, deep_outlier_quantile_threshold)
        deep_mask = z_vals < z_deep_cutoff

    keep_mask = ~(surface_mask | deep_mask)

    surface_outliers = data[surface_mask]
    deep_outliers = data[deep_mask]
    data = data[keep_mask]

    num_surface_outliers = len(surface_outliers)
    num_deep_outliers = len(deep_outliers)
    total_points = len(data)

    print(f"Discarded {num_surface_outliers} surface outliers({100 * num_surface_outliers / total_points:.2f}%)")
    print(f"Discarded {num_deep_outliers} deep outliers({100 * num_deep_outliers / total_points:.2f}%)")

    return data, surface_outliers, deep_outliers


def subsample_data(data, max_points, random_seed):
    """
    Subsamples the data by randomly keeping max_points and discarding the rest
    """

    np.random.seed(random_seed)
    indices = np.random.choice(len(data), max_points, replace=False)
    data = data[indices]

    return data

def calculate_slices(data, slice_thickness):
    """
    Calculates the data slices
    """
    z = data[:, 2]

    # ==== Create Slices Variable ====
    z_min, z_max = z.min(), z.max()
    print(f"z_max={float(z_max)}, z_min={float(z_min)}")
    num_slices = ceil((z_max - z_min)/slice_thickness)

    # Generate slice boundaries from top (z_max) to bottom (z_min)
    slice_bounds, step = np.linspace(z_max, z_max - num_slices*slice_thickness, num_slices+1, retstep=True)

    
    # Initialize list to store slices
    slices = []

    # Create slices
    for i in range(num_slices):
        z_start, z_end = slice_bounds[i], slice_bounds[i + 1]
        
        # Mask for points in the current slice
        mask = (z <= z_start) & (z > z_end)
        
        slices.append(data[mask])

    # Optional: Remove empty slices (unlikely, but safe)
    for s in slices:
        if len(s) <= 0:
            raise ValueError("There are empty slices")

    slices = [s for s in slices if len(s) > 0]


    print(f"{slice_bounds=}")
    print(f"{len(slice_bounds)=}")
    print(f"{num_slices=}")

    return slices, slice_bounds

def compute_centroids(slices, slice_bounds):
    """
    For each slice, finds its centroid
    """
    num_slices = len(slices)

    centroids = [None] * num_slices
    centroid_zs = [None] * num_slices
    centroid_ids = [None] * num_slices

    for i in range(num_slices):
        print(f"Computing centroid for Slice {i + 1} (z ∈ [{slice_bounds[i]:.2f}, {slice_bounds[i + 1]:.2f}] um)")
        slice_points = slices[i]

        if len(slice_points) == 0:
            raise ValueError(
                f"Slice {i + 1} (z ∈ [{slice_bounds[i]:.2f}, {slice_bounds[i + 1]:.2f}] um) contains no points"
            )

        centroid = np.mean(slice_points, axis=0)

        centroids[i] = centroid
        centroid_ids[i] = i + 1

    # Filter out None values
    valid_centroids = [c for c in centroids if c is not None]
    centroids = np.array(valid_centroids)
    centroid_zs = [z for z in centroid_zs if z is not None]
    centroid_ids = [id for id in centroid_ids if id is not None]

    return centroids, centroid_ids


def compute_ellipses(slices, slice_bounds):
    """
    Find the ellipse that best fits each spline using least squares
    """
    num_slices = len(slices)

    ellipses = [None] * num_slices

    for i in range(num_slices):
        print(f"Fitting ellipse for Slice {i + 1} (z ∈ [{slice_bounds[i]:.2f}, {slice_bounds[i + 1]:.2f}] um)")
        slice_points = slices[i]

        if len(slice_points) == 0:
            raise ValueError("Empty slice")

        points_2d = slice_points[:, :2]  # [x, y]  

        center_x, center_y, a, b, eccentricity, angle = fit_ellipse_least_squares(points_2d)
        center_z = np.mean(slice_points[:, 2]) # Set the ellipse's depth as the "center of mass" in depth of the points in the slice

        ellipse_center = np.array([center_x, center_y, center_z])

        ellipse_dict = {"center": ellipse_center, "a": a, "b": b, "eccentricity": eccentricity, "angle": angle}
        ellipses[i] = ellipse_dict            
            

    return ellipses


if __name__ == "__main__":
    # Read filepath as argument
    parser = argparse.ArgumentParser(description="Make plots and extract metrics of the measurements of a bore.")
    parser.add_argument("--filepath", required=True, help="Path to file containing the data file.")

    args = parser.parse_args()

    filepath = args.filepath

    filename = Path(filepath).stem
    


    # ==== Configurable Parameters ====
    sample_x_min = 0
    sample_x_max = 50
    sample_y_min = 0
    sample_y_max = 50
    sample_limits = (sample_x_min, sample_x_max, sample_y_min, sample_y_max)

    max_points = 500000  # Max number of points to use for entire dataset
    subsample_points = False  # Toggle for random subsampling of entire dataset

    slice_thickness = 0.8 # Thickness of each slice in um
    deep_outlier_quantile_threshold = (
        0.0005  # Lower quantile for discarding deep outliers. Set it to None to not discard any.
    )

    shallow_z_threshold = (
        -1
    )  # Discard points shallower than this z (closer to surface). Set it to None to not discard any.
    random_seed = 42  # Seed for reproducibility
    contour_levels = 10  # Number of contour levels for depth plot

    spline_order = 3  # Spline interpolation order
    n_spline_points = 500  # Number of points in the splines

    # Plot toggles
    plot_3d_geometry_toggle = True  # Plot of the point cloud with the fitted spline(s) through the center(s)
    plot_planes = True
    plot_outliers_toggle = False
    plot_outlier_types = "deep"  # Options: "all", "surface", "deep"
    plot_xy_center_approximations = True
    plot_xy_centroids = False
    plot_xy_ellipse_centers = True
    plot_z_histogram = False  # Optional histogram to guide shallow_z_threshold
    plot_superposed_slices_xy_toggle = True  # Toggle for XY view of all slices superposed
    plot_centroids_in_superposed_slices_xy = False  # Toggle for centroids in XY view of superposed slices
    plot_ellipses_in_superposed_slices_xy = True  # Toggle for centroids in XY view of superposed slices
    plot_ellipse_centers_in_superposed_slices_xy = True  # Toggle for ellipses' centers in XY view of superposed slices
    plot_ellipses_and_axes_toggle = True  # Toggle for plot the of all ellipses
    ellipses_shift_to_common_center = False  # Toggle to shift all ellipses so that they have the same center
    plot_contour_toggle = False  # Toggle for contour plot
    plot_ellipses_in_slices = True  # Toggle for plotting ellipses in slice views
    plot_ellipse_metrics_toggle = True  # Toggle for eccentricity and theta vs. depth plot
    plot_individual_slices_grid_toggle = True
    plot_slice_views_one_by_one = False  # Toggle for XY view of each slices

    # ==== Read and Clean the Data ====
    data_raw = read_single_bore(filepath)
    data = clean_data(data_raw)
    data, surface_outliers, deep_outliers = remove_surface_and_outliers(
        data, shallow_z_threshold, deep_outlier_quantile_threshold
    )

    # ==== Random subsample after filtering ====
    if subsample_points and max_points < len(data):
        data = subsample_data(data, max_points, random_seed)

    # ==== Compute data slices ====
    slices, slice_bounds = calculate_slices(data, slice_thickness)


    # ==== Compute centroids for slices ====
    centroids, centroid_ids = compute_centroids(slices, slice_bounds)

    centroids_x, centroids_y, centroids_z = centroids[:, 0], centroids[:, 1], centroids[:, 2]

    # ==== Fit ellipses to slices ====
    ellipses = compute_ellipses(slices, slice_bounds)

    ellipse_centers = []
    for ellipse in ellipses:
        ellipse_centers.append(ellipse["center"])

    ellipse_centers = np.array(ellipse_centers)

    ellipses_x, ellipses_y, ellipses_z = ellipse_centers[:, 0], ellipse_centers[:, 1], ellipse_centers[:, 2]

    # ==== Fit splines ====
    t_arr = np.linspace(0, 1, n_spline_points)

    centroid_spline, centroid_spline_np, centroid_tck = fit_spline_3d(
        centroids_x, centroids_y, centroids_z, t_arr, spline_order=spline_order
    )

    ellipse_spline, ellipse_spline_np, ellipse_tck = fit_spline_3d(
        ellipses_x, ellipses_y, ellipses_z, t_arr, spline_order=spline_order
    )

    # ==== Compute Curvature and Torsion (for Centroid Spline) ====
    D1_centroid, D2_centroid, D3_centroid = compute_derivatives_spline(t_arr, centroid_tck)
    curvature_centroid = compute_curvature(D1_centroid, D2_centroid)
    torsion_centroid = compute_torsion(D1_centroid, D2_centroid, D3_centroid)

    avg_curvature = np.mean(curvature_centroid)
    avg_torsion = np.mean(torsion_centroid)
    print(f"Average Curvature: {avg_curvature:.6f} um^-1")
    print(f"Average Torsion: {avg_torsion:.6f} um^-1")

    # ==== Plotting ====
    x, y, z = data[:, 0], data[:, 1], data[:, 2]

    if plot_3d_geometry_toggle:
        plot_3d_geometry(
            x,
            y,
            z,
            centroid_spline,
            centroids_x,
            centroids_y,
            centroids_z,
            ellipse_spline,
            ellipses_x,
            ellipses_y,
            ellipses_z,
            slice_bounds,
            plot_planes,
            sample_limits
        )

    if plot_outliers_toggle:
        plot_outliers(x, y, z, surface_outliers, deep_outliers, plot_outlier_types)

    if plot_xy_center_approximations:
        plot_xy_centers_path(
            plot_xy_centroids,
            plot_xy_ellipse_centers,
            centroids_x,
            centroids_y,
            centroids_z,
            centroid_ids,
            ellipse_centers,
            sample_limits
        )

    if plot_superposed_slices_xy_toggle:
        plot_superposed_slices_xy(
            slices,
            centroids,
            centroid_ids,
            ellipses,
            sample_limits,
            plot_centroids_in_superposed_slices_xy,
            plot_ellipses_in_all_slices_xy=True,
            plot_ellipse_centers=True
        )

    if plot_ellipses_and_axes_toggle:
        plot_ellipses_and_axes(
            ellipses,
            sample_limits,
            shift_to_common_center=ellipses_shift_to_common_center
        )

    if plot_ellipses_and_axes_toggle:
        plot_ellipses_and_axes(
            ellipses,
            sample_limits,
            shift_to_common_center=not ellipses_shift_to_common_center
        )

    if plot_contour_toggle:
        plot_contour(x, y, z, contour_levels, sample_limits)

    if plot_ellipse_metrics_toggle:
        plot_ellipse_metrics(ellipses)

    if plot_individual_slices_grid_toggle:
        plot_individual_slices_grid(slices, slice_bounds, centroids, ellipses, sample_limits, plot_ellipses_in_slices)

    if plot_slice_views_one_by_one:
        plot_individual_slices_separately(
            slices,
            slice_bounds,
            centroids,
            centroid_ids,
            ellipses,
            sample_limits,
            plot_ellipses_in_slices,
        )
