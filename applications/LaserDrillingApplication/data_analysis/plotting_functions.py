import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata


# ==== Plotting Functions ====
def plot_3d_geometry(x, y, z, spline, cx, cy, cz, ellipse_spline, ecx, ecy, ecz, slice_bounds, plot_planes, sample_limits, filename):
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


def plot_outliers(x, y, z, surface_outliers, deep_outliers, plot_outlier_types, filename):
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
    plot_xy_centroids, plot_xy_ellipse_centers, cx, cy, cz, centroid_ids, ellipse_centers, sample_limits, filename):
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
    plot_ellipse_centers, 
    filename
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
    slice_bounds,
    filename,
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
        z_start, z_end = slice_bounds[i], slice_bounds[i + 1]
        ax.plot(x_ellipse, y_ellipse, "-", color=colors[i], linewidth=2, label=f"Slice {idx} (z ∈ [{z_start:.2f}, {z_end:.2f}] um)")

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


def plot_contour(x, y, z, contour_levels, sample_limits, filename):
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


def plot_ellipse_metrics(ellipses, filename):
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


def plot_individual_slices_grid(slices, slice_bounds, centroids, ellipses, sample_limits, plot_ellipses_in_slices, filename):
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
    n_cols = int(np.ceil(np.sqrt(n_plots)))
    n_rows = int(np.ceil(n_plots / n_cols))

    
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
    plot_ellipses_in_slices,
    filename
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