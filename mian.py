"""
A-BETA: Fully Automated Earthquake Catalog Replenishment Tool
Description:
    This script implements the A-BETA (Auto-Biscale Empirical Transformation Application) method 
    for replenishing incomplete earthquake catalogs, as described in the paper submitted to 
    Seismological Research Letters (under review). It integrates grid-based density analysis 
    and flexible missing region identification to improve catalog completeness.
Authors:
    Chengxiang Zhan, Jiancang Zhuang, Stephen Wu
Dependencies:
    - Required libraries listed in requirements.txt
Citation:
    Zhan, C., Zhuang, J., & Wu, S. (2025). A-BETA: Fully Automated Earthquake Catalog 
    Replenishment Tool. Seismological Research Letters (under review).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from shapely.ops import unary_union
from shapely.geometry import box, LineString, Point, Polygon, MultiPolygon
from scipy.stats import nbinom
from scipy.spatial import cKDTree
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import matplotlib
import sys

matplotlib.use('TkAgg')
plt.rcParams['figure.autolayout'] = False

class TextRedirector:
    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, string):
        self.text_widget.insert(tk.END, string)
        self.text_widget.see(tk.END)
        self.text_widget.update_idletasks()

    def flush(self):
        pass


def load_data(file_path):
    times, magnitudes = [], []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                try:
                    # 分割行并只取前两个值
                    values = line.strip().split()
                    if len(values) < 2:
                        print(f"Skipping invalid line (less than 2 values): {line.strip()}")
                        continue
                    t, m = map(float, values[:2])
                    times.append(t)
                    magnitudes.append(m)
                except ValueError:
                    print(f"Skipping invalid line (non-numeric values): {line.strip()}")
        if not times:
            raise ValueError("File is empty or contains no valid data")
        return np.array(times), np.array(magnitudes)
    except FileNotFoundError:
        raise FileNotFoundError(f"File {file_path} does not exist")
    

def generate_polygon(points, grid_size, min_cells, percentile):
    cell_size = 1.0 / grid_size
    grid_counts = np.zeros((grid_size, grid_size), dtype=int)
    for x, y in points:
        i = min(int(x / cell_size), grid_size - 1)
        j = min(int(y / cell_size), grid_size - 1)
        grid_counts[j, i] += 1
    flattened_counts = grid_counts.flatten()
    # max_points = np.percentile(flattened_counts[flattened_counts > 0], percentile)
    max_points = np.percentile(flattened_counts, percentile)
    print(f"Max points per cell: {max_points}")
    low_density_cells = [(i, j) for i in range(grid_size) for j in range(grid_size) if grid_counts[i, j] <= max_points]
    print(f"Number of low-density cells: {len(low_density_cells)}")

    
    def find_connected_components(cells, grid_size):
        cells_set = set(cells)
        components = []
        visited = set()
        for i, j in cells:
            if (i, j) not in visited:
                component = []
                stack = [(i, j)]
                visited.add((i, j))
                while stack:
                    ci, cj = stack.pop()
                    component.append((ci, cj))
                    for ni, nj in [
                        (ci-1, cj), (ci+1, cj), (ci, cj-1), (ci, cj+1),  # 四连通
                        (ci-1, cj-1), (ci-1, cj+1), (ci+1, cj-1), (ci+1, cj+1)  # 对角方向
                    ]:
                        if (0 <= ni < grid_size and 0 <= nj < grid_size and
                            (ni, nj) in cells_set and (ni, nj) not in visited):
                            stack.append((ni, nj))
                            visited.add((ni, nj))
                components.append(component)
        return components


    connected_components = find_connected_components(low_density_cells, grid_size)

    filtered_components = [comp for comp in connected_components if len(comp) >= min_cells]

    print(f"Initial connected components: {len(connected_components)}, After filtering: {len(filtered_components)}")
    def get_boundary_cells(component, grid_size):
        component_set = set(component)
        boundary_cells = set()
        for i, j in component:
            for ni, nj in [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]:
                if not (0 <= ni < grid_size and 0 <= nj < grid_size and (ni, nj) in component_set):
                    boundary_cells.add((i, j))
                    break
        return boundary_cells
    
    new_components = []
    for component in filtered_components:
        component_set = set(component)
        new_component = component.copy()
        boundary_cells = get_boundary_cells(component, grid_size)
        i_coords = [i for i, _ in component]
        j_coords = [j for _, j in component]
        i_min, i_max = min(i_coords), max(i_coords)
        j_min, j_max = min(j_coords), max(j_coords)
        for i, j in boundary_cells:
            for ni, nj in [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]:
                if (0 <= ni < grid_size and 0 <= nj < grid_size and
                    (ni, nj) not in component_set and
                    i_min <= ni <= i_max and j_min <= nj <= j_max):
                    new_component.append((ni, nj))
                    component_set.add((ni, nj))
        new_components.append(new_component)
    polygons = []
    for component in new_components:
        component_polys = []
        for i, j in component:
            x1 = j * cell_size
            y1 = i * cell_size
            x2 = x1 + cell_size
            y2 = y1 + cell_size
            poly = Polygon([(x1, y1), (x2, y1), (x2, y2), (x1, y2)])
            component_polys.append(poly)
        if component_polys:
            union_poly = unary_union(component_polys)
            union_poly = union_poly.buffer(1e-8, join_style=2, mitre_limit=1.0)
            union_poly = union_poly.simplify(1e-5, preserve_topology=True)
            if isinstance(union_poly, Polygon) and union_poly.is_valid and not union_poly.is_empty:
                polygons.append(union_poly)
            else:
                print(f"Warning: Component polygon invalid, type: {type(union_poly)}, skipping")
    final_poly = unary_union(polygons) if polygons else None
    if final_poly is None or final_poly.is_empty:
        print("Warning: Could not generate a valid polygon, returning empty polygon")
        return Polygon([])
    if isinstance(final_poly, LineString):
        print("Warning: Generated LineString, converting to Polygon")
        final_poly = final_poly.buffer(cell_size * 0.01)
    final_poly = final_poly.intersection(box(0, 0, 1, 1))
    if isinstance(final_poly, LineString):
        print("Warning: Intersection resulted in LineString, converting to Polygon")
        final_poly = final_poly.buffer(cell_size * 0.01)
    if not isinstance(final_poly, (Polygon, MultiPolygon)) or final_poly.is_empty or not final_poly.is_valid:
        print("Warning: Invalid final polygon after intersection, returning empty polygon")
        return Polygon([])
    max_mag_idx = np.argmax(points[:, 1])
    max_mag_point = points[max_mag_idx]
    print(f"Maximum magnitude point: {max_mag_point}")
    minx, miny, maxx, maxy = final_poly.bounds
    right_bottom_point = np.array([maxx, miny])
    print(f"Polygon bottom-right point: {right_bottom_point}")
    x1, y1 = max_mag_point

    y1 -= 0.01     ## set boundary


    x2, y2 = right_bottom_point
    rect = box(min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2))
    print(f"Rectangle bounds: ({min(x1, x2):.4f}, {min(y1, y2):.4f}) to ({max(x1, x2):.4f}, {max(y1, y2):.4f})")
    final_poly = final_poly.intersection(rect)
    if isinstance(final_poly, LineString):
        print("Warning: Rectangle intersection resulted in LineString, converting to Polygon")
        final_poly = final_poly.buffer(cell_size * 0.01)
    if not isinstance(final_poly, (Polygon, MultiPolygon)) or final_poly.is_empty or not final_poly.is_valid or final_poly.area < cell_size**2:
        print("Warning: Rectangle intersection resulted in invalid or too small polygon, returning empty polygon")
        return Polygon([])
    print(f"Final polygon type: {type(final_poly)}, Area: {final_poly.area:.4f}")
    return final_poly

def prepare_data_and_polygon(times, magnitudes, grid_size=20, min_cells=20, percentile=10):
    n = len(times)
    time_order = np.argsort(times)
    mag_order = np.argsort(magnitudes)
    t_transformed = np.zeros(n)
    m_transformed = np.zeros(n)
    t_transformed[time_order] = np.arange(n) / n
    m_transformed[mag_order] = np.arange(n) / n
    points = np.column_stack((t_transformed, m_transformed))
    final_poly = generate_polygon(points, grid_size, min_cells, percentile)
    return points, final_poly

def compute_intersection_lengths(points, poly):
    if not isinstance(poly, (Polygon, MultiPolygon)):
        print(f"Warning: compute_intersection_lengths received invalid geometry type {type(poly)}, returning zero lengths")
        return np.zeros(len(points)), np.zeros(len(points))
    x_lengths = []
    y_lengths = []
    for point in points:
        p = Point(point[0], point[1])
        if poly.contains(p) or poly.boundary.contains(p):
            x_lengths.append(0.0)
            y_lengths.append(0.0)
            continue
        x_line = LineString([(-1, point[1]), (2, point[1])])
        y_line = LineString([(point[0], -1), (point[0], 2)])
        x_length = 0.0
        y_length = 0.0
        # Handle both Polygon and MultiPolygon
        if isinstance(poly, Polygon):
            polygons = [poly]
        else:  # MultiPolygon
            polygons = poly.geoms
        for subpoly in polygons:
            if not subpoly.is_valid or subpoly.is_empty:
                continue
            # Compute x intersection
            x_intersection = subpoly.exterior.intersection(x_line)
            if x_intersection and not x_intersection.is_empty:
                if x_intersection.geom_type == 'MultiPoint':
                    x_coords = sorted([pt.x for pt in x_intersection.geoms])
                    for i in range(0, len(x_coords) - 1, 2):
                        x_length += x_coords[i + 1] - x_coords[i]
                elif x_intersection.geom_type == 'Point':
                    # Single point intersection, no length contribution
                    pass
            # Compute y intersection
            y_intersection = subpoly.exterior.intersection(y_line)
            if y_intersection and not y_intersection.is_empty:
                if y_intersection.geom_type == 'MultiPoint':
                    y_coords = sorted([pt.y for pt in y_intersection.geoms])
                    for i in range(0, len(y_coords) - 1, 2):
                        y_length += y_coords[i + 1] - y_coords[i]
                elif y_intersection.geom_type == 'Point':
                    # Single point intersection, no length contribution
                    pass
        x_lengths.append(x_length)
        y_lengths.append(y_length)
    return np.array(x_lengths), np.array(y_lengths)

def empirical_cdf(values, weights, x):
    indicators = values < x
    weighted_sum = np.sum(weights * indicators)
    total_weight = np.sum(weights)
    if total_weight == 0:
        return 0.0
    return weighted_sum / total_weight

def points_in_polygon(poly, points):
    if not isinstance(poly, (Polygon, MultiPolygon)):
        print(f"Warning: points_in_polygon received invalid geometry type {type(poly)}, returning False for all points")
        return np.zeros(len(points), dtype=bool)
    inside = np.zeros(len(points), dtype=bool)
    for i, point in enumerate(points):
        p = Point(point[0], point[1])
        if isinstance(poly, MultiPolygon):
            inside[i] = any(subpoly.contains(p) or subpoly.boundary.contains(p) for subpoly in poly.geoms)
        elif isinstance(poly, Polygon):
            inside[i] = poly.contains(p) or poly.boundary.contains(p)
    return inside

def estimate_missing_count(transformed_poly, transformed_points):
    if not isinstance(transformed_poly, (Polygon, MultiPolygon)):
        print(f"Warning: Invalid transformed_poly type {type(transformed_poly)}, cannot estimate missing count")
        return 0
    area_s_star = transformed_poly.area
    if area_s_star >= 1 or area_s_star <= 0:
        raise ValueError(f"S* area is invalid: {area_s_star}")
    outside = ~points_in_polygon(transformed_poly, transformed_points)
    n_outside = np.sum(outside)
    print(f"S* area: {area_s_star:.4f}, Number of points outside S*: {n_outside}")
    r = n_outside
    p = 1 - area_s_star
    if r <= 0 or p <= 0 or p >= 1:
        print("Warning: Invalid parameters for negative binomial distribution, returning 0")
        return 0
    n_mis = nbinom.rvs(n=r, p=p)
    print(f"Estimated number of missing events: {n_mis}")
    return n_mis

def generate_missing_points(transformed_poly, n_mis):
    if n_mis == 0 or not isinstance(transformed_poly, (Polygon, MultiPolygon)):
        return np.array([])
    missing_points = []
    minx, miny, maxx, maxy = transformed_poly.bounds
    while len(missing_points) < n_mis:
        u = np.random.uniform(max(minx, 0), min(maxx, 1))
        v = np.random.uniform(max(miny, 0), min(maxy, 1))
        point = Point(u, v)
        if transformed_poly.contains(point):
            missing_points.append([u, v])
    return np.array(missing_points)

def inverse_transform_points(missing_points, original_times, original_magnitudes, final_poly, transformed_points):
    if len(missing_points) == 0:
        return np.array([])
    t_min, t_max = np.min(original_times), np.max(original_times)
    m_min, m_max = np.min(original_magnitudes), np.max(original_magnitudes)
    if t_max == t_min or m_max == m_min:
        print("Warning: Time or magnitude range is invalid, returning empty point set")
        return np.array([])
    points = transformed_points
    inside = points_in_polygon(final_poly, points)
    inside_points = points[inside]
    n_inside = len(inside_points)
    print(f"Number of original points inside polygon: {n_inside}")
    if n_inside > 0 and len(missing_points) > 0:
        print(f"Number of replenished points: {len(missing_points)}, Number of original points inside polygon: {n_inside}")
        points_to_remove = []
        available_indices = np.arange(len(missing_points))
        available_points = missing_points.copy()
        if len(missing_points) < n_inside:
            print(f"Warning: Not enough replenished points, removing only {len(missing_points)}")
            missing_points = np.array([])
        else:
            tree = cKDTree(available_points)
            for _ in range(n_inside):
                if len(available_points) == 0:
                    break
                distances, idx = tree.query(inside_points[len(points_to_remove)], k=1)
                points_to_remove.append(available_indices[idx])
                available_indices = np.delete(available_indices, idx)
                available_points = np.delete(available_points, idx, axis=0)
                tree = cKDTree(available_points) if len(available_points) > 0 else None
            keep_mask = np.ones(len(missing_points), dtype=bool)
            keep_mask[points_to_remove] = False
            missing_points = missing_points[keep_mask]
        print(f"Removed {len(points_to_remove)} replenished points closest to observed points inside polygon, {len(missing_points)} remain")
    if len(missing_points) == 0:
        return np.array([])
    transformed = np.zeros_like(missing_points)
    sorted_t_indices = np.argsort(points[:, 0])
    sorted_t = points[sorted_t_indices, 0]
    sorted_times = original_times[sorted_t_indices]
    sorted_m_indices = np.argsort(points[:, 1])
    sorted_m = points[sorted_m_indices, 1]
    sorted_magnitudes = original_magnitudes[sorted_m_indices]
    for i, (t_star, m_star) in enumerate(missing_points):
        idx = np.searchsorted(sorted_t, t_star, side='right')
        if idx == 0:
            transformed[i, 0] = sorted_times[0]
        elif idx == len(sorted_t):
            transformed[i, 0] = sorted_times[-1]
        else:
            t_left, t_right = sorted_t[idx-1], sorted_t[idx]
            time_left, time_right = sorted_times[idx-1], sorted_times[idx]
            if t_right == t_left:
                transformed[i, 0] = time_left
            else:
                weight = (t_star - t_left) / (t_right - t_left)
                transformed[i, 0] = time_left + weight * (time_right - time_left)
        idx = np.searchsorted(sorted_m, m_star, side='right')
        if idx == 0:
            transformed[i, 1] = sorted_magnitudes[0]
        elif idx == len(sorted_m):
            transformed[i, 1] = sorted_magnitudes[-1]
        else:
            m_left, m_right = sorted_m[idx-1], sorted_m[idx]
            mag_left, mag_right = sorted_magnitudes[idx-1], sorted_magnitudes[idx]
            if m_right == m_left:
                transformed[i, 1] = mag_left
            else:
                weight = (m_star - m_left) / (m_right - m_left)
                transformed[i, 1] = mag_left + weight * (mag_right - mag_left)
    print(f"Replenished points time range: [{transformed[:,0].min():.2f}, {transformed[:,0].max():.2f}]")
    print(f"Replenished points magnitude range: [{transformed[:,1].min():.2f}, {transformed[:,1].max():.2f}]")
    return transformed

def merge_points(original_times, original_magnitudes, missing_points):
    if len(missing_points) == 0:
        return np.column_stack((original_times, original_magnitudes))
    original_points = np.column_stack((original_times, original_magnitudes))
    return np.vstack([original_points, missing_points])

def plot_final_result(ax4, ax5, ax6, original_times, original_magnitudes, missing_points, final_poly, points, transformed_points, missing_points_star, transformed_poly):
    ax4.clear()
    ax4.scatter(transformed_points[:, 0], transformed_points[:, 1], s=1, c='blue', alpha=0.5)
    if len(missing_points_star) > 0:
        ax4.scatter(missing_points_star[:, 0], missing_points_star[:, 1], s=1, c='green', alpha=0.5)
    if isinstance(transformed_poly, MultiPolygon):
        for subpoly in transformed_poly.geoms:
            x, y = subpoly.exterior.xy
            ax4.plot(x, y, 'r-', linewidth=2)
    elif isinstance(transformed_poly, Polygon):
        x, y = transformed_poly.exterior.xy
        ax4.plot(x, y, 'r-', linewidth=2)
    else:
        print(f"Warning: transformed_poly type {type(transformed_poly)} not plotted")
    ax4.set_xlim(-0.05, 1.05)
    ax4.set_ylim(-0.05, 1.05)
    ax4.set_box_aspect(1)
    ax4.set_title('BEPIT with Missing', fontsize=10)
    ax4.set_xlabel('Empirical Time', fontsize=8)
    ax4.set_ylabel('Empirical Magnitude', fontsize=8)
    ax4.tick_params(labelsize=6)
    ax4.grid(True)

    ax5.clear()
    ax5.scatter(original_times, original_magnitudes, s=3, c='blue', alpha=0.5)
    if len(missing_points) > 0:
        ax5.scatter(missing_points[:, 0], missing_points[:, 1], s=3, c='green', alpha=0.5)
    ax5.set_title('Replenished Points', fontsize=10)
    ax5.set_xlabel('Time', fontsize=8)
    ax5.set_ylabel('Magnitude', fontsize=8)
    ax5.tick_params(labelsize=6)
    ax5.grid(True)
    ax5.set_box_aspect(1)

    ax6.clear()
    original_sorted_times = np.sort(original_times)
    original_counts = np.arange(1, len(original_times) + 1)
    all_times = np.concatenate([original_times, missing_points[:, 0] if len(missing_points) > 0 else []])
    all_sorted_times = np.sort(all_times)
    all_counts = np.arange(1, len(all_times) + 1)
    ax6.step(original_sorted_times, original_counts, 'b-', label='Original', where='post')
    ax6.step(all_sorted_times, all_counts, 'r-', label='Replenished', where='post')
    ax6.set_box_aspect(1)
    ax6.set_title('Cumulative Events', fontsize=10)
    ax6.set_xlabel('Time', fontsize=8)
    ax6.set_ylabel('Event Count', fontsize=8)
    ax6.tick_params(labelsize=6)
    ax6.grid(True)
    ax6.legend(fontsize=6)

def combine_polygons(grid_poly, intersection_polygons, union_polygons):
    if not isinstance(grid_poly, (Polygon, MultiPolygon)):
        print("Warning: Grid polygon is invalid, cannot combine")
        return None
    final_poly = grid_poly
    # Apply intersection with intersection-type polygons
    if intersection_polygons:
        try:
            for poly in intersection_polygons:
                if isinstance(poly, (Polygon, MultiPolygon)) and poly.is_valid and not poly.is_empty:
                    final_poly = final_poly.intersection(poly)
                    if isinstance(final_poly, LineString):
                        print("Warning: Intersection resulted in LineString, converting to Polygon")
                        final_poly = final_poly.buffer(1e-5)
                    if not isinstance(final_poly, (Polygon, MultiPolygon)) or final_poly.is_empty or not final_poly.is_valid:
                        print("Warning: Intersection resulted in invalid or empty polygon")
                        return None
        except Exception as e:
            print(f"Warning: Failed to compute intersection with user polygons: {e}")
            return None
    # Apply union with union-type polygons
    if union_polygons:
        try:
            union_polys = [final_poly] + [poly for poly in union_polygons if isinstance(poly, (Polygon, MultiPolygon)) and poly.is_valid and not poly.is_empty]
            if union_polys:
                final_poly = unary_union(union_polys)
                if isinstance(final_poly, LineString):
                    print("Warning: Union resulted in LineString, converting to Polygon")
                    final_poly = final_poly.buffer(1e-5)
                if not isinstance(final_poly, (Polygon, MultiPolygon)) or final_poly.is_empty or not final_poly.is_valid:
                    print("Warning: Union resulted in invalid or empty polygon")
                    return None
        except Exception as e:
            print(f"Warning: Failed to compute union with user polygons: {e}")
            return None
    return final_poly

def iterate_bepit(points, final_poly, grid_size, min_cells, percentile, epsilon=1e-6, max_iter=100, ax=None, canvas=None, option="direct", intersection_polygons=None, union_polygons=None, app=None):
    if not isinstance(final_poly, (Polygon, MultiPolygon)):
        raise ValueError(f"final_poly must be Polygon or MultiPolygon, got {type(final_poly)}")
    points = np.array(points)
    if points.shape[1] != 2:
        raise ValueError("points must be an array of shape [[x_1, y_1], [x_2, y_2], ...]")
    n = len(points)
    current_points = points.copy()
    current_grid_poly = final_poly
    intersection_vertices = [np.array(poly.exterior.coords)[:-1] for poly in intersection_polygons] if intersection_polygons else []
    union_vertices = [np.array(poly.exterior.coords)[:-1] for poly in union_polygons] if union_polygons else []

    def plot_iteration(points, poly, iteration):
        ax.clear()
        ax.scatter(points[:, 0], points[:, 1], s=1, c='blue', alpha=0.5)
        if isinstance(poly, MultiPolygon):
            for subpoly in poly.geoms:
                x, y = subpoly.exterior.xy
                ax.plot(x, y, 'r-', linewidth=2)
        elif isinstance(poly, Polygon):
            x, y = poly.exterior.xy
            ax.plot(x, y, 'r-', linewidth=2)
        else:
            print(f"Warning: Iteration {iteration} skipping plot for geometry type {type(poly)}")
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.set_box_aspect(1)
        ax.set_title(f'Iteration {iteration}', fontsize=10)
        ax.set_xlabel('Empirical Time', fontsize=8)
        ax.set_ylabel('Empirical Magnitude', fontsize=8)
        ax.tick_params(labelsize=6)
        ax.grid(True)
        canvas.draw()
        canvas.flush_events()

    # Initial combination of polygons
    current_poly = current_grid_poly
    if option in ["user", "combine"]:
        current_poly = combine_polygons(current_grid_poly, intersection_polygons, union_polygons)
        if current_poly is None or not isinstance(current_poly, (Polygon, MultiPolygon)):
            print("Warning: Initial polygon combination invalid, using grid polygon")
            current_poly = current_grid_poly
    plot_iteration(current_points, current_poly, 0)
    x_std, y_std = np.std(current_points[:, 0]), np.std(current_points[:, 1])
    print(f"Initial uniformity: x_std={x_std:.4f}, y_std={y_std:.4f}")
    last_w1, last_w2 = None, None

    for iteration in range(max_iter):
        if app is not None and app.stop_iteration:
            print(f"Iteration stopped by user at iteration {iteration + 1}")
            break
        print(f"\nStarting iteration {iteration + 1}-{option}")
        prev_points = current_points.copy()
        print(f"Current combined polygon type: {type(current_poly)}, Area: {current_poly.area:.4f}")
        
        x_lengths, y_lengths = compute_intersection_lengths(current_points, current_poly)
        x_lengths = np.clip(x_lengths, 0, 0.999)
        y_lengths = np.clip(y_lengths, 0, 0.999)
        inside = points_in_polygon(current_poly, current_points)
        w_1 = np.where(inside, 0.0, 1.0 / (1.0 - y_lengths))
        w_2 = np.where(inside, 0.0, 1.0 / (1.0 - x_lengths))
        w_1 = np.nan_to_num(w_1, nan=0.0, posinf=0.0, neginf=0.0)
        w_2 = np.nan_to_num(w_2, nan=0.0, posinf=0.0, neginf=0.0)
        last_w1, last_w2 = w_1.copy(), w_2.copy()
        print(f"x_lengths range: [{x_lengths.min():.4f}, {x_lengths.max():.4f}]")
        print(f"y_lengths range: [{y_lengths.min():.4f}, {y_lengths.max():.4f}]")
        print(f"w_1 range: [{w_1.min():.4f}, {w_1.max():.4f}]")
        print(f"w_2 range: [{w_2.min():.4f}, {w_2.max():.4f}]")

        new_points = np.zeros_like(current_points)
        for i in range(n):
            t_i = current_points[i, 0]
            m_i = current_points[i, 1]
            new_points[i, 0] = empirical_cdf(current_points[:, 0], w_1, t_i)
            new_points[i, 1] = empirical_cdf(current_points[:, 1], w_2, m_i)

        if option == "user":
            # Transform user polygons only
            new_intersection_vertices = []
            new_intersection_polygons = []
            for vertices in intersection_vertices:
                new_vertices = np.zeros_like(vertices)
                for i, (x, y) in enumerate(vertices):
                    new_x = empirical_cdf(current_points[:, 0], w_1, x)
                    new_y = empirical_cdf(current_points[:, 1], w_2, y)
                    new_vertices[i, 0] = np.clip(new_x, 0, 1)
                    new_vertices[i, 1] = np.clip(new_y, 0, 1)
                try:
                    poly = Polygon(new_vertices)
                    if poly.is_valid and not poly.is_empty:
                        new_intersection_polygons.append(poly)
                        new_intersection_vertices.append(new_vertices)
                except Exception as e:
                    print(f"Warning: Iteration {iteration + 1} intersection polygon creation failed: {e}")
            new_union_vertices = []
            new_union_polygons = []
            for vertices in union_vertices:
                new_vertices = np.zeros_like(vertices)
                for i, (x, y) in enumerate(vertices):
                    new_x = empirical_cdf(current_points[:, 0], w_1, x)
                    new_y = empirical_cdf(current_points[:, 1], w_2, y)
                    new_vertices[i, 0] = np.clip(new_x, 0, 1)
                    new_vertices[i, 1] = np.clip(new_y, 0, 1)
                try:
                    poly = Polygon(new_vertices)
                    if poly.is_valid and not poly.is_empty:
                        new_union_polygons.append(poly)
                        new_union_vertices.append(new_vertices)
                except Exception as e:
                    print(f"Warning: Iteration {iteration + 1} union polygon creation failed: {e}")
            # Combine user polygons
            current_poly = combine_polygons(Polygon([(0,0), (1,0), (1,1), (0,1)]), new_intersection_polygons, new_union_polygons)
            if current_poly is None or not isinstance(current_poly, (Polygon, MultiPolygon)):
                print(f"Warning: Iteration {iteration + 1} user polygon combination invalid, retaining previous polygon")
                current_poly = current_poly
            else:
                intersection_vertices = new_intersection_vertices
                union_vertices = new_union_vertices
        elif option == "combine":
            # Generate new grid polygon
            new_grid_poly = generate_polygon(new_points, grid_size, min_cells, percentile)
            if not isinstance(new_grid_poly, (Polygon, MultiPolygon)) or new_grid_poly.is_empty:
                print(f"Warning: Iteration {iteration + 1} generated invalid grid polygon, retaining previous polygon")
                current_poly = current_poly
            else:
                current_grid_poly = new_grid_poly
                # Transform user polygons
                new_intersection_vertices = []
                new_intersection_polygons = []
                for vertices in intersection_vertices:
                    new_vertices = np.zeros_like(vertices)
                    for i, (x, y) in enumerate(vertices):
                        new_x = empirical_cdf(current_points[:, 0], w_1, x)
                        new_y = empirical_cdf(current_points[:, 1], w_2, y)
                        new_vertices[i, 0] = np.clip(new_x, 0, 1)
                        new_vertices[i, 1] = np.clip(new_y, 0, 1)
                    try:
                        poly = Polygon(new_vertices)
                        if poly.is_valid and not poly.is_empty:
                            new_intersection_polygons.append(poly)
                            new_intersection_vertices.append(new_vertices)
                    except Exception as e:
                        print(f"Warning: Iteration {iteration + 1} intersection polygon creation failed: {e}")
                new_union_vertices = []
                new_union_polygons = []
                for vertices in union_vertices:
                    new_vertices = np.zeros_like(vertices)
                    for i, (x, y) in enumerate(vertices):
                        new_x = empirical_cdf(current_points[:, 0], w_1, x)
                        new_y = empirical_cdf(current_points[:, 1], w_2, y)
                        new_vertices[i, 0] = np.clip(new_x, 0, 1)
                        new_vertices[i, 1] = np.clip(new_y, 0, 1)
                    try:
                        poly = Polygon(new_vertices)
                        if poly.is_valid and not poly.is_empty:
                            new_union_polygons.append(poly)
                            new_union_vertices.append(new_vertices)
                    except Exception as e:
                        print(f"Warning: Iteration {iteration + 1} union polygon creation failed: {e}")
                # Combine grid and user polygons
                current_poly = combine_polygons(current_grid_poly, new_intersection_polygons, new_union_polygons)
                if current_poly is None or not isinstance(current_poly, (Polygon, MultiPolygon)):
                    print(f"Warning: Iteration {iteration + 1} polygon combination invalid, using grid polygon")
                    current_poly = current_grid_poly
                else:
                    intersection_vertices = new_intersection_vertices
                    union_vertices = new_union_vertices
        else:
            # Direct iteration with grid-generated polygon
            current_grid_poly = generate_polygon(new_points, grid_size, min_cells, percentile)
            if not isinstance(current_grid_poly, (Polygon, MultiPolygon)):
                print(f"Warning: Iteration {iteration + 1} generated invalid polygon, retaining previous polygon")
                current_poly = current_poly
            else:
                current_poly = current_grid_poly
        
        current_points = new_points
        plot_iteration(current_points, current_poly, iteration + 1)
        x_std, y_std = np.std(current_points[:, 0]), np.std(current_points[:, 1])
        print(f"Uniformity: x_std={x_std:.4f}, y_std={y_std:.4f}")
        max_diff = np.max(np.abs(current_points - prev_points))
        print(f"Maximum coordinate change: {max_diff:.6f}")
        if max_diff < epsilon:
            print(f"Converged after {iteration + 1} iterations")
            break
        elif iteration == max_iter - 1:
            print("Reached maximum iterations without full convergence")
    
    return current_points, current_poly, last_w1, last_w2

class BEPITApp:
    def __init__(self, root):
        self.root = root
        self.root.title("A-BETA")
        self.root.geometry("1200x1000")

        self.reset_state()

        # Control frame (Buttons, Parameters, Options, Output)
        self.control_frame = ttk.Frame(self.root, padding=5)
        self.control_frame.pack(fill=tk.X, padx=5, pady=10)

        # Button frame (6 rows, 1 column, in a LabelFrame)
        self.button_frame = ttk.LabelFrame(self.control_frame, text="Controls", padding=5)
        self.button_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5)
        self.load_data_button = tk.Button(self.button_frame, text="Load Data", command=self.load_data_only, width=12)
        self.load_data_button.grid(row=0, column=0, pady=2)
        self.generate_bepit_button = tk.Button(self.button_frame, text="Transformation", command=self.generate_initial_bepit, width=12, state=tk.DISABLED)
        self.generate_bepit_button.grid(row=1, column=0, pady=2)
        self.select_button = tk.Button(self.button_frame, text="Select Polygon", command=self.start_polygon_selection, width=12, state=tk.DISABLED)
        self.select_button.grid(row=2, column=0, pady=2)
        self.run_button = tk.Button(self.button_frame, text="Run Iteration", command=self.run_iteration, width=12, state=tk.DISABLED)
        self.run_button.grid(row=3, column=0, pady=2)
        self.save_button = tk.Button(self.button_frame, text="Save Data", command=self.save_replenished_data, width=12, state=tk.DISABLED)
        self.save_button.grid(row=4, column=0, pady=2)
        self.reset_button = tk.Button(self.button_frame, text="Reset", command=self.reset_app, width=12)
        self.reset_button.grid(row=5, column=0, pady=2)

        # Left frame (Parameters, Options)
        self.left_frame = ttk.Frame(self.control_frame)
        self.left_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5)

        # Parameter frame
        self.param_frame = ttk.LabelFrame(self.left_frame, text="Parameters", padding=5)
        self.param_frame.pack(fill=tk.X, padx=5, pady=2)
        defaults = {'grid_size': '70', 'percentile': '10', 'min_cells': '200', 'epsilon': '1e-6', 'max_iter': '100'}
        self.param_entries = {}
        for i, (label, var, default) in enumerate([
            ("Grid Size", "grid_size", defaults['grid_size']),
            ("Percentile", "percentile", defaults['percentile']),
            ("Minimum Cells", "min_cells", defaults['min_cells']),
            ("Convergence Threshold", "epsilon", defaults['epsilon']),
            ("Maximum Iterations", "max_iter", defaults['max_iter'])
        ]):
            frame = ttk.Frame(self.param_frame)
            frame.pack(fill=tk.X, pady=1)
            ttk.Label(frame, text=label + ":").pack(side=tk.LEFT)
            entry = ttk.Entry(frame, width=10)
            entry.insert(0, default)
            entry.pack(side=tk.LEFT, padx=5)
            self.param_entries[var] = entry

        # Option frame
        self.option_frame = ttk.LabelFrame(self.left_frame, text="Iteration Option", padding=5)
        self.option_frame.pack(fill=tk.X, padx=5, pady=2)
        self.option_var = tk.StringVar(value="direct")
        self.radio_direct = ttk.Radiobutton(self.option_frame, text="Direct Iteration", value="direct", variable=self.option_var)
        self.radio_direct.pack(anchor=tk.W)
        self.radio_user = ttk.Radiobutton(self.option_frame, text="Use Selected Polygons", value="user", variable=self.option_var, state=tk.DISABLED)
        self.radio_user.pack(anchor=tk.W)
        self.radio_combine = ttk.Radiobutton(self.option_frame, text="Combine with Selected Polygons", value="combine", variable=self.option_var, state=tk.DISABLED)
        self.radio_combine.pack(anchor=tk.W)
        

        # Right frame (Output)
        self.right_frame = ttk.Frame(self.control_frame)
        self.right_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(10, 5))
        self.output_frame = ttk.LabelFrame(self.right_frame, text="Output", padding=5)
        self.output_frame.pack(fill=tk.Y)
        self.output_text = tk.Text(self.output_frame, height=15, width=100, wrap=tk.WORD)
        self.output_text.pack(fill=tk.BOTH, padx=5, pady=5)
        sys.stdout = TextRedirector(self.output_text)

        # Plot frame
        self.fig = plt.Figure(figsize=(12, 6))
        self.ax1 = self.fig.add_subplot(231)
        self.ax2 = self.fig.add_subplot(232)
        self.ax3 = self.fig.add_subplot(233)
        self.ax4 = self.fig.add_subplot(234)
        self.ax5 = self.fig.add_subplot(235)
        self.ax6 = self.fig.add_subplot(236)
        for ax in [self.ax1, self.ax2, self.ax3, self.ax4, self.ax5, self.ax6]:
            ax.set_visible(False)
            ax.set_box_aspect(1)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        self.canvas_has_focus = False

    def reset_state(self):
        self.file_path = None
        self.times = None
        self.magnitudes = None
        self.points = None
        self.final_poly = None
        self.intersection_polygons = []
        self.union_polygons = []
        self.selected_points = []
        self.polygon_line = None
        self.cid = None
        self.is_selecting_polygon = False
        self.is_running = False
        self.stop_iteration = False
        self.transformed_points = None
        self.transformed_poly = None
        self.missing_points = None
        self.missing_points_star = None

    def reset_app(self):
        print("Reset button clicked")
        self.disconnect_polygon_handlers()
        self.reset_state()
        for ax in [self.ax1, self.ax2, self.ax3, self.ax4, self.ax5, self.ax6]:
            ax.clear()
            ax.set_visible(False)
        self.canvas.draw()
        self.load_data_button.config(state=tk.NORMAL)
        self.generate_bepit_button.config(state=tk.DISABLED)
        self.select_button.config(state=tk.DISABLED)
        self.run_button.config(state=tk.DISABLED, text="Run Iteration")
        self.save_button.config(state=tk.DISABLED)
        self.radio_direct.config(state=tk.NORMAL)
        self.radio_user.config(state=tk.DISABLED)
        self.radio_combine.config(state=tk.DISABLED)
        self.output_text.delete(1.0, tk.END)
        for var, entry in self.param_entries.items():
            entry.delete(0, tk.END)
            entry.insert(0, {'grid_size': '70', 'percentile': '10', 'min_cells': '200', 'epsilon': '1e-6', 'max_iter': '100'}[var])
        self.option_var.set("direct")

    def focus_canvas(self, event):
        if not self.canvas_has_focus:
            self.canvas.get_tk_widget().focus_set()
            self.canvas_has_focus = True
            print("Canvas focused")

    def load_data_only(self):
        print("Load Data button clicked")
        self.file_path = filedialog.askopenfilename(filetypes=[("Text files", "*.txt")])
        if not self.file_path:
            return
        try:
            self.times, self.magnitudes = load_data(self.file_path)
            self.points = None
            self.final_poly = None
            self.intersection_polygons = []
            self.union_polygons = []
            self.selected_points = []
            self.transformed_points = None
            self.transformed_poly = None
            self.missing_points = None
            self.missing_points_star = None

            for ax in [self.ax1, self.ax2, self.ax3, self.ax4, self.ax5, self.ax6]:
                ax.clear()
                ax.set_visible(False)
            self.ax1.set_visible(True)
            self.ax1.scatter(self.times, self.magnitudes, s=1, c='blue', alpha=0.5)
            self.ax1.set_title('Original Points', fontsize=10)
            self.ax1.set_xlabel('Time', fontsize=8)
            self.ax1.set_ylabel('Magnitude', fontsize=8)
            self.ax1.tick_params(labelsize=6)
            self.ax1.grid(True)
            self.ax1.set_box_aspect(1)
            self.canvas.draw()
            self.generate_bepit_button.config(state=tk.NORMAL)
            self.select_button.config(state=tk.DISABLED)
            self.run_button.config(state=tk.DISABLED, text="Run Iteration")
            self.save_button.config(state=tk.DISABLED)
            self.radio_direct.config(state=tk.NORMAL)
            self.radio_user.config(state=tk.DISABLED)
            self.radio_combine.config(state=tk.DISABLED)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file: {e}")

    def generate_initial_bepit(self):
        print("Generate BEPIT button clicked")
        try:
            grid_size = int(self.param_entries['grid_size'].get())
            percentile = float(self.param_entries['percentile'].get())
            min_cells = int(self.param_entries['min_cells'].get())
            if grid_size <= 0 or min_cells <= 0 or not (0 <= percentile <= 100):
                raise ValueError("Grid Size and Minimum Cells must be positive integers, Percentile must be between 0 and 100")
            self.points, self.final_poly = prepare_data_and_polygon(self.times, self.magnitudes, grid_size, min_cells, percentile)
            self.intersection_polygons = []
            self.union_polygons = []
            self.selected_points = []
            self.ax2.clear()
            self.ax2.set_visible(True)
            self.ax2.scatter(self.points[:, 0], self.points[:, 1], s=1, c='blue', alpha=0.5)
            if isinstance(self.final_poly, MultiPolygon):
                for subpoly in self.final_poly.geoms:
                    x, y = subpoly.exterior.xy
                    self.ax2.plot(x, y, 'r-', linewidth=2)
            elif isinstance(self.final_poly, Polygon):
                x, y = self.final_poly.exterior.xy
                self.ax2.plot(x, y, 'r-', linewidth=2)
            self.ax2.set_xlim(-0.05, 1.05)
            self.ax2.set_ylim(-0.05, 1.05)
            self.ax2.set_box_aspect(1)
            self.ax2.set_title('Initial BEPIT', fontsize=10)
            self.ax2.set_xlabel('Empirical Time', fontsize=8)
            self.ax2.set_ylabel('Empirical Magnitude', fontsize=8)
            self.ax2.tick_params(labelsize=6)
            self.ax2.grid(True)
            self.canvas.draw()
            self.select_button.config(state=tk.NORMAL)
            self.run_button.config(state=tk.NORMAL)
        except ValueError as e:
            messagebox.showerror("Error", f"Invalid parameter values: {e}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate initial BEPIT: {e}")

    def start_polygon_selection(self):
        print("Select Polygon button clicked")
        try:
            if self.points is None or self.final_poly is None:
                raise ValueError("No initial BEPIT data available")
            self.is_selecting_polygon = True
            self.selected_points = []

            # Prepare ax2 for polygon drawing
            self.ax2.clear()
            self.ax2.scatter(self.points[:, 0], self.points[:, 1], s=1, c='blue', alpha=0.5)
            if isinstance(self.final_poly, MultiPolygon):
                for subpoly in self.final_poly.geoms:
                    x, y = subpoly.exterior.xy
                    self.ax2.plot(x, y, 'r-', linewidth=2)
            elif isinstance(self.final_poly, Polygon):
                x, y = self.final_poly.exterior.xy
                self.ax2.plot(x, y, 'r-', linewidth=2)
            # Plot existing user polygons
            for poly in self.intersection_polygons:
                x, y = poly.exterior.xy
                self.ax2.plot(x, y, 'g-', linewidth=2)  # Green for intersection
            for poly in self.union_polygons:
                x, y = poly.exterior.xy
                self.ax2.plot(x, y, 'm-', linewidth=2)  # Magenta for union
            self.ax2.set_xlim(-0.05, 1.05)
            self.ax2.set_ylim(-0.05, 1.05)
            self.ax2.set_box_aspect(1)
            self.ax2.set_title('Draw Polygon (Press c to close)', fontsize=10)
            self.ax2.set_xlabel('Empirical Time', fontsize=8)
            self.ax2.set_ylabel('Empirical Magnitude', fontsize=8)
            self.ax2.tick_params(labelsize=6)
            self.ax2.grid(True)

            # Initialize line for drawing
            self.polygon_line, = self.ax2.plot([], [], 'k-', linewidth=2)

            # Connect event handlers
            self.cid = self.canvas.mpl_connect('button_press_event', self.on_polygon_click)
            self.root.bind('<Key-c>', self.close_polygon)

            self.select_button.config(state=tk.DISABLED)
            self.radio_user.config(state=tk.DISABLED)
            self.radio_combine.config(state=tk.DISABLED)
            self.load_data_button.config(state=tk.DISABLED)
            self.generate_bepit_button.config(state=tk.DISABLED)
            self.run_button.config(state=tk.DISABLED)
            self.save_button.config(state=tk.DISABLED)
            messagebox.showinfo("Instruction", "Click on ax2 to form a polygon. Press 'c' to close the polygon.")
            self.canvas.draw()
        except Exception as e:
            print(f"Error starting polygon selection: {e}")
            messagebox.showerror("Error", f"Failed to start polygon selection: {e}")
            self.cancel_polygon_selection()

    def on_polygon_click(self, event):
        if not self.is_selecting_polygon or event.inaxes != self.ax2:
            print("Click ignored: not selecting or outside ax2")
            return
        x, y = event.xdata, event.ydata
        if x is None or y is None:
            print("Invalid click coordinates, ignored")
            return
        x = max(0, min(1, x))
        y = max(0, min(1, y))
        print(f"Click registered at ({x:.3f}, {y:.3f})")
        self.selected_points.append((x, y))
        if len(self.selected_points) > 1:
            x_data = [p[0] for p in self.selected_points]
            y_data = [p[1] for p in self.selected_points]
            self.polygon_line.set_data(x_data, y_data)
            self.canvas.draw()

    def close_polygon(self, event):
        print("Close polygon triggered")
        if len(self.selected_points) < 3:
            messagebox.showwarning("Warning", "At least 3 points are required to form a polygon.")
            return
        self.selected_points.append(self.selected_points[0])
        x_data = [p[0] for p in self.selected_points]
        y_data = [p[1] for p in self.selected_points]
        self.polygon_line.set_data(x_data, y_data)
        self.canvas.draw()
        try:
            new_polygon = Polygon(self.selected_points)
            if not new_polygon.is_valid:
                messagebox.showerror("Error", "Invalid polygon created.")
                self.cancel_polygon_selection()
                return
        except Exception as e:
            messagebox.showerror("Error", f"Failed to create polygon: {e}")
            self.cancel_polygon_selection()
            return

        # Ask user for intersection or union
        dialog = tk.Toplevel(self.root)
        dialog.title("Polygon Operation")
        dialog.geometry("300x150")
        dialog.transient(self.root)
        dialog.grab_set()
        ttk.Label(dialog, text="Choose how to use this polygon:").pack(pady=10)
        operation_var = tk.StringVar(value="intersection")
        ttk.Radiobutton(dialog, text="Intersection", value="intersection", variable=operation_var).pack()
        ttk.Radiobutton(dialog, text="Union", value="union", variable=operation_var).pack()
        def confirm():
            operation = operation_var.get()
            if operation == "intersection":
                self.intersection_polygons.append(new_polygon)
                print("Polygon added to intersection polygons")
            else:
                self.union_polygons.append(new_polygon)
                print("Polygon added to union polygons")
            dialog.destroy()
            self.complete_polygon_selection()
        def cancel():
            dialog.destroy()
            self.cancel_polygon_selection()
        ttk.Button(dialog, text="Confirm", command=confirm).pack(pady=5)
        ttk.Button(dialog, text="Cancel", command=cancel).pack(pady=5)
        self.root.wait_window(dialog)

    def complete_polygon_selection(self):
        # Update ax2 with all polygons
        self.ax2.clear()
        self.ax2.scatter(self.points[:, 0], self.points[:, 1], s=1, c='blue', alpha=0.5)
        if isinstance(self.final_poly, MultiPolygon):
            for subpoly in self.final_poly.geoms:
                x, y = subpoly.exterior.xy
                self.ax2.plot(x, y, 'r-', linewidth=2)
        elif isinstance(self.final_poly, Polygon):
            x, y = self.final_poly.exterior.xy
            self.ax2.plot(x, y, 'r-', linewidth=2)
        for poly in self.intersection_polygons:
            x, y = poly.exterior.xy
            self.ax2.plot(x, y, 'g-', linewidth=2)  # Green for intersection
        for poly in self.union_polygons:
            x, y = poly.exterior.xy
            self.ax2.plot(x, y, 'm-', linewidth=2)  # Magenta for union
        self.ax2.set_xlim(-0.05, 1.05)
        self.ax2.set_ylim(-0.05, 1.05)
        self.ax2.set_box_aspect(1)
        self.ax2.set_title('Initial BEPIT with User Polygons', fontsize=10)
        self.ax2.set_xlabel('Empirical Time', fontsize=8)
        self.ax2.set_ylabel('Empirical Magnitude', fontsize=8)
        self.ax2.tick_params(labelsize=6)
        self.ax2.grid(True)
        self.canvas.draw()

        # Clean up
        self.disconnect_polygon_handlers()
        self.is_selecting_polygon = False
        self.select_button.config(state=tk.NORMAL)
        self.radio_user.config(state=tk.NORMAL)
        self.radio_combine.config(state=tk.NORMAL)
        self.load_data_button.config(state=tk.NORMAL)
        self.generate_bepit_button.config(state=tk.NORMAL)
        self.run_button.config(state=tk.NORMAL)
        self.save_button.config(state=tk.NORMAL if self.missing_points is not None else tk.DISABLED)
        messagebox.showinfo("Success", "Polygon created successfully. You can select another polygon or proceed.")

    def disconnect_polygon_handlers(self):
        if self.cid:
            self.canvas.mpl_disconnect(self.cid)
            self.cid = None
        self.root.unbind('<Key-c>')
        if self.polygon_line:
            self.polygon_line.remove()
            self.polygon_line = None
        self.canvas.draw()

    def cancel_polygon_selection(self):
        if not self.is_selecting_polygon:
            return
        print("Polygon selection cancelled")
        self.disconnect_polygon_handlers()
        self.selected_points = []
        self.is_selecting_polygon = False
        # Redraw ax2 to restore original state
        self.ax2.clear()
        self.ax2.set_visible(True)
        self.ax2.scatter(self.points[:, 0], self.points[:, 1], s=1, c='blue', alpha=0.5)
        if isinstance(self.final_poly, MultiPolygon):
            for subpoly in self.final_poly.geoms:
                x, y = subpoly.exterior.xy
                self.ax2.plot(x, y, 'r-', linewidth=2)
        elif isinstance(self.final_poly, Polygon):
            x, y = self.final_poly.exterior.xy
            self.ax2.plot(x, y, 'r-', linewidth=2)
        for poly in self.intersection_polygons:
            x, y = poly.exterior.xy
            self.ax2.plot(x, y, 'g-', linewidth=2)
        for poly in self.union_polygons:
            x, y = poly.exterior.xy
            self.ax2.plot(x, y, 'm-', linewidth=2)
        self.ax2.set_xlim(-0.05, 1.05)
        self.ax2.set_ylim(-0.05, 1.05)
        self.ax2.set_box_aspect(1)
        self.ax2.set_title('Initial BEPIT', fontsize=10)
        self.ax2.set_xlabel('Empirical Time', fontsize=8)
        self.ax2.set_ylabel('Empirical Magnitude', fontsize=8)
        self.ax2.tick_params(labelsize=6)
        self.ax2.grid(True)
        self.canvas.draw()
        self.select_button.config(state=tk.NORMAL)
        self.radio_user.config(state=tk.NORMAL if self.intersection_polygons or self.union_polygons else tk.DISABLED)
        self.radio_combine.config(state=tk.NORMAL if self.intersection_polygons or self.union_polygons else tk.DISABLED)
        self.load_data_button.config(state=tk.NORMAL)
        self.generate_bepit_button.config(state=tk.NORMAL)
        self.run_button.config(state=tk.NORMAL)
        self.save_button.config(state=tk.NORMAL if self.missing_points is not None else tk.DISABLED)

    def run_iteration(self):
        if not self.is_running:
            print("Run Iteration button clicked")
            self.is_running = True
            self.stop_iteration = False
            self.run_button.config(text="Stop Iteration", command=self.stop_iteration_action)
            self.select_button.config(state=tk.DISABLED)
            self.load_data_button.config(state=tk.DISABLED)
            self.generate_bepit_button.config(state=tk.DISABLED)
            self.reset_button.config(state=tk.DISABLED)
            try:
                grid_size = int(self.param_entries['grid_size'].get())
                percentile = float(self.param_entries['percentile'].get())
                min_cells = int(self.param_entries['min_cells'].get())
                epsilon = float(self.param_entries['epsilon'].get())
                max_iter = int(self.param_entries['max_iter'].get())
                if grid_size <= 0 or min_cells <= 0 or epsilon <= 0 or max_iter <= 0 or not (0 <= percentile <= 100):
                    raise ValueError("Parameters must be positive, Percentile must be between 0 and 100")
            except ValueError as e:
                messagebox.showerror("Error", f"Invalid parameter values: {e}")
                self.reset_run_button()
                return
            if self.times is None or self.points is None:
                messagebox.showerror("Error", "Please generate initial BEPIT first.")
                self.reset_run_button()
                return
            option = self.option_var.get()
            print(f"Selected option: {option}")
            if option in ["user", "combine"] and not (self.intersection_polygons or self.union_polygons):
                messagebox.showerror("Error", "Please select at least one polygon first.")
                self.reset_run_button()
                return
            try:
                final_poly = self.final_poly
                if option == "user":
                    final_poly = combine_polygons(Polygon([(0,0), (1,0), (1,1), (0,1)]), self.intersection_polygons, self.union_polygons)
                    if final_poly is None or not isinstance(final_poly, (Polygon, MultiPolygon)):
                        messagebox.showerror("Error", "Combined user polygons are invalid or empty.")
                        self.reset_run_button()
                        return
                elif option == "combine":
                    combined_poly = combine_polygons(self.final_poly, self.intersection_polygons, self.union_polygons)
                    if combined_poly is None or not isinstance(combined_poly, (Polygon, MultiPolygon)):
                        messagebox.showerror("Error", "Initial combination of grid and user polygons is invalid, using grid polygon.")
                        final_poly = self.final_poly
                    else:
                        final_poly = combined_poly
                self.ax3.clear()
                self.ax3.set_visible(True)
                self.ax3.set_title('Iteration Process', fontsize=10)
                self.ax3.set_xlabel('Empirical Time', fontsize=8)
                self.ax3.set_ylabel('Empirical Magnitude', fontsize=8)
                self.ax3.tick_params(labelsize=6)
                self.ax3.grid(True)
                self.ax3.set_xlim(-0.05, 1.05)
                self.ax3.set_ylim(-0.05, 1.05)
                self.ax3.set_box_aspect(1)
                self.canvas.draw()
                self.transformed_points, self.transformed_poly, w_1, w_2 = iterate_bepit(
                    self.points, final_poly, grid_size, min_cells, percentile, epsilon, max_iter,
                    ax=self.ax3, canvas=self.canvas, option=option,
                    intersection_polygons=self.intersection_polygons, union_polygons=self.union_polygons, app=self
                )
                self.complete_iteration()
            except Exception as e:
                messagebox.showerror("Error", f"Failed to run iteration: {e}")
                self.reset_run_button()
        else:
            self.stop_iteration_action()

    def stop_iteration_action(self):
        print("Stop Iteration button clicked")
        self.stop_iteration = True

    def complete_iteration(self):
        try:
            n_mis = estimate_missing_count(self.transformed_poly, self.transformed_points)
            self.missing_points_star = generate_missing_points(self.transformed_poly, n_mis)
            print(f"Generated missing points: {len(self.missing_points_star)}")
            self.missing_points = inverse_transform_points(
                self.missing_points_star, self.times, self.magnitudes, self.transformed_poly, self.transformed_points
            )
            complete_points = merge_points(self.times, self.magnitudes, self.missing_points)
            for ax in [self.ax4, self.ax5, self.ax6]:
                ax.clear()
                ax.set_visible(True)
            plot_final_result(
                self.ax4, self.ax5, self.ax6, self.times, self.magnitudes, self.missing_points,
                self.final_poly, self.points, self.transformed_points, self.missing_points_star, self.transformed_poly
            )
            self.fig.tight_layout()
            self.canvas.draw()
            print("\nFinal point set size:", len(complete_points))
            print("Number of replenished points:", len(self.missing_points))
            self.reset_run_button()
            self.save_button.config(state=tk.NORMAL)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to complete iteration: {e}")
            self.reset_run_button()

    def reset_run_button(self):
        self.is_running = False
        self.stop_iteration = False
        self.run_button.config(text="Run Iteration", command=self.run_iteration)
        self.select_button.config(state=tk.NORMAL if self.points is not None else tk.DISABLED)
        self.load_data_button.config(state=tk.NORMAL)
        self.generate_bepit_button.config(state=tk.NORMAL if self.times is not None else tk.DISABLED)
        self.reset_button.config(state=tk.NORMAL)

    def save_replenished_data(self):
        if self.missing_points is None or self.times is None or self.magnitudes is None:
            messagebox.showerror("Error", "No replenished data to save.")
            return
        save_path = filedialog.asksaveasfilename(
            defaultextension=".txt", filetypes=[("Text files", "*.txt")], title="Save Replenished Data"
        )
        if not save_path:
            return
        try:
            original_points = np.column_stack((self.times, self.magnitudes, np.zeros(len(self.times))))
            replenished_points = np.column_stack((self.missing_points[:, 0], self.missing_points[:, 1], np.ones(len(self.missing_points))))
            all_points = np.vstack([original_points, replenished_points])
            sorted_indices = np.argsort(all_points[:, 0])
            all_points = all_points[sorted_indices]
            np.savetxt(save_path, all_points, fmt='%.6f %.6f %d', header="Time Magnitude Flag")
            print(f"Replenished data saved to {save_path}")
            messagebox.showinfo("Success", f"Data saved successfully to {save_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save data: {e}")

if __name__ == "__main__":
    root = tk.Tk()
    app = BEPITApp(root)
    root.mainloop()