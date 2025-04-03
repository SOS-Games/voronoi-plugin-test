@tool
class_name VoronoiGenerator
extends Node2D

## Generates Voronoi diagrams and applies Lloyd's relaxation using Godot 4 features.

# --- Signals ---
## Emitted when the generation process (including relaxation) is complete.
signal generation_finished(sites: Array[Vector2], cells: Array[Dictionary])

# --- Exports ---
## Seed for random point generation. 0 uses a random seed.
@export var seed : int = 0: set = set_seed
## Number of initial random points (sites) to generate.
@export var num_points : int = 50: set = set_num_points
## The rectangular bounds for point generation and diagram clipping.
@export var bounds := Rect2(0, 0, 500, 500): set = set_bounds
## Number of Lloyd's relaxation iterations to apply.
@export var lloyd_iterations : int = 3: set = set_lloyd_iterations
## If true, automatically generate on _ready() or when parameters change in editor.
@export var generate_on_change := true
@export var boundary_margin_factor : float = 2.0 # How much larger the dummy bounding box is than the site bounds

## --- Visualization ---
@export_group("Visualization")
## Draw the site points.
@export var draw_sites := false: set = _set_redraw
## Draw the edges of the Voronoi cells.
@export var draw_edges := true: set = _set_redraw
## Draw the calculated centroids (useful for debugging Lloyd's).
@export var draw_centroids := false: set = _set_redraw
## Fill the Voronoi cells with color.
@export var fill_cells := true: set = _set_redraw
@export var debug_draw_pre_clip_polygons := true: set = _set_redraw
## Color for drawing sites.
@export var site_color := Color.WHITE
## Color for drawing edges.
@export var edge_color := Color.GRAY
## Color for drawing centroids.
@export var centroid_color := Color.RED
## List of colors to cycle through for filling cells.
@export var cell_colors : Array[Color] = [
	Color.from_hsv(0.0, 0.6, 0.8), Color.from_hsv(0.1, 0.6, 0.8),
	Color.from_hsv(0.2, 0.6, 0.8), Color.from_hsv(0.3, 0.6, 0.8),
	Color.from_hsv(0.4, 0.6, 0.8), Color.from_hsv(0.5, 0.6, 0.8),
	Color.from_hsv(0.6, 0.6, 0.8), Color.from_hsv(0.7, 0.6, 0.8),
	Color.from_hsv(0.8, 0.6, 0.8), Color.from_hsv(0.9, 0.6, 0.8)
]

# --- Internal State ---
var _sites: Array[Vector2] = [] # The input points (using standard Array for flexibility during relaxation)
# Cells store the site and its vertices (PackedVector2Array for performance)
# Structure: [{ "site": Vector2, "vertices": PackedVector2Array }, ...]
var _cells: Array[Dictionary] = []
var _centroids: Array[Vector2] = [] # Calculated centroids for Lloyd's
var _needs_regeneration := true
var _bounds_polygon: PackedVector2Array # Cache bounds as polygon (Godot 4 uses PackedVector2Array)

func _init():
	_update_bounds_polygon()

# --- Setters ---
func set_seed(value: int):
	if seed != value:
		seed = value
		_needs_regeneration = true
		_trigger_generation_and_redraw()

func set_num_points(value: int):
	if num_points != value:
		num_points = max(0, value) # Ensure non-negative
		_needs_regeneration = true
		_trigger_generation_and_redraw()

# Modify set_bounds to ensure position is always (0,0) internally
func set_bounds(value: Rect2):
	# Ensure bounds position is always zero for the parent clipping approach
	var new_bounds = Rect2(Vector2.ZERO, value.size)
	if bounds != new_bounds:
		bounds = new_bounds # Store only size, position is implicitly (0,0)
		_update_bounds_polygon()
		_needs_regeneration = true
		_trigger_generation_and_redraw()
	# If only size changed but position was already 0,0
	elif bounds.size != value.size:
		bounds = new_bounds
		_update_bounds_polygon()
		_needs_regeneration = true
		_trigger_generation_and_redraw()

func set_lloyd_iterations(value: int):
	if lloyd_iterations != value:
		lloyd_iterations = max(0, value) # Ensure non-negative
		_needs_regeneration = true
		_trigger_generation_and_redraw()

func _set_redraw(_value):
	# Generic setter for visualization flags to trigger redraw
	queue_redraw() # Use queue_redraw in Godot 4

func _update_bounds_polygon():
	# Define bounds polygon relative to node's origin (0,0), using bounds.size
	# This represents the rectangle we *attempt* to clip against internally
	_bounds_polygon = PackedVector2Array([
		Vector2.ZERO,
		Vector2(bounds.size.x, 0.0),
		bounds.size, # Equivalent to Vector2(bounds.size.x, bounds.size.y)
		Vector2(0.0, bounds.size.y)
	])

func _trigger_generation_and_redraw():
	if generate_on_change and _needs_regeneration:
		generate() # Regenerate data
	else:
		queue_redraw() # Just redraw if only visualization changed

# --- Godot Methods ---
func _ready():
	# Use queue_redraw instead of update in Godot 4
	if Engine.is_editor_hint():
		if generate_on_change and _needs_regeneration:
			generate()
		else:
			queue_redraw()
	elif generate_on_change: # At runtime
		generate()

func _draw():
	# --- Debug Draw Pre-Clip Polygons (Optional) ---
	if debug_draw_pre_clip_polygons:
		# ... (debug drawing code as before) ...
		# Draw the bounds rectangle clearly relative to (0,0)
		draw_rect(Rect2(Vector2.ZERO, bounds.size), Color.RED.lightened(0.5), false, 2.0)
		# Optional: return here if you only want to see debug shapes

	# --- End Debug Draw ---
	# Only draw the bounds outline if nothing else will be drawn AND debug isn't active
	if _cells.is_empty() and _sites.is_empty():
		if not debug_draw_pre_clip_polygons:
			draw_rect(Rect2(Vector2.ZERO, bounds.size), Color.DARK_GRAY, false, 1.0) # Draw bounds relative to 0,0
		return

	# --- Draw Final Cells (Fill and Edges) using Fallback Logic ---
	# This part uses the fallback logic from the previous step.
	# It draws the pre_clip_vertices if clipping failed,
	# and the parent Control's clip_contents handles the visual bounding.
	if fill_cells or draw_edges:
		for i in range(_cells.size()):
			if not _cells[i] or not _cells[i].has("vertices"): continue

			var cell = _cells[i]
			var vertices_to_draw: PackedVector2Array = cell.vertices

			# Fallback to pre_clip_vertices if 'vertices' is invalid
			if vertices_to_draw.size() < 3 and cell.has("pre_clip_vertices"):
				var pre_clip_verts = cell.pre_clip_vertices
				if pre_clip_verts and pre_clip_verts.size() >= 3:
					vertices_to_draw = pre_clip_verts

			if vertices_to_draw.size() < 3: continue # Skip if still invalid

			# --- Fill Cell ---
			if fill_cells:
				var color_index = i % cell_colors.size()
				var color = cell_colors[color_index] if not cell_colors.is_empty() else Color.WHITE
				draw_colored_polygon(vertices_to_draw, color)

			# --- Draw Edges ---
			if draw_edges:
				for j in range(vertices_to_draw.size()):
					var p1 = vertices_to_draw[j]
					var p2 = vertices_to_draw[(j + 1) % vertices_to_draw.size()]
					draw_line(p1, p2, edge_color, 1.0)

	# --- Draw Centroids ---
	# (Code remains the same - uses _centroids calculated relative to (0,0))
	if draw_centroids:
		if _centroids.size() == _sites.size():
			for centroid in _centroids:
				if is_finite(centroid.x) and is_finite(centroid.y):
					draw_circle(centroid, 3, centroid_color)

	# --- Draw Sites ---
	# (Code remains the same - uses _sites generated relative to (0,0))
	if draw_sites:
		for site in _sites:
			draw_circle(site, 4, site_color)
			draw_circle(site, 2, Color.BLACK)

# --- Core Logic ---

## Clears the current diagram and sites.
func clear():
	_sites.clear()
	_cells.clear()
	_centroids.clear()
	_needs_regeneration = true
	queue_redraw()

## Generates the initial points, Voronoi diagram, and applies Lloyd's relaxation.
func generate():
	print("Generating Voronoi diagram (Godot 4)...")
	if num_points < 3:
		printerr("Voronoi generation requires at least 3 points.")
		clear()
		return

	_generate_initial_points()
	_centroids.clear() # Clear previous centroids

	# --- Lloyd's Relaxation Loop ---
	for i in range(lloyd_iterations + 1): # +1 because we need initial Voronoi
		print("Iteration %d / %d" % [i, lloyd_iterations])
		# 1. Calculate Voronoi Diagram for current sites
		var success = _calculate_voronoi()
		if not success:
			printerr("Failed to calculate Voronoi diagram.")
			break # Stop iteration if Voronoi fails

		# 2. If it's not the last iteration, calculate centroids and update sites
		if i < lloyd_iterations:
			_calculate_centroids()
			if _centroids.size() != _sites.size():
				printerr("Mismatch between number of centroids (%d) and sites (%d). Stopping relaxation." % [_centroids.size(), _sites.size()])
				break # Avoid errors if centroid calculation failed for some cells

			# Update sites to be the centroids for the next iteration
			_sites = _centroids.duplicate() # Use duplicate to avoid reference issues
		else:
			# On the final iteration, centroids might still be useful for drawing
			if draw_centroids:
				_calculate_centroids()


	_needs_regeneration = false
	queue_redraw()
	print("Generation finished.")
	# Use new signal emission syntax in Godot 4
	generation_finished.emit(_sites, _cells)


func _generate_initial_points():
	_sites.clear()
	if seed == 0: randomize()
	else: seed(seed)

	var rng = RandomNumberGenerator.new()
	if seed != 0: rng.seed = seed

	for _i in range(num_points):
		# Generate points within the bounds size, relative to node's (0,0)
		var px = rng.randf_range(0.0, bounds.size.x)
		var py = rng.randf_range(0.0, bounds.size.y)
		_sites.append(Vector2(px, py))

## Calculates the Voronoi diagram for the current `_sites`.
## Includes dummy boundary points strategy and stores pre-clip polygons for debugging.
## Returns true on success (even if some cells are empty), false on major failure.
func _calculate_voronoi() -> bool:
	_cells.clear()
	if _sites.size() < 3:
		printerr("Voronoi requires at least 3 sites.")
		return false

	# --- Strategy: Add Dummy Boundary Points ---
	var max_bound_size = max(bounds.size.x, bounds.size.y)
	if max_bound_size <= 0: max_bound_size = 1000.0 # Avoid zero size
	var margin = max_bound_size * max(1.0, boundary_margin_factor)
	# Define local bounds rect starting at 0,0 for effective bounds calculation
	var local_bounds_rect = Rect2(Vector2.ZERO, bounds.size)
	var effective_bounds = local_bounds_rect.grow(margin) # Grow from local bounds

	# Define dummy points relative to effective_bounds corners
	var dummy_points = [
		effective_bounds.position,
		Vector2(effective_bounds.end.x, effective_bounds.position.y),
		effective_bounds.end,
		Vector2(effective_bounds.position.x, effective_bounds.end.y)
	]

	# Combine original sites and dummy points for triangulation
	var all_points_list : Array[Vector2] = []
	all_points_list.assign(_sites) # Copy original sites
	all_points_list.append_array(dummy_points) # Add dummy points

	var num_original_sites = _sites.size()
	var num_total_points = all_points_list.size()

	if num_total_points < 3: return false # Should not happen now

	var all_points_packed = PackedVector2Array(all_points_list)

	# Perform Delaunay triangulation on the combined set of points
	var delaunay_triangles: PackedInt32Array = Geometry2D.triangulate_delaunay(all_points_packed)
	# --- End Strategy ---

	# Check if triangulation produced any result
	if delaunay_triangles.is_empty() and all_points_packed.size() >= 3:
		printerr("Delaunay triangulation failed (with dummy points), returned empty.")
		return false # Consider this a failure if we expected triangles
	elif delaunay_triangles.is_empty():
		print("Delaunay triangulation resulted in no triangles (num points < 3?).")
		return true # No triangles likely means no sites or fewer than 3 total points

	# --- Construct Voronoi Cells from Delaunay Triangles ---
	var circumcenters = {} # Store circumcenter by triangle index (key = i / 3)

	# Precompute circumcenters using all points
	for i in range(0, delaunay_triangles.size(), 3):
		var p1_idx = delaunay_triangles[i]
		var p2_idx = delaunay_triangles[i+1]
		var p3_idx = delaunay_triangles[i+2]

		# Check for valid indices (safeguard against potential issues in triangulate_delaunay)
		if not (p1_idx >= 0 and p1_idx < num_total_points and \
				p2_idx >= 0 and p2_idx < num_total_points and \
				p3_idx >= 0 and p3_idx < num_total_points):
			printerr("Invalid index from Delaunay triangulation: T:%d, Idxs:(%d, %d, %d), TotalPts:%d" % [i/3, p1_idx, p2_idx, p3_idx, num_total_points])
			continue # Skip this potentially corrupt triangle

		# Indices refer to the combined array (original + dummy)
		var p1 = all_points_packed[p1_idx]
		var p2 = all_points_packed[p2_idx]
		var p3 = all_points_packed[p3_idx]

		var circumcenter = _calculate_circumcenter(p1, p2, p3)

		# Skip INF circumcenters (degenerate triangles)
		if is_inf(circumcenter.x) or is_inf(circumcenter.y):
			print("Skipping INF circumcenter for triangle %d" % (i/3))
			continue # Do not add this circumcenter to the list for any cell

		# If circumcenter is valid, store it
		var triangle_index = i / 3
		circumcenters[triangle_index] = circumcenter

	# Build mapping from original site index to relevant triangle indices
	var site_to_triangles = {} # Map site index (original index) to list of triangle indices
	for i in range(0, delaunay_triangles.size(), 3):
		var triangle_index = i / 3
		# Only consider triangles for which we successfully calculated a circumcenter
		if not circumcenters.has(triangle_index):
			continue

		# Check which *original* sites this triangle involves
		for j in range(3):
			var site_index_in_all = delaunay_triangles[i+j] # Index in the combined array

			# IMPORTANT: Only map triangles to ORIGINAL sites
			if site_index_in_all < num_original_sites:
				var original_site_index = site_index_in_all
				# Initialize array if first time seeing this site
				if not site_to_triangles.has(original_site_index):
					site_to_triangles[original_site_index] = []
				# Add triangle index to this site's list
				site_to_triangles[original_site_index].append(triangle_index)

	# Temporary array to build the final _cells list correctly ordered
	var temp_cells : Array[Dictionary] = []
	temp_cells.resize(num_original_sites) # Pre-allocate space for all original sites

	# Create polygons for each ORIGINAL site
	for original_site_index in range(num_original_sites): # Loop only up to original site count
		var site = _sites[original_site_index] # Get original site Vector2

		# Prepare the dictionary to store cell data, including pre-clip vertices
		var current_cell_data = {
			"site": site,
			"vertices": PackedVector2Array(), # Placeholder for final clipped vertices
			"pre_clip_vertices": PackedVector2Array() # Placeholder for polygon before clipping
		}

		# Check if this site was associated with any valid triangles
		if not site_to_triangles.has(original_site_index):
			print("Warning: Original site %d not associated with any valid triangles." % original_site_index)
			temp_cells[original_site_index] = current_cell_data # Store default empty cell data
			continue

		# Collect the circumcenters for this site's associated triangles
		var associated_triangles = site_to_triangles[original_site_index]
		var cell_points_unordered : Array[Vector2] = []
		for triangle_index in associated_triangles:
			# Circumcenter existence already checked when building site_to_triangles map
			cell_points_unordered.append(circumcenters[triangle_index])

		# --- Order the vertices ---
		# Check if enough points before sorting
		if cell_points_unordered.size() < 3:
			print("Site %d (Original): Not enough points (%d) to form polygon before sorting." % [original_site_index, cell_points_unordered.size()])
			temp_cells[original_site_index] = current_cell_data # Store default empty cell data
			continue

		# Sort points angularly around the site (using the improved robust function)
		var cell_points_ordered: Array[Vector2] = _sort_points_angularly(site, cell_points_unordered)

		# --- Clipping ---
		# Convert ordered points to PackedVector2Array for geometry operations
		var cell_polygon_packed = PackedVector2Array(cell_points_ordered)

		# Store the pre-clip polygon for potential debugging
		current_cell_data["pre_clip_vertices"] = cell_polygon_packed

		# Optional Debugging:
		# print("Site %d (Original): Pre-clipping vertices: %d points" % [original_site_index, cell_polygon_packed.size()])

		# Check if the polygon is valid *before* clipping
		if cell_polygon_packed.size() < 3:
			print("Site %d (Original): Skipping clipping (degenerate polygon pre-clip with %d points)" % [original_site_index, cell_polygon_packed.size()])
			temp_cells[original_site_index] = current_cell_data # Store data even if skipping clip
			continue

		# Perform the clipping: Intersect the raw cell polygon with the bounds polygon
		# _bounds_polygon should be a PackedVector2Array of the original 'bounds' Rect2 corners
		var clipped_polygons: Array[PackedVector2Array] = Geometry2D.clip_polygons(cell_polygon_packed, _bounds_polygon)

		# Check the result of the clipping operation
		if not clipped_polygons.is_empty() and clipped_polygons[0].size() >= 3:
			# Clipping successful, store the resulting polygon
			var final_polygon = clipped_polygons[0]
			# print("Site %d (Original): Post-clipping vertices: %d points" % [original_site_index, final_polygon.size()])
			current_cell_data["vertices"] = final_polygon # Update cell data with clipped vertices
		else:
			# Clipping failed or resulted in a degenerate polygon (line/point)
			# Keep current_cell_data["vertices"] as the default empty PackedVector2Array
			# Optional Debugging: Print failure reason
			if clipped_polygons.is_empty():
				pass
				#print("Site %d (Original): Clipping resulted in empty polygon array." % original_site_index) # Keep this print enabled
			else:
				print("Site %d (Original): Clipping resulted in polygon with < 3 vertices (%d)." % [original_site_index, clipped_polygons[0].size()])
				pass # This case is less common than empty array

		# Store the processed cell data (including pre_clip_vertices and final vertices)
		temp_cells[original_site_index] = current_cell_data

	# Assign the correctly ordered temporary cells to the final _cells array
	_cells = temp_cells

	# Final validation check (optional but good practice)
	if _cells.size() != num_original_sites:
		printerr("CRITICAL Error: Final cell count (%d) mismatch after processing. Expected %d." % [_cells.size(), num_original_sites])
		# Consider how to handle this state - maybe return false?
		return false

	# Return true indicating the process completed, even if some cells are empty.
	return true


## Calculates the geometric centroid of a polygon.
## Returns the original site position as fallback for degenerate polygons.
func _calculate_polygon_centroid(polygon: PackedVector2Array, fallback_point: Vector2) -> Vector2:
	var n = polygon.size()
	if n < 3:
		# print("Centroid calc helper: Polygon has < 3 vertices. Using fallback.")
		return fallback_point # Cannot calculate centroid for line/point

	var area : float = 0.0
	var centroid_x : float = 0.0
	var centroid_y : float = 0.0

	for i in range(n):
		var p1 = polygon[i]
		var p2 = polygon[(i + 1) % n] # Wrap around for the last vertex

		# Cross product (for area and centroid calculation)
		var cross_product = p1.cross(p2) # Equivalent to p1.x * p2.y - p1.y * p2.x

		area += cross_product
		centroid_x += (p1.x + p2.x) * cross_product
		centroid_y += (p1.y + p2.y) * cross_product

	# Finalize area calculation
	area *= 0.5

	# Avoid division by zero for degenerate polygons (collinear points, zero area)
	if abs(area) < 1e-9: # Use a small epsilon for float comparison
		# print("Centroid calc helper: Degenerate polygon (area close to zero). Using fallback.")
		# Fallback: Calculate average of vertices, or just use the site position
		# Averaging vertices is simple:
		# var avg_pos = Vector2.ZERO
		# for p in polygon: avg_pos += p
		# return avg_pos / n
		# Using fallback_point (site position) is often preferred in Lloyd's context:
		return fallback_point
	else:
		# Finalize centroid calculation
		var factor = 1.0 / (6.0 * area)
		centroid_x *= factor
		centroid_y *= factor
		return Vector2(centroid_x, centroid_y)

## Calculates the circumcenter of a triangle defined by p1, p2, p3.
## Returns Vector2.INF if points are collinear. (No changes needed here)
func _calculate_circumcenter(p1: Vector2, p2: Vector2, p3: Vector2) -> Vector2:
	var d = 2 * (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y))
	if abs(d) < 1e-9:
		return Vector2.INF # Indicate degeneracy

	var p1_sq = p1.length_squared()
	var p2_sq = p2.length_squared()
	var p3_sq = p3.length_squared()

	var center_x = (p1_sq * (p2.y - p3.y) + p2_sq * (p3.y - p1.y) + p3_sq * (p1.y - p2.y)) / d
	var center_y = (p1_sq * (p3.x - p2.x) + p2_sq * (p1.x - p3.x) + p3_sq * (p2.x - p1.x)) / d

	return Vector2(center_x, center_y)


func _sort_points_angularly(p_center: Vector2, p_points: Array[Vector2]) -> Array[Vector2]:
	if p_points.size() < 2:
		return p_points # No sorting needed

	# Filter out points that are virtually identical to the center
	var points_to_sort : Array[Vector2] = []
	for p in p_points:
		# Use a small tolerance for floating point comparison
		if (p - p_center).length_squared() > 1e-9: # Avoid sqrt for efficiency
			points_to_sort.append(p)
		#else:
			# Optional: print("Filtered out point coincident with center during sort: %s" % str(p))


	if points_to_sort.size() < 2:
		# If filtering left too few points, return them
		# print("Angular sort: < 2 points remaining after filtering coincident points.")
		return points_to_sort # Result might be empty or single point

	# Use a custom sorter based on angle
	points_to_sort.sort_custom(func(a, b):
		var angle_a = (a - p_center).angle()
		var angle_b = (b - p_center).angle()
		return angle_a < angle_b
	)
	return points_to_sort


## Calculates the geometric centroid for each Voronoi cell,
## using pre-clip vertices as fallback and clamping the result to bounds.
func _calculate_centroids():
	var num_original_sites = _sites.size()
	_centroids.clear()
	_centroids.resize(num_original_sites) # Resize based on original sites

	var clamp_min = Vector2.ZERO
	var clamp_max = bounds.size

	# Iterate based on the expected number of sites, accessing _cells safely
	for i in range(num_original_sites):
		# --- Determine which polygon to use for centroid calculation ---
		var polygon_for_centroid := PackedVector2Array()
		var fallback_pos : Vector2 = Vector2.ZERO # Default fallback

		if i < _sites.size():
			fallback_pos = _sites[i] # Use current site position as primary fallback

		# Check if cell data exists and is valid
		if i < _cells.size() and _cells[i] and _cells[i].has("site"):
			var cell = _cells[i]
			fallback_pos = cell.site # More specific fallback from cell data if available

			# Prioritize successfully clipped vertices if available and valid
			if cell.has("vertices") and cell.vertices and cell.vertices.size() >= 3:
				polygon_for_centroid = cell.vertices
				# print("Site %d: Using clipped vertices for centroid." % i) # Optional debug
			# Otherwise, try falling back to pre-clip vertices if available and valid
			elif cell.has("pre_clip_vertices") and cell.pre_clip_vertices and cell.pre_clip_vertices.size() >= 3:
				polygon_for_centroid = cell.pre_clip_vertices
				# print("Site %d: Using PRE-CLIP vertices for centroid." % i) # Optional debug
			# else: print("Site %d: No valid polygon found for centroid." % i) # Optional debug
		else:
			# print("Centroid calc: Missing or invalid cell data for site index %d." % i)
			pass # Will use fallback_pos below

		# --- Calculate Centroid ---
		var calculated_centroid : Vector2

		if polygon_for_centroid != null:
			# Calculate centroid using the chosen polygon
			calculated_centroid = _calculate_polygon_centroid(polygon_for_centroid, fallback_pos)
		else:
			# If no valid polygon was found, use the fallback position directly
			calculated_centroid = fallback_pos

		# --- Clamp Centroid to Bounds ---
		# Ensure the next site position stays within the defined bounds rectangle.
		# Assumes bounds.position is (0,0) due to parent clipping setup.
		var final_centroid = calculated_centroid.clamp(clamp_min, clamp_max)

		# --- Store Result ---
		# Check if centroid is valid before storing
		if not is_finite(final_centroid.x) or not is_finite(final_centroid.y):
			# print("Centroid calc: Clamped centroid is invalid for site %d. Using fallback pos." % i)
			_centroids[i] = fallback_pos # Ultimate fallback
		else:
			_centroids[i] = final_centroid


# --- Public Methods ---

## Get the current list of site points.
func get_sites() -> Array[Vector2]:
	return _sites.duplicate() # Return a copy

## Get the calculated Voronoi cells.
## Each cell is a Dictionary: {site: Vector2, vertices: PackedVector2Array}
func get_cells() -> Array[Dictionary]:
	var cells_copy = []
	for cell in _cells:
		cells_copy.append({
			"site": cell.site,
			# Duplicate PackedVector2Array
			"vertices": cell.vertices.duplicate() if cell.vertices else PackedVector2Array()
		})
	return cells_copy
