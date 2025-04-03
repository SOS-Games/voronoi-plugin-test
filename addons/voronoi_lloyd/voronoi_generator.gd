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
	# --- Debug Draw (Optional - can remove if not needed) ---
	if debug_draw_pre_clip_polygons:
		# This will now draw the same as the final vertices if calc works
		var pre_clip_color = Color.YELLOW
		var vertex_color = Color.ORANGE
		for i in range(_cells.size()):
			if _cells[i] and _cells[i].has("pre_clip_vertices"):
				var verts: PackedVector2Array = _cells[i].pre_clip_vertices
				if verts.size() >= 2:
					for j in range(verts.size()):
						draw_line(verts[j], verts[(j + 1) % verts.size()], pre_clip_color, 1.0, false)
					for p in verts: draw_circle(p, 3, vertex_color)
		# Draw bounds for reference
		draw_rect(Rect2(Vector2.ZERO, bounds.size), Color.RED.lightened(0.5), false, 2.0)
		# Optional: return here if only debugging

	# --- End Debug Draw ---

	# Draw background bounds if no cells exist
	if _cells.is_empty() and _sites.is_empty():
		if not debug_draw_pre_clip_polygons:
			draw_rect(Rect2(Vector2.ZERO, bounds.size), Color.DARK_GRAY, false, 1.0)
		return

	# --- Draw Final Cells (Fill and Edges) ---
	if fill_cells or draw_edges:
		for i in range(_cells.size()):
			# Check if cell data and vertices exist and are valid
			if not _cells.get(i) or not _cells[i].has("vertices"): continue # Safe check using .get()
			var cell = _cells[i]
			# *** Use 'vertices' directly - no fallback or snapping needed ***
			var vertices_to_draw : PackedVector2Array = cell.vertices

			if vertices_to_draw.size() < 3: continue # Skip degenerate polygons

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
	# _calculate_centroids should now work correctly using the accurate 'vertices'
	if draw_centroids:
		if _centroids.size() == _cells.size(): # Compare with actual cell count
			for centroid in _centroids:
				if is_finite(centroid.x) and is_finite(centroid.y):
					draw_circle(centroid, 3, centroid_color)

	# --- Draw Sites ---
	# Draw the sites that were actually processed (unique sites)
	if draw_sites:
		for cell in _cells: # Iterate through the generated cells
			if cell and cell.has("site"):
				var site = cell.site
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

## Calculates the Voronoi diagram using iterative Sutherland-Hodgman clipping.
## Returns true on success, false on major failure.
func _calculate_voronoi() -> bool:
	_cells.clear()

	if _sites.is_empty(): return true
	if _sites.size() == 1:
		var bounds_poly = PackedVector2Array([ Vector2.ZERO, Vector2(bounds.size.x, 0), bounds.size, Vector2(0, bounds.size.y) ])
		_cells.append({"site": _sites[0], "vertices": bounds_poly, "pre_clip_vertices": bounds_poly.duplicate()})
		return true

	# --- Sutherland-Hodgman Approach ---
	var num_sites = _sites.size()
	var temp_cells : Array[Dictionary] = []
	temp_cells.resize(num_sites)

	# Initial bounds polygon
	var initial_bounds_polygon = PackedVector2Array([
		Vector2.ZERO, Vector2(bounds.size.x, 0),
		bounds.size, Vector2(0, bounds.size.y)
	])

	# --- Pre-filter duplicate sites (Important for this method) ---
	var sites_to_process : Array[Vector2] = []
	var site_positions_seen = {}
	for site in _sites:
		var rounded_pos_str = "%.4f,%.4f" % [site.x, site.y]
		if not site_positions_seen.has(rounded_pos_str):
			sites_to_process.append(site)
			site_positions_seen[rounded_pos_str] = true

	var num_unique_sites = sites_to_process.size()
	if num_unique_sites < _sites.size():
		printerr("Warning: Duplicate site positions detected and removed.")
	if num_unique_sites < 1: return true

	temp_cells.clear()
	temp_cells.resize(num_unique_sites)
	# --- End filtering ---


	# Process each unique site
	for i in range(num_unique_sites):
		var site_s = sites_to_process[i]
		# Start with the full bounds polygon
		var current_polygon = initial_bounds_polygon.duplicate()

		# Clip against all *other* unique sites using S-H
		for j in range(num_unique_sites):
			if i == j: continue

			var site_o = sites_to_process[j]

			# Calculate perpendicular bisector components
			var mid = (site_s + site_o) * 0.5
			var dir_s_to_o = site_o - site_s
			var len_sq = dir_s_to_o.length_squared()
			if len_sq < 1e-9: continue # Skip coincident points

			var dir_s_to_o_normalized = dir_s_to_o / sqrt(len_sq)
			# Normal vector pointing from bisector towards site_s (the half-plane to keep)
			var normal_to_s = -dir_s_to_o_normalized

			# Clip current_polygon against the line defined by (mid, normal_to_s)
			current_polygon = _clip_polygon_sutherland_hodgman(current_polygon, mid, normal_to_s)

			# If polygon becomes degenerate after a clip, stop clipping for this cell
			if current_polygon.size() < 3:
				# print("Site %d: Polygon became degenerate after clip against Site %d." % [i, j])
				current_polygon = PackedVector2Array() # Ensure it's empty
				break

		# Store final result for site i
		if current_polygon.size() < 3:
			# print("Site %d (%s): Final polygon degenerate." % [i, str(site_s)])
			current_polygon = PackedVector2Array()

		# Store in both vertices and pre_clip_vertices for consistency
		temp_cells[i] = {"site": site_s, "vertices": current_polygon, "pre_clip_vertices": current_polygon.duplicate()}


	_cells = temp_cells

	# Rebuild _sites based on the sites actually stored in _cells
	var final_sites : Array[Vector2] = []
	for cell_data in _cells:
		if cell_data and cell_data.has("site"):
			final_sites.append(cell_data.site)
	_sites = final_sites

	return true


## Clips a polygon using the Sutherland-Hodgman algorithm against a single line (half-plane).
## Keeps the side where dot(point - point_on_line, normal) >= 0.
func _clip_polygon_sutherland_hodgman(subject_polygon: PackedVector2Array,
									 point_on_line: Vector2, normal: Vector2) -> PackedVector2Array:

	if subject_polygon.size() < 3: return PackedVector2Array() # Cannot clip degenerate polygon

	var output_list : Array[Vector2] = []
	var input_polygon = subject_polygon # Use subject directly

	# Ensure normal is normalized (important for distance checks)
	var norm = normal.normalized()
	if not is_finite(norm.x): norm = Vector2(0,1) # Fallback if normal was zero

	# Iterate through each edge of the subject polygon (p1 -> p2)
	var p1 = input_polygon[input_polygon.size() - 1] # Start with the last vertex
	for p2_idx in range(input_polygon.size()):
		var p2 = input_polygon[p2_idx]

		# Calculate signed distance from line for p1 and p2
		# Positive distance means "inside" the half-plane we want to keep
		var dist_p1 = (p1 - point_on_line).dot(norm)
		var dist_p2 = (p2 - point_on_line).dot(norm)

		var p1_is_inside = dist_p1 >= -1e-7 # Use tolerance for floating point checks
		var p2_is_inside = dist_p2 >= -1e-7

		if p1_is_inside and p2_is_inside:
			# Case 1: Both points inside - Output p2
			output_list.append(p2)
		elif p1_is_inside and not p2_is_inside:
			# Case 2: Edge goes from inside to outside - Output intersection
			if abs(dist_p1 - dist_p2) > 1e-9: # Avoid division by zero if points are very close
				var intersect_ratio = dist_p1 / (dist_p1 - dist_p2)
				var intersection_point = p1.lerp(p2, intersect_ratio)
				output_list.append(intersection_point)
			# else: print("SH Clip: Degenerate edge detected (inside->outside)") # Optional debug
		elif not p1_is_inside and p2_is_inside:
			# Case 3: Edge goes from outside to inside - Output intersection, then p2
			if abs(dist_p1 - dist_p2) > 1e-9: # Avoid division by zero
				var intersect_ratio = dist_p1 / (dist_p1 - dist_p2)
				var intersection_point = p1.lerp(p2, intersect_ratio)
				output_list.append(intersection_point)
			# else: print("SH Clip: Degenerate edge detected (outside->inside)") # Optional debug
			output_list.append(p2)
		# Case 4: Both points outside - Output nothing

		# Move to next edge
		p1 = p2

	return PackedVector2Array(output_list) # Return the clipped vertices


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
## using the accurately calculated 'vertices' and clamping the result.
func _calculate_centroids():
	# Calculate based on the actual number of cells generated
	var num_cells = _cells.size()
	_centroids.clear()
	_centroids.resize(num_cells) # Resize based on actual cell count

	# Define the min and max corners for clamping based on bounds.size
	var clamp_min = Vector2.ZERO
	var clamp_max = bounds.size

	# Iterate through the generated cells
	for i in range(num_cells):
		# --- Determine fallback position ---
		var fallback_pos : Vector2 = Vector2.ZERO
		if _cells[i] and _cells[i].has("site"):
			fallback_pos = _cells[i].site
		# elif i < _sites.size(): fallback_pos = _sites[i] # Less reliable if duplicates removed

		# --- Get Polygon ---
		var polygon_for_centroid := PackedVector2Array()
		if _cells[i] and _cells[i].has("vertices") and _cells[i].vertices and _cells[i].vertices.size() >= 3:
			polygon_for_centroid = _cells[i].vertices
		# else: print("Centroid Calc: Cell %d has no valid vertices." % i)

		# --- Calculate Centroid ---
		var calculated_centroid : Vector2
		if polygon_for_centroid != null:
			calculated_centroid = _calculate_polygon_centroid(polygon_for_centroid, fallback_pos)
		else:
			calculated_centroid = fallback_pos # Use fallback if no valid polygon

		# --- Clamp Centroid to Bounds using Vector2.clamp ---
		var final_centroid = calculated_centroid.clamp(clamp_min, clamp_max)

		# --- Store Result ---
		if not is_finite(final_centroid.x) or not is_finite(final_centroid.y):
			# print("Centroid calc: Clamped centroid is invalid for cell %d. Using fallback pos." % i)
			_centroids[i] = fallback_pos # Ultimate fallback
		else:
			_centroids[i] = final_centroid

		# --- Update Original _sites Array for Next Iteration ---
		# IMPORTANT: We need to map the centroid back to the correct original site index
		# if duplicates were removed. This is tricky.
		# For simplicity NOW, let's update the internal _sites array directly with unique sites
		# after the first iteration. This breaks if the user expects _sites to remain unmodified.
		# A better approach requires mapping unique back to original indices.

	# --- Simple (but potentially problematic) Site Update ---
	# After calculating all centroids for the unique sites, update _sites array
	# This assumes the order in _centroids matches the order of unique sites in _cells
	if _centroids.size() == num_cells: # Only update if centroid calc completed for all cells
		# --- Explicitly type new_sites ---
		var new_sites : Array[Vector2] = []
		# ---------------------------------
		for i in range(num_cells):
			# Check if centroid is valid before adding
			if _centroids[i] != null and is_finite(_centroids[i].x) and is_finite(_centroids[i].y):
				new_sites.append(_centroids[i])
			else:
				# If centroid invalid, try to keep the corresponding site from _cells?
				if i < _cells.size() and _cells[i] and _cells[i].has("site"):
					new_sites.append(_cells[i].site) # Keep old site position
				# else: # If we can't find old site, maybe skip? Or add Vector2.ZERO?
				#   print("Warning: Could not preserve site position during centroid update for index %d" % i)


		# If fewer unique sites than original, the _sites array shrinks.
		# Assigning Array[Vector2] to Array[Vector2]
		_sites = new_sites
	# else: print("Warning: Centroid count did not match cell count. _sites not updated.")

	# --- End Simple Site Update ---


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
