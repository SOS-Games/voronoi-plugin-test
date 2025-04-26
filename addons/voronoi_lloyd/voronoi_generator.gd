@tool
class_name VoronoiGenerator
extends Node2D

## Generates Voronoi diagrams using iterative clipping and applies Lloyd's relaxation.

# --- Signals ---
signal generation_finished(sites: Array[Vector2], cells: Array[Dictionary])

# --- Exports ---
@export var seed : int = 0: set = set_seed
@export var num_points : int = 50: set = set_num_points
## Bounds are relative to this node's origin (0,0). Use a parent Control node with Clip Contents for visual bounding.
@export var bounds := Rect2(0, 0, 500, 500): set = set_bounds
@export var lloyd_iterations : int = 3: set = set_lloyd_iterations
@export var generate_on_change := true

## --- Visualization ---
@export_group("Visualization")
@export var draw_sites := false: set = _set_redraw
@export var draw_edges := true: set = _set_redraw
@export var draw_centroids := false: set = _set_redraw
@export var fill_cells := true: set = _set_redraw
## Optional debug view: draws the polygon stored in 'pre_clip_vertices' (should match 'vertices' now)
@export var debug_draw_pre_clip_polygons := false: set = _set_redraw
@export var site_color := Color.WHITE
@export var edge_color := Color.GRAY
@export var centroid_color := Color.RED
@export var cell_colors : Array[Color] = [
	Color.from_hsv(0.0, 0.6, 0.8), Color.from_hsv(0.1, 0.6, 0.8),
	Color.from_hsv(0.2, 0.6, 0.8), Color.from_hsv(0.3, 0.6, 0.8),
	Color.from_hsv(0.4, 0.6, 0.8), Color.from_hsv(0.5, 0.6, 0.8),
	Color.from_hsv(0.6, 0.6, 0.8), Color.from_hsv(0.7, 0.6, 0.8),
	Color.from_hsv(0.8, 0.6, 0.8), Color.from_hsv(0.9, 0.6, 0.8)
]

# --- Internal State ---
var _sites: Array[Vector2] = [] # Current site positions (updated by Lloyd's)
# Structure: [{ "site": Vector2, "vertices": PackedVector2Array, "pre_clip_vertices": PackedVector2Array}, ...]
# Note: 'vertices' and 'pre_clip_vertices' should be identical with S-H method
var _cells: Array[Dictionary] = []
var _centroids: Array[Vector2] = []
var _needs_regeneration := true
# No _bounds_polygon needed for S-H method, bounds rect used directly

func _init():
	# No _update_bounds_polygon needed
	pass

# --- Setters ---
func set_seed(value: int):
	if seed != value:
		seed = value
		_needs_regeneration = true
		_trigger_generation_and_redraw()

func set_num_points(value: int):
	if num_points != value:
		num_points = max(0, value)
		_needs_regeneration = true
		_trigger_generation_and_redraw()

func set_bounds(value: Rect2):
	# Ensure bounds position is always zero for the parent clipping approach
	var new_bounds = Rect2(Vector2.ZERO, value.size)
	if bounds != new_bounds:
		bounds = new_bounds
		_needs_regeneration = true
		_trigger_generation_and_redraw()
	elif bounds.size != value.size: # Handle size-only change
		bounds = new_bounds
		_needs_regeneration = true
		_trigger_generation_and_redraw()


func set_lloyd_iterations(value: int):
	if lloyd_iterations != value:
		lloyd_iterations = max(0, value)
		_needs_regeneration = true
		_trigger_generation_and_redraw()

func _set_redraw(_value):
	queue_redraw()

func _trigger_generation_and_redraw():
	if generate_on_change and _needs_regeneration:
		generate()
	else:
		queue_redraw()

# --- Godot Methods ---
func _ready():
	if Engine.is_editor_hint():
		if generate_on_change and _needs_regeneration:
			generate()
		else:
			queue_redraw()
	#elif generate_on_change:
	#	generate()

func _draw():
	# --- Debug Draw (Optional) ---
	if debug_draw_pre_clip_polygons:
		var pre_clip_color = Color.YELLOW
		var vertex_color = Color.ORANGE
		for i in range(_cells.size()):
			if _cells[i] and _cells[i].has("pre_clip_vertices"): # Check key existence
				var verts: PackedVector2Array = _cells[i].pre_clip_vertices
				if verts and verts.size() >= 2: # Check if array is valid
					for j in range(verts.size()):
						draw_line(verts[j], verts[(j + 1) % verts.size()], pre_clip_color, 1.0, false)
					for p in verts: draw_circle(p, 3, vertex_color)
		draw_rect(Rect2(Vector2.ZERO, bounds.size), Color.RED.lightened(0.5), false, 2.0)
		# Optional: return # Uncomment to only show debug draw

	# --- End Debug Draw ---

	if _cells.is_empty() and _sites.is_empty():
		if not debug_draw_pre_clip_polygons:
			draw_rect(Rect2(Vector2.ZERO, bounds.size), Color.DARK_GRAY, false, 1.0)
		return

	# --- Draw Final Cells (Fill and Edges) ---
	if fill_cells or draw_edges:
		for i in range(_cells.size()):
			# Use .get() for safer access, check key and validity
			if not _cells.get(i) or not _cells[i].has("vertices"): continue
			var cell = _cells[i]
			var vertices_to_draw : PackedVector2Array = cell.vertices
			if not vertices_to_draw or vertices_to_draw.size() < 3: continue # Skip invalid/degenerate

			if fill_cells:
				var color_index = i % cell_colors.size()
				var color = cell_colors[color_index] if not cell_colors.is_empty() else Color.WHITE
				draw_colored_polygon(vertices_to_draw, color)

			if draw_edges:
				for j in range(vertices_to_draw.size()):
					var p1 = vertices_to_draw[j]
					var p2 = vertices_to_draw[(j + 1) % vertices_to_draw.size()]
					draw_line(p1, p2, edge_color, 1.0)

	# --- Draw Centroids ---
	if draw_centroids:
		if _centroids.size() == _cells.size():
			for centroid in _centroids:
				if is_finite(centroid.x) and is_finite(centroid.y):
					draw_circle(centroid, 3, centroid_color)

	# --- Draw Sites ---
	if draw_sites:
		# Draw sites based on the potentially modified _sites array
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
	print("Generating Voronoi diagram (Sutherland-Hodgman)...")
	if num_points == 0 :
		clear()
		return # Allow zero points without error

	_generate_initial_points()
	if _sites.is_empty(): # Handle case where num_points was 0 or generation failed
		clear()
		return

	# --- Lloyd's Relaxation Loop ---
	for i in range(lloyd_iterations + 1):
		print("Iteration %d / %d" % [i, lloyd_iterations])
		# 1. Calculate Voronoi Diagram for current sites
		var success = _calculate_voronoi()
		if not success:
			printerr("Failed to calculate Voronoi diagram.")
			break # Stop iteration if Voronoi fails

		# 2. If it's not the last iteration, calculate centroids and update sites
		if i < lloyd_iterations:
			# Check if there are cells to process before calculating centroids
			if _cells.is_empty():
				print("No cells generated, stopping relaxation.")
				break
			_calculate_centroids()
			# Site update is now handled within _calculate_centroids
			if _sites.is_empty(): # Check if site update resulted in empty array
				print("Site array became empty after centroid update, stopping.")
				break
		else:
			# On the final iteration, centroids might still be useful for drawing
			if draw_centroids and not _cells.is_empty():
				_calculate_centroids()

	_needs_regeneration = false
	queue_redraw()
	print("Generation finished.")
	generation_finished.emit(get_sites(), get_cells()) # Emit copies


func _generate_initial_points():
	_sites.clear()
	if num_points <= 0: return # Don't generate if count is zero or less

	if seed == 0: randomize()
	else: seed(seed)

	var rng = RandomNumberGenerator.new()
	if seed != 0: rng.seed = seed

	for _i in range(num_points):
		var px = rng.randf_range(0.0, bounds.size.x)
		var py = rng.randf_range(0.0, bounds.size.y)
		_sites.append(Vector2(px, py))


## Clips a polygon using the Sutherland-Hodgman algorithm against a single line (half-plane).
## Keeps the side where dot(point - point_on_line, normal) >= 0.
func _clip_polygon_sutherland_hodgman(subject_polygon: PackedVector2Array,
									 point_on_line: Vector2, normal: Vector2) -> PackedVector2Array:

	if subject_polygon.size() < 3: return PackedVector2Array()

	var output_list : Array[Vector2] = []
	var input_polygon = subject_polygon

	# Normalize the normal vector for dot product checks
	var norm = normal.normalized()
	# Fallback for zero normal (shouldn't happen with distinct sites)
	if not is_finite(norm.x) or norm.length_squared() < 1e-9: return subject_polygon

	# Tolerance for floating point comparisons near zero
	var tolerance : float = 1e-6

	var p1 = input_polygon[input_polygon.size() - 1]
	var dist_p1 = (p1 - point_on_line).dot(norm)

	for p2_idx in range(input_polygon.size()):
		var p2 = input_polygon[p2_idx]
		var dist_p2 = (p2 - point_on_line).dot(norm)

		var p1_is_inside = dist_p1 >= -tolerance
		var p2_is_inside = dist_p2 >= -tolerance

		if p1_is_inside and p2_is_inside: # Both inside
			output_list.append(p2)
		elif p1_is_inside and not p2_is_inside: # Going out
			var dist_diff = dist_p1 - dist_p2
			if abs(dist_diff) > tolerance:
				var intersect_ratio = dist_p1 / dist_diff
				output_list.append(p1.lerp(p2, intersect_ratio))
		elif not p1_is_inside and p2_is_inside: # Coming in
			var dist_diff = dist_p1 - dist_p2
			if abs(dist_diff) > tolerance:
				var intersect_ratio = dist_p1 / dist_diff
				output_list.append(p1.lerp(p2, intersect_ratio))
			output_list.append(p2)
		# else: # Both outside - output nothing

		p1 = p2
		dist_p1 = dist_p2

	return PackedVector2Array(output_list)


## Calculates the Voronoi diagram using iterative Sutherland-Hodgman clipping.
func _calculate_voronoi() -> bool:
	_cells.clear()

	if _sites.is_empty(): return true

	# --- Pre-filter duplicate sites ---
	var sites_to_process : Array[Vector2] = []
	var site_positions_seen = {}
	for site in _sites:
		var rounded_pos_str = "%.4f,%.4f" % [site.x, site.y] # Adjust precision as needed
		if not site_positions_seen.has(rounded_pos_str):
			sites_to_process.append(site)
			site_positions_seen[rounded_pos_str] = true

	var num_unique_sites = sites_to_process.size()
	if num_unique_sites < _sites.size():
		printerr("Warning: Duplicate site positions detected and removed.")
	if num_unique_sites < 1: return true # No unique sites left
	if num_unique_sites == 1: # Handle single unique site case
		var bounds_poly = PackedVector2Array([ Vector2.ZERO, Vector2(bounds.size.x, 0), bounds.size, Vector2(0, bounds.size.y) ])
		_cells.append({"site": sites_to_process[0], "vertices": bounds_poly, "pre_clip_vertices": bounds_poly.duplicate()})
		_sites = sites_to_process # Update _sites to only contain the unique one
		return true

	# --- Sutherland-Hodgman Iterative Clipping ---
	var temp_cells : Array[Dictionary] = []
	temp_cells.resize(num_unique_sites)

	var initial_bounds_polygon = PackedVector2Array([
		Vector2.ZERO, Vector2(bounds.size.x, 0),
		bounds.size, Vector2(0, bounds.size.y)
	])

	for i in range(num_unique_sites):
		var site_s = sites_to_process[i]
		var current_polygon = initial_bounds_polygon.duplicate()

		for j in range(num_unique_sites):
			if i == j: continue
			var site_o = sites_to_process[j]

			var mid = (site_s + site_o) * 0.5
			var dir_s_to_o = site_o - site_s
			var len_sq = dir_s_to_o.length_squared()
			if len_sq < 1e-9: continue

			var dir_s_to_o_normalized = dir_s_to_o / sqrt(len_sq)
			var normal_to_s = -dir_s_to_o_normalized # Normal points towards site_s

			current_polygon = _clip_polygon_sutherland_hodgman(current_polygon, mid, normal_to_s)

			if current_polygon.size() < 3:
				current_polygon = PackedVector2Array()
				break # Stop early if polygon becomes degenerate

		# Store final result
		temp_cells[i] = {"site": site_s, "vertices": current_polygon, "pre_clip_vertices": current_polygon.duplicate()}


	_cells = temp_cells
	# Update _sites to contain only the unique sites that were processed
	var final_sites : Array[Vector2] = []
	for cell_data in _cells:
		if cell_data and cell_data.has("site"):
			final_sites.append(cell_data.site)
	_sites = final_sites

	return true


## Calculates the geometric centroid of a polygon.
func _calculate_polygon_centroid(polygon: PackedVector2Array, fallback_point: Vector2) -> Vector2:
	var n = polygon.size()
	if n < 3: return fallback_point

	var area : float = 0.0
	var centroid_x : float = 0.0
	var centroid_y : float = 0.0

	for i in range(n):
		var p1 = polygon[i]
		var p2 = polygon[(i + 1) % n]
		var cross_product = p1.cross(p2)
		area += cross_product
		centroid_x += (p1.x + p2.x) * cross_product
		centroid_y += (p1.y + p2.y) * cross_product

	area *= 0.5
	if abs(area) < 1e-9: return fallback_point
	else:
		var factor = 1.0 / (6.0 * area)
		centroid_x *= factor
		centroid_y *= factor
		# Check for NaN/Inf just in case
		if not is_finite(centroid_x) or not is_finite(centroid_y): return fallback_point
		return Vector2(centroid_x, centroid_y)


## Calculates the geometric centroid for each Voronoi cell and clamps the result.
func _calculate_centroids():
	var num_cells = _cells.size()
	if num_cells == 0:
		_centroids.clear()
		return # No cells to process

	_centroids.clear()
	_centroids.resize(num_cells)

	var clamp_min = Vector2.ZERO
	var clamp_max = bounds.size

	for i in range(num_cells):
		var fallback_pos = _cells[i].site if (_cells.get(i) and _cells[i].has("site")) else Vector2.ZERO
		var polygon_for_centroid := PackedVector2Array()

		if _cells.get(i) and _cells[i].has("vertices") and _cells[i].vertices and _cells[i].vertices.size() >= 3:
			polygon_for_centroid = _cells[i].vertices

		var calculated_centroid : Vector2
		if polygon_for_centroid != null:
			calculated_centroid = _calculate_polygon_centroid(polygon_for_centroid, fallback_pos)
		else:
			calculated_centroid = fallback_pos

		var final_centroid = calculated_centroid.clamp(clamp_min, clamp_max)

		if not is_finite(final_centroid.x) or not is_finite(final_centroid.y):
			_centroids[i] = fallback_pos
		else:
			_centroids[i] = final_centroid

	# --- Site Update Handled in _calculate_voronoi post-processing ---
	# The logic at the end of _calculate_voronoi rebuilds _sites based on cells,
	# and the centroids calculated here will be used in the *next* iteration's _calculate_voronoi.
	# We need to update _sites array based on _centroids *before* next iteration.
	if _centroids.size() == num_cells:
		var next_sites : Array[Vector2] = []
		for i in range(num_cells):
			# Use calculated centroid if valid, otherwise keep old site position
			if _centroids[i] != null and is_finite(_centroids[i].x) and is_finite(_centroids[i].y):
				next_sites.append(_centroids[i])
			elif _cells.get(i) and _cells[i].has("site"): # Fallback to site from cell data
				next_sites.append(_cells[i].site)
			# else: # If even cell data is bad, maybe skip or add ZERO?
			 #   print("Warning: Cannot determine next site position for cell index %d" % i)

		# Check if any sites were added before overwriting
		if not next_sites.is_empty():
			_sites = next_sites
		# else: print("Warning: Failed to determine any site positions for next iteration.")


# --- Public Methods ---
func get_sites() -> Array[Vector2]:
	return _sites.duplicate()

func get_cells() -> Array[Dictionary]:
	var cells_copy : Array[Dictionary] = []
	for cell in _cells:
		# Ensure data exists before trying to copy
		if cell and cell.has("site") and cell.has("vertices"):
			cells_copy.append({
				"site": cell.site,
				"vertices": cell.vertices.duplicate() if cell.vertices else PackedVector2Array()
				# No longer need to copy pre_clip_vertices unless debug is needed externally
			})
	return cells_copy
