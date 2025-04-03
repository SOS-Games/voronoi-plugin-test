# Godot Voronoi & Lloyd's Relaxation Addon

This addon provides a `VoronoiGenerator` Node2D for creating Voronoi diagrams and applying Lloyd's relaxation algorithm within the Godot Engine.

## Features

*   Generates Voronoi diagrams based on a set of points (sites).
*   Applies Lloyd's relaxation algorithm to distribute points more evenly.
*   Generates initial points randomly within specified bounds.
*   Clips Voronoi cells to specified rectangular bounds.
*   Visualizes sites, cell edges, filled cells, and centroids (for debugging).
*   Configurable via the Godot Inspector (`num_points`, `bounds`, `lloyd_iterations`, visualization options).
*   Optionally regenerates automatically when parameters change.
*   Provides `generate()`, `clear()`, `get_sites()`, and `get_cells()` methods for scripting.
*   Emits `generation_finished(sites, cells)` signal.

## Installation

1.  Download or clone this repository/files.
2.  Copy the `addons/voronoi_lloyd` folder into your Godot project's `addons/` directory.
3.  Enable the plugin in Godot: Go to `Project` -> `Project Settings` -> `Plugins` tab, find "VoronoiLloyd", and set Status to "Enabled".

## Usage

1.  Add a `VoronoiGenerator` node to your scene.
2.  Configure its properties in the Inspector:
	*   Set `Bounds` to define the area.
	*   Set `Num Points` for the initial random sites.
	*   Set `Lloyd Iterations` for the number of relaxation steps (0 for none).
	*   Adjust `Seed` for reproducible results (0 for random).
	*   Toggle visualization options (`Draw Sites`, `Draw Edges`, `Fill Cells`, etc.).
	*   Enable/disable `Generate On Change`.
3.  If `Generate On Change` is off, you can call the `generate()` method from a script to trigger the process.
4.  Access the results using `get_sites()` and `get_cells()` or connect to the `generation_finished` signal.

## Limitations

*   **Voronoi Algorithm:** The Voronoi generation uses Godot's `Geometry.triangulate_delaunay_2d` and derives the diagram from its dual. This is simpler to implement in GDScript but may be less robust or performant than dedicated Voronoi algorithms (like Fortune's) especially regarding edge cases (collinear points, points exactly on circumcircles) and handling points near or on the boundary.
*   **Performance:** GDScript calculation can be slow for a very large number of points (`num_points > 1000`). For high performance, consider GDExtension or C#.
*   **Boundary Cells:** Cells for sites near the boundary rely heavily on Godot's `Geometry.clip_polygons_2d`. The construction of the initial (potentially infinite) cell polygon before clipping is simplified in this implementation.

## Example Script

```gdscript
extends Node

@onready var voronoi_gen = $VoronoiGenerator

func _ready():
	# Connect to the signal
	voronoi_gen.generation_finished.connect(_on_generation_finished)

	# Example: Trigger generation manually if 'generate_on_change' is off
	if not voronoi_gen.generate_on_change:
		voronoi_gen.generate()

func _on_generation_finished(sites: Array[Vector2], cells: Array[Dictionary]):
	print("Voronoi generation complete!")
	print("Number of sites: ", sites.size())
	print("Number of cells: ", cells.size())

	# Access cell data (example: print vertices of the first cell)
	if not cells.is_empty():
		print("Vertices of first cell: ", cells[0].vertices)

func _input(event):
	# Example: Regenerate on mouse click
	if event is InputEventMouseButton and event.button_index == MOUSE_BUTTON_LEFT and event.pressed:
		print("Regenerating on click...")
		voronoi_gen.generate()
