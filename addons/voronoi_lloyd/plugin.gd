@tool
extends EditorPlugin

const VoronoiGenerator = preload("res://addons/voronoi_lloyd/voronoi_generator.gd")
# Optional: uncomment and provide path if you created an icon
# const VoronoiGeneratorIcon = preload("res://addons/voronoi_lloyd/icons/VoronoiGenerator.svg")

func _enter_tree():
	# Register the custom node
	# Optional: add ", VoronoiGeneratorIcon" after "Node2D" if you have an icon
	add_custom_type("VoronoiGenerator", "Node2D", VoronoiGenerator, null)


func _exit_tree():
	# Unregister the custom node
	remove_custom_type("VoronoiGenerator")
