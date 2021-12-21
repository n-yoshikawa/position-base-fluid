import bpy
import math

def look_at(obj_camera, point):
    loc_camera = obj_camera.matrix_world.to_translation()

    direction = point - loc_camera
    # point the cameras '-Z' and use its 'Y' as up
    rot_quat = direction.to_track_quat('-Z', 'Y')

    # assume we're using euler rotation
    obj_camera.rotation_euler = rot_quat.to_euler()

# Create a new material
def add_water_material():
    material = bpy.data.materials.new(name="Plane Light Emission Shader")
    material.use_nodes = True

    # Remove default
    material.node_tree.nodes.remove(material.node_tree.nodes.get('Principled BSDF'))
    material_output = material.node_tree.nodes.get('Material Output')
    
    # Add Principaled BSDF
    bsdf = material.node_tree.nodes.new('ShaderNodeBsdfPrincipled')
    bsdf.inputs['Base Color'].default_value = (0.102,0.196,0.799,1)
    bsdf.inputs['Metallic'].default_value = 0.2
    bsdf.inputs['Specular'].default_value = 0.5
    bsdf.inputs['Sheen Tint'].default_value = 0.5
    bsdf.inputs['Transmission'].default_value = 1.0
    bsdf.inputs['Metallic'].default_value = 0.2
    # link emission shader to material
    material.node_tree.links.new(material_output.inputs[0], bsdf.outputs[0])

    # Add Bump
    bump = material.node_tree.nodes.new('ShaderNodeBump')
    bump.inputs['Strength'].default_value = 0.2
    material.node_tree.links.new(bsdf.inputs['Normal'],bump.outputs['Normal'])

    # Musgrave
    musgrave = material.node_tree.nodes.new('ShaderNodeTexMusgrave')
    musgrave.musgrave_dimensions = '3D'
    musgrave.musgrave_type = 'HYBRID_MULTIFRACTAL'
    musgrave.inputs['Scale'].default_value = 6.8
    musgrave.inputs['Detail'].default_value = 10.6
    musgrave.inputs['Dimension'].default_value = 1.2
    musgrave.inputs['Lacunarity'].default_value = 1.72
    musgrave.inputs['Offset'].default_value = 0
    musgrave.inputs['Gain'].default_value = 1
    # link bump to musgrave
    print( musgrave.outputs.keys())
    material.node_tree.links.new(bump.inputs['Height'], musgrave.outputs['Fac'])

    # Mapping
    mapping = material.node_tree.nodes.new('ShaderNodeMapping')
    # link mapping to musgrave
    material.node_tree.links.new(musgrave.inputs['Vector'], mapping.outputs['Vector'])
    
    # Texture Coordinate
    tex = material.node_tree.nodes.new('ShaderNodeTexCoord')
    # link mapping to musgrave
    material.node_tree.links.new(musgrave.inputs['Vector'], tex.outputs['Object'])
    return material

def add_transparent_water_material():
    material = bpy.data.materials.new(name="Transparent Water Shader")
    material.use_nodes = True
    # Remove default
    material.node_tree.nodes.remove(material.node_tree.nodes.get('Principled BSDF'))
    material_output = material.node_tree.nodes.get('Material Output')

    # Add Mixer
    mixer = material.node_tree.nodes.new('ShaderNodeMixShader')
    mixer.inputs['Fac'].default_value = 0.059
    # Link Mixer
    material.node_tree.links.new(material_output.inputs[0], mixer.outputs[0])
    
    # Add Glass
    glass = material.node_tree.nodes.new('ShaderNodeBsdfGlass')
    glass.inputs['Roughness'].default_value = 0
    glass.inputs['IOR'].default_value = 1.45
    glass.inputs['Color'].default_value = (0.540872,0.582622,1.0,1.0)
    # Link Glass
    material.node_tree.links.new(mixer.inputs[1], glass.outputs[0])
    
    # Add Transparent
    transparent = material.node_tree.nodes.new('ShaderNodeBsdfTransparent')
    # Link Transparent
    material.node_tree.links.new(mixer.inputs[2], transparent.outputs[0])
    
    # Set Blend Mode
    material.blend_method = "BLEND"
    return material
    
def add_hdr_to_world():
    C = bpy.context
    scn = C.scene

    # Get the environment node tree of the current scene
    node_tree = scn.world.node_tree
    tree_nodes = node_tree.nodes

    # Clear all nodes
    tree_nodes.clear()

    # Add Background node
    node_background = tree_nodes.new(type='ShaderNodeBackground')

    # Add Environment Texture node
    node_environment = tree_nodes.new('ShaderNodeTexEnvironment')
    # Load and assign the image to the node property
    node_environment.image = bpy.data.images.load("Applications/Blender.app/Contents/Resources/3.0/datafiles/studiolights/world/interior.exr") # Relative path
    node_environment.location = -300,0

    # Add Output node
    node_output = tree_nodes.new(type='ShaderNodeOutputWorld')   
    node_output.location = 200,0

    # Link all nodes
    links = node_tree.links
    link = links.new(node_environment.outputs["Color"], node_background.inputs["Color"])
    link = links.new(node_background.outputs["Background"], node_output.inputs["Surface"])


# Render Out Image
def render(i):
    bpy.context.scene.render.filepath = '/Users/Marta/csc2549/position-base-/fluid_sim_output/hdr_particle_' + str(i) + '.jpg'
    bpy.context.scene.render.resolution_x = 1920 #perhaps set resolution in code
    bpy.context.scene.render.resolution_y = 1080
    bpy.ops.render.render(write_still=True)

waterMaterial = add_transparent_water_material()

file_loc = '/Users/Marta/csc2549/position-base-fluid/marching_cubes_bin20_8000p_padded/bin20_step'
bpy.context.scene.camera.location.y = -49.54
for i in range(0,999):
    obj = bpy.ops.import_scene.obj(filepath=file_loc+str(i)+".obj")
    obj_object = bpy.context.selected_objects[-1]
    obj_object.active_material = waterMaterial
    # set the shading of all polygons to flat 
    for f in obj_object.data.polygons:
        f.use_smooth = True
    obj_object.location = (-8,8,-8)
    bpy.context.scene.camera.location.y += 0.09
    render(i)
    bpy.ops.object.delete()
