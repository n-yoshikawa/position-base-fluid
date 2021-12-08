import bpy
import numpy as np

with open('positions.npy', 'rb') as f:
    positions = np.load(f)

N = positions.shape[1]

for i, p in enumerate(positions[0]):
    bpy.ops.mesh.primitive_uv_sphere_add(location=(p[0], p[1], p[2]), radius=0.02)
    bpy.context.active_object.name = 'obj{}'.format(i)
    mat = bpy.data.materials.new(name= 'colored')
    mat.diffuse_color = (0., 0., 1., 1)
    bpy.context.active_object.active_material = mat
    bpy.context.object.active_material.shadow_method = 'NONE'
    bpy.ops.object.shade_smooth()

for frame in range(positions.shape[0]):
    for i, p in enumerate(positions[frame]):
        obj = bpy.data.objects['obj{}'.format(i)]
        obj.location = (p[0], p[1], p[2])
        obj.keyframe_insert(data_path = "location", index = -1)
    bpy.context.scene.frame_set(frame)
