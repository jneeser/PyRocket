import gmsh
import sys
import numpy as np 

import config 

# list of x and y points describing the inner and outer chamber wall
y_inner = config.geometry[:,1]
y_outer = y_inner + config.t_w

x_inner = config.geometry[:,0]
x_outer = config.geometry[:,0]

# list of x and y points for left and right vertex
n_points_v = 6                                      # number of points on vertical
dx = 0.001
y_left = np.linspace(y_inner[0], y_outer[0], n_points_v)
y_right = np.linspace(y_inner[-1], y_outer[-1], n_points_v)

x_left = np.ones(len(y_left)) * x_inner[0]
x_right = np.ones(len(y_right)) * x_inner[-1]


# initialise gmsh
gmsh.initialize()
gmsh.model.add("chamber")


# create gmsh point object lists for top and bottom lines (x,y,z)
line_outer = []
line_inner = []
for i in range(len(y_inner)):
    line_outer.append(gmsh.model.occ.addPoint(x_outer[i],y_outer[i],0, dx))
    line_inner.append(gmsh.model.occ.addPoint(x_inner[i],y_inner[i],0, dx))


# points of left and right vertices (x,y,z)
line_left = []
line_right = []
for i in range(len(y_left)):
    line_left.append(gmsh.model.occ.addPoint(x_left[i],y_left[i],0, dx))
    line_right.append(gmsh.model.occ.addPoint(x_right[i],y_right[i],0, dx))


# create lines
outer = gmsh.model.occ.addSpline(line_outer)
inner = gmsh.model.occ.addSpline(line_inner)
left = gmsh.model.occ.addSpline(line_left)
right = gmsh.model.occ.addSpline(line_right)
gmsh.model.occ.synchronize()

#curve = gmsh.model.geo.addCurveLoop([right, outer, -left, -inner])          #[4, 1, -3. -2]
#gmsh.model.mesh.setTransfiniteCurve(curve, 100)
#gmsh.model.geo.addPlaneSurface([1], 1)
# extrude surface 
#top_layer = gmsh.model.geo.extrude(surface, 0, 0, 1, numElements = [1],recombine=True)
curve = gmsh.model.occ.addCurveLoop([right, outer, -left, -inner])
surface = gmsh.model.occ.addPlaneSurface([curve])
gmsh.model.occ.synchronize()
top_layer = gmsh.model.occ.extrude([(2,surface)], 0, 0, dx)

left_wall = gmsh.model.addPhysicalGroup(1,[left])
right_wall = gmsh.model.addPhysicalGroup(1,[right])
outer_wall = gmsh.model.addPhysicalGroup(1,[outer])
inner_wall = gmsh.model.addPhysicalGroup(1,[inner])

gmsh.model.setPhysicalName(1, left_wall, "Left")
gmsh.model.setPhysicalName(1, right_wall, "Right")
gmsh.model.setPhysicalName(1, outer_wall, "Outer")
gmsh.model.setPhysicalName(1, inner_wall, "Inner")


gmsh.model.occ.synchronize()


gmsh.model.mesh.generate(3)

# ... and save it to disk
gmsh.write("chamber.msh")

if '-nopopup' not in sys.argv:

    gmsh.fltk.run()

gmsh.finalize()
