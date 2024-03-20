import numpy as np
try:
    import PyDealII.Release as dealii
except:
    print('Unable to import PyDealII.')
    print('Please make sure to add the PyDealII module to your PYTHONPATH by using the following command:')
    print('  export PYTHONPATH=$PYTHONPATH:$(spack location -i dealii)/lib/python3.10/site-packages')
    exit()


# get the material id for each point
def get_cell_indices(vertices, triangles, points):
    indices = []
    for p in points:
        for tri in triangles:
            v0 = vertices[int(tri[0]) - 1]
            v1 = vertices[int(tri[1]) - 1]
            v2 = vertices[int(tri[2]) - 1]

            attr = int(tri[3])

            A1 = (v1[0] - v0[0]) * (p[1] - v0[1]) - (v1[1] - v0[1]) * (p[0] - v0[0])
            A2 = (v2[0] - v1[0]) * (p[1] - v1[1]) - (v2[1] - v1[1]) * (p[0] - v1[0])
            A3 = (v0[0] - v2[0]) * (p[1] - v2[1]) - (v0[1] - v2[1]) * (p[0] - v2[0])

            # whether the point is inside the triangle
            if (A1 < 0 and A2 < 0 and A3 < 0) or (A1 > 0 and A2 > 0 and A3 > 0):
                indices.append(attr)
                break

    return indices


# load triangle mesh from file
vertices = np.loadtxt('../oil-gas-model.node', skiprows=1, usecols=(1, 2))
triangles = np.loadtxt('../oil-gas-model.ele', skiprows=1, usecols=(1, 2, 3, 4))

h_steps = [1000 for i in range(10)]  # horizontal center steps
v_steps = [1000 for i in range(6)]  # vertical center steps

# add padding cells
for i in range(5):
    h_steps.insert(0, 1000 * (1.3**(i+1)))
    h_steps.append(1000 * (1.3**(i+1)))
    v_steps.insert(0, 1000 * (1.3**(i+1)))
    v_steps.append(1000 * (1.3**(i+1)))

p1 = dealii.Point([-sum(h_steps[:10]), -sum(h_steps[:10]), -sum(v_steps[:5])])  # bottom left point
p2 = dealii.Point([sum(h_steps[10:]), sum(h_steps[10:]), sum(v_steps[5:])])  # top right point

# create the mesh
tria = dealii.Triangulation('3D')
tria.generate_subdivided_steps_hyper_rectangle([h_steps, h_steps, v_steps], p1, p2)

# refine the mesh crosses more than one material iteratively
max_level = 5
while True:
    refined = False

    for cell in tria.active_cells():
        c = np.array(cell.center().to_list())

        if not (np.abs(c[0]) < 1000 and np.abs(c[1]) < 10000 and c[2] > 0 and c[2] < 3000):
            continue

        points = []
        for i in range(cell.n_vertices()):
            v = np.array(cell.get_vertex(i).to_list())
            p = v * 0.99 + c * 0.01
            points.append(p[1:])
        points.append(c[1:])

        indices = get_cell_indices(vertices, triangles, points)

        if len(set(indices)) > 1 and cell.level() < max_level:
            cell.refine_flag = 'isotropic'
            refined = True

    if not refined:
        break

    tria.execute_coarsening_and_refinement()

    print(tria.n_active_cells())

# assign material id for each cell
for cell in tria.active_cells():
    c = np.array(cell.center().to_list())

    if c[2] < 0:
        cell.material_id = 0 # air
    elif c[2] < 6000:
        cell.material_id = 1 # layer 1
    else:
        cell.material_id = 2 # layer 2

    if np.abs(c[0]) < 1000 and np.abs(c[1]) < 10000 and c[2] > 0 and c[2] < 3000:
        points = [c[1:]]
        cell.material_id = get_cell_indices(vertices, triangles, points)[0]

# save the mesh
tria.save('oil-gas.tria')
tria.write('oil-gas.vtu', 'vtu')

# save the resistivity file
rho = [1E08, 100, 1000, 500, 5, 5]
with open('oil-gas.rho', 'w') as rhof:
    print('# attributes', file=rhof)
    print('%d' % (len(rho)), file=rhof)
    for i in range(0, len(rho)):
        print('%E ' % (rho[i]), file=rhof)
    print('', file=rhof)
