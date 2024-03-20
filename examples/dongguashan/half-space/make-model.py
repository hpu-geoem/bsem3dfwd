import numpy as np
import scipy.interpolate as si
try:
    import PyDealII.Release as dealii
except:
    print('Unable to import PyDealII.')
    print('Please make sure to add the PyDealII module to your PYTHONPATH by using the following command:')
    print('  export PYTHONPATH=$PYTHONPATH:$(spack location -i dealii)/lib/python3.10/site-packages')
    exit()


# read the model materials from file
model = np.loadtxt('../dongguashan.dat')

# use nearest neighbor interpolation to get the material index
material_interp = si.NearestNDInterpolator(model[:, :3], model[:, 3])

h_steps = [100 for i in range(20)] # horizontal center steps
v_steps = [100 for i in range(15)] # vertical center steps

# add padding cells
for i in range(10):
    h_steps.insert(0, 100 * (1.4**(i+1)))
    h_steps.append(100 * (1.4**(i+1)))
    v_steps.insert(0, 100 * (1.4**(i+1)))
    v_steps.append(100 * (1.4**(i+1)))

p1 = dealii.Point([-sum(h_steps[:20]), -sum(h_steps[:20]),
                  -sum(v_steps[:10])])  # bottom left point
p2 = dealii.Point([sum(h_steps[20:]), sum(h_steps[20:]),
                  sum(v_steps[10:])])  # top right point

# create the mesh
tria = dealii.Triangulation('3D')
tria.generate_subdivided_steps_hyper_rectangle(
    [h_steps, h_steps, v_steps], p1, p2)

# refine the cells cross more than one material
max_level = 3
while True:
    refined = False

    for cell in tria.active_cells():
        c = np.array(cell.center().to_list())

        if c[2] < 0 or c[2] > 2000 or np.abs(c[0]) > 2000 or np.abs(c[1]) > 1500:
            continue

        indices = set()
        for i in range(cell.n_vertices()):
            v = np.array(cell.get_vertex(i).to_list())
            indices.add(int(material_interp(v)))

        if len(indices) > 1 and cell.level() < max_level:
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
        cell.material_id = 3
    elif c[2] > 2000 or np.abs(c[0]) > 2000 or np.abs(c[1]) > 1500:
        cell.material_id = 0
    else:
        cell.material_id = int(material_interp(c))

# save the mesh
tria.save('dongguashan.tria')
tria.write('dongguashan.vtu', 'vtu')

# save the resistivity file
rho = [100, 100, 100, 1E08]
with open('dongguashan.rho', 'w') as rhof:
    print('# attributes', file=rhof)
    print('%d' % (len(rho)), file=rhof)
    for i in range(0, len(rho)):
        print('%E ' % (rho[i]), file=rhof)
    print('', file=rhof)
