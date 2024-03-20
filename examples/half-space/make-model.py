try:
    import PyDealII.Release as dealii
except:
    print('Unable to import PyDealII.')
    print('Please make sure to add the PyDealII module to your PYTHONPATH by using the following command:')
    print('  export PYTHONPATH=$PYTHONPATH:$(spack location -i dealii)/lib/python3.10/site-packages')
    exit()


h_steps = [128 for i in range(10)] # horizontal center steps
v_steps = [100 for i in range(5)] # vertical center steps

# add padding cells
for i in range(7):
    h_steps.insert(0, 128 * (1.3**(i+1)))
    h_steps.append(128 * (1.3**(i+1)))
    v_steps.insert(0, 100 * (1.3**(i+1)))
    v_steps.append(100 * (1.3**(i+1)))

p1 = dealii.Point([-sum(h_steps[:12]), -sum(h_steps[:12]), -sum(v_steps[:7])])  # bottom left point
p2 = dealii.Point([sum(h_steps[12:]), sum(h_steps[12:]), sum(v_steps[7:])])  # top right point

# create the mesh
tria = dealii.Triangulation('3D')
tria.generate_subdivided_steps_hyper_rectangle([h_steps, h_steps, v_steps], p1, p2)

# assign material id for each cell
for cell in tria.active_cells():
    c = cell.center().to_list()
    if c[2] < 0:
        cell.material_id = 0
    else:
        cell.material_id = 1

# save the mesh
tria.save('half-space.tria')
tria.write('half-space.vtu', 'vtu')

# save the resistivity file
rho = [1E08, 100]
with open('half-space.rho', 'w') as rhof:
    print('# attributes', file=rhof)
    print('%d' % (len(rho)), file=rhof)
    for i in range(0, len(rho)):
        print('%E ' % (rho[i]), file=rhof)
    print('', file=rhof)
