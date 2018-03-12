from pyrosetta import *
from pyrosetta.toolbox import *
init()

pose = pose_from_pdb('test.pdb')
scorefxn = get_fa_scorefxn()

filters = rosetta.protocols.simple_filters.PackStatFilter()
read = pyrosetta.rosetta.core.pack.task.operation.ReadResfile('Resfile.resfile')
task = pyrosetta.rosetta.core.pack.task.TaskFactory()
task.push_back(read)

movemap = MoveMap()
movemap.set_bb(False)
movemap.set_chi(True)

mover = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()
mover.set_task_factory(task)
mover.set_movemap(movemap)
mover.set_scorefxn(scorefxn)

MC = pyrosetta.rosetta.protocols.simple_moves.GenericMonteCarloMover()
MC.set_mover(mover)
MC.set_scorefxn(scorefxn) # <--- Problem Here
MC.set_maxtrials(1)
MC.set_temperature(1)
MC.set_preapply(True)
MC.set_drift(True)
MC.set_sampletype('high')
MC.add_filter(filters , False , 1.0 , 'high' , True) # <--- Problem Here
print('\n\n++++++++++++++++++++++++++++++\n' , pose , '\n++++++++++++++++++++++++++++++\n\n')
MC.apply(pose)
print('\n\n++++++++++++++++++++++++++++++\n' , pose , '\n++++++++++++++++++++++++++++++\n\n')
