import os
from pyrosetta import *
from pyrosetta.toolbox import *
init()

pose = pose_from_pdb('Backbone.pdb')

mover = pyrosetta.rosetta.protocols.simple_moves.AddConstraintsToCurrentConformationMover()
mover.CA_only()
mover.generate_constraints(pose)
mover.apply(pose)

scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('fldsgn_cen')
scorefxn.set_weight(rosetta.core.scoring.atom_pair_constraint , 1.0)
#scorefxn.set_weight(rosetta.core.scoring.cart_bonded_angle , 1.0)
#scorefxn.set_weight(rosetta.core.scoring.cart_bonded_improper , 1.0)
#scorefxn.set_weight(rosetta.core.scoring.cart_bonded_length , 1.0)
#scorefxn.set_weight(rosetta.core.scoring.cart_bonded_proper , 1.0)
#scorefxn.set_weight(rosetta.core.scoring.cart_bonded_ring , 1.0)
#scorefxn.set_weight(rosetta.core.scoring.cart_bonded_torsion , 1.0)

relax = pyrosetta.rosetta.protocols.relax.FastRelax()
relax.set_scorefxn(scorefxn)
relax.cartesian()
relax.ramp_down_constraints(False)
relax.apply(pose)


'''
Resfile = open('Resfile.resfile' , 'w')
Resfile.write('NATAA\n')
Resfile.write('start\n')
for line in range(len(pose) + 3):
	Resfile.write(str(line + 1) + ' ' + 'A' + ' V\n')
Resfile.close()
read = pyrosetta.rosetta.core.pack.task.operation.ReadResfile('Resfile.resfile')
task = pyrosetta.rosetta.core.pack.task.TaskFactory()
task.push_back(read)
movemap = MoveMap()
movemap.set_bb(True)
movemap.set_chi(False)

#mover = pyrosetta.rosetta.protocols.forge.remodel.RemodelMover()
#mover.apply(pose)

#mover = pyrosetta.rosetta.protocols.forge.remodel.RemodelDesignMover()
#mover.task()
#mover.show()
#mover.apply(pose)
'''

pose.dump_pdb('test.pdb')
os.system('pymol test.pdb')
os.system('rm test.pdb')
