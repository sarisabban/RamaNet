import os
from pyrosetta import *
from pyrosetta.toolbox import *
init()

pose = pose_from_pdb('Backbone.pdb')

mover = pyrosetta.rosetta.protocols.simple_moves.AddConstraintsToCurrentConformationMover()
mover.CA_only()
mover.generate_constraints(pose)
mover.apply(pose)

scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('talaris2013_cart')
scorefxn.set_weight(rosetta.core.scoring.coordinate_constraint , 1.0)

#cart_bonded energy term to the scorefunction in question, as well as remove the pro_close term (to avoid double counting). 

#scorefxn.set_weight(rosetta.core.scoring.cart_bonded_angle , 1.0)
#scorefxn.set_weight(rosetta.core.scoring.cart_bonded_improper , 1.0)
#scorefxn.set_weight(rosetta.core.scoring.cart_bonded_length , 1.0)
#scorefxn.set_weight(rosetta.core.scoring.cart_bonded_proper , 1.0)
#scorefxn.set_weight(rosetta.core.scoring.cart_bonded_ring , 1.0)
#scorefxn.set_weight(rosetta.core.scoring.cart_bonded_torsion , 1.0)

relax = pyrosetta.rosetta.protocols.relax.FastRelax()
relax.set_scorefxn(scorefxn)
relax.cartesian(True)
relax.ramp_down_constraints(False)
relax.apply(pose)

pose.dump_pdb('test.pdb')
os.system('pymol test.pdb')
os.system('rm test.pdb')
