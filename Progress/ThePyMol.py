import pymol
#Run using pymol -c ThePyMol.py

cmd.wizard('mutagenesis')
cmd.load('Backbone.pdb')
cmd.refresh_wizard()
for res in range(138):
	cmd.get_wizard().do_select('/Backbone//A/GLY`{}'.format(res))
	cmd.get_wizard().set_mode('VAL')
	cmd.get_wizard().apply()
cmd.set_wizard()
cmd.save('123.pdb')
