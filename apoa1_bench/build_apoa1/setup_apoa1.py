def split_pdb():
    import os
    for seg in ['PRO1','PRO2','LIP1','LIP2']:
        os.system(f'convpdb.pl -segnames -readseg -nsel {seg}: apoa1.pdb > apoa1_{seg.lower()}.pdb')
    for i in range(1,10):
        os.system(f'convpdb.pl -segnames -readseg -nsel WAT{str(i)}: -setseg WT0{str(i-1)} -renumber 1 apoa1.pdb > apoa1_wt0{i-1}.pdb')
    return

import pycharmm
from pycharmm import read, write, generate, settings, psf, coor, select_atoms, lingo, minimize, energy
cbl = settings.set_bomb_level(-5)
read.rtf('toppar/top_all36_prot.rtf')
read.rtf('toppar/top_all36_lipid.rtf', append=True)
read.prm('toppar/par_all36m_prot.prm', flex=True)
read.prm('toppar/par_all36_lipid.prm', append=True, flex=True)
read.stream('toppar/toppar_water_ions.str')
settings.set_bomb_level(cbl)

read.sequence_pdb('apoa1_pro1.pdb')
generate.new_segment('PRO1')
read.pdb('apoa1_pro1.pdb', resid=True)

read.sequence_pdb('apoa1_pro2.pdb')
generate.new_segment('PRO2')
read.pdb('apoa1_pro2.pdb', resid=True)

read.sequence_pdb('apoa1_lip1.pdb')
generate.new_segment('LIP1')
generate.rename(to_rename='ATOM',new_name='H9R',
                selection=select_atoms.SelectAtoms(atom_type='H91'))
generate.rename(to_rename='ATOM',new_name='H10R',
                selection=select_atoms.SelectAtoms(atom_type='H101'))
read.pdb('apoa1_lip1.pdb',append=True)

read.sequence_pdb('apoa1_lip2.pdb')
generate.new_segment('LIP2')
generate.rename(to_rename='ATOM',new_name='H9R',
                selection=select_atoms.SelectAtoms(atom_type='H91'))
generate.rename(to_rename='ATOM',new_name='H10R',
                selection=select_atoms.SelectAtoms(atom_type='H101'))
read.pdb('apoa1_lip2.pdb', append=True)
generate.rename(to_rename='ATOM',new_name='H91',
                selection=select_atoms.SelectAtoms(atom_type='H9R'))
generate.rename(to_rename='ATOM',new_name='H101',
                selection=select_atoms.SelectAtoms(atom_type='H10R'))

for i in range(9):
    read.sequence_pdb(f'apoa1_wt0{i}.pdb')
    generate.new_segment(f'WT0{i}', angles=False, dihedrals=False)
    read.pdb(f'apoa1_wt0{i}.pdb', resid=True)


write.psf_card('apoa1_new.psf')
write.coor_card('apoa1_new.crd')

lingo.charmm_script('faster on')
energy.show()
