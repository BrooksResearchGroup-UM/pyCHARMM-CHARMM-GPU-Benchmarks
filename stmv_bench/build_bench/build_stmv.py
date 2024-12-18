def get_coordinates():
    '''1) Read the psf to extract the segnames.
       2) Loop through all segnames and use convpdb.pl to extract ATOM
          records for that segment and renumber all residues from 1.
       3) Note process did not work for the segments with only a single
          letter, see function get_I_and_N.'''
    import os
    psf = open('stmv.psf','r').readlines()
    segs = []
    for l in psf:
        ls = l.strip()
        if len(ls)<1: continue
        ls = ls.split()
        if len(ls)<2: continue
        if ls[1] == '!NATOM' : readit = True
        elif ls[1] == '!NBOND:' : break
        elif len(ls) == 9: segs.append(ls[1])
    segs = sorted(list(set(segs)))
    for s in segs:
        command = f'convpdb.pl -segnames -readseg -nsel {s}: -renumber 1 stmv.pdb > pdb/stmv_{s}.pdb'
        print(command)
        os.system(command)
    return

def get_I_and_N():
    '''Extract ATOM records for segments N (nucleic acid) and I (Mg2+ ions)'''
    cn = open('pdb/stmv_N.pdb','w')
    ci = open('pdb/stmv_I.pdb','w')
    coor = open('stmv.pdb','r').readlines()
    for l in coor:
        ls = l.strip()
        if len(ls)<1: continue
        ls = ls.split()
        if len(ls)<11: continue
        if ls[0] == 'ATOM' and (ls[-1]=='N'):
            print(l.strip(),file=cn)
        elif ls[0] == 'ATOM' and ls[-1]=='I':
            print(l.strip(),file=ci)
    print('END',file=ci)
    print('END',file=cn)
    return
# Called above functions to construct the coordinates for each segment
# Use pycharmm to build system coordinates
def build_complete():
    import pycharmm
    from pycharmm import read, generate, settings, write, psf, lingo, energy
    read.rtf('toppar/top_all36_prot.rtf')
    read.rtf('toppar/top_all36_na.rtf', append=True)
    read.prm('toppar/par_all36m_prot.prm', flex=True)
    read.prm('toppar/par_all36_na.prm', flex=True, append=True)
    read.stream('toppar/toppar_water_ions.str')
    read.psf_card('stmv_virus.psf')
    read.coor_card('stmv_virus.crd')
    read.psf_card('stmv_water.psf', append=True)
    read.coor_card('stmv_water.crd', append=True)
    write.psf_card('stmv_system.psf')
    write.coor_card('stmv_system.crd')
    exit()
    lingo.charmm_script('faster on')
    # Check energy
    energy.show()
    exit()
    return

def build_segments():
    import pycharmm
    from pycharmm import read, generate, settings, write, psf, lingo, energy
    read.rtf('toppar/top_all36_prot.rtf')
    read.rtf('toppar/top_all36_na.rtf', append=True)
    read.prm('toppar/par_all36m_prot.prm', flex=True)
    read.prm('toppar/par_all36_na.prm', flex=True, append=True)
    read.stream('toppar/toppar_water_ions.str')
    # Build protein coat
    for i in range(60):
        read.sequence_pdb(f'pdb/stmv_C{i}.pdb')
        generate.new_segment(f'C{i}')
        read.pdb(f'pdb/stmv_C{i}.pdb', resid=True)

    # Now add RNA
    read.sequence_pdb('pdb/stmv_N.pdb')
    generate.new_segment('N',first_patch='5PHO',last_patch='3TER')
    read.pdb('pdb/stmv_N.pdb', resid=True)

    # Now add the ions - Mg2+
    read.sequence_pdb('pdb/stmv_I.pdb')
    generate.new_segment('I')
    read.pdb('pdb/stmv_I.pdb', resid=True)

    # Now add the ions - Cl-
    read.sequence_pdb('pdb/stmv_CI.pdb')
    generate.new_segment('CI')
    read.pdb('pdb/stmv_CI.pdb', resid=True)

    write.psf_card('stmv_virus.psf')
    write.coor_card('stmv_virus.crd')
    psf.delete_atoms()
    # Build water
    for i in range(1,65):
        # W
        read.sequence_pdb(f'pdb/stmv_W{i}.pdb')
        generate.new_segment(f'W{i}', angles=False, dihedrals=False)
        read.pdb(f'pdb/stmv_W{i}.pdb', resid=True)
    # Build water
    for i in range(1,39):
        # Z
        read.sequence_pdb(f'pdb/stmv_Z{i}.pdb')
        generate.new_segment(f'Z{i}', angles=False, dihedrals=False)
        read.pdb(f'pdb/stmv_Z{i}.pdb', resid=True)
    for i in range(40,65):
        # Z
        read.sequence_pdb(f'pdb/stmv_Z{i}.pdb')
        generate.new_segment(f'Z{i}', angles=False, dihedrals=False)
        read.pdb(f'pdb/stmv_Z{i}.pdb', resid=True)
    write.psf_card('stmv_water.psf')
    write.coor_card('stmv_water.crd')
    return
