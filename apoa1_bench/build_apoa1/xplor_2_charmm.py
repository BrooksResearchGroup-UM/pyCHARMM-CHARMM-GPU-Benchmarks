parameters = {'BOND':[],'ANGLE':[],'DIHEDRAL':[],
              'IMPROPER':[],'NONBONDED':[],'NBFIX':[]}
all_atoms = []
for p_file in ['par_all22_popc.xplor','par_all22_prot_lipid.xplor']:
    pin = open(p_file,'r').readlines()
    for i,l in enumerate(pin):
        if len(l) > 5:
            key = l.split()[0]
            if key in ['BOND','ANGLE','DIHEDRAL',
                                'IMPROPER','NONBONDED','NBFIX']:
                if key == 'NONBONDED': all_atoms.append(l.split()[1])
                parameters[key].append(l.replace(l.split()[0],'').lstrip().strip().replace('UB',''))
            elif key == '{' and l.split()[1] == 'NONBONDED':
                parameters['NONBONDED'].append(l.replace(l.split()[0],'')
                                       .lstrip().strip())
                for j in range(i+1,i+10):
                    if pin[j].strip().split()[-1] == '}':
                        parameters['NONBONDED'].append(pin[j].replace(pin[j].split()[-1],'')
                                               .strip())
                        break
                    else:
                        parameters[key].append(pin[j].strip())
                        

po = open('par_all22_prot_lipid.inp','w')
for i in range(50):
    if pin[i][0] == '{':
        print(f'* {pin[i].strip()}',file=po)
print('*   \n\n',file=po)
for key in parameters.keys():
    if key != 'NONBONDED':
        print(f'\n\n{key}',file=po)
    else: print('\n\n',file=po)
    for item in parameters[key]:
        if key == 'NBFIX':
            terms = item.split()
            a = float(terms[2])
            b = float(terms[3])
            ap = float(terms[4])
            bp = float(terms[5])
            print(f'{terms[0]}  {terms[1]}  {-b*b/a:.4f} {(a/b)**(1/6):.4f}  {-bp*bp/ap:.4f}  {(ap/bp)**(1/6):.4f}',file=po)
        else: print(f'{item}',file=po)
po.close()
        
all_atoms = list(set(all_atoms))
psf = open('apoa1.psf','r').readlines()
atoms = {}
aread = False
for i,l in enumerate(psf):
    if 'NATOM' in l: aread = True
    ln = l.strip().split()
    if 'bonds' in l: break
    if aread:
        if len(ln) <= 7: continue
        elif 'bonds' in ln[2].lower(): break
        else:
            atoms[ln[5]] = ln[7]
for a in all_atoms:
    if not a in atoms.keys(): atoms[a] = 20.0000
rtf = open('apoa1.rtf','w')
print('* title\n*  \n   22 1\n\n', file=rtf)
for a in atoms.keys():
    print(f'mass -1 {a} {atoms[a]}',file=rtf)
print('\n',file=rtf)
rtf.close()
