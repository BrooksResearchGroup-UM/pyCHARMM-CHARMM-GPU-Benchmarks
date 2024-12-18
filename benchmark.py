# Cell -1
import argparse

# Define the parser
parser = argparse.ArgumentParser(description='Specify which benchmark to run')

# Declare an argument (`--algo`), saying that the 
# corresponding value should be stored in the `algo` 
# field, and using a default value if the argument 
# isn't given
parser.add_argument('--benchmark', action="store", dest='benchmark', default='5_newDHFR/openmm')

# Now, parse the command line arguments and store the 
# values in the `args` variable
args = parser.parse_args()
benchdir = args.benchmark
print(benchdir)

# Cell 0
def get_GPUInfo(i=None):
    from crimm.Utils import cuda_info
    from crimm.Utils.cuda_info import CUDAInfo
    device_info = CUDAInfo()
    if i == None:
        for i in range(device_info.n_gpus):
            print(device_info.devices[i].name)
    else:
        return device_info.devices[i].name
    return
# Cell 1
import os
import time
import numpy as np
import yaml
#print(os.environ['CHARMM_LIB_DIR'])
#benchdir = '2_T4L/charmm/template'
#benchdir = '2_T4L/bladelib/template'
#benchdir = '4_HSP90/charmm/template'
#benchdir = '4_HSP90/bladelib/template'
#benchdir = '5_newDHFR/openmm'
#benchdir = '5_newDHFR/bladelib'
#benchdir = '6_dmpg/bladelib'
#benchdir = '6_dmpg/openmm'
#benchdir = 'apoa1_bench/openmm'
#benchdir = 'apoa1_bench/bladelib'
#benchdir = 'stmv_bench/openmm'
#benchdir = 'stmv_bench/bladelib'
home = os.getcwd()
# If working in stmv_bench uncompress gzipped files
if 'stmv_bench' in benchdir:
    if 'openmm' in benchdir: os.system(f'gunzip {benchdir}/*.gz')
    elif 'bladelib' in benchdir: os.system(f'gunzip {benchdir}/*.gz')
    elif 'charmm' in benchdir: os.system(f'gunzip {benchdir}/*.gz')
# Read the benchmark data from yaml file
yi = open(f'{benchdir}/bench.yml','r')
input_ = yaml.safe_load(yi)
#print(input_)
yi.close()

# Cell 2
import pycharmm_init
pycharmm_init.chsize = 1500000
from pycharmm import *
import pycharmm

# Cell 3
def setup_nonbond(cutnb=9.0, ctofnb=9.0, ctonnb=9.0,
                  vswit=True, vfswit=False, kappa=0.32,
                  fftx=64,ffty=64,fftz=64):
    script.CommandScript('faster', on=True).run()
    nbond = pycharmm.NonBondedScript(cutnb=cutnb, cutim=cutnb,
                                       ctonnb=ctonnb, ctofnb=ctofnb,
                                       bycb=True,cdie=True, eps=1,
                                       atom=True, vatom=True,
                                       switch=True, vfswit=vfswit, vswit=vswit,
                                       inbfrq=-1, imgfrq=-1,
                                       ewald=True, pmewald=True, kappa=kappa,
                                       fftx=fftx,ffty=ffty,fftz=fftz,order=4)
    
    return nbond

def setup_PBC(boxsize=[], cutoff=9.0, segments=[],
              xtltyp='cubic', blade=False):
    from pycharmm import crystal, image, coor
    import numpy as np
    """input: boxhalf [0.0]
              cutoff [9.0]
              segments  []
              blade [False]
    defines the periodic boundary conditions for a cubic volume of boxsize. 
    Uses: crystal_define_cubic(), crystal.build(), image.setup_residue,
    image.setup_segment to construct symmetry operations. 

    If global variable openmm is true
    the image centering is at [boxhalf,boxhalf,boxhalf] otherwise at [0,0,0].
    """
    if xtltyp == 'ortho': crystal.define_ortho(boxsize[0],boxsize[1],boxsize[2])
    elif xtltyp == 'cubic': crystal.define_cubic(boxsize[0])
    crystal.build(cutoff)

    if len(segments)>0:
        # This is unnessary for GPU benchmarks
        res_sel = pycharmm.SelectAtoms(select_all=True)
        for segid in segments:
            res_sel = res_sel & ~pycharmm.SelectAtoms(seg_id=segid)
        resnames = list(set(pycharmm.SelectAtoms(res_sel).get_res_names()))

        for segment in segments:
            image.setup_segment(0.0, 0.0, 0.0, segment)
        for residue in resnames:
            image.setup_residue(0.0, 0.0, 0.0, residue)
        if not blade:
            pos = coor.get_positions()
            pos.x += boxsize[0]/2
            pos.y += boxsize[1]/2
            pos.z += boxsize[2]/2
            coor.set_positions(pos)

    return

def run_md(useomm=False,useblade=False,nsteps=50000,nprint=10000,time_step=0.002,
           leap=True,lang=True, cpt=False, temp=298.15, pconstant=False, pmass=False,
           pref=False, pgamma=False, hoover=False, tmass=False, nsavl=0, ntrfrq=0):
    
    dynamics.set_fbetas(np.full(psf.get_natom(),1.0,dtype=float))
    # open unit 3 write form name lambda.lam
    lam_file = pycharmm.CharmmFile(file_name='lambda.lmd', 
                                   file_unit=3,
                                   formatted=False,read_only=False)
   
    my_dyn = pycharmm.DynamicsScript(leap=leap, cpt=cpt, lang=lang, start=True,
                                     nstep=nsteps, timest=time_step,
                                     firstt=temp, finalt=temp, tbath=temp,
                                     tstruc=temp, reft=temp,
                                     teminc=0.0, twindh=0.0, twindl=0.0,
                                     pconstant=pconstant, pmass=pmass, pref=pref,
                                     pgamma=pgamma, hoover=hoover, tmass=tmass,
                                     inbfrq=0, imgfrq=0,  # on GPUs these can be zero
                                     iasors=0, iasvel=1, ichecw=0, iscale=0,
                                     iscvel=0, echeck=-1, nsavc=0, nsavv=0, nsavl=nsavl, ntrfrq=ntrfrq,
                                     iprfrq=2*nprint, nprint=nprint, ihtfrq=0, ieqfrq=0,
                                     ilbfrq=0, ihbfrq=0,iunlam=lam_file.file_unit,
                                     omm=useomm, blade=useblade )
    my_dyn.run()
    lam_file.close()

    return

def setup_system(rtf=[],param=[],pflex=False,msld=False,msld_variables='',stream=[],psf='',crd='',pdb=''):
    # set bomlev to -1 incase error in NBFIX stuff
    obl = settings.set_bomb_level(-1)
    if len(rtf)>0:
        append=False
        for i, top in enumerate(rtf):
            if i > 0: append=True
            read.rtf(top,append=append)
    
    if len(param)>0:
        append=False
        for i, par in enumerate(param):
            if i > 0: append=True
            read.prm(par,append=append,flex=pflex)
    if msld: read.stream(msld_variables) 
    if len(stream)>0:
        for f in stream:
            read.stream(f)
    settings.set_bomb_level(obl)
    if len(psf)>0: read.psf_card(psf)
    if len(crd)>0: read.coor_card(crd)
    if len(pdb)>0: read.pdb(pdb,resid=True)
    return

# Cell 4
def hmr(newpsf=''):
    # HMR revisited
    # Get all the masses from the current atoms
    masses = np.array(psf.get_amass())
    resnames = np.array(atom_info.get_res_names(np.arange(0,psf.get_natom(),1)))
    # Build a logical array of all the atoms which are not 'TIP3' waters
    not_waters = resnames != 'TIP3'
    # Build a logical array of all hydrogen atoms based on criterion m_H <= 2
    hydrogens = masses <= 2
    # Build a logical array of all hydrogen atoms not belonging to 'TIP3' waters
    not_water_hydrogen = hydrogens*not_waters
    # Augment each non-water hydrogen mass by 2x original hydrogen mass
    masses += 2*masses*not_water_hydrogen
    # Process the bond array to find heavy atom - hydrogen bonds, 
    # reduce mass of heavy atom by 2*m_H
    bonds = np.array(psf.get_ib_jb())
    for ibnd in range(bonds.shape[1]):
        ib = bonds[0,ibnd]-1
        jb = bonds[1,ibnd]-1
        if not_water_hydrogen[ib] or not_water_hydrogen[jb]:
            if not_water_hydrogen[ib]: masses[jb] -= 2.016
            else: masses[ib] -= 2.016
    # Reset masses to new values
    scalar.set_masses(masses)
    # Write the new psf if requested
    if newpsf != '': write.psf_card(newpsf)

    # Cell 5
os.chdir(benchdir)
psf_file = input_['files']['psf']

msld = 'msld' in input_.keys()
if (input_['dyn']['time_step'] > 0.002 and\
   os.path.isfile(f'{input_["system"]}_hmr.psf')) and\
   not msld:
    psf_file = f'{input_["system"]}_hmr.psf'

if msld:
    msld_variables = input_['msld']['variables_file']
    lingo.set_charmm_variable('FNEX',5.5)
    lingo.set_charmm_variable('TEMP',298.15)
else: msld_variables = ''
param1 = ''
if len(input_['files']['param']) > 1:
    param0 = input_['files']['param'][0]
    param1 = input_['files']['param'][1]
elif len(input_['files']['param']) == 1:
    param0 = input_['files']['param'][0]
else:
    param0 = ''
    param1 = ''
setup_system(rtf=input_['files']['rtf'],
             param=param0,
             pflex=param1,
             msld=msld,
             msld_variables=msld_variables,
             stream=input_['files']['stream'],
             psf=psf_file,
             crd=input_['files']['crd'],
             pdb=input_['files']['pdb'])

if input_['dyn']['time_step'] > 0.002 and\
  not os.path.isfile(f'{input_["system"]}_hmr.psf'):
    hmr(newpsf=f'{input_["system"]}_hmr.psf')

# Cell 6
locals().update(input_['nnb'])
nnb = setup_nonbond(cutnb=cutnb, ctofnb=ctofnb, ctonnb=ctonnb,
                    vswit=vswit, vfswit=vfswit,
                    kappa=kappa, fftx=fftx, ffty=ffty, fftz=fftz)
nnb.run()

# for benchmarks on GPUs using BLaDE and OpenMM APIs image centering not required
# thus setting empty segments and resnames list acheive this
# This is done in stream file for msld benchmarks
if not 'msld' in input_.keys(): 
    setup_PBC(boxsize=input_['pbc']['boxsize'],xtltyp=input_['pbc']['xtltyp'],
          segments=[], blade=(input_['platform'] == 'blade'))

shake.on(param=True,fast=True,tol=1e-7,bonh=True)
if input_['platform'] == 'openmm' and input_['ngpus'] == 2:
    script.CommandScript('omm',deviceid='0,1').run()
elif input_['platform'] == 'domdec':
    lingo.set_charmm_variable('NNODES',input_['ngpus'])    
    script.CommandScript('domdec', gpu='only', dlb='off', ndir='1 @nnodes 1').run()
timing = []
locals().update(input_['dyn'])
if pmass:
    pmass = psf.get_natom() * 0.12
for run in range(input_['benchmark']['niter']):
    timing.append(time.time())
    run_md(useomm=useomm,useblade=useblade,
           nsteps=nsteps,nprint=nprint,
           time_step=time_step,
           leap=leap,lang=lang,temp=temp, pconstant=pconstant, pmass=pmass,
           pref=pref, pgamma=pgamma, hoover=hoover, tmass=tmass, nsavl=nsavl,
           ntrfrq=ntrfrq)
    timing[-1] = time.time() - timing[-1]

timing = 1.0/np.asarray(timing) # nsteps/sec
timing *= (nsteps*time_step)*(60*60*24)/(1000000.0*0.001)
bench = benchdir.split('/')[0]
code = benchdir.split('/')[1]
os.chdir(home)
bout = open('benchmark.csv','a')
print(f'{bench}, {code}, {np.mean(timing):.2f} +/- {np.std(timing):.2f}, {get_GPUInfo(0)}', file=bout)
print(f'Benchmark {bench} on {code}: {np.mean(timing):.2f} +/- {np.std(timing):.2f} ns/day on {get_GPUInfo(0)}')
bout.close()
