from openmm.app import *
from openmm import *
from openmm.unit import *
#from sys import stdout
from simtk import unit as u


def cala_minimize_energy(pdb_file, log_file, steps):
    pdb = PDBFile(pdb_file)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(forcefield, 
                        model='tip3p', 
                        padding=10*angstroms)
    
    system = forcefield.createSystem(modeller.topology, 
                                    nonbondedMethod=PME,
                                    nonbondedCutoff=10*angstrom, 
                                    constraints=HBonds)

    integrator = LangevinMiddleIntegrator(300*kelvin, 
                                        1/u.picosecond, 
                                        0.004*u.picoseconds)

    platform = Platform.getPlatformByName('OpenCL')
    platformProperties = {'Precision': 'double'}

    simulation = Simulation(modeller.topology, system, integrator, platform, platformProperties)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    simulation.reporters.append(PDBReporter(output_file, 1000))
    simulation.reporters.append(StateDataReporter(log_file, 
                                                1000, 
                                                #step=True,
                                                totalEnergy=True,
                                                kineticEnergy=True,
                                                potentialEnergy=True, 
                                                #temperature=True,
                                                ))

    simulation.step(steps)

if __name__ == '__main__':
    pdb_file = '/mnt/sdc/lanwei/TLR2/TLR2.pdb'
    log_file = '/mnt/sdc/lanwei/TLR2/energy.txt'
    output_file = '/mnt/sdc/lanwei/TLR2/output.pdb'
    steps = 1000
    cala_minimize_energy('output.pdb', 'log.txt')