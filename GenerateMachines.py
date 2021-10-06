# Generation of cavity mode configurations and lattices according to configuration files

import configparser, re, random, os, shutil, math, numpy as np

random.seed()

Ee = 0.5109989461e6 # Electron rest mass, eV
c = 2.99792458e8    # Speed of light, m/s

# Read configuration files
machine = configparser.ConfigParser()
machine.read("machine.ini")
config = configparser.ConfigParser()
config.read(machine['General']['Machine'].strip('"'))

# Create necessary directories
if not os.path.exists('Temp'):
    os.mkdir('Temp')
if int(config['Spread']['PreserveCavities']) == 0:
    if os.path.exists('machines/cavities'):
        shutil.rmtree('machines/cavities')
if os.path.exists('machines/lattices'):
    shutil.rmtree('machines/lattices')
if os.path.exists('machines/matrices'):
    shutil.rmtree('machines/matrices')
if os.path.exists('machines/transition'):
    shutil.rmtree('machines/transition')
if int(config['Spread']['PreserveCavities']) == 0:
    os.makedirs('machines/cavities')
os.makedirs('machines/lattices')
os.makedirs('machines/matrices')
os.makedirs('machines/transition')

# Read HOM parameters
Ncav = int(config['General']['NumberOfCavities'])
Nmodes = int(config['General']['NumberOfModes'])
Npass = int(config['General']['NumberOfCavityPasses']) if int(config['Run']['UsedCavityPasses']) < 0 else int(config['Run']['UsedCavityPasses'])
stream = open(config['General']['CavityParameters'].strip('"'), 'r')
t0 = 1/float(config['General']['FundamentalFrequency'])
n = [0] * Nmodes        # HOM index (in CST Studio output)
ord = [0] * Nmodes      # HOM azimuthal order
f = [0.] * Nmodes       # HOM frequency (Hz)
Q = [0.] * Nmodes       # HOM Q-factor
RoQd = [0.] * Nmodes    # HOM (R/Q)[dipole] (Ohm/m^2)
theta = [0.] * Nmodes   # HOM polarization angle (deg)
df = [0.] * Nmodes      # HOM frequency spread (Hz)
dQQ = [0.] * Nmodes     # HOM Q-factor relative spread
for i in range(Nmodes):
    (n[i], ord[i], f[i], Q[i], RoQd[i], theta[i], df[i], dQQ[i]) = [float(j) for j in stream.readline().split()]
stream.close()
MaxMachineIndex = int(config['Spread']['NumberOfMachines']) if int(config['Run']['UseSweep']) == 0 else int(config['Spread']['NumberOfMachines'])*int(config['Run']['NstepsX'])*int(config['Run']['NstepsY'])

# Randomize HOM parameters
for MachineIndex in range(MaxMachineIndex):
    if (int(config['Spread']['SameCavities']) != 0) and (MachineIndex > 0):
        shutil.copy('machines/cavities/0000.txt', 'machines/cavities/%.4d.txt' % MachineIndex)
    else:
        if MachineIndex < int(config['Spread']['NumberOfMachines']):
            if int(config['Spread']['PreserveCavities']) == 0:
                stream = open('machines/cavities/%.4d.txt' % MachineIndex, 'w')
                for j in range(Ncav):
                    for k in range(Nmodes):
                        if (ord[k] > 0) and (int(config['Spread']['RandomPolarization']) != 0):
                            theta[k] = 90/ord[k]*(2*random.random()-1) if k%2==0 else theta[k-1] + 90/ord[k]
                        df[k] = 1e9
                        while abs(df[k]) > float(config['Spread']['FrequencySpread'])*float(config['Spread']['DistrCutoffSigma']):
                            df[k] = np.random.normal(0., float(config['Spread']['FrequencySpread']))
                        dQQ[k] = 10
                        while abs(dQQ[k]) > 1+float(config['Spread']['DistrCutoffSigma'])/math.sqrt(2.):
                            dQQ[k] = np.random.gamma(2., 0.5)
                        stream.write("%d\t%d\t%e\t%e\t%e\t%e\n" % (n[k], ord[k], f[k]+df[k], Q[k]*dQQ[k], RoQd[k], theta[k]))
                stream.close()
        else:
            shutil.copy('machines/cavities/%.4d.txt' % (MachineIndex % int(config['Spread']['NumberOfMachines'])), 'machines/cavities/%.4d.txt' % MachineIndex)

# Create input file 'run.ele' for elegant
stream = open('Temp/run.ele', 'w')
stream.write('&run_setup\nlattice=\"Temp/lattice.lte\"\ndefault_order=1\nuse_beamline=%s\np_central=%e\nprint_statistics=0\n&end\n\n' %
    (config['General']['Beamline'], float(config['General']['p_central'])))
stream.write('&run_control\n&end\n\n')
stream.write('&twiss_output\nfilename=%%s.twi\nmatched=0\nbeta_x=%e\nalpha_x=%e\nbeta_y=%e\nalpha_y=%e\n&end\n\n' %
    (float(config['General']['beta_x']), float(config['General']['alpha_x']), float(config['General']['beta_y']), float(config['General']['alpha_y'])))
stream.write('&matrix_output\nSDDS_output=%s.mat\n&end\n\n')
stream.write('&bunched_beam\nemit_x=%e\nemit_y=%e\nsigma_s=%e\nsigma_dp=%e\nuse_twiss_command_values=1\n&end\n\n' %
    (float(config['General']['emit_x']), float(config['General']['emit_y']), float(config['General']['sigma_s']), float(config['General']['sigma_dp'])))
stream.write('&track\n&end\n')
stream.close()

# Initialize sigma-matrices
S0 = np.array([[0.,0.,0.,0.,0.,0.] for i in range(6)])
S0[0][0] = float(config['General']['emit_x'])*float(config['General']['beta_x'])
S0[0][1] = -float(config['General']['emit_x'])*float(config['General']['alpha_x'])
S0[1][0] = S0[0][1]
S0[1][1] = float(config['General']['emit_x'])*(1+float(config['General']['alpha_x'])**2)/float(config['General']['beta_x'])
S0[2][2] = float(config['General']['emit_y'])*float(config['General']['beta_y'])
S0[2][3] = -float(config['General']['emit_y'])*float(config['General']['alpha_y'])
S0[3][2] = S0[2][3]
S0[3][3] = float(config['General']['emit_y'])*(1+float(config['General']['alpha_y'])**2)/float(config['General']['beta_y'])
S0[4][4] = float(config['General']['sigma_s'])**2
S0[5][5] = float(config['General']['sigma_dp'])**2
S0 = np.array(S0)
S1 = np.array([[0.,0.,0.,0.,0.,0.] for i in range(6)])

for MachineIndex in range(MaxMachineIndex):
    if (int(config['Spread']['SameLattice']) != 0) and (MachineIndex > 0):
        # Copy the first lattice
        shutil.copy('machines/lattices/0000.txt', 'machines/lattices/%.4d.txt' % MachineIndex)
        shutil.copy('machines/matrices/0000.txt', 'machines/matrices/%.4d.txt' % MachineIndex)
        shutil.copy('machines/transition/0000.txt', 'machines/transition/%.4d.txt' % MachineIndex)
    else:
        if MachineIndex < int(config['Spread']['NumberOfMachines']):
            # Generate lattice files with randomized quadrupole gradient and roll
            stream = open(config['General']['LatticeFile'].strip('"'), 'r')
            streamout = open('machines/lattices/%.4d.lte' % MachineIndex, 'w')
            GradError = 0 if MachineIndex == 0 and int(config['Spread']['IdealLatticeFirst']) != 0 else float(config['Spread']['QuadGradient'])
            RollError = 0 if MachineIndex == 0 and int(config['Spread']['IdealLatticeFirst']) != 0 else float(config['Spread']['QuadRoll'])
            while True:
                str = stream.readline()
                if str == "":
                    break
                A = re.split('[:,=\(\)&\s]+', str)
                if ((A[0].find("!") == -1)and(A[1].upper().find("QUAD") == 0)):
                    L = 0
                    K1 = 0
                    TILT = 0
                    for i in range(1, len(A)):
                        if (A[i-1].upper() == 'L'):
                            L = float(A[i])
                        if (A[i-1].upper() == 'K1'):
                            K1 = float(A[i])
                        if (A[i-1].upper() == 'TILT'):
                            TILT = float(A[i])
                    K1 *= 1 + GradError*(2*random.random() - 1)
                    TILT += RollError*(2*random.random() - 1)
                    str = ("%s: QUAD, L = %e, K1 = %e, TILT = %e\n") % (A[0], L, K1, TILT)
                if (A[0].find("%") == 0) and (A[2].find("sto") == 0) and (A[3].find(config['Run']['NameX'].strip('"')) == 0):
                    str = ("%% %e sto %s\n") % (float(config['Run']['StartX']), config['Run']['NameX'].strip('"'))
                if (A[0].find("%") == 0) and (A[2].find("sto") == 0) and (A[3].find(config['Run']['NameY'].strip('"')) == 0):
                    str = ("%% %e sto %s\n") % (float(config['Run']['StartY']), config['Run']['NameY'].strip('"'))
                streamout.write(str)
            streamout.close()
            stream.close()
        else:
            # Generate lattice files according to parameter sweep
            stream = open('machines/lattices/%.4d.lte' % (MachineIndex % int(config['Spread']['NumberOfMachines'])), 'r')
            streamout = open('machines/lattices/%.4d.lte' % MachineIndex, 'w')
            while True:
                str = stream.readline()
                if str == "":
                    break
                A = re.split('[:,=\(\)&\s]+', str)
                #print(A)
                if (A[0].find("%") == 0) and (A[2].find("sto") == 0) and (A[3].find(config['Run']['NameX'].strip('"')) == 0):
                    ValueX = float(config['Run']['StartX']) + float(config['Run']['StepX'])*((MachineIndex // int(config['Spread']['NumberOfMachines'])) % int(config['Run']['NstepsX']))
                    str = ("%% %e sto %s\n") % (ValueX, config['Run']['NameX'].strip('"'))
                if (A[0].find("%") == 0) and (A[2].find("sto") == 0) and (A[3].find(config['Run']['NameY'].strip('"')) == 0):
                    ValueY = float(config['Run']['StartY']) + float(config['Run']['StepY'])*(MachineIndex // (int(config['Spread']['NumberOfMachines']) * int(config['Run']['NstepsX'])))
                    str = ("%% %e sto %s\n") % (ValueY, config['Run']['NameY'].strip('"'))
                streamout.write(str)
            streamout.close()
            stream.close()
            print("machine %d, %s = %e, %s = %e\n" % (MachineIndex % int(config['Spread']['NumberOfMachines']), config['Run']['NameX'].strip('"'), ValueX, config['Run']['NameY'].strip('"'), ValueY))
        
        # Save transport matrices and cavity positions to 'prematrices.txt' and bunch energies at cavity positions to 'energies.txt'
        shutil.copy('machines/lattices/%.4d.lte' % MachineIndex, 'Temp/lattice.lte')
        os.system('elegant Temp/run.ele')
        os.system('sddsprocess Temp/run.mat Temp/cav.mat -match=column,ElementName=%s' % config['General']['CavityElement'])
        os.system('sddsprintout Temp/cav.mat -column=R* -column=s -width=100 > Temp/prematrices.txt')
        os.system('sddsprocess Temp/run.twi Temp/cav.twi -match=column,ElementName=%s' % config['General']['CavityElement'])
        os.system('sddsprintout Temp/cav.twi -column=pCentral0 > Temp/energies.txt')
        
        # Read transport matrices and cavity positions from 'prematrices.txt'
        E = [0.] * Npass
        beta = [0.] * Npass
        l = [0.] * Npass
        t = [0.] * Npass
        stream = open('Temp/prematrices.txt', 'r')
        for i in range(17):
            stream.readline()
        for i in range(Npass):
            t[i] = []
            for j in range(6):
                t[i].append([float(k) for k in stream.readline().split()])
            t[i] = np.array(t[i])                       # Transform transport matrices into NumPy matrices
            l[i] = float(stream.readline())
        stream.close()
        
        # Read bunch energies at cavity positions from 'energies.txt'
        stream = open('Temp/energies.txt', 'r')
        for i in range(5):
            stream.readline()
        for i in range(Npass):
            E[i] = float(stream.readline()) # beta*gamma
            beta[i] = E[i]/math.sqrt(E[i]*E[i] + 1)
            E[i] = Ee*math.sqrt(E[i]*E[i] + 1)
        stream.close()
        
        # Transform absolute transport matrices and cavity positions to relative transport matrices and distances
        T = [0.] * Npass
        T[0] = t[0]
        L = [0.] * Npass
        L[0] = l[0]
        for i in range(Npass-1):
            T[i+1] = t[i+1].dot(np.linalg.inv(t[i]))
            L[i+1] = l[i+1] - l[i]

        # Save energies, relative transport matrices and distances
        stream = open("machines/matrices/%.4d.txt" % MachineIndex, 'w')
        for i in range(Npass):
            stream.write("%e\n" % E[i])
            stream.write("%e\n" % L[i])
            for j in range(6):
                for k in range(6):
                    stream.write("%e\t" % T[i][j][k])
                stream.write("\n")
        stream.close()

        # Calculate sigma-matrices
        S1 = S0
        S = []
        for i in range(Npass):
            S1 = T[i].dot(S1.dot(T[i].transpose()))
            S.append(S1)

        # Calculate and save transition coefficients between all cavity pairs
        stream = open("machines/transition/%.4d.txt" % MachineIndex, 'w')
        for FirstCavity in range(Ncav):
            for LastCavity in range(Ncav):
                p0 = []     # Particle momentum at FirstCavity for each pass
                tr = []     # Recirculation time for each pass
                trt0 = []   # Ratio tr[i]/t0 (not saved into output file)
                Tmatr = []  # Transport matrix between FirstCavity and LastCavity for each pass
                Smatr = []  # Sigma-matrix at FirstCavity for each pass
                aquad = []  # Quadrupole transition coefficient a
                bquad = []  # Quadrupole transition coefficient b
                cquad = []  # Quadrupole transition coefficient c
                dquad = []  # Quadrupole transition coefficient d
                for i in range(math.ceil(Npass/Ncav)):
                    FirstArc = FirstCavity + Ncav*i
                    for j in range(math.ceil(Npass/Ncav)):
                        LastArc = LastCavity + Ncav*j
                        if (LastArc > FirstArc)and(LastArc < Npass):
                            p0.append((E[FirstArc]/c*beta[FirstArc] + E[FirstArc+1]/c*beta[FirstArc+1])/2)
                            dt = L[FirstArc+1]/(c*beta[FirstArc+1])
                            M = T[FirstArc+1]
                            for k in range(LastArc - FirstArc - 1):
                                dt = dt + L[FirstArc+k+2]/(c*beta[FirstArc+k+2])
                                M = T[FirstArc+k+2].dot(M)
                            S1 = S[FirstArc]
                            tr.append(dt)
                            trt0.append(dt/t0)
                            Tmatr.append(M)
                            Smatr.append(S1)
                            aquad.append(2*((M[2][0]*M[2][1]-M[0][0]*M[0][1])*S1[0][0] + (M[2][1]**2-M[0][1]**2)*S1[0][1] + 
                                (M[0][0]*M[0][3]-M[2][0]*M[2][3]+M[2][1]*M[2][2]-M[0][1]*M[0][2])*S1[0][2] + 
                                (M[2][1]*M[2][3]-M[0][1]*M[0][3])*S1[0][3] + (M[0][1]*M[0][3]-M[2][1]*M[2][3])*S1[1][2] + 
                                (M[0][2]*M[0][3]-M[2][2]*M[2][3])*S1[2][2] + (M[0][3]**2-M[2][3]**2)*S1[2][3]))
                            bquad.append(2*((M[0][0]*M[0][3]-M[2][0]*M[2][3])*S1[0][0] + (M[0][1]*M[0][3]-M[2][1]*M[2][3])*S1[0][1] + 
                                (M[0][0]*M[0][1]-M[2][0]*M[2][1]+M[0][2]*M[0][3]-M[2][2]*M[2][3])*S1[0][2] + 
                                (M[0][3]**2-M[2][3]**2)*S1[0][3] + (M[0][1]**2-M[2][1]**2)*S1[1][2] + 
                                (M[0][1]*M[0][2]-M[2][1]*M[2][2])*S1[2][2] + (M[0][1]*M[0][3]-M[2][1]*M[2][3])*S1[2][3]))
                            cquad.append(2*(-(M[0][0]*M[2][1]+M[0][1]*M[2][0])*S1[0][0] - 2*M[0][1]*M[2][1]*S1[0][1] + 
                                (M[0][3]*M[2][0]-M[0][1]*M[2][2]+M[0][0]*M[2][3]-M[0][2]*M[2][1])*S1[0][2] - 
                                (M[0][1]*M[2][3]+M[0][3]*M[2][1])*S1[0][3] + (M[0][3]*M[2][1]+M[0][1]*M[2][3])*S1[1][2] + 
                                (M[0][2]*M[2][3]+M[0][3]*M[2][2])*S1[2][2] + 2*M[0][3]*M[2][3]*S1[2][3]))
                            dquad.append(2*((M[0][0]*M[2][3]+M[0][3]*M[2][0])*S1[0][0] + (M[0][3]*M[2][1]+M[0][1]*M[2][3])*S1[0][1] + 
                                (M[0][3]*M[2][2]+M[0][0]*M[2][1]+M[0][1]*M[2][0]+M[0][2]*M[2][3])*S1[0][2] + 2*M[0][3]*M[2][3]*S1[0][3] + 
                                2*M[0][1]*M[2][1]*S1[1][2] + (M[0][2]*M[2][1]+M[0][1]*M[2][2])*S1[2][2] + (M[0][3]*M[2][1]+M[0][1]*M[2][3])*S1[2][3]))
                                
                print(*p0, "\n", *tr, file=stream)
                print(*aquad, "\n", *bquad, "\n", *cquad, "\n", *dquad, file=stream)
                for k in range(len(Tmatr)):
                    stream.write('%e\t' % (Tmatr[k][0][1]))
                stream.write('\n')
                for k in range(len(Tmatr)):
                    stream.write('%e\t' % (Tmatr[k][0][3]))
                stream.write('\n')
                for k in range(len(Tmatr)):
                    stream.write('%e\t' % (Tmatr[k][2][1]))
                stream.write('\n')
                for k in range(len(Tmatr)):
                    stream.write('%e\t' % (Tmatr[k][2][3]))
                stream.write('\n')
                for k in range(len(Tmatr)):
                    stream.write('%e\t' % (Tmatr[k][4][5]))
                stream.write('\n')
                
                p0.clear()
                tr.clear()
                trt0.clear()
                Tmatr.clear()
                Smatr.clear()
                aquad.clear()
                bquad.clear()
                cquad.clear()
                dquad.clear()
        stream.close()
