# Stability analysis for monopole modes (approximated in case of multiple recirculations)

import configparser, os, math, cmath, shutil, numpy as np, time

# Read configuration files
machine = configparser.ConfigParser()
machine.read("machine.ini")
config = configparser.ConfigParser()
config.read(machine['General']['Machine'].strip('"'))

# Define constants and calculate the number of modes
c = 2.99792458e8
f0 = float(config['General']['FundamentalFrequency'])
w0 = 2*math.pi*f0
t0 = 1/f0
Nc = int(config['General']['NumberOfCavities'])
NmTotal = int(config['General']['NumberOfModes'])
NmMonopole = [int(i) for i in config['General']['NumberOfModesByOrder'].split(',')][0]
StartMode = 0

# Calculate the number of machines and active modes
ActiveMachine = int(config['Run']['Machine'])
ActiveModes = [int(i) for i in config['Run']['Modes'].split(',')]
InitialMachineIndex = 0 if ActiveMachine < 0 else ActiveMachine
MaxMachineIndex = int(config['Spread']['NumberOfMachines']) if int(config['Run']['UseSweep']) == 0 else int(config['Spread']['NumberOfMachines'])*int(config['Run']['NstepsX'])*int(config['Run']['NstepsY'])
FinalMachineIndex = MaxMachineIndex if ActiveMachine < 0 else ActiveMachine + 1
if not os.path.exists('machines/results'):
    os.makedirs('machines/results')
StreamByMachine = open('machines/results/ThresholdCurrentByMachine.txt', 'w') if ActiveMachine < 0 else 0
NumberOfModeEnsembles = NmMonopole if ActiveModes[0] < 0 else 1
Nm = len(ActiveModes) if NumberOfModeEnsembles == 1 else 1
NcNm = Nc*Nm
StreamByMode = open('machines/results/ThresholdCurrentByMode.txt', 'w') if ActiveModes[0] < 0 else 0
if NumberOfModeEnsembles > 1:
    stream = open(config['General']['CavityParameters'].strip('"'), 'r')
    for i in range(StartMode):
        stream.readline()
    for i in range(NumberOfModeEnsembles):
        f = [float(j) for j in stream.readline().split()][2]
        StreamByMode.write("%e\t" % (f/float(config['CCP']['OutputFrequencyUnits'])))
    StreamByMode.write("\n")
    stream.close()
StreamCCP = open('machines/results/ComplexCurrentPlot.txt', 'w') if ActiveMachine >= 0 and ActiveModes[0] >= 0 else 0

# The main cycle
start_time = time.time()
Imax = float(config['CCP']['MaxCurrent'])
step = float(config['CCP']['Step'])/t0
for MachineIndex in range(InitialMachineIndex, FinalMachineIndex):
    print("machine: %d" % MachineIndex)
    if os.path.exists('machines/cavities/WorstCaseIncluded') and MachineIndex <= 1:
        print('worst case configuration')
    if StreamByMachine != 0:
        StreamByMachine.write("%d\t" % MachineIndex)
    for ModeEnsembleIndex in range(NumberOfModeEnsembles):
        # Read mode data for the chosen machine and mode ensemble
        stream = open('machines/cavities/%.4d.txt' % MachineIndex, 'r')
        f = [0.] * Nc
        w = [0.] * Nc
        Q = [0.] * Nc
        Gamma = [0.] * Nc
        RoQd = [0.] * Nc
        theta = [0.] * Nc
        Theta = [0.] * Nc
        C = [0.] * Nc
        S = [0.] * Nc
        for i in range(Nc):
            f[i] = []
            w[i] = []
            Q[i] = []
            Gamma[i] = []
            RoQd[i] = []
            theta[i] = []
            Theta[i] = []
            C[i] = []
            S[i] = []
        if NumberOfModeEnsembles > 1:
            ActiveModes.clear()
        for i in range(Nc):
            for j in range(NmTotal):
                A = [float(j) for j in stream.readline().split()]
                if A[1] == 0:
                    if (NumberOfModeEnsembles == 1 and A[0] in ActiveModes) or (NumberOfModeEnsembles > 1 and (j - StartMode) == ModeEnsembleIndex):
                        if NumberOfModeEnsembles > 1 and i == 0:
                            ActiveModes.append(int(A[0]))
                        f[i].append(A[2])
                        Q[i].append(A[3])
                        RoQd[i].append(A[4])
                        theta[i].append(A[5])
            w[i] = [2*math.pi*j for j in f[i]]
            Gamma[i] = [w[i][j]/(2*Q[i][j]) for j in range(Nm)]
        if len(ActiveModes) == 0:
            print('ERROR: No active monopole modes.')
            sys.exit(0)
        stream.close()
        print("modes: ", end = "")
        print(*ActiveModes)

        # Read transition coefficients for the chosen machine and mode ensemble
        stream = open('machines/transition/%.4d.txt' % MachineIndex, 'r')
        p0 = [[0.] * Nc for i in range(Nc)]
        tr = [[0.] * Nc for i in range(Nc)]
        deltat0 = [[0.] * Nc for i in range(Nc)]
        T56 = [[0.] * Nc for i in range(Nc)]
        for i in range(Nc):
            for j in range(Nc):
                p0[i][j] = list(map(float, stream.readline().split()))
                tr[i][j] = list(map(float, stream.readline().split()))
                deltat0[i][j] = [math.ceil(tr[i][j][k]/t0)*t0 - tr[i][j][k] for k in range(len(tr[i][j]))]
                for n in range(8):
                    stream.readline()
                T56[i][j] = list(map(float, stream.readline().split()))
        stream.close()

        # Perform stability analysis
        W = [[0.] * NcNm for i in range(NcNm)]
        dww = float(config['CCP']['RelativeRange'])
        weff = (1-dww)*2*math.pi*f[0][0]
        weff_final = (1+dww)*2*math.pi*f[0][0]
        flag = 0
        Icur = [0.] * NcNm
        sortedIcur = [0.] * NcNm
        Iprev = [0.] * NcNm
        Iglobmin = 1e100
        Iglobmin2 = 1e100
        while weff < weff_final:
            # Calculate W(omega) matrix
            for i in range(Nc):
                for im in range(Nm):
                    wp = complex(w[i][im], Gamma[i][im])
                    wm = complex(w[i][im], -Gamma[i][im])
                    wpp = weff + wp
                    wpm = weff + wm
                    wmp = weff - wp
                    wmm = weff - wm
                    tmp1 = 1/(1-cmath.exp(complex(0,wmm*t0)))
                    tmp2 = 1/(1-cmath.exp(complex(0,-wpm*t0)))
                    tmp3 = 1/(1-cmath.exp(complex(0,wpp*t0)))
                    tmp4 = 1/(1-cmath.exp(complex(0,-wmp*t0)))
                    for j in range(Nc):
                        for jm in range(Nm):
                            summ = 0.
                            for k in range(len(p0[i][j])):
                                KoI = -T56[i][j][k]*RoQd[j][im]*w[i][im]*t0/(4*p0[i][j][k]*c**2)
                                tmp5 = cmath.exp(complex(0,weff*tr[i][j][k]+wmm*deltat0[i][j][k]))
                                tmp6 = cmath.exp(complex(0,-weff*tr[i][j][k]-wpm*deltat0[i][j][k]))
                                tmp7 = cmath.exp(complex(0,weff*tr[i][j][k]+wpp*deltat0[i][j][k]))
                                tmp8 = cmath.exp(complex(0,-weff*tr[i][j][k]-wmp*deltat0[i][j][k]))
                                summ = summ + KoI*complex(0,1)/2*(wm*(tmp5*tmp1+tmp6*tmp2)-wp*(tmp7*tmp3+tmp8*tmp4))
                            W[Nc*im+i][Nc*jm+j] = summ
            W = np.array(W)
            Weig = np.linalg.eigvals(W)
            
            # Find minimal positive eigenvalue (because all W(omega) eigenvalues are real)
            for k in range(NcNm):
                Icur[k] = 1/Weig[k]
                if (abs(Icur[k].imag) < 1e-10):
                    if (Icur[k] > 0)and(Icur[k] < Iglobmin):
                        Iglobmin = abs(Icur[k])
                        wglobmin = weff
                    if (Icur[k] < 0)and(abs(Icur[k]) < Iglobmin2):
                        Iglobmin2 = abs(Icur[k])
            if StreamCCP != 0:
                for k in range(NcNm):
                    if abs(Icur[k]) < Imax:
                        StreamCCP.write("%e\t%e\t%e\n" % (weff, Icur[k].real/float(config['CCP']['OutputCurrentUnits']), Icur[k].imag/float(config['CCP']['OutputCurrentUnits'])))
            flag = 1
            weff = weff + step
        print("Ith = %e A" % (Iglobmin))
        if StreamByMachine != 0:
            StreamByMachine.write("%e\t" % (Iglobmin/float(config['CCP']['OutputCurrentUnits'])))
        if StreamByMode != 0:
            StreamByMode.write("%e\t" % (Iglobmin/float(config['CCP']['OutputCurrentUnits'])))
    time_elapsed = time.time() - start_time
    time_remaining = time_elapsed * (FinalMachineIndex - MachineIndex - 1)/(MachineIndex - InitialMachineIndex + 1)
    print("elapsed %.3f h  remaining %.3f h" % (time_elapsed/3600, time_remaining/3600))
    print("")
    if StreamByMachine != 0:
        StreamByMachine.write("\n")
    if StreamByMode != 0:
        StreamByMode.write("\n")
if StreamByMachine != 0:
    StreamByMachine.close()
if StreamByMode != 0:
    StreamByMode.close()
if StreamCCP != 0:
    StreamCCP.close()
print("--- %s seconds ---" % (time.time() - start_time))
