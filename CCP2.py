# Stability analysis for quadrupole modes (approximated in case of multiple recirculations)

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
NmDipole = [int(i) for i in config['General']['NumberOfModesByOrder'].split(',')][1]
NmQuadrupole = [int(i) for i in config['General']['NumberOfModesByOrder'].split(',')][2]
StartMode = NmMonopole + NmDipole

# Calculate the number of machines and active modes
ActiveMachine = int(config['Run']['Machine'])
ActiveModes = [int(i) for i in config['Run']['Modes'].split(',')]
InitialMachineIndex = 0 if ActiveMachine < 0 else ActiveMachine
MaxMachineIndex = int(config['Spread']['NumberOfMachines']) if int(config['Run']['UseSweep']) == 0 else int(config['Spread']['NumberOfMachines'])*int(config['Run']['NstepsX'])*int(config['Run']['NstepsY'])
FinalMachineIndex = MaxMachineIndex if ActiveMachine < 0 else ActiveMachine + 1
if not os.path.exists('machines/results'):
    os.makedirs('machines/results')
StreamByMachine = open('machines/results/ThresholdCurrentByMachine.txt', 'w') if ActiveMachine < 0 else 0
NumberOfModeEnsembles = NmQuadrupole // 2 if ActiveModes[0] < 0 else 1
Nm = len(ActiveModes) if NumberOfModeEnsembles == 1 else 2
NcNm = Nc*Nm
StreamByMode = open('machines/results/ThresholdCurrentByMode.txt', 'w') if ActiveModes[0] < 0 else 0
if NumberOfModeEnsembles > 1:
    stream = open(config['General']['CavityParameters'].strip('"'), 'r')
    for i in range(StartMode):
        stream.readline()
    for i in range(NumberOfModeEnsembles):
        f1 = [float(j) for j in stream.readline().split()][2]
        f2 = [float(j) for j in stream.readline().split()][2]
        StreamByMode.write("%e\t" % (0.5*(f1+f2)/float(config['CCP']['OutputFrequencyUnits'])))
    StreamByMode.write("\n")
    stream.close()
StreamCCP = open('machines/results/ComplexCurrentPlot.txt', 'w') if ActiveMachine >= 0 and ActiveModes[0] >= 0 else 0
StreamCCPsorted = open('machines/results/ComplexCurrentPlotSorted.txt', 'w') if ActiveMachine >= 0 and ActiveModes[0] >= 0 else 0
StreamCCPunsorted = open('machines/results/ComplexCurrentPlotUnsorted.txt', 'w') if ActiveMachine >= 0 and ActiveModes[0] >= 0 else 0

# The main cycle
start_time = time.time()
Imax = float(config['CCP']['MaxCurrent'])
step = float(config['CCP']['Step'])/t0
for MachineIndex in range(InitialMachineIndex, FinalMachineIndex):
    print("machine: %d" % MachineIndex)
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
                if A[1] == 2:
                    if (NumberOfModeEnsembles == 1 and A[0] in ActiveModes) or (NumberOfModeEnsembles > 1 and (j - StartMode) // 2 == ModeEnsembleIndex):
                        if NumberOfModeEnsembles > 1 and i == 0:
                            ActiveModes.append(int(A[0]))
                        f[i].append(A[2])
                        Q[i].append(A[3])
                        RoQd[i].append(A[4])
                        theta[i].append(A[5])
            w[i] = [2*math.pi*j for j in f[i]]
            Gamma[i] = [w[i][j]/(2*Q[i][j]) for j in range(Nm)]
            Theta[i] = [j/180*math.pi for j in theta[i]]
            C[i] = [math.cos(2*j) for j in Theta[i]]
            S[i] = [math.sin(2*j) for j in Theta[i]]
        if len(ActiveModes) == 0:
            print('ERROR: No active quadrupole modes.')
            sys.exit(0)
        stream.close()
        print("modes: ", end = "")
        print(*ActiveModes)

        # Read transition coefficients for the chosen machine and mode ensemble
        stream = open('machines/transition/%.4d.txt' % MachineIndex, 'r')
        p0 = [[0.] * Nc for i in range(Nc)]
        tr = [[0.] * Nc for i in range(Nc)]
        deltat0 = [[0.] * Nc for i in range(Nc)]
        aquad = [[0.] * Nc for i in range(Nc)]
        bquad = [[0.] * Nc for i in range(Nc)]
        cquad = [[0.] * Nc for i in range(Nc)]
        dquad = [[0.] * Nc for i in range(Nc)]
        for i in range(Nc):
            for j in range(Nc):
                p0[i][j] = list(map(float, stream.readline().split()))
                tr[i][j] = list(map(float, stream.readline().split()))
                deltat0[i][j] = [math.ceil(tr[i][j][k]/t0)*t0 - tr[i][j][k] for k in range(len(tr[i][j]))]
                aquad[i][j] = list(map(float, stream.readline().split()))
                bquad[i][j] = list(map(float, stream.readline().split()))
                cquad[i][j] = list(map(float, stream.readline().split()))
                dquad[i][j] = list(map(float, stream.readline().split()))
                for n in range(5):
                    stream.readline()
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
        Ith = 1e100
        TotalIntersections = [0.] * NcNm
        MinIntersection = [0.] * NcNm
        IthByMode = [1e100] * NcNm
        while weff < weff_final:
            # Calculate W(omega) matrix
            for i in range(Nc):
                for im in range(Nm):
                    ci = C[i][im]
                    si = S[i][im]
                    tmp1 = 1/(1-cmath.exp(complex(-Gamma[i][im], weff + w[i][im])*t0))
                    tmp2 = 1/(1-cmath.exp(complex(-Gamma[i][im], weff - w[i][im])*t0))
                    for j in range(Nc):
                        for jm in range(Nm):
                            summ = 0.
                            cj = C[j][jm]
                            sj = S[j][jm]
                            for k in range(len(p0[i][j])):
                                D = aquad[i][j][k]*ci*cj + bquad[i][j][k]*si*cj + cquad[i][j][k]*ci*sj + dquad[i][j][k]*si*sj
                                KoI = D*RoQd[j][im]*t0/(2*p0[i][j][k])
                                tmp3 = cmath.exp(complex(-Gamma[i][im], weff + w[i][im])*deltat0[i][j][k])
                                tmp4 = cmath.exp(complex(-Gamma[i][im], weff - w[i][im])*deltat0[i][j][k])
                                summ = summ + KoI/(2*complex(0,1))*cmath.exp(complex(0,weff*tr[i][j][k]))*(tmp3*tmp1 - tmp4*tmp2)
                            W[Nc*im+i][Nc*jm+j] = summ
            W = np.array(W)
            Weig = np.linalg.eigvals(W)
            
            # Sort W(omega) eigenvalues 
            for k in range(NcNm):
                Icur[k] = 1/Weig[k]
            if flag == 0:
                for k in range(NcNm):
                    sortedIcur[k] = Icur[k]
            if flag == 1:
                usedIcur = [0] * NcNm
                usedIprev = [0] * NcNm
                for i in range(NcNm):
                    Imin = 1e100
                    for j in range(NcNm):
                        if (abs(Icur[j]) < Imin)and(usedIcur[j] == 0):
                            Imin = abs(Icur[j])
                            indImin = j
                    usedIcur[indImin] = 1
                    Dmin = 1e100
                    for j in range(Nc*Nm):
                        if (abs(Icur[indImin] - Iprev[j]) < Dmin)and(usedIprev[j] == 0):
                            Dmin = abs(Icur[indImin] - Iprev[j])
                            indDmin = j
                    usedIprev[indDmin] = 1
                    sortedIcur[indDmin] = Icur[indImin]
                for i in range(NcNm):
                    if ((sortedIcur[i].imag>=0)and(Iprev[i].imag<0))or((sortedIcur[i].imag<0)and(Iprev[i].imag>=0)):
                        Istar = (Iprev[i].real*sortedIcur[i].imag - sortedIcur[i].real*Iprev[i].imag)/(sortedIcur[i].imag - Iprev[i].imag)
                        if (Istar > 0)and(Istar < Ith):
                            Ith = Istar
                            MinIndex = i
                        if (Istar > 0)and(Istar < IthByMode[i]):
                            IthByMode[i] = Istar
                            MinIntersection[i] = TotalIntersections[i]
                        TotalIntersections[i] += 1
            if StreamCCP != 0:
                for k in range(NcNm):
                    if abs(sortedIcur[k]) < Imax:
                        StreamCCP.write("%e\t%e\t%e\n" % (weff, sortedIcur[k].real/float(config['CCP']['OutputCurrentUnits']), sortedIcur[k].imag/float(config['CCP']['OutputCurrentUnits'])))
                    StreamCCPsorted.write("%e\t%e\t" % (sortedIcur[k].real/float(config['CCP']['OutputCurrentUnits']), sortedIcur[k].imag/float(config['CCP']['OutputCurrentUnits'])))
                    StreamCCPunsorted.write("%e\t%e\t" % (Icur[k].real/float(config['CCP']['OutputCurrentUnits']), Icur[k].imag/float(config['CCP']['OutputCurrentUnits'])))
                StreamCCPsorted.write("\n")
                StreamCCPunsorted.write("\n")
            for k in range(NcNm):
                Iprev[k] = sortedIcur[k]
            flag = 1
            weff = weff + step
        print("Ith = %e A" % Ith)
        if StreamCCP != 0:
            for i in range(NcNm):
                print('      Curve %2d : minimal positive intersection No. %2d of %2d at %e A%s' % (i+1, MinIntersection[i], TotalIntersections[i], IthByMode[i], (' is I[th]' if i == MinIndex else '')))
        for i in range(NcNm):
            if (MinIntersection[i] == TotalIntersections[i]) or (MinIntersection[i] == 0):
                print('                 WARNING: Frequency range is too small')
        if StreamByMachine != 0:
            StreamByMachine.write("%e\t" % (Ith/float(config['CCP']['OutputCurrentUnits'])))
        if StreamByMode != 0:
            StreamByMode.write("%e\t" % (Ith/float(config['CCP']['OutputCurrentUnits'])))
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
    StreamCCPsorted.close()
    StreamCCPunsorted.close()
print("--- %s seconds ---" % (time.time() - start_time))
