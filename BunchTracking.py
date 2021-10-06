# Bunch tracking (exact)

import configparser, os, math, cmath, random, numpy as np, time

Ee = 0.5109989461e6
c = 2.99792458e8

# Read configuration files
machine = configparser.ConfigParser()
machine.read("machine.ini")
config = configparser.ConfigParser()
config.read(machine['General']['Machine'].strip('"'))

f0 = float(config['General']['FundamentalFrequency'])
t0 = 1/f0
w0 = 2*math.pi*f0

I = float(config['Tracking']['Current'])
Na = int(config['General']['NumberOfCavityPasses']) if int(config['Run']['UsedCavityPasses']) < 0 else int(config['Run']['UsedCavityPasses'])
Nb = int(config['Tracking']['NumberOfBunches'])
Nc = int(config['General']['NumberOfCavities'])

stream = open('machines/matrices/%.4d.txt' % int(config['Run']['Machine']), 'r')
Ea = [0.] * Na
La = [0.] * Na
Ta = [0.] * Na
for i in range(Na):
    Ea[i] = float(stream.readline())
    La[i] = float(stream.readline())
    Ta[i] = []
    for j in range(6):
        Ta[i].append([float(k) for k in stream.readline().split()])
stream.close()

Va = [c*math.sqrt(1-(Ee/i)**2) for i in Ea]
Zcav = [La[0]] * Na
for i in range(1,Na):
    Zcav[i] = Zcav[i-1] + La[i]
        
stream = open('machines/cavities/%.4d.txt' % int(config['Run']['Machine']), 'r')
if os.path.exists('machines/cavities/WorstCaseIncluded') and int(config['Run']['Machine']) <= 1:
    print('worst case configuration')
ActiveModes = [int(i) for i in config['Run']['Modes'].split(',')]
Nm = len(ActiveModes)
NmTotal = int(config['General']['NumberOfModes'])
ord = [0.] * Nc
f = [0.] * Nc
w = [0.] * Nc
Q = [0.] * Nc
Gamma = [0.] * Nc
RoQd = [0.] * Nc
theta = [0.] * Nc
Theta = [0.] * Nc
for i in range(Nc):
    ord[i] = []
    f[i] = []
    w[i] = []
    Q[i] = []
    Gamma[i] = []
    RoQd[i] = []
    theta[i] = []
    Theta[i] = []
QuadrupoleModesActive = False
for i in range(Nc):
    for j in range(NmTotal):
        A = [float(j) for j in stream.readline().split()]
        if A[0] in ActiveModes:
            ord[i].append(int(A[1]))
            if int(A[1]) == 2:
                QuadrupoleModesActive = True
            f[i].append(A[2])
            Q[i].append(A[3])
            RoQd[i].append(A[4])
            theta[i].append(A[5])
    w[i] = [2*math.pi*j for j in f[i]]
    Gamma[i] = [w[i][j]/(2*Q[i][j]) for j in range(Nm)]
    Theta[i] = [j/180*math.pi for j in theta[i]]
stream.close()

q = I*t0
BunchSpacing = Va[0]*t0
X = [[0.] * 6 for i in range(Nb)]
Z = [-i*BunchSpacing for i in range(Nb)]
BunchUpdateTime = [0.] * Nb
BunchesToCheck = [Nb] * Na
BunchesToCheck[0] = 0
InitialModeAmplitudes = [int(i) for i in config['Tracking']['InitialModeAmplitudes'].split(',')]
V = [[0.] * Nm for i in range(Nc)]
for ic in range(Nc):
    for im in range(Nm):
        V[ic][im] = InitialModeAmplitudes[ord[ic][im]]
CavityUpdateTime = [0.] * Nc
t = 0.
LinesCount = [0] * Nc

sumt = [[0.] * Nm for i in range(Nc)]
sumt2 = [[0.] * Nm for i in range(Nc)]
sumlnV = [[0.] * Nm for i in range(Nc)]
sumtlnV = [[0.] * Nm for i in range(Nc)]
countN = [[0] * Nm for i in range(Nc)]
decr = [[0] * Nm for i in range(Nc)]

sumV_1 = [[0.] * Nm for i in range(Nc)]
sumV2_1 = [[0.] * Nm for i in range(Nc)]
countN_1 = [[0] * Nm for i in range(Nc)]
sumV_2 = [[0.] * Nm for i in range(Nc)]
sumV2_2 = [[0.] * Nm for i in range(Nc)]
countN_2 = [[0] * Nm for i in range(Nc)]

if QuadrupoleModesActive:
    for i in range(Nb):
        X[i] = np.array(X[i])
    for i in range(Na):
        Ta[i] = np.array(Ta[i])
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
    #print(S0)
    S = []
    for i in range(Nb):
        S.append(S0.copy())

if not os.path.exists('Temp'):
    os.mkdir('Temp')
streamdecr = []
streamvolt = []
for i in range(Nc):
    streamdecr.append(open('Temp/decrement%.4d.txt'%i, 'w'))
    streamvolt.append(open('Temp/voltage%.4d.txt'%i, 'w'))
DataPointsPerWrite = max(1, int(int(config['Tracking']['NumberOfBunches']) * \
    int(config['General']['NumberOfCavityPasses']) / int(config['General']['NumberOfCavities']) / int(config['Tracking']['MaxVoltLines'])))
print('%d bunches, I = %lf A     Performing tracking...' % (int(config['Tracking']['NumberOfBunches']), float(config['Tracking']['Current'])))
start_time = time.time()
LinearityCheck = False
if LinearityCheck:
    dummy1 = [[] for i in range(Nb)]
    dummy2 = [[] for i in range(Nb)]
for turn in range(Na*Nb):
    dtmin = 1e100
    for i in range(Na):
        if BunchesToCheck[i] < Nb:
            dt = (Zcav[i] - Z[BunchesToCheck[i]])/Va[i] - (t - BunchUpdateTime[BunchesToCheck[i]])
            if dt < dtmin:
                dtmin = dt
                ia = i                  # index of the arc which loses the bunch
                ib = BunchesToCheck[i]  # index of the bunch

    t = t + dtmin
    Xtmp = [0.] * 6
    for i in range(6):
        for j in range(6):
            Xtmp[i] = Xtmp[i] + Ta[ia][i][j]*X[ib][j]
    X[ib] = Xtmp
    if QuadrupoleModesActive:
        S[ib] = Ta[ia].dot(S[ib].dot(Ta[ia].transpose()))
    Z[ib] = Zcav[ia]
    BunchUpdateTime[ib] = t
    BunchesToCheck[ia] += 1
    if (ia < Na-1)and(ib == 0):
        BunchesToCheck[ia+1] = 0
    
    ic = ia % Nc                        # index of the cavity
    dt = t - CavityUpdateTime[ic] if int(config['Tracking']['LongitudinalMotion']) == 0 else t + X[ib][4]/c - CavityUpdateTime[ic]
    p = Ea[ia]*Va[ia]/c**2
    if ia < Na-1:
        p = (Ea[ia]*Va[ia]/c**2 + Ea[ia+1]*Va[ia+1]/c**2)/2
    for im in range(Nm):
        if ord[ic][im] == 0:
            V[ic][im] = V[ic][im]*cmath.rect(math.exp(-Gamma[ic][im]*dt), w[ic][im]*dt) - 0.5*RoQd[ic][im]*q*w[ic][im] / 2
            if LinearityCheck:
                dummy1[ib].append(V[ic][im].real)
                dummy2[ib].append(- 0.5*RoQd[ic][im]*q*w[ic][im] / 2)
            X[ib][5] += V[ic][im].real/Ea[ia]
            V[ic][im] += - 0.5*RoQd[ic][im]*q*w[ic][im] / 2
        elif ord[ic][im] == 1:
            sintheta = math.sin(Theta[ic][im])
            costheta = math.cos(Theta[ic][im])
            V[ic][im] = V[ic][im]*cmath.rect(math.exp(-Gamma[ic][im]*dt), w[ic][im]*dt) + RoQd[ic][im]*q*c*(X[ib][0]*costheta + X[ib][2]*sintheta)
            if LinearityCheck:
                dummy1[ib].append(V[ic][im].imag)
                dummy2[ib].append(RoQd[ic][im]*q*c*(X[ib][0]*costheta + X[ib][2]*sintheta))
            X[ib][1] += V[ic][im].imag/(2*p*c)*costheta
            X[ib][3] += V[ic][im].imag/(2*p*c)*sintheta
        elif ord[ic][im] == 2:
            V[ic][im] = (V[ic][im]*cmath.rect(math.exp(-Gamma[ic][im]*dt), w[ic][im]*dt) + 
                RoQd[ic][im]*2*q*c**2/w[ic][im]*((S[ib][0][0]-S[ib][2][2])*math.cos(2*Theta[ic][im]) + 2*S[ib][0][2]*math.sin(2*Theta[ic][im])))
            if LinearityCheck:
                dummy1[ib].append(V[ic][im].imag)
                dummy2[ib].append(RoQd[ic][im]*2*q*c**2/w[ic][im]*((S[ib][0][0]-S[ib][2][2])*math.cos(2*Theta[ic][im]) + 2*S[ib][0][2]*math.sin(2*Theta[ic][im])))
            P = w[ic][im]/(4*p*c**2)*V[ic][im].imag
            s2 = P*math.sin(2*Theta[ic][im])
            c2 = P*math.cos(2*Theta[ic][im])
            Tq = np.array([[1.,0.,0.,0.,0.,0.],[-c2,1.,s2,0.,0.,0.],[0.,0.,1.,0.,0.,0.],[s2,0.,c2,1.,0.,0.],[0.,0.,0.,0.,1.,0.],[0.,0.,0.,0.,0.,1.]])
            X[ib] = Tq.dot(X[ib])
            S[ib] = Tq.dot(S[ib].dot(Tq.transpose()))
    CavityUpdateTime[ic] = t if int(config['Tracking']['LongitudinalMotion']) == 0 else t + X[ib][4]/c

    if (t > 0.5*Nb*t0) and (t < 0.75*Nb*t0):
        for i in range(len(ActiveModes)):
            sumV_1[ic][i] += abs(V[ic][i])
            sumV2_1[ic][i] += abs(V[ic][i])**2
            countN_1[ic][i] += 1
    if (t > 0.75*Nb*t0) and (t < 1*Nb*t0):
        for i in range(len(ActiveModes)):
            sumV_2[ic][i] += abs(V[ic][i])
            sumV2_2[ic][i] += abs(V[ic][i])**2
            countN_2[ic][i] += 1
        
    for i in range(len(ActiveModes)):
        if (t > 0.5*Nb*t0) and (t < 0.75*Nb*t0):
            sumt[ic][i] += t
            sumt2[ic][i] += t**2
            sumlnV[ic][i] += math.log(abs(V[ic][i]))
            sumtlnV[ic][i] += t*math.log(abs(V[ic][i]))
            countN[ic][i] += 1
    if LinesCount[ic] % DataPointsPerWrite == 0:
        if countN[ic][0] > 1:
            streamdecr[ic].write("%e    " % (t/float(config['Tracking']['OutputTimeUnits'])))
        for i in range(len(ActiveModes)):
            if countN[ic][i] > 1:
                decr[ic][i] = (countN[ic][i]*sumtlnV[ic][i] - sumlnV[ic][i]*sumt[ic][i])/(countN[ic][i]*sumt2[ic][i] - sumt[ic][i]**2)
                streamdecr[ic].write("%e    " % (decr[ic][i]))
        if countN[ic][0] > 1:
            streamdecr[ic].write("\n");
    
    if LinesCount[ic] % DataPointsPerWrite == 0:
        streamvolt[ic].write("%e    " % (t/float(config['Tracking']['OutputTimeUnits'])))
        for i in range(len(ActiveModes)):
            streamvolt[ic].write("%e    " % abs(V[ic][i]));
        streamvolt[ic].write("\n");
    LinesCount[ic] += 1
for i in range(Nc):
    streamdecr[i].close()
    streamvolt[i].close()

if LinearityCheck:
    for i in range(Nb):
        print('%e' % (dummy2[i][1]/dummy1[i][0] if dummy1[i][0] != 0 else 0))

status = 'Decrement: stable'
for ic in range(Nc):
    for i in range(len(ActiveModes)):
        if decr[ic][i] > 0:
            status = 'Decrement: unstable'
print(status)

status = 'Dispersion: stable'
for ic in range(Nc):
    for i in range(len(ActiveModes)):
        #print(sumV2_1[ic][i]/countN_1[ic][i] - (sumV_1[ic][i]/countN_1[ic][i])**2, sumV2_2[ic][i]/countN_2[ic][i] - (sumV_2[ic][i]/countN_2[ic][i])**2)
        if (sumV2_1[ic][i]/countN_1[ic][i] - (sumV_1[ic][i]/countN_1[ic][i])**2 < sumV2_2[ic][i]/countN_2[ic][i] - (sumV_2[ic][i]/countN_2[ic][i])**2):
            status = 'Dispersion: unstable'
print(status)

for i in range(Nc):
    streamdecr[i] = open('Temp/decrement%.4d.txt'%i, 'r')
stream = open('machines/results/Decrement.txt', 'w')
flag = 0
while True:
    outstr = ''
    for i in range(Nc):
        str = streamdecr[i].readline()
        if str == "":
            flag = 1
            break
        outstr += str.rstrip() + '    '
    stream.write(outstr + '\n')
    if flag == 1:
        break
stream.close()
for i in range(Nc):
    streamdecr[i].close()
   
for i in range(Nc):
    streamvolt[i] = open('Temp/voltage%.4d.txt'%i, 'r')
stream = open('machines/results/Voltage.txt', 'w')
flag = 0
while True:
    outstr = ''
    for i in range(Nc):
        str = streamvolt[i].readline()
        if str == "":
            flag = 1
            break
        outstr += str.rstrip() + '    '
    stream.write(outstr + '\n')
    if flag == 1:
        break
stream.close()
for i in range(Nc):
    streamvolt[i].close()
   

DataPointsPerWrite = int(int(config['Tracking']['NumberOfBunches']) / int(config['Tracking']['MaxCoordLines']))
stream = open('machines/results/FinalCoords.txt', 'w')
for i in range(Nb):
    if i % DataPointsPerWrite == 0:
        stream.write("%d" % i)
        for j in range(6):
            stream.write("\t%e" % X[i][j])
        stream.write("\n")
stream.close()
print("--- %s seconds ---" % (time.time() - start_time))
