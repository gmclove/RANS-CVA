from scipy.signal import periodogram
import numpy as np
import matplotlib.pyplot as plt

with open('dump/ics_temporals.txt', 'r') as f:
    tempData = np.loadtxt(f, skiprows=1, usecols=(5,6,7,8))  # [::5]
print('temporal data')
print(tempData)
print(len(tempData), len(tempData[0]))

# P_OVER_RHO_INTGRL_1 = tempData[:, 5]
# print(P_OVER_RHO_INTGRL_1)
# print(len(P_OVER_RHO_INTGRL_1))

# TAU_INTGRL_NORM = tempData[:, 9]
# print(TAU_INTGRL_NORM)

# x = tempData[:,1]
# freq, psd=periodogram(tempData, fs=4)
# freq, psd=periodogram(P_OVER_RHO_INTGRL_1, fs=4)
# freq, psd=periodogram(x, fs=4)

# psd=np.transpose(psd)

# print('Frequencies')
# print(freq)
#
# print('PSD')
# print(psd)
# print('PSD Transposed')
# print(np.transpose(psd))
# print('Square Root of PSD Max')
# print(np.sqrt(psd.max()))
# print('shape of psd matix')
# print(len(psd), len(psd[0]))
col = ['P_OVER_RHO_INTGRL_(1)', 'P_OVER_RHO_INTGRL_(2)', 'TAU_INTGRL_(1)', 'TAU_INTGRL_(2)']

def plotPSD(freq, psd, meas):
    plt.plot(freq, psd)
    plt.xlabel('Frequencies')
    plt.ylabel('PSD')
    plt.title(meas)
    plt.savefig(f'plots/freq-psd_{meas.replace("_", "-")}.png')
    plt.close()

    plt.plot(freq, psd)
    plt.xlabel('Frequencies')
    plt.ylabel('PSD')
    plt.title(meas)
    plt.xlim((0,0.25))
    plt.savefig(f'plots/freq-psd_{meas.replace("_", "-")}-zoom.png')
    plt.close()
    # plt.show()

tempData = np.transpose(tempData)
for i, meas in enumerate(tempData):
    print(col[i])
    print(meas)
    freq, psd=periodogram(meas, fs=10000/100) #, scaling='spectrum')
    plotPSD(freq, psd, col[i])
    print('Dominant frequency')
    print(freq[psd.argmax()])

# plt.semilogy(freq, psd)
# plt.xlabel('Frequencies')
# plt.ylabel('PSD')
# plt.show()

# plt.psd(tempData)
# plt.show()
