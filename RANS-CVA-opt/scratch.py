# n_sims = 10
# n = 0
#
#
# while n < n_sims:
#
#     if os.path.exists(f'{caseDir}/obj.txt'): #solver01_rank00.log'):
#         n +=1
#         continue

# import numpy as np
#
# meanC_dns = np.load('flowdata/meanC-dns.npy')
# meanC_cva = np.load('dump/meanC-cva.npy')
# print(meanC_dns[30,:,0].shape)
# print(meanC_cva.shape)

# for bin in meanC_dns:
#     print(bin.shape)
#     print(bin)
#
# print(meanC_dns[0])
# print(meanC_cva[0])
# print(meanC_dns)
#
# print(meanC_cva)
import numpy as np

x = np.array([[1,2],[3,5],[5,6]])
indices = [[i] for i in range(len(x))]
x = np.append(x, indices, axis=1)

print(x)

x = np.array([[1,2],[3,5],[5,6]])
indices = [*range(len(x))]
x = np.insert(x, 0, indices, axis=1)

print(x)
