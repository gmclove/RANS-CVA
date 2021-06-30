# from pymooCFD.runPreProc import runPreProc
# runPreProc()

#from pymooCFD.optPreProc.gen0Map import gen0Map
#gen0Map()
class PtsList(object):
    pass
        
from pymooCFD.runOpt import runOpt
runOpt(restart=False)

# from pymooCFD import runPostProc
# runPostProc()

# from pymooCFD.util.handleData import compressDir
# compressDir('/dump')

# import base_case.rans_cva_sim