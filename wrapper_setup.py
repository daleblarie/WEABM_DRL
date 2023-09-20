# compile all .cpp files with the collowing command
# g++ -fPIC -Wall -shared -o Simulation.so ./CPP_Base/*.cpp -std=c++11 -Ofast
# g++ -fPIC -Wall -shared -o Simulation_estim.so ./CPP_Estim/*.cpp -std=c++11 -Ofast
# g++ -fPIC -Wall -shared -o Simulation_estim_ecm.so ./CPP_Estim+ECM/*.cpp -std=c++11 -Ofast
# g++ -fPIC -Wall -shared -o Simulation_F1.so ./CPP_F1/*.cpp -std=c++11 -Ofast
# g++ -fPIC -Wall -shared -o Simulation_F2.so ./CPP_F2/*.cpp -std=c++11 -Ofast


import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer

def setUpWrapper():
    xDim=30;
    yDim=30;
    zDim=200;
    MAX_STEPS = 1400;
    # Change which CDLL is loaded to change which CPP simulation is actually run
    sim = ctypes.CDLL('./Simulation.so')

    # sim.CreateInstance.argtypes = (ctypes.POINTER(ctypes.c_float),)
    sim.CreateInstance.restype = ctypes.POINTER(ctypes.c_void_p)

    sim.InitializeSimulation.argtypes = (ctypes.POINTER(ctypes.c_void_p),)

    sim.DoInjury.argtypes = (ctypes.POINTER(ctypes.c_void_p),)

    sim.StepSimulation.argtypes = (ctypes.POINTER(ctypes.c_void_p),)

    sim.getSurfaceCytokineLevels.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.getSurfaceCytokineLevels.restype = ndpointer(dtype=ctypes.c_float, shape=(28))

    sim.getTotalCytokineLevels.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.getTotalCytokineLevels.restype = ndpointer(dtype=ctypes.c_float, shape=(28))

    sim.getNumPMN.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.getNumPMN.restype = ctypes.c_int
    sim.getNumMacrophage.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.getNumMacrophage.restype = ctypes.c_int
    sim.getNumFibroblast.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.getNumFibroblast.restype = ctypes.c_int
    sim.getNumTH0.argtypes =  (ctypes.POINTER(ctypes.c_void_p),)
    sim.getNumTH0.restype = ctypes.c_int
    sim.getNumTH1.argtypes =  (ctypes.POINTER(ctypes.c_void_p),)
    sim.getNumTH1.restype = ctypes.c_int
    sim.getNumTH2.argtypes =  (ctypes.POINTER(ctypes.c_void_p),)
    sim.getNumTH2.restype = ctypes.c_int
    sim.getNumTH17.argtypes =  (ctypes.POINTER(ctypes.c_void_p),)
    sim.getNumTH17.restype = ctypes.c_int


    sim.get_xDim.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_xDim.restype = ctypes.c_int
    sim.get_yDim.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_yDim.restype = ctypes.c_int
    sim.get_current_step.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_current_step.restype = ctypes.c_int

    sim.getInnervation.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.getInnervation.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))

    sim.getCollagen1.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.getCollagen1.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))

    sim.getCollagen3.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.getCollagen3.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))

    sim.getViableMuscle.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.getViableMuscle.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))

    sim.getVascularization.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.getVascularization.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))

    sim.getLife.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.getLife.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))

    sim.getNecrosis.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.getNecrosis.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))

    sim.getContamination.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.getContamination.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_IL8.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_IL8.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_MIP1b.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_MIP1b.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_TNF.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_TNF.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_HGF.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_HGF.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_PAF.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_PAF.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_IL1.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_IL1.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_PDGF.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_PDGF.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_IL13.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_IL13.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_TGFb.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_TGFb.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_VEGF.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_VEGF.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_GCSF.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_GCSF.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_IL10.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_IL10.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_IL4.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_IL4.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_IL17.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_IL17.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_Myf5.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_Myf5.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_IL6.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_IL6.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_IL1ra.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_IL1ra.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_sIL1r.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_sIL1r.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_sTNFr.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_sTNFr.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_endotoxin.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_endotoxin.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_IFNg.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_IFNg.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_IL12.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_IL12.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_MRF4.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_MRF4.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_cytotox.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_cytotox.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_IGF.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_IGF.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_DAMP.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_DAMP.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_MCP1.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_MCP1.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))
    sim.get_FNE.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.get_FNE.restype = ndpointer(dtype=ctypes.c_float, shape=(zDim,yDim,xDim))

    sim.apply_IL8.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_MIP1b.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_TNF.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_HGF.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_PAF.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_IL1.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_PDGF.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_IL13.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_TGFb.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_VEGF.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_GCSF.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_IL10.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_IL4.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_IL17.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_Myf5.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_IL6.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_IL1ra.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_sIL1r.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_sTNFr.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_endotoxin.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_IFNg.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_IL12.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_MRF4.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_cytotox.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_IGF.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_DAMP.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_MCP1.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)
    sim.apply_FNE.argtypes = (ctypes.POINTER(ctypes.c_void_p), ctypes.c_float)

    sim.endSimulation.argtypes = (ctypes.POINTER(ctypes.c_void_p),)
    sim.endSimulation.restype = ndpointer(dtype=ctypes.c_float, shape=(43,MAX_STEPS))

    sim.setSeed.argtypes = (ctypes.POINTER(ctypes.c_void_p),ctypes.c_int)


    return sim
if __name__ == '__main__':
    setUpWrapper(simObj)
