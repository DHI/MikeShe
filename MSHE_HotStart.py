'''
Script to initiate a MIKESHE WM simulation using unsat zone water content from a previous run as a hot start condition.
'''

print('Importing Required Libraries...')
# Import MIKESHE
import importlib.util
spec = importlib.util.spec_from_file_location('MShePy311', r'C:\Program Files (x86)\DHI\MIKE Zero\2025\bin\x64\MShePy311.pyd')
MShePy = importlib.util.module_from_spec(spec)

# Import other required libraries including MIKEIO
import mikeio
import numpy as np

# Load DFS3 Result File of Unsat Zone and get Water Content Array
print('Importing Water Content Data From Previous Model Run...')
uz_data = mikeio.read(r'D:\Projects\MIKESHE_Hotstart\TestModel\TestModel_3DUZ.dfs3')
wc_arr = uz_data['water content in unsaturated zone']
wc_arr_finalts = wc_arr[-1,:,:,:].values # Access array of 3D water content values at final timestep
#Transpose Z dimensions from first dim of array to last to align with MShePy dataset structure
wc_arr_finalts = np.transpose(wc_arr_finalts, (1, 2, 0))

# Initiate MIKE SHE Model
print('Initializing MIKESHE Model...')
model = r'D:\Projects\MIKESHE_Hotstart\TestModel\TestModel.she'
MShePy.wm.initialize(model)
MShePy

simperiod = MShePy.wm.simPeriod()
start_time = simperiod[0]
end_time = simperiod[-1]
assign_timesteps = 10 # Number of timesteps to execute water content overwrite

for ts in range(assign_timesteps):
    print(ts)
    (t0, t1, wc) = MShePy.wm.getValues(MShePy.paramTypes.UZ_WC)
    wc[:] = wc_arr_finalts.tolist()
    MShePy.wm.setValues(wc)
    MShePy.wm.performTimeStep()

MShePy.wm.runToTime(end_time)
MShePy.wm.terminate(True)