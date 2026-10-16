import ctypes
import numpy as np
import os
import sys

# assume we are running mshe_watermovement(d).exe. If we are running from python.exe we are screwed!
exe_dir = os.path.dirname(sys.executable)
if "d.exe" in sys.executable:
  lib = ctypes.CDLL(os.path.join(exe_dir, "MShe_wmd.dll"), winmode=0)
else:
  lib = ctypes.CDLL(os.path.join(exe_dir, "MShe_wm.dll"), winmode=0)


#---------------------------------------------------

lib.MSHEWM_GETELEMENTCOUNT.argtypes = [
  ctypes.c_char_p,
  ctypes.POINTER(ctypes.c_int),
  ctypes.c_size_t
]
lib.MSHEWM_GETELEMENTCOUNT.restype = ctypes.c_int

def get_element_count(element_set_id):
  elementSetId = element_set_id.encode()  
  str_len = ctypes.c_size_t(len(elementSetId))
  elementCount = ctypes.c_int(0)
  res = lib.MSHEWM_GETELEMENTCOUNT(elementSetId, ctypes.byref(elementCount), str_len)
  if res == 0:
    raise Exception("Error getting element counts from MIKE-She WM")
  return elementCount.value

#---------------------------------------------------

lib.MSHEWM_ISINTERNALCELL.argtypes = [
  ctypes.POINTER(ctypes.c_int)
]
lib.MSHEWM_ISINTERNALCELL.restype = ctypes.c_int

def is_internal(i_cell):
  iCell = ctypes.c_int(i_cell + 1) # 0-based => 1-based
  int_res = lib.MSHEWM_ISINTERNALCELL(ctypes.byref(iCell))  
  return int_res != 0

#---------------------------------------------------

lib.MSHEWM_GETUZDIMENSIONS.argtypes = [
  ctypes.POINTER(ctypes.c_int),     # nuz      (scalar)
  ctypes.POINTER(ctypes.c_int)      # lluz     (scalar)
]
lib.MSHEWM_GETUZDIMENSIONS.restype = None

def get_uz_dimensions():  
  nuz  = ctypes.c_int(0)
  lluz = ctypes.c_int(0)
  res = lib.MSHEWM_GETUZDIMENSIONS(ctypes.byref(nuz), ctypes.byref(lluz))
  return nuz.value, lluz.value

#---------------------------------------------------

lib.MSHEWM_GETNUMBEROFSZLAYERS.argtypes = [
  ctypes.POINTER(ctypes.c_int)     # nlay      (scalar)
]
lib.MSHEWM_GETNUMBEROFSZLAYERS.restype = None

def get_sz_layers():
  nlay  = ctypes.c_int(0)
  res = lib.MSHEWM_GETNUMBEROFSZLAYERS(ctypes.byref(nlay))
  return nlay.value

#---------------------------------------------------

lib.MSHEWM_GET3DGRIDINFORMATION.argtypes = [
  ctypes.POINTER(ctypes.c_int),     # nuz       (scalar)
  ctypes.POINTER(ctypes.c_int),     # lluz      (scalar)
  ctypes.POINTER(ctypes.c_double),  # dGridX    (nref)
  ctypes.POINTER(ctypes.c_double),  # dGridY    (nref)
  ctypes.POINTER(ctypes.c_double),  # z_top     (nref*nlay)
  ctypes.POINTER(ctypes.c_double),  # zbottom   (nref*nlay)
  ctypes.POINTER(ctypes.c_double),  # dsize     (scalar)
  ctypes.POINTER(ctypes.c_double)   # dorient   (scalar)
]
lib.MSHEWM_GET3DGRIDINFORMATION.restype = ctypes.c_int

def get_3d_gridinfo():
  nref = get_element_count("BaseGrid")
  nlay = get_sz_layers()
  nref_c    = ctypes.c_int(nref)
  nlay_c    = ctypes.c_int(nlay)
  nref_nlay = nref * nlay
  dGridX_c  = (ctypes.c_double * nref)()
  dGridY_c  = (ctypes.c_double * nref)()
  z_top_c   = (ctypes.c_double * nref_nlay)()
  z_bot_c   = (ctypes.c_double * nref_nlay)()
  dsize_c   = ctypes.c_double(0.0)
  dOrnt_c   = ctypes.c_double(0.0)
  res = lib.MSHEWM_GET3DGRIDINFORMATION(ctypes.byref(nref_c), ctypes.byref(nlay_c), dGridX_c, dGridY_c, z_top_c, z_bot_c, ctypes.byref(dsize_c), ctypes.byref(dOrnt_c))
  if res == 0:
    raise Exception("Error getting 3D grid info from MIKE-She WM")
  z_top = np.array(z_top_c).reshape((nlay, nref)) # C order: cell fastest, layer slowest, top layer first
  z_bot = np.array(z_bot_c).reshape((nlay, nref)) # C order: cell fastest, layer slowest, top layer first
  return z_top, z_bot

#---------------------------------------------------

lib.MSHEWM_GEWMUZ3DGRIDINFORMATION.argtypes = [
  ctypes.POINTER(ctypes.c_int),     # nuz       (scalar)
  ctypes.POINTER(ctypes.c_int),     # lluz      (scalar)
  ctypes.POINTER(ctypes.c_double),  # dGridX    (nuz)
  ctypes.POINTER(ctypes.c_double),  # dGridY    (nuz)
  ctypes.POINTER(ctypes.c_double),  # z_top      (nuz*lluz)
  ctypes.POINTER(ctypes.c_double),  # zbottom   (nuz*lluz)
  ctypes.POINTER(ctypes.c_double),  # dsize     (scalar)
  ctypes.POINTER(ctypes.c_double)   # dorient   (scalar)
]
lib.MSHEWM_GEWMUZ3DGRIDINFORMATION.restype = ctypes.c_int

def get_wmuz_3d_gridinfo():
  nuz, lluz = get_uz_dimensions()
  nuz_c     = ctypes.c_int(nuz)
  lluz_c    = ctypes.c_int(lluz)
  nuzluz    = nuz * lluz
  dGridX_c  = (ctypes.c_double * nuz)()
  dGridY_c  = (ctypes.c_double * nuz)()
  z_top_c   = (ctypes.c_double * nuzluz)() # not sure why, but using nuz * lluz gives an error
  z_bot_c   = (ctypes.c_double * nuzluz)() # not sure why, but using nuz * lluz gives an error
  dsize_c   = ctypes.c_double(0.0)
  dOrnt_c   = ctypes.c_double(0.0)
  res = lib.MSHEWM_GEWMUZ3DGRIDINFORMATION(ctypes.byref(nuz_c), ctypes.byref(lluz_c), dGridX_c, dGridY_c, z_top_c, z_bot_c, ctypes.byref(dsize_c), ctypes.byref(dOrnt_c))
  if res == 0:
    raise Exception("Error getting UZ grid info from MIKE-She WM")
  z_top = np.array(z_top_c).reshape((lluz, nuz), order='F') # 'F': Fortran order (z, y, x)
  z_bot = np.array(z_bot_c).reshape((lluz, nuz), order='F') # 'F': Fortran order (z, y, x)
  return z_top, z_bot

#---------------------------------------------------

lib.MSHEWM_GETVALUES.argtypes = [
  ctypes.c_char_p,
  ctypes.c_char_p,
  ctypes.POINTER(ctypes.c_int),
  ctypes.POINTER(ctypes.c_double),
  ctypes.POINTER(ctypes.c_double),
  ctypes.POINTER(ctypes.c_double),
  ctypes.c_size_t,
  ctypes.c_size_t
]
lib.MSHEWM_GETVALUES.restype = ctypes.c_int

def get_values(quantity_ID, element_set_id):
  quantityID = quantity_ID.encode()
  q_len = ctypes.c_size_t(len(quantityID))
  elementSetID = element_set_id.encode()
  e_len = ctypes.c_size_t(len(elementSetID))
  n = get_element_count(element_set_id)
  n_c = ctypes.c_int(n)
  values_c = (ctypes.c_double * n)()
  startHrs = ctypes.c_double(0.0)
  endHrs = ctypes.c_double(0.0)
  res = lib.MSHEWM_GETVALUES(quantityID, elementSetID, ctypes.byref(n_c), ctypes.byref(startHrs), ctypes.byref(endHrs), values_c, q_len, e_len)
  if res == 0:
    raise Exception("Error getting values from MIKE-She WM")
  return startHrs.value, endHrs.value, values_c[:]

#---------------------------------------------------