import ctypes
import numpy as np
import os
import sys

# assume we are running mshe_watermovement(d).exe. If we are running from python.exe we are screwed!
exe_dir = os.path.dirname(sys.executable)
if "d.exe" in sys.executable:
  lib = ctypes.CDLL(os.path.join(exe_dir, "MShe_wqd.dll"), winmode=0)
else:
  lib = ctypes.CDLL(os.path.join(exe_dir, "MShe_wq.dll"), winmode=0)

LeID =  50 # Copy from WQ MOpenMI.f90, keep up to date!

#---------------------------------------------------

lib.MSHEWQ_GETNUMBEROFSHENODES.argtypes = [
  ctypes.POINTER(ctypes.c_int)
]
lib.MSHEWQ_GETNUMBEROFSHENODES.restype = None

def get_nref():
  nref_c = ctypes.c_int(0)
  res = lib.MSHEWQ_GETNUMBEROFSHENODES(ctypes.byref(nref_c))
  return nref_c.value

#---------------------------------------------------

lib.MSHEWQ_GETNUMBEROFSZLAYERS.argtypes = [
  ctypes.POINTER(ctypes.c_int)     # nlay      (scalar)
]
lib.MSHEWQ_GETNUMBEROFSZLAYERS.restype = None

def get_sz_layers():
  nlay_c  = ctypes.c_int(0)
  res = lib.MSHEWQ_GETNUMBEROFSZLAYERS(ctypes.byref(nlay_c))
  return nlay_c.value

#---------------------------------------------------
lib.MSHEWQ_INIT.argtypes = [
  ctypes.c_char_p,
  ctypes.POINTER(ctypes.c_int),
  ctypes.POINTER(ctypes.c_int),
  ctypes.c_size_t
]
lib.MSHEWQ_INIT.restype = ctypes.c_int

def init(she_path, wm_coup):
  path_encoded = she_path.encode()
  batch = ctypes.c_int(0)
  wmwq = ctypes.c_int(wm_coup)
  str_len = ctypes.c_size_t(len(she_path)) # * 20) # use * 20 to test/break it...
  res = lib.MSHEWQ_INIT(path_encoded, ctypes.byref(batch), ctypes.byref(wmwq), str_len)
  return res != 0

#---------------------------------------------------

lib.MSHEWQ_OPENMI_INIT.argtypes = []
lib.MSHEWQ_OPENMI_INIT.restype = ctypes.c_int

def openmi_init():
  res = 0
  res = lib.MSHEWQ_OPENMI_INIT()
  return res != 0

#---------------------------------------------------

lib.MSHEWQ_GETINPUTEXCHANGEITEMCOUNT.argtypes = []
lib.MSHEWQ_GETINPUTEXCHANGEITEMCOUNT.restype = ctypes.c_int

def get_input_exchange_item_count():
  return lib.MSHEWQ_GETINPUTEXCHANGEITEMCOUNT()

#---------------------------------------------------

lib.MSHEWQ_GETINPUTEXCHANGEITEM.argtypes = [
  ctypes.POINTER(ctypes.c_int),     # (in)  item index
  ctypes.c_char_p,                  # (out) item ID
  ctypes.c_char_p,                  # (out) element set ID
  ctypes.c_char_p,                  # (out) item description
  ctypes.POINTER(ctypes.c_int),     # (out) eumType
  ctypes.POINTER(ctypes.c_int),     # (out) eumUnit
  ctypes.POINTER(ctypes.c_int),     # (out, bool) hidden?
  ctypes.POINTER(ctypes.c_int),     # (out, bool) neededForWmCoupling?
  ctypes.c_size_t,                  # string length 1
  ctypes.c_size_t,                  # string length 2
  ctypes.c_size_t                   # string length 3
]
lib.MSHEWQ_GETINPUTEXCHANGEITEM.restype = ctypes.c_int

def get_input_exchange_item(item):  
  i = ctypes.c_int(item + 1)
  buf1 = ctypes.create_string_buffer(LeID)
  buf2 = ctypes.create_string_buffer(LeID)
  buf3 = ctypes.create_string_buffer(LeID)
  l = ctypes.c_size_t(LeID)
  eum_type   = ctypes.c_int(0)
  eum_unit   = ctypes.c_int(0)
  hidden     = ctypes.c_int(0)
  needed_wm  = ctypes.c_int(0)
  res = lib.MSHEWQ_GETINPUTEXCHANGEITEM(
    ctypes.byref(i),
    buf1, buf2, buf3,
    ctypes.byref(eum_type), ctypes.byref(eum_unit), ctypes.byref(hidden), ctypes.byref(needed_wm),
    l, l, l
    )
  if res == 0:
    raise Exception("Error getting input exchange item from MIKE-She WQ")
  return (
    buf1.value.decode().rstrip(), buf2.value.decode().rstrip(), buf3.value.decode(),
    eum_type.value, eum_unit.value, hidden.value != 0, needed_wm.value != 0
  )

#---------------------------------------------------

lib.MSHEWQ_GETOUTPUTEXCHANGEITEMCOUNT.argtypes = []
lib.MSHEWQ_GETOUTPUTEXCHANGEITEMCOUNT.restype = ctypes.c_int

def get_output_exchange_item_count():
  return lib.MSHEWQ_GETOUTPUTEXCHANGEITEMCOUNT()

#---------------------------------------------------

lib.MSHEWQ_GETOUTPUTEXCHANGEITEM.argtypes = [
  ctypes.POINTER(ctypes.c_int),     # (in)  item index
  ctypes.c_char_p,                  # (out) item ID
  ctypes.c_char_p,                  # (out) element set ID
  ctypes.c_char_p,                  # (out) item description
  ctypes.POINTER(ctypes.c_int),     # (out) eumType
  ctypes.POINTER(ctypes.c_int),     # (out) eumUnit
  ctypes.POINTER(ctypes.c_int),     # (out, bool) hidden?
  ctypes.c_size_t,                  # string length 1
  ctypes.c_size_t,                  # string length 2
  ctypes.c_size_t                   # string length 3
]
lib.MSHEWQ_GETOUTPUTEXCHANGEITEM.restype = ctypes.c_int

def get_output_exchange_item(item):
  i = ctypes.c_int(item + 1)
  buf1 = ctypes.create_string_buffer(LeID)
  buf2 = ctypes.create_string_buffer(LeID)
  buf3 = ctypes.create_string_buffer(LeID)
  l = ctypes.c_size_t(LeID)
  eum_type   = ctypes.c_int(0)
  eum_unit   = ctypes.c_int(0)
  needed_wm  = ctypes.c_int(0)
  res = lib.MSHEWQ_GETOUTPUTEXCHANGEITEM(
    ctypes.byref(i),
    buf1, buf2, buf3,
    ctypes.byref(eum_type), ctypes.byref(eum_unit), ctypes.byref(needed_wm),
    l, l, l
    )
  if res == 0:
    raise Exception("Error getting output exchange item from MIKE-She WQ")
  return (
    buf1.value.decode().rstrip(), buf2.value.decode().rstrip(), buf3.value.decode(),
    eum_type.value, eum_unit.value, needed_wm.value != 0
  )

#---------------------------------------------------

lib.MSHEWQ_PREPARE.argtypes = [
  ctypes.POINTER(ctypes.c_int),
  ctypes.c_char_p,
  ctypes.c_size_t
]
lib.MSHEWQ_PREPARE.restype = ctypes.c_int

def openmi_prepare(used_quantity_ids):  
  nexquantities = ctypes.c_int(len(used_quantity_ids))
  listOfUsedquantityIDs = used_quantity_ids.encode()
  str_len = ctypes.c_size_t(len(listOfUsedquantityIDs))
  res = lib.MSHEWQ_PREPARE(nexquantities, listOfUsedquantityIDs, str_len)
  if res == 0:
    raise Exception("Error preparing MIKE-She WQ exchange interface")

#---------------------------------------------------

lib.MSHEWQ_GETELEMENTCOUNT.argtypes = [
  ctypes.c_char_p,
  ctypes.POINTER(ctypes.c_int),
  ctypes.c_size_t
]
lib.MSHEWQ_GETELEMENTCOUNT.restype = ctypes.c_int

def get_element_count(element_id):
  elementID = element_id.encode()
  str_len = ctypes.c_size_t(len(element_id))
  n = ctypes.c_int(0)
  res = lib.MSHEWQ_GETELEMENTCOUNT(elementID, ctypes.byref(n), str_len)
  if res == 0:
    return None
  return n.value

#---------------------------------------------------

lib.MSHEWQ_GETUZDIMENSIONS.argtypes = [
  ctypes.POINTER(ctypes.c_int),
  ctypes.POINTER(ctypes.c_int)
]
lib.MSHEWQ_GETUZDIMENSIONS.restype = ctypes.c_int

def get_uz_dimensions():
  nUz = ctypes.c_int(0)
  llUz = ctypes.c_int(0)
  res = lib.MSHEWQ_GETUZDIMENSIONS(ctypes.byref(nUz), ctypes.byref(llUz))
  if res == 0:
    raise Exception('Error retrieving UZ dimensions from WQ')
  return nUz.value, llUz.value  

#---------------------------------------------------

lib.MSHEWQ_GETVALUES.argtypes = [
  ctypes.c_char_p,
  ctypes.c_char_p,
  ctypes.POINTER(ctypes.c_int),
  ctypes.POINTER(ctypes.c_double),
  ctypes.POINTER(ctypes.c_double),
  ctypes.POINTER(ctypes.c_double),
  ctypes.c_size_t,
  ctypes.c_size_t
]
lib.MSHEWQ_GETVALUES.restype = ctypes.c_int

def get_values(quantity_ID, element_set_id):
  quantityID = quantity_ID.encode()
  q_len = ctypes.c_size_t(len(quantityID))
  elementSetID = element_set_id.encode()
  e_len = ctypes.c_size_t(len(elementSetID))
  n = get_element_count(element_set_id)
  n_c = ctypes.c_int(n)
  values_c = (ctypes.c_double * n)()
  startHrs = ctypes.c_double(0.0)
  endHrs   = ctypes.c_double(0.0)
  res = lib.MSHEWQ_GETVALUES(quantityID, elementSetID, ctypes.byref(n_c), ctypes.byref(startHrs), ctypes.byref(endHrs), values_c, q_len, e_len)
  if res == 0:
    raise Exception("Error getting values from MIKE-She WQ")
  if(element_set_id == "WMUZ3DGrid"):
    nuz = get_element_count("WMUZ2DGrid")
    lluz = int(n / nuz)
    values = np.array(values_c).reshape((lluz, nuz), order='F') # 'F': Fortran order (z, y, x)
  elif(element_set_id == "BaseGrid"):
    values = np.array(values_c)
  elif(element_set_id == "SZ3DGrid"):    
    nref = get_nref()
    nlay = get_sz_layers()
    values = np.array(values_c).reshape((nlay, nref)) # C order: cell fastest, layer slowest, top layer first
  else:
    raise Exception(f'Not implemented for element_set_id {element_set_id}')
  return startHrs.value, endHrs.value, values


#---------------------------------------------------

lib.MSHEWQ_SETVALUES.argtypes = [
  ctypes.c_char_p,
  ctypes.c_char_p,
  ctypes.POINTER(ctypes.c_int),
  ctypes.POINTER(ctypes.c_double),
  ctypes.c_size_t,
  ctypes.c_size_t
]
lib.MSHEWQ_SETVALUES.restype = ctypes.c_int

def set_values(quantity_ID, element_set_id, values_2d):
  ''' 3d values, flattened 1d/depth'''
  quantityID = quantity_ID.encode()
  q_len = ctypes.c_size_t(len(quantityID))
  elementSetID = element_set_id.encode()
  e_len = ctypes.c_size_t(len(elementSetID))
  n = ctypes.c_int(np.size(values_2d))
  values_1d = np.array(values_2d).reshape(np.size(values_2d), order='F') # 'F': Fortran order (z, y, x)
  values_c = (ctypes.c_double * len(values_1d))(*values_1d)
  res = lib.MSHEWQ_SETVALUES(quantityID, elementSetID, ctypes.byref(n), values_c, q_len, e_len)
  if res == 0:
    raise Exception("Error setting values to MIKE-She WQ")

#---------------------------------------------------

lib.MSHEWQ_PERFORMTIMESTEP.argtypes = [
  ctypes.POINTER(ctypes.c_double),
  ctypes.POINTER(ctypes.c_double)
]
lib.MSHEWQ_PERFORMTIMESTEP.restype = ctypes.c_int

def perform_time_step():
  dtUsedSeconds = ctypes.c_double(0.0)
  relativeSecsFromStart = ctypes.c_double(0.0)
  res = lib.MSHEWQ_PERFORMTIMESTEP(ctypes.byref(dtUsedSeconds), ctypes.byref(relativeSecsFromStart)) != 0 # r8 out, r8 out
  if res == 0:
    raise Exception("Error performing MIKE-She WQ time step")
  return dtUsedSeconds.value / 3600.0, relativeSecsFromStart.value / 3600.0

#---------------------------------------------------

lib.MSHEWQ_SETWMCOUPLINGTIME.argtypes = [
  ctypes.POINTER(ctypes.c_double)
]
lib.MSHEWQ_SETWMCOUPLINGTIME.restype = ctypes.c_int

def set_wm_coupling_time_h(wm_time_h):
  wmTimeSec = ctypes.c_double(wm_time_h * 3600.0)
  res = lib.MSHEWQ_SETWMCOUPLINGTIME(ctypes.byref(wmTimeSec))
  if res == 0:
    raise Exception("Error setting wm-coupling-time to MIKE-She WQ")

#---------------------------------------------------

lib.MSHEWQ_TERMINATE.argtypes = [
  ctypes.POINTER(ctypes.c_int),
  ctypes.POINTER(ctypes.c_int)
]
lib.MSHEWQ_TERMINATE.restype = ctypes.c_int

def terminate(success):
  lNormalTermination = ctypes.c_int(1 if success else 0)
  lNoReturn = ctypes.c_int(0)
  res = lib.MSHEWQ_TERMINATE(ctypes.byref(lNormalTermination), ctypes.byref(lNoReturn)) # in, in
  if success and res == 0: # if not success the function always returns false, it's not another error then!
    raise Exception("Error shutting down MIKE-She WQ")