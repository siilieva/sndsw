 ################################################################################
 #    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
 #                                                                              #
 #              This software is distributed under the terms of the             # 
 #         GNU Lesser General Public Licence version 3 (LGPL) version 3,        #  
 #                  copied verbatim in the file "LICENSE"                       #
 ################################################################################
# - Try to find EvtGen instalation
# Once done this will define
#

MESSAGE(STATUS "Looking for EvtGen ...")

FIND_PATH(EVTGEN_INCLUDE_DIR NAMES EvtGen/EvtGen.hh PATHS
  ${SIMPATH}/include/EvtGen
  NO_DEFAULT_PATH
)

FIND_PATH(EVTGEN_LIB_DIR  NAMES libEvtGen.so PATHS
  ${SIMPATH}/lib
  NO_DEFAULT_PATH
)

Find_Path(EVTGENDATA NAMES evt.pdl PATHS
  ${SIMPATH}/share/EvtGen/
)

If(NOT EVTGENDATA)
  Message(STATUS "Could not find EvtGen data files")
EndIf()

if (PYTHIA8_INCLUDE_DIR AND PYTHIA8_LIB_DIR)
   set(PYTHIA8_FOUND TRUE)
endif (PYTHIA8_INCLUDE_DIR AND PYTHIA8_LIB_DIR)

if (PYTHIA8_FOUND)
  if (NOT PYTHIA8_FOUND_QUIETLY)
    MESSAGE(STATUS "Looking for PYTHIA8... - found ${PYTHIA8_LIB_DIR}")
    SET(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${PYTHIA8_LIB_DIR})
  endif (NOT PYTHIA8_FOUND_QUIETLY)
else (PYTHIA8_FOUND)
  if (PYTHIA8_FOUND_REQUIRED)
    message(FATAL_ERROR "Looking for PYTHIA8... - Not found")
  endif (PYTHIA8_FOUND_REQUIRED)
endif (PYTHIA8_FOUND)
