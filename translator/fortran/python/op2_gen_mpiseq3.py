##########################################################################
#
# OpenMP code generator
#
# This routine is called by op2 which parses the input files
#
# It produces a file xxx_kernel.F90 for each kernel,
# plus a master kernel file
#
##########################################################################

import re
import datetime
import os
import glob
import util

def comm(line):
  global file_text, FORTRAN, CPP
  global depth
  if len(line) == 0:
    prefix = ''
  else:
    prefix = ' '*depth
  if len(line) == 0:
    file_text +='\n'
  elif FORTRAN:
    file_text +='! '+line+'\n'
  elif CPP:
    file_text +=prefix+'//'+line+'\n'

def rep(line,m):
  global dims, idxs, typs, indtyps, inddims

  if FORTRAN:
    if m < len(inddims):
      line = re.sub('INDDIM',str(inddims[m]),line)
      line = re.sub('INDTYP',str(indtyps[m]),line)

    line = re.sub('INDARG','ind_arg'+str(m+1),line)
    line = re.sub('DIMS',str(dims[m]),line)
    line = re.sub('ARG','arg'+str(m+1),line)
    line = re.sub('TYP',typs[m],line)
    line = re.sub('IDX',str(int(idxs[m])),line)
  elif CPP:
    line = re.sub('INDDIM',str(inddims[m]),line)
    line = re.sub('INDTYP',str(indtyps[m]),line)

    line = re.sub('INDARG','ind_arg'+str(m),line)
    line = re.sub('DIM',str(dims[m]),line)
    line = re.sub('ARG','arg'+str(m),line)
    line = re.sub('TYP',typs[m],line)
    line = re.sub('IDX',str(int(idxs[m])),line)
  return line

def code(text):
  global file_text, FORTRAN, CPP, g_m
  global depth
  if len(text) == 0:
    file_text += '\n'
    return
  if len(text) == 0:
    prefix = ''
  else:
    prefix = ' '*depth
  if FORTRAN:
    file_text += prefix+rep(text,g_m)+'\n'
  elif CPP:
    file_text += prefix+rep(text,g_m)+'\n'

def code_pre(text):
  global file_text, FORTRAN, CPP, g_m
  if FORTRAN:
    file_text += rep(text,g_m)+'\n'
  elif CPP:
    file_text += rep(text,g_m)+'\n'

def DO(i,start,finish):
  global file_text, FORTRAN, CPP, g_m
  global depth
  if FORTRAN:
    code('DO '+i+' = '+start+', '+finish+'-1, 1')
  elif CPP:
    code('for ( int '+i+'='+start+'; '+i+'<'+finish+'; '+i+'++ ){')
  depth += 2

def FOR(i,start,finish):
  global file_text, FORTRAN, CPP, g_m
  global depth
  if FORTRAN:
    code('DO '+i+' = '+start+', '+finish+'-1')
  elif CPP:
    code('for ( int '+i+'='+start+'; '+i+'<'+finish+'; '+i+'++ ){')
  depth += 2

def ENDDO():
  global file_text, FORTRAN, CPP, g_m
  global depth
  depth -= 2
  if FORTRAN:
    code('END DO')
  elif CPP:
    code('}')

def ENDFOR():
  global file_text, FORTRAN, CPP, g_m
  global depth
  depth -= 2
  if FORTRAN:
    code('END DO')
  elif CPP:
    code('}')

def IF(line):
  global file_text, FORTRAN, CPP, g_m
  global depth
  if FORTRAN:
    code('IF ('+line+') THEN')
  elif CPP:
    code('if ('+ line + ') {')
  depth += 2

def ELSE():
  global file_text, FORTRAN, CPP, g_m
  global depth
  depth -= 2
  if FORTRAN:
    code('ELSE')
  elif CPP:
    code('else {')
  depth += 2

def ENDIF():
  global file_text, FORTRAN, CPP, g_m
  global depth
  depth -= 2
  if FORTRAN:
    code('END IF')
  elif CPP:
    code('}')


def op2_gen_mpiseq3(master, date, consts, kernels, hydra, bookleaf):

  global dims, idxs, typs, indtyps, inddims
  global FORTRAN, CPP, g_m, file_text, depth

  OP_ID   = 1;  OP_GBL   = 2;  OP_MAP = 3;

  OP_READ = 1;  OP_WRITE = 2;  OP_RW  = 3;
  OP_INC  = 4;  OP_MAX   = 5;  OP_MIN = 6;

  accsstring = ['OP_READ','OP_WRITE','OP_RW','OP_INC','OP_MAX','OP_MIN' ]

  grouped = 0

  any_soa = 0
  for nk in range (0,len(kernels)):
    any_soa = any_soa or sum(kernels[nk]['soaflags'])

##########################################################################
#  create new kernel file
##########################################################################

  for nk in range (0,len(kernels)):
    name  = kernels[nk]['name']
    nargs = kernels[nk]['nargs']
    dims  = kernels[nk]['dims']
    maps  = kernels[nk]['maps']
    var   = kernels[nk]['var']
    typs  = kernels[nk]['typs']
    accs  = kernels[nk]['accs']
    idxs  = kernels[nk]['idxs']
    inds  = kernels[nk]['inds']
    soaflags = kernels[nk]['soaflags']
    optflags = kernels[nk]['optflags']
    ninds   = kernels[nk]['ninds']
    inddims = kernels[nk]['inddims']
    indaccs = kernels[nk]['indaccs']
    indtyps = kernels[nk]['indtyps']
    invinds = kernels[nk]['invinds']
    mapnames = kernels[nk]['mapnames']
    invmapinds = kernels[nk]['invmapinds']
    mapinds = kernels[nk]['mapinds']
    nmaps = 0
    if ninds > 0:
      nmaps = max(mapinds)+1

    optidxs = [0]*nargs
    indopts = [-1]*nargs
    nopts = 0
    for i in range(0,nargs):
      if optflags[i] == 1 and maps[i] == OP_ID:
        optidxs[i] = nopts
        nopts = nopts+1
      elif optflags[i] == 1 and maps[i] == OP_MAP:
        if i == invinds[inds[i]-1]: #i.e. I am the first occurence of this dat+map combination
          optidxs[i] = nopts
          indopts[inds[i]-1] = i
          nopts = nopts+1
        else:
          optidxs[i] = optidxs[invinds[inds[i]-1]]

#
# set two logicals
#
    j = -1
    for i in range(0,nargs):
      if maps[i] == OP_MAP and (accs[i] == OP_INC or accs[i] == OP_RW):
        j = i
    ind_inc = j >= 0

    j = 0
    for i in range(0,nargs):
      if maps[i] == OP_MAP and accs[i] == OP_RW:
        j = i
    ind_rw = j > 0

    j = -1
    for i in range(0,nargs):
      if maps[i] == OP_GBL and accs[i] != OP_READ:
        j = i
    reduct = j >= 0

    reproducible=util.reproducible
    repr_temp_array=util.repr_temp_array
    repr_coloring=util.repr_coloring
    trivial_coloring=util.trivial_coloring
    repro_reduction = -1
    repro_loop = reproducible
    if reproducible:
      for i in range(0,nargs):
        if maps[i] == OP_GBL and accs[i] == OP_INC:
          if typs[i] == 'real(8)' or typs[i] == 'REAL(kind=8)' or typs[i] == 'real*8' or typs[i] == 'r8' or typs[i] == 'real(4)' or typs[i] == 'REAL(kind=4)' or typs[i] == 'real*4' or typs[i] == 'r4':
            repro_reduction = i
      repro_loop = ind_inc or ind_rw

    #if we have an indirect RW, we can't do temp array
    if (reproducible==1 and repr_temp_array==1 and ind_rw==1) or (repr_temp_array==1 and repr_coloring==1):
      repr_temp_array=0
      repr_coloring=1

    FORTRAN = 1;
    CPP     = 0;
    g_m = 0;
    file_text = ''
    depth = 0

    needDimList = []
    for g_m in range(0,nargs):
      if (not dims[g_m].isdigit()):# and not (dims[g_m] in ['NPDE','DNTQMU','DNFCROW','1*1']):
        needDimList = needDimList + [g_m]

##########################################################################
#  Generate Header
##########################################################################
    if hydra:
      code('MODULE '+kernels[nk]['master_file']+'_'+kernels[nk]['mod_file'][9:]+'_module_MODULE')
    else:
      code('MODULE '+name.upper()+'_MODULE')
    code('USE OP2_FORTRAN_DECLARATIONS')
    code('USE OP2_FORTRAN_RT_SUPPORT')
    code('USE ISO_C_BINDING')
    if hydra == 0 and bookleaf == 0:
      code('USE OP2_CONSTANTS')
    if bookleaf:
      code('USE kinds_mod,    ONLY: ink,rlk')
      code('USE parameters_mod,ONLY: LI')

    code('')

##########################################################################
#  Inline user kernel function
##########################################################################
    code('')
    code('CONTAINS')
    code('')

    if hydra == 1:
      file_text += '!DEC$ ATTRIBUTES FORCEINLINE :: ' + name + '\n'
      modfile = kernels[nk]['mod_file'][9:]+'_module'
      filename = 'kernels/'+kernels[nk]['master_file']+'_'+name+'.inc'
      if not os.path.isfile(filename):
        files = [f for f in glob.glob('kernels/*'+name+'.inc')]
        if len(files)>0:
          filename = files[0]
        else:
          print('kernel for '+name+' not found')
#      modfile = modfile.replace('INIT_INIT','INIT')
#      name2 = name.replace('INIT_INIT','INIT')
#      filename = modfile.split('_')[1].lower() + '/' + modfile.split('_')[0].lower() + '/' + name2 + '.F95'
#      if not os.path.isfile(filename):
#        filename = modfile.split('_')[1].lower() + '/' + modfile.split('_')[0].lower() + '/' + name + '.F95'
#      if not os.path.isfile(filename):
#        filename = modfile.split('_')[1].lower() + '/' + modfile.split('_')[0].lower() + '/' + name2[:-1] + '.F95'
      fid = open(filename, 'r')
      text = fid.read()
      fid.close()
      text = text.replace('recursive subroutine','subroutine')
      #text = text.replace('module','!module')
      text = text.replace('contains','!contains')
      text = text.replace('end !module','!end module')

      #
      # substitute npdes with DNPDE
      #
#      using_npdes = 0
#      for g_m in range(0,nargs):
#        if var[g_m] == 'npdes':
#          using_npdes = 1
#      if using_npdes:
#        i = re.search('\\bnpdes\\b',text)
#        j = i.start()
#        i = re.search('\\bnpdes\\b',text[j:])
#        j = j + i.start()+5
#        i = re.search('\\bnpdes\\b',text[j:])
#        j = j + i.start()+5
#        text = text[1:j] + re.sub('\\bnpdes\\b','NPDE',text[j:])
#
      file_text += text
      file_text += '\n#undef MIN\n'
      #code(kernels[nk]['mod_file'])
    elif bookleaf == 1:
      file_text += '!DEC$ ATTRIBUTES FORCEINLINE :: ' + name + '\n'
      modfile = kernels[nk]['mod_file']
      prefixes=['./','ale/','utils/','io/','eos/','hydro/','mods/']
      prefix_i=0
      while (prefix_i<7 and (not os.path.exists(prefixes[prefix_i]+modfile))):
        prefix_i=prefix_i+1
      fid = open(prefixes[prefix_i]+modfile, 'r')
      text = fid.read()
      i = re.search('SUBROUTINE '+name+'\\b',text).start() #text.find('SUBROUTINE '+name)
      j = i + 10 + text[i+10:].find('SUBROUTINE '+name) + 11 + len(name)
      file_text += text[i:j]+'\n\n'
    else:
      comm('user function')
      code('#include "'+name+'.inc"')
      code('')

    code('')

##########################################################################
#  Generate wrapper to iterate over set
##########################################################################

    code('SUBROUTINE op_wrap_'+name+'( &')
    depth = depth + 2
    for g_m in range(0,ninds):
      if invinds[g_m] in needDimList:
        code('& opDat'+str(invinds[g_m]+1)+'Dim, &')
      code('& opDat'+str(invinds[g_m]+1)+'Local, &')
    for g_m in range(0,nargs):
      if maps[g_m] != OP_MAP:
        if g_m in needDimList:
          code('& opDat'+str(g_m+1)+'Dim, &')
        code('& opDat'+str(g_m+1)+'Local, &')
    if nmaps > 0:
      k = []
      for g_m in range(0,nargs):
        if maps[g_m] == OP_MAP and (not mapnames[g_m] in k):
          k = k + [mapnames[g_m]]
          code('& opDat'+str(invinds[inds[g_m]-1]+1)+'Map, &')
          code('& opDat'+str(invinds[inds[g_m]-1]+1)+'MapDim, &')
    if repro_reduction>=0:
      code('& exec_size, &') 
    if ind_inc and repro_loop and repr_coloring:
      code('& color_based_exec_row_starts, color_based_exec, num_colors)')
    elif ind_inc and repro_loop and repr_temp_array:
      for g_m in range(0,ninds):
        if indaccs[g_m] == OP_INC:
          code('& opDat'+str(invinds[g_m]+1)+'Tmp, &')
      code('& row_start_idx, reversed_map, to_size, bottom, top)')
    else:
      code('& bottom,top)')
    code('implicit none')
    for g_m in range(0,ninds):
      if invinds[g_m] in needDimList:
        code('INTEGER(kind=4) opDat'+str(invinds[g_m]+1)+'Dim')
        code(typs[invinds[g_m]]+' opDat'+str(invinds[g_m]+1)+'Local(opDat'+str(invinds[g_m]+1)+'Dim,*)')
      else:
        code(typs[invinds[g_m]]+' opDat'+str(invinds[g_m]+1)+'Local('+str(dims[invinds[g_m]])+',*)')
    if repro_reduction>=0:
      code('INTEGER(kind=4) exec_size') 
    for g_m in range(0,nargs):
      if maps[g_m] != OP_MAP:
        if g_m in needDimList:
          code('INTEGER(kind=4) opDat'+str(g_m+1)+'Dim')
      if maps[g_m] == OP_ID:
        if g_m in needDimList:
          code(typs[g_m]+' opDat'+str(g_m+1)+'Local(opDat'+str(g_m+1)+'Dim,*)')
        else:
          code(typs[g_m]+' opDat'+str(g_m+1)+'Local('+str(dims[g_m])+',*)')
      elif maps[g_m] == OP_GBL:
        repro_mult = ''
        if repro_reduction>=0 and accs[g_m] == OP_INC and (typs[g_m] == 'real(8)' or typs[g_m] == 'REAL(kind=8)' or typs[g_m] == 'real*8' or typs[g_m] == 'r8' or typs[g_m] == 'real(4)' or typs[g_m] == 'REAL(kind=4)' or typs[g_m] == 'real*4' or typs[g_m] == 'r4'):
          repro_mult = ',exec_size'
        if g_m in needDimList:
          code(typs[g_m]+' opDat'+str(g_m+1)+'Local(opDat'+str(g_m+1)+'Dim'+repro_mult+')')
        else:
          code(typs[g_m]+' opDat'+str(g_m+1)+'Local('+str(dims[g_m])+repro_mult+')')
    if nmaps > 0:
      k = []
      for g_m in range(0,nargs):
        if maps[g_m] == OP_MAP and (not mapnames[g_m] in k):
          k = k + [mapnames[g_m]]
          code('INTEGER(kind=4) opDat'+str(invinds[inds[g_m]-1]+1)+'Map(*)')
          code('INTEGER(kind=4) opDat'+str(invinds[inds[g_m]-1]+1)+'MapDim')

    if ind_inc and repro_loop and repr_coloring:
      code('INTEGER(kind=4) num_colors,i1,i2,c')
      code('integer(kind=4) color_based_exec_row_starts(*), color_based_exec(*)')
    elif ind_inc and repro_loop and repr_temp_array:
      code('integer(kind=4) :: row_start_idx(*), reversed_map(*), to_size')
      code('INTEGER(kind=4) bottom,top,i1,i2,d')
      for g_m in range(0, ninds):
        if indaccs[g_m] == OP_INC:
          code(typs[invinds[g_m]]+' opDat'+str(invinds[g_m]+1)+'Tmp(*)')
    else:
      code('INTEGER(kind=4) bottom,top,i1')
    if nmaps > 0:
      k = []
      line = 'INTEGER(kind=4) '
      for g_m in range(0,nargs):
        if maps[g_m] == OP_MAP and (not mapinds[g_m] in k):
          k = k + [mapinds[g_m]]
          line += 'map'+str(mapinds[g_m]+1)+'idx, '
      code(line[:-2])
    code('')
#    if ind_inc == 0 and reduct == 0:
#      code('!DIR$ simd')
    if ind_inc and repro_loop and repr_coloring:
      if trivial_coloring:
        DO('i2','0','num_colors')
      else:
        DO('c','0','num_colors')
        DO('i2','color_based_exec_row_starts(c+1)','color_based_exec_row_starts(c+2)')
      code('i1 = color_based_exec(i2+1)')
    else:
      DO('i1','bottom','top')
    k = []
    for g_m in range(0,nargs):
      if maps[g_m] == OP_MAP and (not mapinds[g_m] in k):
        k = k + [mapinds[g_m]]
        code('map'+str(mapinds[g_m]+1)+'idx = opDat'+str(invmapinds[inds[g_m]-1]+1)+'Map(1 + i1 * opDat'+str(invmapinds[inds[g_m]-1]+1)+'MapDim + '+str(int(idxs[g_m])-1)+')+1')
    comm('kernel call')
    for g_m in range(0,nargs):
      if repro_reduction>=0 and maps[g_m] == OP_GBL and accs[g_m] == OP_INC and (typs[g_m] == 'real(8)' or typs[g_m] == 'REAL(kind=8)' or typs[g_m] == 'real*8' or typs[g_m] == 'r8' or typs[g_m] == 'real(4)' or typs[g_m] == 'REAL(kind=4)' or typs[g_m] == 'real*4' or typs[g_m] == 'r4'):
        code('opDat'+str(g_m+1)+'Local(:,i1+1) = 0.0')
    if ind_inc and repro_loop and repr_temp_array:
      for g_m in range(0, ninds):
        if invinds[g_m] in needDimList:
          datdim = 'opDat'+str(invinds[g_m]+1)+'Dim'
        else:
          datdim = str(inddims[g_m])
        if indaccs[g_m] == OP_INC:
          code('opDat'+str(invinds[g_m]+1)+'Tmp(i1*'+datdim+'*opDat'+str(invmapinds[g_m]+1)+'MapDim+1:(i1+1)*'+datdim+'*opDat'+str(invmapinds[g_m]+1)+'MapDim) = 0.0_8')
    line = 'CALL '+name+'( &'
    indent = '\n'+' '*depth
    for g_m in range(0,nargs):
      if maps[g_m] == OP_ID:
        line = line + indent + '& opDat'+str(g_m+1)+'Local(1,i1+1)'
      if maps[g_m] == OP_MAP:
        if ind_inc and repro_loop and repr_temp_array and accs[g_m] == OP_INC:
          if g_m in needDimList:
            datdim = 'opDat'+str(g_m+1)+'Dim' 
          else:
            datdim = str(dims[g_m])
          line = line + indent + '& opDat'+str(invinds[inds[g_m]-1]+1)+'Tmp((i1*opDat'+str(invmapinds[inds[g_m]-1]+1)+'MapDim + '+str(int(idxs[g_m])-1)+')*'+datdim+'+1)'
        else:
          line = line +indent + '& opDat'+str(invinds[inds[g_m]-1]+1)+'Local(1,map'+str(mapinds[g_m]+1)+'idx)'
      if maps[g_m] == OP_GBL:
        if repro_reduction>=0 and accs[g_m] == OP_INC and (typs[g_m] == 'real(8)' or typs[g_m] == 'REAL(kind=8)' or typs[g_m] == 'real*8' or typs[g_m] == 'r8' or typs[g_m] == 'real(4)' or typs[g_m] == 'REAL(kind=4)' or typs[g_m] == 'real*4' or typs[g_m] == 'r4'):
          line = line + indent +'& opDat'+str(g_m+1)+'Local(1,i1+1)'
        else:
          line = line + indent +'& opDat'+str(g_m+1)+'Local(1)'
      if g_m < nargs-1:
        line = line +', &'
      else:
         line = line +' &'
    depth = depth - 2
    code(line + indent + '& )')
    depth = depth + 2

    ENDDO()
    if ind_inc and repro_loop and repr_coloring and not trivial_coloring:
      ENDDO()

    if ind_inc and repro_loop and repr_temp_array:
      for g_m in range(0,ninds):
        if indaccs[g_m] == OP_INC:
          DO('i1','0','to_size')
          DO('i2','row_start_idx(i1+1)','row_start_idx(i1+2)')
          if g_m in needDimList:
            d = 'opDat'+str(invinds[g_m]+1)+'Dim'
          else:
            d = str(inddims[g_m])
          DO('d','0',d)
          code('opDat'+str(invinds[g_m]+1)+'Local(d+1,i1+1) = opDat'+str(invinds[g_m]+1)+'Local(d+1,i1+1) + &')
          code('    & opDat'+str(invinds[g_m]+1)+'Tmp(reversed_map(i2+1) * '+d+' + 1 + d)')
          ENDDO()
          ENDDO()
          ENDDO()
    depth = depth - 2
    code('END SUBROUTINE')

##########################################################################
#  Generate SEQ host stub
##########################################################################
    code('SUBROUTINE '+name+'_host( userSubroutine, set, &'); depth = depth + 2
    for g_m in range(0,nargs):
      if g_m == nargs-1:
        code('& opArg'+str(g_m+1)+' )')
      else:
        code('& opArg'+str(g_m+1)+', &')

    code('')
    code('IMPLICIT NONE')
    code('character(kind=c_char,len=*), INTENT(IN) :: userSubroutine')
    code('type ( op_set ) , INTENT(IN) :: set')
    code('')

    for g_m in range(0,nargs):
      code('type ( op_arg ) , INTENT(IN) :: opArg'+str(g_m+1))
    code('')

    code('type ( op_arg ) , DIMENSION('+str(nargs)+') :: opArgArray')
    code('INTEGER(kind=4) :: numberOfOpDats')
    code('INTEGER(kind=4), DIMENSION(1:8) :: timeArrayStart')
    code('INTEGER(kind=4), DIMENSION(1:8) :: timeArrayEnd')
    code('REAL(kind=8) :: startTime')
    code('REAL(kind=8) :: endTime')
    code('INTEGER(kind=4) :: returnSetKernelTiming')
    code('INTEGER(kind=4) :: n_upper')
    code('type ( op_set_core ) , POINTER :: opSetCore')
    code('')
    for g_m in range(0,ninds):
      code('INTEGER(kind=4), POINTER, DIMENSION(:) :: opDat'+str(invinds[g_m]+1)+'Map')
      code('INTEGER(kind=4) :: opDat'+str(invinds[g_m]+1)+'MapDim')
      code(typs[invinds[g_m]]+', POINTER, DIMENSION(:) :: opDat'+str(invinds[g_m]+1)+'Local')
      code('INTEGER(kind=4) :: opDat'+str(invinds[g_m]+1)+'Cardinality')
      code('')
    for g_m in range(0,nargs):
      if maps[g_m] == OP_ID:
        code(typs[g_m]+', POINTER, DIMENSION(:) :: opDat'+str(g_m+1)+'Local')
        code('INTEGER(kind=4) :: opDat'+str(g_m+1)+'Cardinality')
        code('')
      if maps[g_m] == OP_GBL:
        code(typs[g_m]+', POINTER, DIMENSION(:) :: opDat'+str(g_m+1)+'Local')

    for g_m in range(0,nargs):
      if maps[g_m] == OP_MAP and optflags[g_m]==1:
        code(typs[g_m]+', POINTER, DIMENSION(:) :: opDat'+str(g_m+1)+'OptPtr')


    if ind_inc and repro_loop and repr_coloring:
      code('type(op_reversed_map_core) :: rev_map')
      code('INTEGER(kind=4), POINTER, DIMENSION(:) :: color_based_exec_row_starts')
      code('INTEGER(kind=4), POINTER, DIMENSION(:) :: color_based_exec')
    elif ind_inc and repro_loop and repr_temp_array:
      code('type(op_reversed_map_core) :: rev_map')
      code('INTEGER(kind=4), POINTER, DIMENSION(:) :: row_start_idx')
      code('INTEGER(kind=4), POINTER, DIMENSION(:) :: reversed_map')
      code('INTEGER(kind=4) :: repr_from_size, repr_to_size')
      for i in range(0,nargs):
        if maps[i] == OP_MAP and accs[i] == OP_INC:
          code(typs[i]+', POINTER, DIMENSION(:) :: opDat'+str(i+1)+'TmpPtr')
    code('')
    code('INTEGER(kind=4) :: i1')
    code('REAL(kind=4) :: dataTransfer')

    code('')
    code('numberOfOpDats = '+str(nargs))
    code('')

    for g_m in range(0,nargs):
      code('opArgArray('+str(g_m+1)+') = opArg'+str(g_m+1))
    code('')

    code('returnSetKernelTiming = setKernelTime('+str(nk)+' , userSubroutine//C_NULL_CHAR, &')
    code('& 0.0_8, 0.00000_4,0.00000_4, 0)')

    code('call op_timers_core(startTime)')
    code('')
    #mpi halo exchange call
    if grouped:
      code('n_upper = op_mpi_halo_exchanges_grouped(set%setCPtr,numberOfOpDats,opArgArray,1)')
    else:
      code('n_upper = op_mpi_halo_exchanges(set%setCPtr,numberOfOpDats,opArgArray)')

    if repro_reduction>=0:
      code('call prepareScratch(opArgArray, numberOfOpDats, n_upper)')

    code('')
    code('opSetCore => set%setPtr')
    code('')
    if ind_inc and repro_loop and repr_coloring:
      j = 0
      for i in range(0,nargs):
        if maps[i] == OP_MAP and (accs[i] == OP_INC or accs[i] == OP_RW):
          j=i
          break
      code('rev_map = op_get_reversed_map(opArg'+str(j+1)+')')
      code('call c_f_pointer(rev_map%color_based_exec_row_starts, color_based_exec_row_starts, (/rev_map%number_of_colors+1/))')
      code('call c_f_pointer(rev_map%color_based_exec, color_based_exec, (/opSetCore%size+opSetCore%exec_size/))')
    code('')

    if ind_inc and repro_loop and repr_temp_array:
      j = 0
      for i in range(0,nargs):
        if maps[i] == OP_MAP and (accs[i] == OP_INC or accs[i] == OP_RW):
          j=i
          break
      code('rev_map = op_get_reversed_map(opArg'+str(j+1)+')')
      code('call c_f_pointer(rev_map%row_start_idx, row_start_idx, (/getSetSizeFromOpArg(opArg'+str(j+1)+')+1/))')
      code('call c_f_pointer(rev_map%reversed_map, reversed_map, (/row_start_idx(getSetSizeFromOpArg(opArg'+str(j+1)+'))/))')
      for i in range(0,nargs):
        if maps[i] == OP_MAP and accs[i] == OP_INC:
          code('call c_f_pointer(op_get_reproducible_tmparray(opArg'+str(i+1)+',repr_from_size, repr_to_size),opDat'+str(i+1)+'TmpPtr, (/n_upper * getMapDimFromOpArg(opArg'+str(i+1)+') * opArg'+str(i+1)+'%dim/))')
      code('repr_to_size = getSetSizeFromOpArg(opArg'+str(j+1)+')')
          

    for g_m in range(0,ninds):
      code('opDat'+str(invinds[g_m]+1)+'Cardinality = opArg'+str(invinds[g_m]+1)+'%dim * getSetSizeFromOpArg(opArg'+str(invinds[g_m]+1)+')')
      code('opDat'+str(invinds[g_m]+1)+'MapDim = getMapDimFromOpArg(opArg'+str(invinds[g_m]+1)+')')
    for g_m in range(0,nargs):
      if maps[g_m] == OP_ID:
        code('opDat'+str(g_m+1)+'Cardinality = opArg'+str(g_m+1)+'%dim * getSetSizeFromOpArg(opArg'+str(g_m+1)+')')

    for g_m in range(0,ninds):
      code('CALL c_f_pointer(opArg'+str(invinds[g_m]+1)+'%data,opDat'+str(invinds[g_m]+1)+'Local,(/opDat'+str(invinds[g_m]+1)+'Cardinality/))')
      code('CALL c_f_pointer(opArg'+str(invinds[g_m]+1)+'%map_data,opDat'+str(invinds[g_m]+1)+'Map,(/opSetCore%size*opDat'+str(invinds[g_m]+1)+'MapDim/))')
    for g_m in range(0,nargs):
      if maps[g_m] == OP_ID:
        code('CALL c_f_pointer(opArg'+str(g_m+1)+'%data,opDat'+str(g_m+1)+'Local,(/opDat'+str(g_m+1)+'Cardinality/))')
      elif maps[g_m] == OP_GBL:
        if repro_reduction>=0 and accs[g_m] == OP_INC and (typs[g_m] == 'real(8)' or typs[g_m] == 'REAL(kind=8)' or typs[g_m] == 'real*8' or typs[g_m] == 'r8' or typs[g_m] == 'real(4)' or typs[g_m] == 'REAL(kind=4)' or typs[g_m] == 'real*4' or typs[g_m] == 'r4'):
          code('CALL c_f_pointer(opArgArray('+str(g_m+1)+')%data_d,opDat'+str(g_m+1)+'Local, (/opArg'+str(g_m+1)+'%dim*n_upper/))')
        else:
          code('CALL c_f_pointer(opArg'+str(g_m+1)+'%data,opDat'+str(g_m+1)+'Local, (/opArg'+str(g_m+1)+'%dim/))')
    code('')

    code('')
    if not (repro_loop and (repr_coloring or repr_temp_array)):
      code('CALL op_wrap_'+name+'( &')
      for g_m in range(0,ninds):
        if invinds[g_m] in needDimList:
          code('& opArg'+str(invinds[g_m]+1)+'%dim, &')
        code('& opDat'+str(invinds[g_m]+1)+'Local, &')
      for g_m in range(0,nargs):
        if maps[g_m] != OP_MAP:
          if g_m in needDimList:
            code('& opArg'+str(g_m+1)+'%dim, &')
          code('& opDat'+str(g_m+1)+'Local, &')

      if nmaps > 0:
        k = []
        for g_m in range(0,nargs):
          if maps[g_m] == OP_MAP and (not mapnames[g_m] in k):
            k = k + [mapnames[g_m]]
            code('& opDat'+str(invinds[inds[g_m]-1]+1)+'Map, &')
            code('& opDat'+str(invinds[inds[g_m]-1]+1)+'MapDim, &')
      if repro_reduction>=0:
        code('& n_upper, &') 
      code('& 0, opSetCore%core_size)')
    if grouped:
      code('CALL op_mpi_wait_all_grouped(numberOfOpDats,opArgArray,1)')
    else:
      code('CALL op_mpi_wait_all(numberOfOpDats,opArgArray)')
    code('CALL op_wrap_'+name+'( &')
    for g_m in range(0,ninds):
      if invinds[g_m] in needDimList:
          code('& opArg'+str(invinds[g_m]+1)+'%dim, &')
      code('& opDat'+str(invinds[g_m]+1)+'Local, &')
    for g_m in range(0,nargs):
      if maps[g_m] != OP_MAP:
        if g_m in needDimList:
            code('& opArg'+str(g_m+1)+'%dim, &')
        code('& opDat'+str(g_m+1)+'Local, &')

    if nmaps > 0:
      k = []
      for g_m in range(0,nargs):
        if maps[g_m] == OP_MAP and (not mapnames[g_m] in k):
          k = k + [mapnames[g_m]]
          code('& opDat'+str(invinds[inds[g_m]-1]+1)+'Map, &')
          code('& opDat'+str(invinds[inds[g_m]-1]+1)+'MapDim, &')
    #code('& 0, n_upper)')
    if repro_reduction>=0:
      code('& n_upper, &') 
    if repro_loop and repr_coloring:
      if trivial_coloring:
        code('& color_based_exec_row_starts, color_based_exec, n_upper)')
      else:
        code('& color_based_exec_row_starts, color_based_exec, rev_map%number_of_colors)')
    elif repro_loop and repr_temp_array:
      for g_m in range(0,ninds):
        if indaccs[g_m] == OP_INC:
          code('& opDat'+str(invinds[g_m]+1)+'TmpPtr, &')
      code('& row_start_idx, reversed_map, repr_to_size, 0, n_upper)')
    else:
      code('& opSetCore%core_size, n_upper)')


    IF('(n_upper .EQ. 0) .OR. (n_upper .EQ. opSetCore%core_size)')
    if grouped:
      code('CALL op_mpi_wait_all_grouped(numberOfOpDats,opArgArray,1)')
    else:
      code('CALL op_mpi_wait_all(numberOfOpDats,opArgArray)')
    ENDIF()
    code('')



    code('')
    code('CALL op_mpi_set_dirtybit(numberOfOpDats,opArgArray)')
    code('')

    #reductions
    for g_m in range(0,nargs):
        if repro_reduction>=0 and maps[g_m] == OP_GBL and accs[g_m] == OP_INC and (typs[g_m] == 'real(8)' or typs[g_m] == 'REAL(kind=8)' or typs[g_m] == 'real*8' or typs[g_m] == 'r8' or typs[g_m] == 'real(4)' or typs[g_m] == 'REAL(kind=4)' or typs[g_m] == 'real*4' or typs[g_m] == 'r4'):
          code('call op_reproducible_sum(opArg'+str(g_m+1)+', opSetCore%size, opDat'+str(g_m+1)+'Local)')
        elif maps[g_m] == OP_GBL and (accs[g_m] == OP_INC or accs[g_m] == OP_MIN or accs[g_m] == OP_MAX or accs[g_m] == OP_WRITE):
          if optflags[g_m] == 1:
            IF('opArg'+str(g_m+1)+'%opt == 1')
          if typs[g_m] == 'real(8)' or typs[g_m] == 'REAL(kind=8)' or typs[g_m] == 'real*8' or typs[g_m] == 'r8':
            code('CALL op_mpi_reduce_double(opArg'+str(g_m+1)+',opArg'+str(g_m+1)+'%data)')
          elif typs[g_m] == 'real(4)' or typs[g_m] == 'REAL(kind=4)' or typs[g_m] == 'real*4' or typs[g_m] == 'r4':
            code('CALL op_mpi_reduce_float(opArg'+str(g_m+1)+',opArg'+str(g_m+1)+'%data)')
          elif typs[g_m] == 'integer(4)' or typs[g_m] == 'INTEGER(kind=4)' or typs[g_m] == 'integer*4' or typs[g_m] == 'i4':
            code('CALL op_mpi_reduce_int(opArg'+str(g_m+1)+',opArg'+str(g_m+1)+'%data)')
          elif typs[g_m] == 'logical' or typs[g_m] == 'logical*1':
            code('CALL op_mpi_reduce_bool(opArg'+str(g_m+1)+',opArg'+str(g_m+1)+'%data)')
          else:
            print('Error, reduction type '+typs[g_m]+' unrecognised')
          code('')
          if optflags[g_m] == 1:
            ENDIF()

    code('call op_timers_core(endTime)')
    code('')
    code('dataTransfer = 0.0')
    if ninds == 0:
      for g_m in range(0,nargs):
        if optflags[g_m] == 1:
          IF('opArg'+str(g_m+1)+'%opt == 1')
        if accs[g_m] == OP_READ or accs[g_m] == OP_WRITE:
          if maps[g_m] == OP_GBL:
            code('dataTransfer = dataTransfer + opArg'+str(g_m+1)+'%size')
          else:
            code('dataTransfer = dataTransfer + opArg'+str(g_m+1)+'%size * opSetCore%size')
        else:
          if maps[g_m] == OP_GBL:
            code('dataTransfer = dataTransfer + opArg'+str(g_m+1)+'%size * 2.d0')
          else:
            code('dataTransfer = dataTransfer + opArg'+str(g_m+1)+'%size * opSetCore%size * 2.d0')
        if optflags[g_m] == 1:
          ENDIF()
    else:
      names = []
      for g_m in range(0,ninds):
        mult=''
        if indaccs[g_m] != OP_WRITE and indaccs[g_m] != OP_READ:
          mult = ' * 2.d0'
        if not var[invinds[g_m]] in names:
          if optflags[invinds[g_m]] == 1:
            IF('opArg'+str(g_m+1)+'%opt == 1')
          code('dataTransfer = dataTransfer + opArg'+str(invinds[g_m]+1)+'%size * MIN(n_upper,getSetSizeFromOpArg(opArg'+str(invinds[g_m]+1)+'))'+mult)
          names = names + [var[invinds[g_m]]]
          if optflags[invinds[g_m]] == 1:
            ENDIF()
      for g_m in range(0,nargs):
        mult=''
        if accs[g_m] != OP_WRITE and accs[g_m] != OP_READ:
          mult = ' * 2.d0'
        if not var[g_m] in names:
          if optflags[g_m] == 1:
            IF('opArg'+str(g_m+1)+'%opt == 1')
          names = names + [var[invinds[g_m]]]
          if maps[g_m] == OP_ID:
            code('dataTransfer = dataTransfer + opArg'+str(g_m+1)+'%size * MIN(n_upper,getSetSizeFromOpArg(opArg'+str(g_m+1)+'))'+mult)
          elif maps[g_m] == OP_GBL:
            code('dataTransfer = dataTransfer + opArg'+str(g_m+1)+'%size'+mult)
          if optflags[g_m] == 1:
            ENDIF()
      if nmaps > 0:
        k = []
        for g_m in range(0,nargs):
          if maps[g_m] == OP_MAP and (not mapnames[g_m] in k):
            k = k + [mapnames[g_m]]
            code('dataTransfer = dataTransfer + n_upper * opDat'+str(invinds[inds[g_m]-1]+1)+'MapDim * 4.d0')

    code('returnSetKernelTiming = setKernelTime('+str(nk)+' , userSubroutine//C_NULL_CHAR, &')
    code('& endTime-startTime, dataTransfer, 0.00000_4, 1)')
    #code('returnSetKernelTiming = setKernelTime('+str(nk)+' , userSubroutine//C_NULL_CHAR, &')
    #code('& endTime-startTime,0.00000,0.00000, 1)')
    depth = depth - 2
    code('END SUBROUTINE')
    code('END MODULE')
    code('')

##########################################################################
#  output individual kernel file
##########################################################################
    if hydra:
      name = 'kernels/'+kernels[nk]['master_file']+'_'+name
      fid = open(name+'_seqkernel.F90','w')
    elif bookleaf:
      fid = open(prefixes[prefix_i]+name+'_seqkernel.f90','w')
    else:
      fid = open(name+'_seqkernel.F90','w')
    date = datetime.datetime.now()
    fid.write('!\n! auto-generated by op2.py\n!\n\n')
    fid.write(file_text.strip())
    fid.close()
