##########################################################################
#
# CUDA code generator
#
# This routine is called by op2 which parses the input files
#
# It produces a file xxx_kernel.cu for each kernel,
# plus a master kernel file
#
##########################################################################

import re
import datetime
import glob
import os
import op2_gen_common

def comm(line):
  global file_text, FORTRAN, CPP
  global depth
  prefix = ' '*depth
  if len(line) == 0:
    file_text +='\n'
  elif FORTRAN:
    file_text +='!  '+line+'\n'
  elif CPP:
    file_text +=prefix+'//'+line.rstrip()+'\n'

def rep(line,m):
  global dims, idxs, typs, indtyps, inddims
  if m < len(inddims):
    line = re.sub('<INDDIM>',str(inddims[m]),line)
    line = re.sub('<INDTYP>',str(indtyps[m]),line)

  line = re.sub('<INDARG>','ind_arg'+str(m),line)
  line = re.sub('<DIM>',str(dims[m]),line)
  line = re.sub('<ARG>','arg'+str(m),line)
  line = re.sub('<TYP>',typs[m],line)
  line = re.sub('<IDX>',str(int(idxs[m])),line)
  return line

def code(text):
  global file_text, FORTRAN, CPP, g_m
  global depth
  if text == '':
    prefix = ''
  else:
    prefix = ' '*depth
  file_text += prefix+rep(text,g_m).rstrip()+'\n'


def FOR(i,start,finish):
  global file_text, FORTRAN, CPP, g_m
  global depth
  if FORTRAN:
    code('do '+i+' = '+start+', '+finish+'-1')
  elif CPP:
    code('for ( int '+i+'='+start+'; '+i+'<'+finish+'; '+i+'++ ){')
  depth += 2

def FOR_INC(i,start,finish,inc):
  global file_text, FORTRAN, CPP, g_m
  global depth
  if FORTRAN:
    code('do '+i+' = '+start+', '+finish+'-1')
  elif CPP:
    code('for ( int '+i+'='+start+'; '+i+'<'+finish+'; '+i+'+='+inc+' ){')
  depth += 2

def ENDFOR():
  global file_text, FORTRAN, CPP, g_m
  global depth
  depth -= 2
  if FORTRAN:
    code('enddo')
  elif CPP:
    code('}')

def IF(line):
  global file_text, FORTRAN, CPP, g_m
  global depth
  if FORTRAN:
    code('if ('+line+') then')
  elif CPP:
    code('if ('+ line + ') {')
  depth += 2

def ENDIF():
  global file_text, FORTRAN, CPP, g_m
  global depth
  depth -= 2
  if FORTRAN:
    code('endif')
  elif CPP:
    code('}')

def op2_gen_cuda_simple(master, date, consts, kernels,sets, macro_defs):

  global dims, idxs, typs, indtyps, inddims
  global FORTRAN, CPP, g_m, file_text, depth

  OP_ID   = 1;  OP_GBL   = 2;  OP_MAP = 3;

  OP_READ = 1;  OP_WRITE = 2;  OP_RW  = 3;
  OP_INC  = 4;  OP_MAX   = 5;  OP_MIN = 6;

  accsstring = ['OP_READ','OP_WRITE','OP_RW','OP_INC','OP_MAX','OP_MIN' ]



  reproducible=op2_gen_common.reproducible
  repr_temp_array=op2_gen_common.repr_temp_array
  repr_coloring=op2_gen_common.repr_coloring

  if reproducible and repr_temp_array and repr_coloring:
    repr_coloring = 0

  inc_stage=0
  op_color2=0
  op_color2_force=1
##########################################################################
#  create new kernel file
##########################################################################

  for nk in range (0,len(kernels)):

#Optimization settings
    inc_stage=0
    op_color2_force=0
    atomics=1

    if (reproducible and repr_coloring and atomics):
      atomics=0


    name, nargs, dims, maps, var, typs, accs, idxs, inds, soaflags, optflags, decl_filepath, \
           ninds, inddims, indaccs, indtyps, invinds, mapnames, invmapinds, mapinds, nmaps, nargs_novec, \
           unique_args, vectorised, cumulative_indirect_index = op2_gen_common.create_kernel_info(kernels[nk], inc_stage)

    varrrrrr=var

    if reproducible and repr_temp_array:
      for i in range(0,nargs):
        if maps[i] == OP_MAP:
          for g_m in range(0,ninds):
            if indaccs[g_m] == OP_RW:
              print('Warning - OP_RW reproducibility is not supported with temporary array method. Changing to reproducible coloring.')
              repr_temp_array = 0
              repr_coloring = 1
              break
        if repr_coloring:
          break

    mapnames2=[]
    repro_if=0
    repro_prime_map=''
    if reproducible:
      mapnames2 = mapnames[:]
      for i in range(0,len(mapnames)):
        if mapnames[i].find('[')>=0:
          mapnames2[i] = mapnames[i][:mapnames[i].find('[')]

      if reproducible:
        if ninds>0:
          if nmaps > 0:
            k = []
            for g_m in range(0,nargs):
              if accs[g_m] == OP_INC and maps[g_m] == OP_MAP and (not mapnames2[g_m] in k):
                k = k + [mapnames2[g_m]]
                repro_if=1
                repro_prime_map=mapnames2[g_m]
    any_soa = 0
    any_soa = any_soa or sum(soaflags)
    op_color2=0

    
#
# set logicals
#
    j = -1
    for i in range(0,nargs):
      if maps[i] == OP_MAP and accs[i] == OP_INC:
        j = i
    ind_inc = j >= 0

    j = -1
    for i in range(0,nargs):
      if maps[i] == OP_MAP and accs[i] == OP_RW:
        j = i
    ind_rw = j >= 0

    if atomics and ind_rw:
      atomics = 0

    if ind_rw or op_color2_force:
        op_color2 = 1
    else:
        op_color2 = 0

    #no staging with 2 level colouring
    if op_color2:
      inc_stage=0

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

    j = -1
    for i in range(0,nargs):
      if maps[i] == OP_GBL and accs[i] != OP_READ and accs[i] != OP_WRITE:
        j = i
    reduct = j >= 0

    if inc_stage:
      ninds_staged = 0
      inds_staged = [-1]*nargs
      for i in range(0,nargs):
        if maps[i]==OP_MAP and accs[i]==OP_INC:
          if inds_staged[invinds[inds[i]-1]] == -1:
            inds_staged[i] = ninds_staged
            ninds_staged = ninds_staged + 1
          else:
            inds_staged[i] = inds_staged[invinds[inds[i]-1]]
      invinds_staged = [-1]*ninds_staged
      inddims_staged = [-1]*ninds_staged
      indopts_staged = [-1]*ninds_staged
      for i in range(0,nargs):
        if inds_staged[i] >= 0 and invinds_staged[inds_staged[i]] == -1:
          invinds_staged[inds_staged[i]] = i
          inddims_staged[inds_staged[i]] = dims[i]
          if optflags[i] == 1:
            indopts_staged[inds_staged[i]] = i
      for i in range(0,nargs):
        inds_staged[i] = inds_staged[i] + 1

##########################################################################
#  start with CUDA kernel function
##########################################################################

    FORTRAN = 0;
    CPP     = 1;
    g_m = 0;
    file_text = ''
    depth = 0


    #strides for SoA
    if any_soa:
      if nmaps > 0:
        k = []
        for g_m in range(0,nargs):
          if repr_temp_array and repro_if and maps[g_m] == OP_MAP and accs[g_m] == OP_INC and (not (mapnames[g_m]+'tmpa') in k):
            k = k + [(mapnames[g_m]+'tmpa')]
            code('__constant__ int opMap_'+str(mapnames2[g_m])+'_'+name+'_stride_temp_inc_OP2CONSTANT;')
            code('int opMap_'+str(mapnames2[g_m])+'_'+name+'_stride_temp_inc_OP2HOST;')
          if maps[g_m] == OP_MAP and (not mapnames[g_m] in k):
            k = k + [mapnames[g_m]]
            code('__constant__ int opDat'+str(invinds[inds[g_m]-1])+'_'+name+'_stride_OP2CONSTANT;')
            code('int opDat'+str(invinds[inds[g_m]-1])+'_'+name+'_stride_OP2HOST=-1;')
      dir_soa = -1
      for g_m in range(0,nargs):
        if maps[g_m] == OP_ID and ((not dims[g_m].isdigit()) or int(dims[g_m]) > 1):
          code('__constant__ int direct_'+name+'_stride_OP2CONSTANT;')
          code('int direct_'+name+'_stride_OP2HOST=-1;')
          dir_soa = g_m
          break


    file_name = decl_filepath

    f = open(file_name, 'r')
    kernel_text = f.read()
    f.close()

    if CPP:
      includes = op2_gen_common.extract_includes(kernel_text)
      if len(includes) > 0:
        for include in includes:
          code(include)
        code("")

    comm('user function')

    kernel_text = op2_gen_common.comment_remover(kernel_text)
    kernel_text = op2_gen_common.remove_trailing_w_space(kernel_text)

    p = re.compile('void\\s+\\b'+name+'\\b')
    i = p.search(kernel_text).start()

    if(i < 0):
      print("\n********")
      print("Error: cannot locate user kernel function name: "+name+" - Aborting code generation")
      exit(2)
    i2 = i

    #i = kernel_text[0:i].rfind('\n') #reverse find
    j = kernel_text[i:].find('{')
    k = op2_gen_common.para_parse(kernel_text, i+j, '{', '}')
    signature_text = kernel_text[i:i+j]
    l = signature_text[0:].find('(')
    head_text = signature_text[0:l].strip() #save function name
    m = op2_gen_common.para_parse(signature_text, 0, '(', ')')
    signature_text = signature_text[l+1:m]
    body_text = kernel_text[i+j+1:k]

    ## Replace occurrences of '#include "<FILE>"' within loop with the contents of <FILE>:
    body_text = op2_gen_common.replace_local_includes_with_file_contents(body_text, os.path.dirname(master))

    # check for number of arguments
    if len(signature_text.split(',')) != nargs_novec:
        print('Error parsing user kernel('+name+'): must have '+str(nargs)+' arguments')
        return

    for i in range(0,nargs_novec):
        var = signature_text.split(',')[i].strip()
        if (kernels[nk]['soaflags'][i] and ((op_color2 or  not (kernels[nk]['maps'][i] == OP_MAP and kernels[nk]['accs'][i] == OP_INC) ) or ( kernels[nk]['maps'][i] == OP_MAP and kernels[nk]['accs'][i] == OP_INC and repr_temp_array and repro_if) ) ): 
          var = var.replace('*','')
          #locate var in body and replace by adding [idx]
          length = len(re.compile('\\s+\\b').split(var))
          var2 = re.compile('\\s+\\b').split(var)[length-1].strip()

          if int(kernels[nk]['idxs'][i]) < 0 and kernels[nk]['maps'][i] == OP_MAP:
            if (kernels[nk]['maps'][i] == OP_MAP and kernels[nk]['accs'][i] == OP_INC and repr_temp_array and repro_if):
              body_text = re.sub(r'\b'+var2+'(\[[^\]]\])\[([\\s\+\*A-Za-z0-9_]*)\]'+'', var2+r'\1[(\2)*'+op2_gen_common.get_stride_string(unique_args[i]-1,maps,mapnames2,name,1)+']', body_text)
              
            else:
              body_text = re.sub(r'\b'+var2+'(\[[^\]]\])\[([\\s\+\*A-Za-z0-9_]*)\]'+'', var2+r'\1[(\2)*'+op2_gen_common.get_stride_string(unique_args[i]-1,maps,mapnames,name)+']', body_text)
          else:
            body_text = re.sub('\*\\b'+var2+'\\b\\s*(?!\[)', var2+'[0]', body_text)
            if (kernels[nk]['maps'][i] == OP_MAP and kernels[nk]['accs'][i] == OP_INC and repr_temp_array and repro_if):
              body_text = re.sub(r'\b'+var2+'\[([\\s\+\*A-Za-z0-9_]*)\]'+'', var2+r'[(\1)*'+ \
                                 op2_gen_common.get_stride_string(unique_args[i]-1,maps,mapnames2,name,1)+']', body_text)
            else:
              body_text = re.sub(r'\b'+var2+'\[([\\s\+\*A-Za-z0-9_]*)\]'+'', var2+r'[(\1)*'+ \
                                 op2_gen_common.get_stride_string(unique_args[i]-1,maps,mapnames,name)+']', body_text)

    for nc in range(0,len(consts)):
      varname = consts[nc]['name']
      body_text = re.sub('\\b'+varname+'\\b', varname+'_cuda',body_text)

    signature_text = '__device__ '+head_text + '_gpu( '+signature_text + ') {'
    file_text += signature_text + body_text + '}\n'

    comm('')
    comm(' CUDA kernel function')

    if FORTRAN:
      code('subroutine op_cuda_'+name+'(')
    elif CPP:
      code('__global__ void op_cuda_'+name+'(')

    depth = 2

    if nopts > 0:
      code('int optflags,')

    for g_m in range(0,ninds):
      if (indaccs[g_m]==OP_READ):
        code('const <INDTYP> *__restrict <INDARG>,')
      else:
        if reproducible and repr_temp_array and indaccs[g_m]==OP_INC and (indtyps[g_m] == 'double' or indtyps[g_m] == 'float'):
          code('<TYP> *__restrict tmp_incs'+str(g_m)+'_d,')
        else:
          code('<INDTYP> *__restrict <INDARG>,')

    if nmaps > 0:
      k = []
      for g_m in range(0,nargs):
        if maps[g_m] == OP_MAP and (not mapnames[g_m] in k):
          k = k + [mapnames[g_m]]
          code('const int *__restrict opDat'+str(invinds[inds[g_m]-1])+'Map, ')



    for g_m in range(0,nargs):
      if maps[g_m] == OP_ID:
        if accs[g_m] == OP_READ:
          code('const <TYP> *__restrict <ARG>,')
        else:
          code('<TYP> *<ARG>,')
      elif maps[g_m] == OP_GBL:
        if accs[g_m] == OP_INC or accs[g_m] == OP_MIN or accs[g_m] == OP_MAX or accs[g_m] == OP_WRITE:
          code('<TYP> *<ARG>,')
        elif accs[g_m] == OP_READ:
          code('const <TYP> *<ARG>,')

    if ind_inc and inc_stage==1:
      code('int   *ind_map,')
      code('short *arg_map,')
      code('int   *ind_arg_sizes,')
      code('int   *ind_arg_offs, ')

    if ninds>0:
      if op_color2 or (repr_coloring and repro_if):
        code('int start,           ')
        code('int end,             ')
        code('int *col_reord,      ')
      elif not atomics:
        code('int    block_offset, ')
        code('int   *blkmap,       ')
        code('int   *offset,       ')
        code('int   *nelems,       ')
        code('int   *ncolors,      ')
        code('int   *colors,       ')
        code('int   nblocks,       ')
      else:
        code('int start,           ')
        code('int end,             ')
      if repr_temp_array and repro_if:
        code('int   set_size,')
        code('int   prime_map_dim) {    ')
      else:
        code('int   set_size) {    ')
    else:
      code('int   set_size ) {')
      code('')


    for g_m in range(0,nargs):
      if maps[g_m]==OP_GBL and accs[g_m]!=OP_READ and accs[g_m] != OP_WRITE:
        if not (reproducible and repr_coloring and accs[g_m]==OP_INC and typs[g_m]=="double"):
          code('<TYP> <ARG>_l[<DIM>];')
          if accs[g_m] == OP_INC:
            FOR('d','0','<DIM>')
            code('<ARG>_l[d]=ZERO_<TYP>;')
            ENDFOR()
          else:
            FOR('d','0','<DIM>')
            code('<ARG>_l[d]=<ARG>[d+blockIdx.x*<DIM>];')
            ENDFOR()
      elif maps[g_m]==OP_MAP and accs[g_m]==OP_INC and not op_color2 and not atomics and not (repr_coloring and repro_if):
        code('<TYP> <ARG>_l[<DIM>];')

    if not op_color2 and not atomics and not (repr_coloring and repro_if):
      for m in range (1,ninds+1):
        g_m = m -1
        v = [int(inds[i]==m) for i in range(len(inds))]
        v_i = [vectorised[i] for i in range(len(inds)) if inds[i] == m]
        if sum(v)>1 and sum(v_i)>0: #check this sum(v_i)
          if indaccs[m-1] == OP_INC:
            ind = int(max([idxs[i] for i in range(len(inds)) if inds[i]==m])) + 1
            code('<INDTYP> *arg'+str(invinds[m-1])+'_vec['+str(ind)+'] = {'); depth += 2;
            for n in range(0,nargs):
              if inds[n] == m:
                g_m = n
                code('<ARG>_l,')
            depth -= 2
            code('};')
#
# lengthy code for general case with indirection
#
    if ninds>0 and not op_color2 and not atomics and not (repr_coloring and repro_if):
      code('')
      if inc_stage==1:
        for g_m in range (0,ninds):
          if indaccs[g_m] == OP_INC:
            code('__shared__  int  *<INDARG>_map, <INDARG>_size;')
            code('__shared__  <INDTYP> *<INDARG>_s;')
        code('')
      if ind_inc:
        code('__shared__ int    nelems2, ncolor;')

      code('__shared__ int    nelem, offset_b;')
      code('')
      code('extern __shared__ char shared[];')
      code('')
      IF('blockIdx.x+blockIdx.y*gridDim.x >= nblocks')
      code('return;')
      ENDIF()
      IF('threadIdx.x==0')
      code('')
      comm('get sizes and shift pointers and direct-mapped data')
      code('')
      code('int blockId = blkmap[blockIdx.x + blockIdx.y*gridDim.x  + block_offset];')
      code('')
      code('nelem    = nelems[blockId];')
      code('offset_b = offset[blockId];')
      code('')

      if ind_inc:
        code('nelems2  = blockDim.x*(1+(nelem-1)/blockDim.x);')
        code('ncolor   = ncolors[blockId];')
        code('')

      if inc_stage==1 and ind_inc:
        for g_m in range (0,ninds_staged):
          if indopts_staged[g_m-1] > 0:
            IF('optflags & 1<<'+str(optidxs[indopts_staged[g_m-1]]))
          code('ind_arg'+str(inds[invinds_staged[g_m]]-1)+'_size = ind_arg_sizes['+str(g_m)+'+blockId*'+ str(ninds_staged)+'];')
          if indopts_staged[g_m-1] > 0:
            ENDIF()

        code('')
        for m in range (1,ninds_staged+1):
          g_m = m - 1
          c = [i for i in range(nargs) if inds_staged[i]==m]
          code('ind_arg'+str(inds[invinds_staged[g_m]]-1)+'_map = &ind_map['+str(cumulative_indirect_index[c[0]])+\
          '*set_size] + ind_arg_offs['+str(m-1)+'+blockId*'+str(ninds_staged)+'];')

        code('')
        comm('set shared memory pointers')
        code('int nbytes = 0;')

        for g_m in range(0,ninds_staged):
          code('ind_arg'+str(inds[invinds_staged[g_m]]-1)+'_s = ('+typs[invinds_staged[g_m]]+' *) &shared[nbytes];')
          if g_m < ninds_staged-1:
            if indopts_staged[g_m-1] > 0:
              IF('optflags & 1<<'+str(optidxs[indopts_staged[g_m-1]]))
            code('nbytes    += ROUND_UP(ind_arg'+str(inds[invinds_staged[g_m]]-1)+'_size*sizeof('+typs[invinds_staged[g_m]]+')*'+dims[invinds_staged[g_m]]+');')
            if indopts_staged[g_m-1] > 0:
              ENDIF()


      ENDIF()
      code('__syncthreads(); // make sure all of above completed')
      code('')

      if inc_stage==1:
        for g_m in range(0,ninds):
          if indaccs[g_m] == OP_INC:
            FOR_INC('n','threadIdx.x','<INDARG>_size*<INDDIM>','blockDim.x')
            code('<INDARG>_s[n] = ZERO_<INDTYP>;')
            ENDFOR()
        if ind_inc:
          code('')
          code('__syncthreads();')
          code('')

      if ind_inc:
        FOR_INC('n','threadIdx.x','nelems2','blockDim.x')
        code('int col2 = -1;')
        k = []
        for g_m in range(0,nargs):
          if maps[g_m] == OP_MAP and (not mapinds[g_m] in k):
            k = k + [mapinds[g_m]]
            code('int map'+str(mapinds[g_m])+'idx;')
        IF('n<nelem')
        comm('initialise local variables')

        for g_m in range(0,nargs):
          if maps[g_m]==OP_MAP and accs[g_m]==OP_INC:
            FOR('d','0','<DIM>')
            code('<ARG>_l[d] = ZERO_<TYP>;')
            ENDFOR()
      else:
        FOR_INC('n','threadIdx.x','nelem','blockDim.x')
        k = []
        for g_m in range(0,nargs):
          if maps[g_m] == OP_MAP and (not mapinds[g_m] in k):
            k = k + [mapinds[g_m]]
            code('int map'+str(mapinds[g_m])+'idx;')

      #non-optional maps
      k = []
      for g_m in range(0,nargs):
        if maps[g_m] == OP_MAP and (not optflags[g_m]) and (not mapinds[g_m] in k):
          k = k + [mapinds[g_m]]
          code('map'+str(mapinds[g_m])+'idx = opDat'+str(invmapinds[inds[g_m]-1])+'Map[n + offset_b + set_size * '+str(int(idxs[g_m]))+'];')

      #whatever didn't come up and is opt
      for g_m in range(0,nargs):
        if maps[g_m] == OP_MAP and (not mapinds[g_m] in k):
          if optflags[g_m]==1:
            IF('optflags & 1<<'+str(optidxs[g_m]))
          else:
            k = k + [mapinds[g_m]]

          code('map'+str(mapinds[g_m])+'idx = opDat'+str(invmapinds[inds[g_m]-1])+'Map[n + offset_b + set_size * '+str(int(idxs[g_m]))+'];')
          if optflags[g_m]==1:
            ENDIF()

      code('')
      for g_m in range (0,nargs):
        if accs[g_m] != OP_INC: #TODO: add opt handling here
          u = [i for i in range(0,len(unique_args)) if unique_args[i]-1 == g_m]
          if len(u) > 0 and vectorised[g_m] > 0:
            if accs[g_m] == OP_READ:
              line = 'const <TYP>* <ARG>_vec[] = {\n'
            else:
              line = '<TYP>* <ARG>_vec[] = {\n'

            v = [int(vectorised[i] == vectorised[g_m]) for i in range(0,len(vectorised))]
            first = [i for i in range(0,len(v)) if v[i] == 1]
            first = first[0]

            indent = ' '*(depth+2)
            for k in range(0,sum(v)):
              if soaflags[g_m]:
                line = line + indent + ' &ind_arg'+str(inds[first]-1)+'[map'+str(mapinds[g_m+k])+'idx],\n'
              else:
                line = line + indent + ' &ind_arg'+str(inds[first]-1)+'[<DIM> * map'+str(mapinds[g_m+k])+'idx],\n'
            line = line[:-2]+'};'
            code(line)
#
# simple version for atomics/global coloring
#
    elif ninds>0:
      code('int tid = threadIdx.x + blockIdx.x * blockDim.x;')
      IF('tid + start < end')
      if atomics:
        code('int n = tid + start;')
      else:
        code('int n = col_reord[tid + start];')
      comm('initialise local variables')

      for g_m in range(0,nargs):
        if not (reproducible and accs[g_m]==OP_INC and (typs[g_m]=="double" or typs[g_m]=="float") and repro_if):
          if maps[g_m]==OP_MAP and accs[g_m]==OP_INC:
            code('<TYP> <ARG>_l[<DIM>];')
            FOR('d','0','<DIM>')
            code('<ARG>_l[d] = ZERO_<TYP>;')
            ENDFOR()

      #mapidx declarations
      k = []
      for g_m in range(0,nargs):
        if maps[g_m] == OP_MAP and (not mapinds[g_m] in k):
          k = k + [mapinds[g_m]]
          code('int map'+str(mapinds[g_m])+'idx;')

      #non-optional maps
      k = []
      for g_m in range(0,nargs):
        if maps[g_m] == OP_MAP and (not optflags[g_m]) and (not mapinds[g_m] in k):
          k = k + [mapinds[g_m]] #non-opt
          code('map'+str(mapinds[g_m])+'idx = opDat'+str(invmapinds[inds[g_m]-1])+'Map[n + set_size * '+str(int(idxs[g_m]))+'];')

      #whatever didn't come up and is opt
      for g_m in range(0,nargs):
        if maps[g_m] == OP_MAP and (not mapinds[g_m] in k):
          if optflags[g_m]==1:
            IF('optflags & 1<<'+str(optidxs[g_m]))
          else:
            k = k + [mapinds[g_m]]

          code('map'+str(mapinds[g_m])+'idx = opDat'+str(invmapinds[inds[g_m]-1])+'Map[n + set_size * '+str(int(idxs[g_m]))+'];')
          if optflags[g_m]==1:
            ENDIF()


      for g_m in range(0,ninds):
        if reproducible and repr_temp_array and indaccs[g_m]==OP_INC and (indtyps[g_m] == 'double' or indtyps[g_m] == 'float') and repro_if:
          if optflags[invinds[g_m]]==1:
            IF('optflags & 1<<'+str(optidxs[ invinds[g_m] ]))
          if not (any_soa and soaflags[int(invinds[g_m])]):
            FOR('i','0',' prime_map_dim * '+ str(inddims[g_m]))
            code('tmp_incs'+str(g_m)+'_d[i+n*prime_map_dim * '+str(inddims[g_m])+']=(<TYP>)0.0;\n')
            ENDFOR()          
          if optflags[invinds[g_m]]==1:
            ENDIF()

      for g_m in range (0,nargs):
          u = [i for i in range(0,len(unique_args)) if unique_args[i]-1 == g_m]
          if len(u) > 0 and vectorised[g_m] > 0:
            if accs[g_m] == OP_READ:
              line = 'const <TYP>* <ARG>_vec[] = {\n'
            else:
              line = '<TYP>* <ARG>_vec[] = {\n'
            if atomics and accs[g_m] == OP_INC:
              indent = ' '*(depth+2)
              if repr_temp_array and repro_if and accs[g_m]==OP_INC and (typs[g_m] == 'double' or typs[g_m] == 'float'):
                k=0
                for n in range(0,nargs):
                  if vectorised[n] == vectorised[g_m]:
                    line += indent+'&tmp_incs'+str(inds[g_m]-1)+'_d[(n*prime_map_dim+'+str(k)+')*'+str(dims[g_m])+'],\n'
                    k+=1
              else:
                for n in range(0,nargs):
                  if vectorised[n] == vectorised[g_m]:
                    line = line + indent + 'arg'+str(n)+'_l,\n'
              line = line[:-2]+'};'
              code(line)
            else:
              v = [int(vectorised[i] == vectorised[g_m]) for i in range(0,len(vectorised))]
              first = [i for i in range(0,len(v)) if v[i] == 1]
              first = first[0]

              indent = ' '*(depth+2)
              for k in range(0,sum(v)):
                if soaflags[g_m]:
                  line = line + indent + ' &ind_arg'+str(inds[first]-1)+'[map'+str(mapinds[g_m+k])+'idx],\n'
                else:
                  line = line + indent + ' &ind_arg'+str(inds[first]-1)+'[<DIM> * map'+str(mapinds[g_m+k])+'idx],\n'
              line = line[:-2]+'};'
              code(line)




#
# simple alternative when no indirection
#
    else:
      code('')
      comm('process set elements')
      FOR_INC('n','threadIdx.x+blockIdx.x*blockDim.x','set_size','blockDim.x*gridDim.x')
      for g_m in range (0,nargs):
        if (reproducible and accs[g_m]==OP_INC and typs[g_m]=="double") and reduct and maps[g_m]==OP_GBL and accs[g_m]!=OP_READ and accs[g_m]!=OP_WRITE:
          FOR('d','0','<DIM>')
          code('<ARG>[n+d]=ZERO_<TYP>;')
          ENDFOR()

#
# kernel call
#
    code('')
    comm('user-supplied kernel call')
    line = name+'_gpu('
    prefix = ' '*len(name)
    a = 0 #only apply indentation if its not the 0th argument
    indent =''
    for m in range (0, nargs):
      if a > 0:
        indent = '     '+' '*len(name)

      if maps[m] == OP_GBL:
        if accs[m] == OP_READ or accs[m] == OP_WRITE:
          line += rep(indent+'<ARG>,\n',m)
        else:
          if reproducible and accs[m]==OP_INC and typs[m]=="double":
            line += rep(indent+'<ARG>+n*<DIM>,\n',m)
          else:
            line += rep(indent+'<ARG>_l,\n',m);
        a =a+1
      elif maps[m]==OP_MAP and  accs[m]==OP_INC and not op_color2:
        if vectorised[m]:
          if m+1 in unique_args:
            line += rep(indent+'<ARG>_vec,\n',m)
        else:
          if reproducible and repro_if and accs[m]==OP_INC and (typs[m] == 'double' or typs[m] == 'float'):
            if repr_temp_array:
              if soaflags[m]:
                #line += indent+'tmp_incs'+str(inds[m]-1)+'_d+map'+str(mapinds[m])+'idx,'+'\n'
                line += indent+'tmp_incs'+str(inds[m]-1)+'_d+offset'+str(int(idxs[m]))+',\n'
                if optflags[invinds[m]]==1:
                  IF('optflags & 1<<'+str(optidxs[ invinds[g_m] ]))
                code('int offset'+str(int(idxs[m]))+'=n+set_size*'+str(int(idxs[m]))+';')
                FOR('d','0',str(dims[m]))
                code('tmp_incs'+str(inds[m]-1)+'_d[offset'+str(int(idxs[m]))+'+d*opMap_'+str(mapnames2[m])+'_'+name+'_stride_temp_inc_OP2CONSTANT' +']=('+str(typs[m])+')0.0;\n')
                ENDFOR()
                if optflags[invinds[m]]==1:
                  ENDIF()
              else:
                k=0
                for i in range(0,m):
                  if inds[m]==inds[i]:
                    k+=1
                line += indent+'&tmp_incs'+str(inds[m]-1)+'_d[(n*prime_map_dim+'+str(k)+')*'+str(dims[m])+'],\n'
            elif repr_coloring:
              line += rep(indent+'ind_arg'+str(inds[m]-1)+'+map'+str(mapinds[m])+'idx*<DIM>,'+'\n',m)
          else:
            line += rep(indent+'<ARG>_l,\n',m)
        a =a+1
      elif maps[m]==OP_MAP:
        if vectorised[m]:
          if m+1 in unique_args:
            line += rep(indent+'<ARG>_vec,\n',m)
        else:
          if soaflags[m]:
            line += rep(indent+'ind_arg'+str(inds[m]-1)+'+map'+str(mapinds[m])+'idx,'+'\n',m)
          else:
            line += rep(indent+'ind_arg'+str(inds[m]-1)+'+map'+str(mapinds[m])+'idx*<DIM>,'+'\n',m)
        a =a+1
      elif maps[m]==OP_ID:
        if ninds>0 and not op_color2 and not atomics and not (repr_coloring and repro_if):
          if soaflags[m]:
            line += rep(indent+'<ARG>+(n+offset_b),\n',m)
          else:
            line += rep(indent+'<ARG>+(n+offset_b)*<DIM>,\n',m)
          a =a+1
        else:
          if soaflags[m]:
            line += rep(indent+'<ARG>+n,\n',m)
          else:
            line += rep(indent+'<ARG>+n*<DIM>,\n',m)
          a =a+1
      else:
        print('internal error 1 ')

    code(line[0:-2]+');') #remove final ',' and \n

#
# updating for indirect kernels ...
#
    if ninds>0 and not op_color2 and not atomics and not repr_coloring:
      if ind_inc:
        code('col2 = colors[n+offset_b];')
        ENDIF()
        code('')
        comm('store local variables')
        code('')
        if inc_stage==1:
          for g_m in range(0,nargs):
            if maps[g_m]==OP_MAP and accs[g_m]==OP_INC:
              code('int <ARG>_map;')
          IF('col2>=0')
          for g_m in range(0,nargs):
            if maps[g_m] == OP_MAP and accs[g_m] == OP_INC:
              code('<ARG>_map = arg_map['+str(cumulative_indirect_index[g_m])+'*set_size+n+offset_b];')
          ENDIF()
          code('')

        FOR('col','0','ncolor')
        IF('col2==col')

        if inc_stage==1:
          for g_m in range(0,nargs):
            if maps[g_m] == OP_MAP and accs[g_m] == OP_INC:
              if optflags[g_m]==1:
                IF('optflags & 1<<'+str(optidxs[g_m]))
              for d in range(0,int(dims[g_m])):
                if soaflags[g_m]:
                  code('<ARG>_l['+str(d)+'] += ind_arg'+str(inds[g_m]-1)+'_s[<ARG>_map+'+str(d)+'*ind_arg'+str(inds[g_m]-1)+'_size];')
                else:
                  code('<ARG>_l['+str(d)+'] += ind_arg'+str(inds[g_m]-1)+'_s['+str(d)+'+<ARG>_map*<DIM>];')
#          for g_m in range(0,nargs):
#            if maps[g_m] == OP_MAP and accs[g_m] == OP_INC:
              for d in range(0,int(dims[g_m])):
                if soaflags[g_m]:
                  code('ind_arg'+str(inds[g_m]-1)+'_s[<ARG>_map+'+str(d)+'*ind_arg'+str(inds[g_m]-1)+'_size] = <ARG>_l['+str(d)+'];')
                else:
                  code('ind_arg'+str(inds[g_m]-1)+'_s['+str(d)+'+<ARG>_map*<DIM>] = <ARG>_l['+str(d)+'];')

              if optflags[g_m]==1:
                ENDIF()
        else:
          for g_m in range(0,nargs):
            if maps[g_m] == OP_MAP and accs[g_m] == OP_INC:
              if optflags[g_m]==1:
                IF('optflags & 1<<'+str(optidxs[g_m]))
              for d in range(0,int(dims[g_m])):
                if soaflags[g_m]:
                  code('<ARG>_l['+str(d)+'] += ind_arg'+str(inds[g_m]-1)+'['+str(d)+'*'+op2_gen_common.get_stride_string(g_m,maps,mapnames,name)+'+map'+str(mapinds[g_m])+'idx];')
                else:
                  code('<ARG>_l['+str(d)+'] += ind_arg'+str(inds[g_m]-1)+'['+str(d)+'+map'+str(mapinds[g_m])+'idx*<DIM>];')
#          for g_m in range(0,nargs):
#            if maps[g_m] == OP_MAP and accs[g_m] == OP_INC:
              for d in range(0,int(dims[g_m])):
                if soaflags[g_m]:
                  code('ind_arg'+str(inds[g_m]-1)+'['+str(d)+'*'+op2_gen_common.get_stride_string(g_m,maps,mapnames,name)+'+map'+str(mapinds[g_m])+'idx] = <ARG>_l['+str(d)+'];')
                else:
                  code('ind_arg'+str(inds[g_m]-1)+'['+str(d)+'+map'+str(mapinds[g_m])+'idx*<DIM>] = <ARG>_l['+str(d)+'];')
              if optflags[g_m]==1:
                ENDIF()

        ENDFOR()
        code('__syncthreads();')
        ENDFOR()
    if ninds>0 and atomics:
          for g_m in range(0,nargs):
            if maps[g_m] == OP_MAP and accs[g_m] == OP_INC and not (repr_temp_array and repro_if and accs[g_m]==OP_INC and (typs[g_m] == 'double' or typs[g_m] == 'float')):
              if optflags[g_m]==1:
                IF('optflags & 1<<'+str(optidxs[g_m]))
              for d in range(0,int(dims[g_m])):
                if soaflags[g_m]:
                  code('atomicAdd(&ind_arg'+str(inds[g_m]-1)+'['+str(d)+'*'+op2_gen_common.get_stride_string(g_m,maps,mapnames,name)+'+map'+str(mapinds[g_m])+'idx],<ARG>_l['+str(d)+']);')
                else:
                  code('atomicAdd(&ind_arg'+str(inds[g_m]-1)+'['+str(d)+'+map'+str(mapinds[g_m])+'idx*<DIM>],<ARG>_l['+str(d)+']);')
              if optflags[g_m]==1:
                ENDIF()


    ENDFOR()

    if inc_stage:
      for g_m in range(0,ninds):
        if indaccs[g_m]==OP_INC:
          if indopts[g_m] > 0:
            IF('optflags & 1<<'+str(optidxs[indopts[g_m-1]]))
          if soaflags[invinds[g_m]]:
            FOR_INC('n','threadIdx.x','<INDARG>_size','blockDim.x')
            for d in range(0,int(dims[invinds[g_m]])):
              code('arg'+str(invinds[g_m])+'_l['+str(d)+'] = <INDARG>_s[n+'+str(d)+'*<INDARG>_size] + <INDARG>[<INDARG>_map[n]+'+str(d)+'*'+op2_gen_common.get_stride_string(invinds[g_m],maps,mapnames,name)+'];')
            for d in range(0,int(dims[invinds[g_m]])):
              code('<INDARG>[<INDARG>_map[n]+'+str(d)+'*'+op2_gen_common.get_stride_string(invinds[g_m],maps,mapnames,name)+'] = arg'+str(invinds[g_m])+'_l['+str(d)+'];')
            ENDFOR()
          else:
            FOR_INC('n','threadIdx.x','<INDARG>_size*<INDDIM>','blockDim.x')
            code('<INDARG>[n%<INDDIM>+<INDARG>_map[n/<INDDIM>]*<INDDIM>] += <INDARG>_s[n];')
            ENDFOR()
          if indopts[g_m] > 0:
            ENDIF()

#
# global reduction
#
    if reduct:
       code('')
       comm('global reductions')
       code('')
       for m in range (0,nargs):
         g_m = m
         if maps[m]==OP_GBL and accs[m]!=OP_READ and accs[m] != OP_WRITE and not (reproducible and accs[m]==OP_INC and typs[m]=="double"):
           FOR('d','0','<DIM>')
           if accs[m]==OP_INC:
             code('op_reduction<OP_INC>(&<ARG>[d+blockIdx.x*<DIM>],<ARG>_l[d]);')
           elif accs[m]==OP_MIN:
             code('op_reduction<OP_MIN>(&<ARG>[d+blockIdx.x*<DIM>],<ARG>_l[d]);')
           elif accs[m]==OP_MAX:
             code('op_reduction<OP_MAX>(&<ARG>[d+blockIdx.x*<DIM>],<ARG>_l[d]);')
           else:
             print('internal error: invalid reduction option')
             sys.exit(2);
           ENDFOR()
    depth -= 2
    code('}')
    code('')

##########################################################################
# then C++ stub function
##########################################################################

    code('')
    comm('host stub function')
    code('void op_par_loop_'+name+'(char const *name, op_set set,')
    depth += 2

    for m in unique_args:
      g_m = m - 1
      if m == unique_args[len(unique_args)-1]:
        code('op_arg <ARG>){')
        code('')
      else:
        code('op_arg <ARG>,')

    for g_m in range (0,nargs):
      if maps[g_m]==OP_GBL:
        code('<TYP>*<ARG>h = (<TYP> *)<ARG>.data;')

    code('int nargs = '+str(nargs)+';')
    code('op_arg args['+str(nargs)+'];')
    code('')


    for g_m in range (0,nargs):
      u = [i for i in range(0,len(unique_args)) if unique_args[i]-1 == g_m]
      if len(u) > 0 and vectorised[g_m] > 0:
        code('<ARG>.idx = 0;')
        code('args['+str(g_m)+'] = <ARG>;')

        v = [int(vectorised[i] == vectorised[g_m]) for i in range(0,len(vectorised))]
        first = [i for i in range(0,len(v)) if v[i] == 1]
        first = first[0]
        if (optflags[g_m] == 1):
          argtyp = 'op_opt_arg_dat(arg'+str(first)+'.opt, '
        else:
          argtyp = 'op_arg_dat('

        FOR('v','1',str(sum(v)))
        code('args['+str(g_m)+' + v] = '+argtyp+'arg'+str(first)+'.dat, v, arg'+\
        str(first)+'.map, <DIM>, "<TYP>", '+accsstring[accs[g_m]-1]+');')
        ENDFOR()
        code('')

      elif vectorised[g_m]>0:
        pass
      else:
        code('args['+str(g_m)+'] = <ARG>;')

    if nopts>0:
      code('int optflags = 0;')
      for i in range(0,nargs):
        if optflags[i] == 1:
          IF('args['+str(i)+'].opt')
          code('optflags |= 1<<'+str(optidxs[i])+';')
          ENDIF()
    if nopts > 30:
      print('ERROR: too many optional arguments to store flags in an integer')
#
# start timing
#
    code('')
    comm(' initialise timers')
    code('double cpu_t1, cpu_t2, wall_t1, wall_t2;')
    code('op_timing_realloc('+str(nk)+');')
    code('op_timers_core(&cpu_t1, &wall_t1);')
    code('OP_kernels[' +str(nk)+ '].name      = name;')
    code('OP_kernels[' +str(nk)+ '].count    += 1;')
    code('')

#
#   indirect bits
#
    if ninds>0:
      code('')
      code('int    ninds   = '+str(ninds)+';')
      line = 'int    inds['+str(nargs)+'] = {'
      for m in range(0,nargs):
        line += str(inds[m]-1)+','
      code(line[:-1]+'};')
      code('')

      IF('OP_diags>2')
      code('printf(" kernel routine with indirection: '+name+'\\n");')
      ENDIF()

      if not atomics:
        code('')
        comm('get plan')
        code('#ifdef OP_PART_SIZE_'+ str(nk))
        code('  int part_size = OP_PART_SIZE_'+str(nk)+';')
        code('#else')
        code('  int part_size = OP_part_size;')
        code('#endif')
        code('')
      #code('int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);')
      code('int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);')

#
# direct bit
#
    else:
      code('')
      IF('OP_diags>2')
      code('printf(" kernel routine w/o indirection:  '+ name + '");')
      ENDIF()
      code('')
      code('int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);')
      #code('op_mpi_halo_exchanges_cuda(set, nargs, args);')

#
# for reproducible incs method
#
    
    repro_prime_map=''
    if reproducible:
      if repr_temp_array and repro_if:
        repro_if=0
        if ninds>0:
          if nmaps > 0:
            k = []
            line='set->size > 0 && '
            for g_m in range(0,nargs):
              if accs[g_m] == OP_INC and maps[g_m] == OP_MAP and (not mapnames2[g_m] in k):
                k = k + [mapnames2[g_m]]
                code('op_map prime_map_'+str(mapnames2[g_m])+' = <ARG>.map;\n')
                code('op_reversed_map rev_map_'+str(mapnames2[g_m])+' = OP_reversed_map_list[prime_map_'+str(mapnames2[g_m])+'->index];\n')
                line = line + 'rev_map_'+str(mapnames2[g_m])+' != NULL && '
                code('')
                repro_if=1
                repro_prime_map=mapnames2[g_m]
            
            if repro_if:
              IF(line[:-3])
            
            for g_map in k:
              code('int prime_map_'+str(g_map)+'_dim = prime_map_'+str(g_map)+'->dim;\n')
              code('int set_from_size_'+str(g_map)+' = prime_map_'+str(g_map)+'->from->size + prime_map_'+str(g_map)+'->from->exec_size;\n')
              code('int set_to_size_'+str(g_map)+' = prime_map_'+str(g_map)+'->to->size + prime_map_'+str(g_map)+'->to->exec_size + prime_map_'+str(g_map)+'->to->nonexec_size;\n')
              code('')
            
            k=[]
            for g_m in range(0,nargs):              
              if accs[g_m] == OP_INC and maps[g_m] == OP_MAP:
                first=0
                for i in range(0,g_m+1):
                  if maps[g_m]==maps[i] and varrrrrr[g_m]==varrrrrr[i]:
                    first=i
                    break
                if not first in k:
                  k = k + [first] 
                  code('<TYP> *tmp_incs'+str(first)+'_d = NULL;\n')
                  if optflags[g_m]==1:
                    IF('<ARG>.opt')
                  code('int required_tmp_incs_size'+str(first)+' = set_from_size_'+str(mapnames2[first])+' * prime_map_'+str(mapnames2[first])+'_dim * arg'+str(first)+'.dat->size;\n')
                  

                  code('reallocTempArrays(arg'+str(first)+'.dat->index, required_tmp_incs_size'+str(first)+');')
                  code('tmp_incs'+str(first)+'_d = (<TYP> *)op_repr_incs[arg'+str(first)+'.dat->index].tmp_incs_d;\n')
                  
                  if optflags[g_m]==1:
                    ENDIF()
                  code('')
      elif repr_coloring and repro_if:
        if ninds>0:
          if nmaps > 0:
            line='set->size > 0 && '
            for g_m in range(0,nargs):
              if (accs[g_m] == OP_INC or accs[g_m] == OP_RW) and maps[g_m] == OP_MAP:
                code('op_map prime_map = <ARG>.map;\n')
                code('op_reversed_map rev_map = OP_reversed_map_list[prime_map->index];\n')
                line = line + 'rev_map != NULL && '
                code('')
                repro_if=1
                break
            
            if repro_if:
              IF(line[:-3])
        if not repro_if:
          IF('set->size > 0')
      else:
        IF('set->size > 0')
    else:
      IF('set->size > 0')
    code('')

#
# kernel call for indirect version
#
    if ninds>0 and not atomics:
      if inc_stage==1 and ind_inc:
        code('op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_STAGE_INC);')
      elif op_color2:
        code('op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_COLOR2);')
      else:
        code('op_plan *Plan = op_plan_get(name,set,part_size,nargs,args,ninds,inds);')
      code('')


#
# transfer constants
#
    g = [i for i in range(0,nargs) if maps[i] == OP_GBL and (accs[i] == OP_READ or accs[i] == OP_WRITE)]
    if len(g)>0:
      comm('transfer constants to GPU')
      code('int consts_bytes = 0;')
      for m in range(0,nargs):
        g_m = m
        if maps[m]==OP_GBL and (accs[m]==OP_READ or accs[m] == OP_WRITE):
          code('consts_bytes += ROUND_UP(<DIM>*sizeof(<TYP>));')

      code('reallocConstArrays(consts_bytes);')
      code('consts_bytes = 0;')

      for m in range(0,nargs):
        if maps[m]==OP_GBL and (accs[m] == OP_READ or accs[m] == OP_WRITE):
          g_m = m
          code('<ARG>.data   = OP_consts_h + consts_bytes;')
          code('<ARG>.data_d = OP_consts_d + consts_bytes;')
          FOR('d','0','<DIM>')
          code('((<TYP> *)<ARG>.data)[d] = <ARG>h[d];')
          ENDFOR()
          code('consts_bytes += ROUND_UP(<DIM>*sizeof(<TYP>));')
      code('mvConstArraysToDevice(consts_bytes);')
      code('')

      #managing constants
    if any_soa:
      if nmaps > 0:
        k = []
        for g_m in range(0,nargs):
          if repr_temp_array and repro_if and maps[g_m] == OP_MAP and accs[g_m] == OP_INC and (not (mapnames[g_m]+'tmpa') in k):
            k = k + [(mapnames[g_m]+'tmpa')]
            IF('(OP_kernels[' +str(nk)+ '].count==1) || (opMap_'+str(mapnames2[g_m])+'_'+name+'_stride_temp_inc_OP2HOST != set_from_size_'+str(mapnames2[g_m])+'* prime_map_'+str(mapnames2[g_m])+'_dim)')
            code('opMap_'+str(mapnames2[g_m])+'_'+name+'_stride_temp_inc_OP2HOST = set_from_size_'+str(mapnames2[g_m])+'* prime_map_'+str(mapnames2[g_m])+'_dim;')
            code('cudaMemcpyToSymbol(opMap_'+str(mapnames2[g_m])+'_'+name+'_stride_temp_inc_OP2CONSTANT, &opMap_'+str(mapnames2[g_m])+'_'+name+'_stride_temp_inc_OP2HOST,sizeof(int));')
            ENDIF()
          if maps[g_m] == OP_MAP and (not mapnames[g_m] in k):
            k = k + [mapnames[g_m]]
            IF('(OP_kernels[' +str(nk)+ '].count==1) || (opDat'+str(invinds[inds[g_m]-1])+'_'+name+'_stride_OP2HOST != getSetSizeFromOpArg(&arg'+str(g_m)+'))')
            code('opDat'+str(invinds[inds[g_m]-1])+'_'+name+'_stride_OP2HOST = getSetSizeFromOpArg(&arg'+str(g_m)+');')
            code('cudaMemcpyToSymbol(opDat'+str(invinds[inds[g_m]-1])+'_'+name+'_stride_OP2CONSTANT, &opDat'+str(invinds[inds[g_m]-1])+'_'+name+'_stride_OP2HOST,sizeof(int));')
            ENDIF()
      if dir_soa!=-1:
          IF('(OP_kernels[' +str(nk)+ '].count==1) || (direct_'+name+'_stride_OP2HOST != getSetSizeFromOpArg(&arg'+str(dir_soa)+'))')
          code('direct_'+name+'_stride_OP2HOST = getSetSizeFromOpArg(&arg'+str(dir_soa)+');')
          code('cudaMemcpyToSymbol(direct_'+name+'_stride_OP2CONSTANT,&direct_'+name+'_stride_OP2HOST,sizeof(int));')
          ENDIF()

#
# transfer global reduction initial data
#

    if ninds == 0 or atomics:
      comm('set CUDA execution parameters')
      code('#ifdef OP_BLOCK_SIZE_'+str(nk))
      code('  int nthread = OP_BLOCK_SIZE_'+str(nk)+';')
      code('#else')
      code('  int nthread = OP_block_size;')
      code('#endif')
      code('')
      if ninds==0:
        code('int nblocks = 200;')
        code('')

    if reduct:
      comm('transfer global reduction data to GPU')
      if not (reproducible and accs[g_m]==OP_INC and typs[g_m]=="double"):
        if ninds>0 and not atomics:
          code('int maxblocks = 0;')
          FOR('col','0','Plan->ncolors')
          code('maxblocks = MAX(maxblocks,Plan->ncolblk[col]);')
          ENDFOR()
        elif atomics and ninds>0:
          code('int maxblocks = (MAX(set->core_size, set->size+set->exec_size-set->core_size)-1)/nthread+1;')
        else:
          code('int maxblocks = nblocks;')

      code('int reduct_bytes = 0;')
      code('int reduct_size  = 0;')

      for g_m in range(0,nargs):
        if maps[g_m]==OP_GBL and accs[g_m]!=OP_READ and accs[g_m]!=OP_WRITE:
          if reproducible and accs[g_m]==OP_INC and typs[g_m]=="double":
            code('reduct_bytes += ROUND_UP(set_size*<ARG>.size);')
          else:
            code('reduct_bytes += ROUND_UP(maxblocks*<DIM>*sizeof(<TYP>));')
          code('reduct_size   = MAX(reduct_size,sizeof(<TYP>));')

      code('reallocReductArrays(reduct_bytes);')
      code('reduct_bytes = 0;')

      for g_m in range(0,nargs):
        if maps[g_m]==OP_GBL and accs[g_m]!=OP_READ and accs[g_m]!=OP_WRITE:
          code('<ARG>.data   = OP_reduct_h + reduct_bytes;')
          code('<ARG>.data_d = OP_reduct_d + reduct_bytes;')
          if reproducible and accs[g_m]==OP_INC and typs[g_m]=="double":
            code('reduct_bytes += ROUND_UP(set_size*<ARG>.size);')
          else:
            FOR('b','0','maxblocks')
            FOR('d','0','<DIM>')
            if accs[g_m]==OP_INC:
              code('((<TYP> *)<ARG>.data)[d+b*<DIM>] = ZERO_<TYP>;')
            else:
              code('((<TYP> *)<ARG>.data)[d+b*<DIM>] = <ARG>h[d];')
            ENDFOR()
            ENDFOR()
            code('reduct_bytes += ROUND_UP(maxblocks*<DIM>*sizeof(<TYP>));')
      
      if not reproducible:
        code('mvReductArraysToDevice(reduct_bytes);')
      code('')

    if repro_if:
    #if repr_coloring:
      #code('op_mpi_wait_all_cuda(nargs, args);')
      code('op_mpi_wait_all_grouped(nargs, args, 2);')

#
# kernel call for indirect version
#
    if ninds>0 and not atomics:
      comm('execute plan')
      if not op_color2:
        code('')
        code('int block_offset = 0;')
      if repro_if and repr_coloring:
        code('op_mpi_wait_all_grouped(nargs, args, 2);')
        #code('op_mpi_wait_all_cuda(nargs, args);')
        FOR('col','0','rev_map->number_of_colors')
      else:
        FOR('col','0','Plan->ncolors')
      if not (reproducible and repr_coloring and repro_if):
        IF('col==Plan->ncolors_core')
        code('op_mpi_wait_all_grouped(nargs, args, 2);')
        #code('op_mpi_wait_all_cuda(nargs, args);')
        ENDIF()
      code('#ifdef OP_BLOCK_SIZE_'+str(nk))
      code('int nthread = OP_BLOCK_SIZE_'+str(nk)+';')
      code('#else')
      code('int nthread = OP_block_size;')
      code('#endif')
      code('')
      if repro_if and repr_coloring:
        code('int start = rev_map->color_based_exec_row_starts[col];')
        code('int end = rev_map->color_based_exec_row_starts[col+1];')
        code('int nblocks = (end - start - 1)/nthread + 1;')
      else:
        if op_color2:
          code('int start = Plan->col_offsets[0][col];')
          code('int end = Plan->col_offsets[0][col+1];')
          code('int nblocks = (end - start - 1)/nthread + 1;')
        else:
          code('dim3 nblocks = dim3(Plan->ncolblk[col] >= (1<<16) ? 65535 : Plan->ncolblk[col],')
          code('Plan->ncolblk[col] >= (1<<16) ? (Plan->ncolblk[col]-1)/65535+1: 1, 1);')
          IF('Plan->ncolblk[col] > 0')

      if reduct or (inc_stage==1 and ind_inc):
        if reduct and inc_stage==1:
          code('int nshared = MAX(Plan->nshared,reduct_size*nthread);')
        elif reduct:
          code('int nshared = reduct_size*nthread;')
        else:
          code('int nshared = Plan->nsharedCol[col];')
        code('op_cuda_'+name+'<<<nblocks,nthread,nshared>>>(')
      else:
        code('op_cuda_'+name+'<<<nblocks,nthread>>>(')

      if nopts > 0:
        code('optflags,')
      for m in range(1,ninds+1):
        g_m = invinds[m-1]
        code('(<TYP> *)<ARG>.data_d,')
      if nmaps > 0:
        k = []
        for g_m in range(0,nargs):
          if maps[g_m] == OP_MAP and (not mapnames[g_m] in k):
            k = k + [mapnames[g_m]]
            code('arg'+str(g_m)+'.map_data_d, ')
      for g_m in range(0,nargs):
        if inds[g_m]==0:
          code('(<TYP>*)<ARG>.data_d,')

      if inc_stage==1 and ind_inc:
        code('Plan->ind_map,')
        code('Plan->loc_map,')
        code('Plan->ind_sizes,')
        code('Plan->ind_offs,')
      if repro_if and repr_coloring:
        code('start,')
        code('end,')
        code('rev_map->color_based_exec_d,')
      else:
        if op_color2:
          code('start,')
          code('end,')
          code('Plan->col_reord,')
        else:
          code('block_offset,')
          code('Plan->blkmap,')
          code('Plan->offset,')
          code('Plan->nelems,')
          code('Plan->nthrcol,')
          code('Plan->thrcol,')
          code('Plan->ncolblk[col],')
      code('set->size+set->exec_size);')
      code('')
      if reduct:
        comm('transfer global reduction data back to CPU')
        IF('col == Plan->ncolors_owned-1')
        code('mvReductArraysToHost(reduct_bytes);')
        ENDIF()
      if not op_color2 and not (repr_coloring and repro_if):
        ENDIF()
        code('block_offset += Plan->ncolblk[col];')
      ENDFOR()

#
#
#
    elif ninds>0 and atomics:
      if reduct:
        FOR('round','0','3')
      else:
        FOR('round','0','2')
      IF('round==1')
      code('op_mpi_wait_all_grouped(nargs, args, 2);')
      #code('op_mpi_wait_all_cuda(nargs, args);')
      ENDIF()
      if reduct:
        code('int start = round==0 ? 0 : (round==1 ? set->core_size : set->size);')
        code('int end = round==0 ? set->core_size : (round==1? set->size :  set->size + set->exec_size);')
      else:
        code('int start = round==0 ? 0 : set->core_size;')
        code('int end = round==0 ? set->core_size : set->size + set->exec_size;')
      IF('end-start>0')
      code('int nblocks = (end-start-1)/nthread+1;')
      if reduct:
        code('int nshared = reduct_size*nthread;')
        code('op_cuda_'+name+'<<<nblocks,nthread,nshared>>>(')
      else:
        code('op_cuda_'+name+'<<<nblocks,nthread>>>(')
      if nopts > 0:
        code('optflags,')
      for m in range(1,ninds+1):
        g_m = invinds[m-1]
        if reproducible and repr_temp_array and accs[g_m]==OP_INC and (typs[g_m] == 'double' or typs[g_m] == 'float'):
          code('tmp_incs'+str(g_m)+'_d,')
        else:
          code('(<TYP> *)<ARG>.data_d,')
      if nmaps > 0:
        k = []
        for g_m in range(0,nargs):
          if maps[g_m] == OP_MAP and (not mapnames[g_m] in k):
            k = k + [mapnames[g_m]]
            code('arg'+str(g_m)+'.map_data_d, ')
      for g_m in range(0,nargs):
        if inds[g_m]==0:
          code('(<TYP>*)<ARG>.data_d,')
      if repro_if and repr_temp_array:
        code('start,end,set->size+set->exec_size,prime_map_'+repro_prime_map+'_dim);')
      else:
        code('start,end,set->size+set->exec_size);')
      ENDIF()
      if reduct:
        code('if (round==1) mvReductArraysToHost(reduct_bytes);')

      ENDFOR()
#
# kernel call for direct version
#
    else:
      if reduct:
        code('int nshared = reduct_size*nthread;')
        code('op_cuda_'+name+'<<<nblocks,nthread,nshared>>>(')
      else:
        code('op_cuda_'+name+'<<<nblocks,nthread>>>(')

      indent = '  '#*(len(name)+42)
      if nopts > 0:
        code(indent+'optflags,')
      for g_m in range(0,nargs):
        if g_m > 0:
          code(indent+'(<TYP> *) <ARG>.data_d,')
        else:
          code(indent+'(<TYP> *) <ARG>.data_d,')

      code(indent+'set->size );')



    #apply increments to actual data
    if reproducible and repr_temp_array:    
      if ninds>0:
        if nmaps > 0:
          code('int nblocks;')
          k=[]
          for g_m in range(0,nargs):              
            if accs[g_m] == OP_INC and maps[g_m] == OP_MAP:
              first=0
              for i in range(0,g_m+1):
                if maps[g_m]==maps[i] and varrrrrr[g_m]==varrrrrr[i]:
                  first=i
                  break
              
              if not first in k:   
                code('')
                k = k + [first]
                if optflags[g_m]==1:
                  IF('optflags & 1<<'+str(optidxs[g_m]))
                code('nblocks = (set_to_size_'+str(mapnames2[first])+'-1)/nthread+1;')
                if any_soa and soaflags[g_m] :
                  code('apply_tmp_incs_soa<<<nblocks,nthread>>>(tmp_incs'+str(g_m)+'_d, <ARG>.data_d, rev_map_'+str(mapnames2[first])+'->reversed_map_d, rev_map_'+str(mapnames2[first])+ \
                  '->row_start_idx_d, <ARG>.dim, prime_map_'+str(mapnames2[first])+'_dim,set_to_size_'+str(mapnames2[first])+',set_from_size_'+str(mapnames2[first])+');')
                else:
                  code('apply_tmp_incs<<<nblocks,nthread>>>(tmp_incs'+str(g_m)+'_d, <ARG>.data_d, rev_map_'+str(mapnames2[first])+'->reversed_map_d, rev_map_'+str(mapnames2[first])+'->row_start_idx_d, <ARG>.dim, set_to_size_'+str(mapnames2[first])+');')
                if optflags[g_m]==1:
                  ENDIF()

    if ninds>0 and not atomics:
      code('OP_kernels['+str(nk)+'].transfer  += Plan->transfer;')
      code('OP_kernels['+str(nk)+'].transfer2 += Plan->transfer2;')


#
# transfer global reduction initial data
#
    if reduct:
      if ninds == 0:
        comm('transfer global reduction data back to CPU')
        code('mvReductArraysToHost(reduct_bytes);')

      for m in range(0,nargs):
        g_m = m
        if maps[m]==OP_GBL and accs[m]!=OP_READ and accs[m] != OP_WRITE:
          if reproducible and accs[m]==OP_INC and typs[m]=="double":
            code('reprLocalSum(&<ARG>,set_size,(double*)<ARG>.data);')
            code('<ARG>.data = (char *)<ARG>h;')
            code('op_mpi_repr_inc_reduce_double(&<ARG>,(double*)<ARG>.data);')
          else:
            FOR('b','0','maxblocks')
            FOR('d','0','<DIM>')
            if accs[m]==OP_INC:
              code('<ARG>h[d] = <ARG>h[d] + ((<TYP> *)<ARG>.data)[d+b*<DIM>];')
            elif accs[m]==OP_MIN:
              code('<ARG>h[d] = MIN(<ARG>h[d],((<TYP> *)<ARG>.data)[d+b*<DIM>]);')
            elif accs[m]==OP_MAX:
              code('<ARG>h[d] = MAX(<ARG>h[d],((<TYP> *)<ARG>.data)[d+b*<DIM>]);')
            ENDFOR()
            ENDFOR()

            code('<ARG>.data = (char *)<ARG>h;')
            code('op_mpi_reduce(&<ARG>,<ARG>h);')
          
    for g_m in range(0,nargs):
      if maps[g_m] == OP_GBL and accs[g_m] == OP_WRITE:
        code('')
        code('mvConstArraysToHost(consts_bytes);')
        break

    for g_m in range(0,nargs):
      if maps[g_m] == OP_GBL and accs[g_m] == OP_WRITE:
        FOR('d','0','<DIM>')
        code('<ARG>h[d] = ((<TYP> *)<ARG>.data)[d];')
        ENDFOR()
        code('<ARG>.data = (char *)<ARG>h;')
        code('op_mpi_reduce(&<ARG>,<ARG>h);')

    ENDIF()
    code('op_mpi_set_dirtybit_cuda(nargs, args);')

#
# update kernel record
#

    code('cutilSafeCall(cudaDeviceSynchronize());')
    comm('update kernel record')
    code('op_timers_core(&cpu_t2, &wall_t2);')
    code('OP_kernels[' +str(nk)+ '].time     += wall_t2 - wall_t1;')

    if ninds == 0:
      line = 'OP_kernels['+str(nk)+'].transfer += (float)set->size *'

      for g_m in range (0,nargs):
        if optflags[g_m]==1:
          IF('<ARG>.opt')
        if maps[g_m]!=OP_GBL:
          if accs[g_m]==OP_READ:
            code(line+' <ARG>.size;')
          else:
            code(line+' <ARG>.size * 2.0f;')
        if optflags[g_m]==1:
          ENDIF()
    depth = depth - 2
    code('}')


##########################################################################
#  output individual kernel file
##########################################################################
    if not os.path.exists('cuda'):
        os.makedirs('cuda')
    fid = open('cuda/'+name+'_kernel.cu','w')
    date = datetime.datetime.now()
    fid.write('//\n// auto-generated by op2.py\n//\n\n')
    fid.write(file_text)
    fid.close()

# end of main kernel call loop


##########################################################################
#  output one master kernel file
##########################################################################

  file_text = ''

  comm('global constants')

  code('#ifndef MAX_CONST_SIZE')
  code('#define MAX_CONST_SIZE 128')
  code('#endif')
  code('')

  for nc in range (0,len(consts)):
    if consts[nc]['dim']==1:
      code('__constant__ '+consts[nc]['type'][1:-1]+' '+consts[nc]['name']+'_cuda;')
    else:
      if consts[nc]['dim'].isdigit() and int(consts[nc]['dim']) > 0:
        num = str(consts[nc]['dim'])
      else:
        num = 'MAX_CONST_SIZE'

      code('__constant__ '+consts[nc]['type'][1:-1]+' '+consts[nc]['name']+'_cuda['+num+'];')
  code('')

  comm('header')

  if os.path.exists('./user_types.h'):
    code('#ifndef OP_FUN_PREFIX\n#define OP_FUN_PREFIX __host__ __device__\n#endif')
    code('#include "../user_types.h"')
  code('#include "op_lib_cpp.h"')
  code('#include "op_cuda_rt_support.h"')
  code('#include "op_cuda_reduction.h"')

  code('')
  code('void op_decl_const_char(int dim, char const *type,')
  code('int size, char *dat, char const *name){')
  depth = depth + 2

  code('if (!OP_hybrid_gpu) return;')
  for nc in range(0,len(consts)):
    IF('!strcmp(name,"'+consts[nc]['name']+'")')
    if not consts[nc]['dim'] or int(consts[nc]['dim']) > 1:
      IF('!strcmp(name,"'+consts[nc]['name']+'") && size>MAX_CONST_SIZE')
      code('printf("error: MAX_CONST_SIZE not big enough\\n"); exit(1);')
      ENDIF()
    code('cutilSafeCall(cudaMemcpyToSymbol('+consts[nc]['name']+'_cuda, dat, dim*size));')
    ENDIF()
    code('else ')

  code('{')
  depth = depth + 2
  code('printf("error: unknown const name\\n"); exit(1);')
  ENDIF()


  depth = depth - 2
  code('}')
  code('')
  comm('user kernel files')

  for nk in range(0,len(kernels)):
    file_text = file_text +\
    '#include "'+kernels[nk]['name']+'_kernel.cu"\n'

  master = master.split('.')[0]
  fid = open('cuda/'+master.split('.')[0]+'_kernels.cu','w')
  fid.write('//\n// auto-generated by op2.py\n//\n\n')
  fid.write(file_text)
  fid.close()

