#*******************************************************************************
# PROGRAMME	: IO
# AUTEUR	: C. BACOUR
# CREATION	: 11/2007
# COMPILATEUR	: PYTHON
#
# ORCHidee Assimilation System
#
# Functions to read/write NetCDF or FORTRAN files
#
#*******************************************************************************

import os, sys, copy, glob
from Scientific.IO import NetCDF
import numpy as np

#from TOOLS import *
#import *

# =============================================================================
# Creation recursive des repertoires
# =============================================================================
def mkdir(dirname):

    print
    
    for i in range(len(dirname.split('/'))):
        el = dirname.split('/')[i]
        if i == 0:
            dir = os.path.join('/',el)
        else:
            dir = os.path.join(dir,el)

        if len(glob.glob(dir)) == 0:
            print " + creation " + dir
            os.system('mkdir '+dir)
# END
# ------------------------------------------------------------------------------


# =============================================================================
# If batch mode : copy the results to the output directory
# =============================================================================
def batch_cp(Config,Site,logfile):


    ###, filecmp

    # - Create the output directory if it does not exist
    if len(glob.glob( os.path.join(Config.PATH_MAIN,Site.base_name) )) == 0:
        mkdir(os.path.join(Config.PATH_MAIN,Site.base_name))
        
    # - Copy
    print '\n### BATCH_MODE : cp -pR ' + Config.PATH_EXEC + ' ' + os.path.join(Config.PATH_MAIN, Site.base_name)
    logfile.write('\n### BATCH_MODE : cp -pR ' + Config.PATH_EXEC + '/* ' + os.path.join(Config.PATH_MAIN, Site.base_name)+'\n')
    os.system('cp -pR ' + Config.PATH_EXEC + '/* ' + os.path.join(Config.PATH_MAIN, Site.base_name))
        
    # - Compare the 2 directories : if everything's similar => erase PATH_EXEC
    ## cmpfiles = filecmp.dircmp(Config.PATH_EXEC, os.path.join(Config.PATH_MAIN, Site.base_name))
##     if len(cmpfiles.diff_files)==0:
##         print '### BATCH_MODE : NO difference between \n  +'+Config.PATH_EXEC+ \
##               '\n    and\n  +'+ os.path.join(Config.PATH_MAIN, Site.base_name)
##         print '### BATCH_MODE : => erasing of '+Config.PATH_EXEC
##         logfile.write('### BATCH_MODE : NO difference between \n  +'+Config.PATH_EXEC+\
##                       '\n    and\n  +'+ os.path.join(Config.PATH_MAIN, Site.base_name)+'\n')
##         logfile.write('### BATCH_MODE : => erasing of '+Config.PATH_EXEC+'\n')
##         os.system('rm -fR '+Config.PATH_EXEC)
        
##     else:
##         print '### BATCH_MODE : !!!!! WARNING !!!!!'
##         print '### BATCH_MODE : SOME differences detected between\n  +'+Config.PATH_EXEC+\
##               '\n    and\n  +'+ os.path.join(Config.PATH_MAIN, Site.base_name)
##         print '### BATCH_MODE : we do not erase '+Config.PATH_EXEC
##         logfile.write('### BATCH_MODE : !!!!! WARNING !!!!!\n')
##         logfile.write('### BATCH_MODE : SOME differences detected between \n  +'+Config.PATH_EXEC+\
##                       '\n    and\n  +'+ os.path.join(Config.PATH_MAIN, Site.base_name)+'\n')
##         logfile.write('### BATCH_MODE : we do not erase '+Config.PATH_EXEC+'\n')
        
            
# END
# ------------------------------------------------------------------------------




# =============================================================================
# Lecture fichiers NetCDF
# Entrees : filename, nom du fichier d'entree
#           ovar (0|1), retourne ou pas le dictionnaire des infos sur les
#           variables
#           ogattr (0|1), retourne ou pas le dictionnaire des infos sur les
#           attributs globaux
#
# Sortie  : [vars, attr, dims] TUPLE
#           ou vars : LISTE de dictionnaires (1 par variable)
#           ou attr : LISTE de dictionnaires (1 par attribut global)
#           ou dims : LISTE de dictionnaires (1 par dimension)
# =============================================================================

def readnc(filename):
    """[vars, gattr, dims] = readnc(filename)

    Lecture des donnees contenues dans un fichier NetCDF
    Retourne les donnees sous forme de deux tableaux de dictionnaires,
    1 pour les variables et 1 pour les attributs globaux:

       vars = {
               'Varname1' :{ 
                          'datatype':,
                          'ndims':,
                          'dim_name':,   (TUPLE)
                          'dim_size':,
                          'value':,
                          'attr_name':,
                          'attr_value':,'
                          }
               'Varname2' :{ 
                          'datatype':,
                          'ndims':,
                          'dim_name':,   (TUPLE)
                          'dim_size':,
                          'value':,
                          'attr_name':,
                          'attr_value':,'
                          }
               ...
               }

       gattr = {'name':,
                'value':,
               }

       dims = {'name':,
               'size':,
               }
    
    """


    # -- Constantes --
    # attributs NetCDF predefinis a ne pas lire
    attributes_var_except = ['assignValue','getValue','typecode'] # pour les variables
    attributes_except = ['close', 'createDimension', 'createVariable','flush', 'sync'] # pour les attributs globaux


    vars=[]
    gattr=[]
    dims=[]
    
    
    # -- Ouverture du fichier --
    try:
        file = NetCDF.NetCDFFile(filename,'r')
    except:
        sys.exit('* READNC : Le fichier %s n\'est pas trouve' %(filename))
        
        
    # -- Attributs globaux --
    attributes = dir(file)
    gattr = []
    
    for elem in attributes:
        if elem not in attributes_except:
            dico = {}
            dico['name']  = elem
            dico['value'] = getattr(file,elem)
            
            gattr.append(dico)


    # -- Dimensions --
    #  Remarque: pour recuperer la dimensions d'une dimension UNLIMITED,
    #  var.dimensions ne fonctionne pas => utiliser var.variable.shape
    file_dim_keys  = file.dimensions.keys()
    
    dimensions = []
    for elem in file_dim_keys:
        dico = {}
        dico['name'] = elem
        if file.dimensions[elem] != None:
            dico['size'] = file.dimensions[elem]
        else:
            dico['size'] = 'UNLIMITED'

        dims.append(dico)
    

    # -- Variables ---
    
    keys = file.variables.keys()
    vars = {}
    for elem in keys:

        varObj = file.variables[elem]
        
        vars[elem]= {}
        vars[elem]['datatype']   = varObj.typecode()
        vars[elem]['dim_name']   = varObj.dimensions
        vars[elem]['dim_size']   = varObj.shape
        vars[elem]['value']      = varObj.getValue()
        
        #print 'lecture ' + elem
        
        # attributs
        attr_name = dir(varObj)         
        for elem_remove in attributes_var_except:
            attr_name.remove(elem_remove)
            
        attr_value = []
        for elem_attr in attr_name:
            attr_value.append(getattr(varObj,elem_attr))
            
        vars[elem]['attr_name']  = attr_name
        vars[elem]['attr_value'] = attr_value
            
            
    file.close()

    # -- Tuple de sortie ---
    return vars, gattr,dims
   
# END
# ------------------------------------------------------------------------------



# ===============================================================================
# Ecriture fichiers NetCDF
# Entrees : filename, nom du fichier d'entree
#           vars : LISTE de dictionnaires (1 par variable)
#           attr : LISTE de dictionnaires (1 par attribut global)
#           dims : LISTE de dictionnaires (1 par dimension)
#           append : ajout ou non a 1 fichier existant
# ===============================================================================
def writenc(filename, vars =[], gattr=[], dims = [], append = 0, dim_order = None, var_order = None):
    """def writenc(filename, vars =(-1), gattr=(-1), dims = [-1], append = 0):

    Ecriture de donnees dans un fichier NetCDF
    
       vars = {
               'Varname1' :  {
                             'datatype':,
                             'ndims':,
                             'dim_name':,   (TUPLE)
                             'dim_size':,
                             'value':,
                             'attr_name' : [''],
                             'attr_value': [''],'
                             }
                ......

       dims = [{'name':,
               'size':  },
              {'name':,
               'size':  }]

       gattr = [{'name':,
                'value': },
               {'name':,
                'value':   }]
    
    """


    # -- Ouverture du fichier --
    if append:
        mode = 'a'
    else:
        mode = 'w'
    
        
    try:
        file = NetCDF.NetCDFFile(filename,mode)
    except:
        sys.exit('* WRITENC : Le fichier %s n\'a pu etre cree' %(filename))
        
        
    # -- Attributs globaux --
    if len(gattr) > 0:
                
        if type(gattr) == type({}):   # dictionnaire
            gattr = [gattr]
        
        for i in range(len(gattr)):
            elem = gattr[i]
            setattr(file,elem['name'] , elem['value'])

    # -- Dimensions --
    if len(dims) > 0:

        dim_unlimited = -999.
        
        if type(dims) == type({}):   # dictionnaire
            dims = [dims]

        # if a dimension order is prescribed, change the order of the dimensions
        if dim_order != None:

            if len(dim_order) != len(dims): sys.exit('* WRITENC : you have provided DIM_ORDER with a length that does not match that of DIMS')

            buf = copy.copy(dims)
            bufname=[];ind=[]
            for i in range(len(dim_order)): bufname.append(buf[i]['name'])
            for i in range(len(dim_order)): ind.append( bufname.index(dim_order[i]) )

            for i in range(len(dim_order)):
                buf[i] = dims[ind[i]]
            
            dims = buf
        

        # create the dimension in the NetCDF file
        for i in range(len(dims)):
            elem = dims[i]

            if elem['size'] == 'UNLIMITED': # on recherche nombre d'elements de cette dimension
                dim_name = elem['name']
                var_name = vars.keys()
                icnt = 0
                while dim_unlimited == -999.: 
                    for vdim in list(vars[var_name[icnt]]['dim_name']):
                        if vdim == dim_name:
                            indice=list(vars[var_name[icnt]]['dim_name']).index(vdim)
                            dim_unlimited=int(vars[var_name[icnt]]['value'].shape[indice])
                            break
                    icnt = icnt+1

                elem['size'] = dim_unlimited

                
            else:
                elem['size'] = int(elem['size'])
           
            file.createDimension(elem['name'], elem['size'])

            
    file_dim_keys  = file.dimensions.keys()
    dim_size = []
    for elt in file_dim_keys:  dim_size.append(file.dimensions[elt])
            
    
    # -- Variables ---
    if len(vars) > 0:

        
         # if a variable order is prescribed, change the order of the variables
        if var_order != None:
            if len(var_order) != len(vars.keys()): sys.exit('* WRITENC : you have provided VAR_ORDER with a length that does not match that of VAR')
            var_name = var_order
        else:
            var_name = vars.keys() # nom des elements

            
        # write the variables

        for elname in var_name:          

            #si type 'list'
            if type(vars[elname]['value']) == type([]): 
                vars[elname]['value'] = np.array(vars[elname]['value'])
           
            # Definition de la variable
            if vars[elname].has_key('datatype') == True:
                vtype = vars[elname]['datatype']
            else:
                vtype = vars[elname]['value'].dtype

            if vtype != 'c' and  type(vars[elname]['value']) != type(np.array(1)):
                vars[elname]['value'] = np.array(vars[elname]['value'])

            if vtype == 'i' :
                var_val = vars[elname]['value'].astype(np.int32) #np.array(vars[elname]['value'], np.int32)
            if vtype == 'l' :
                var_val = vars[elname]['value'].astype(np.int64) #np.array(vars[elname]['value'], np.int64)
            if vtype == 'f':
                var_val = vars[elname]['value'].astype(np.float32) #np.array(vars[elname]['value'], np.float64)
            if vtype == 'd':
                var_val = vars[elname]['value'].astype(np.float64) #np.array(vars[elname]['value'], np.float64)


            netCDFVar = file.createVariable(elname, vtype,
                                            (vars[elname]['dim_name'])) 
            # Attributs de la variable
            attnames = vars[elname]['attr_name']
            attvalues = vars[elname]['attr_value']
            for j in range(len(attnames)):
                setattr(netCDFVar, attnames[j], attvalues[j])
            
            # Redimensionnement pour etre en accord avec les dimensions                     
            var_dim_name = vars[elname]['dim_name']
            dim_id = []
            shape_resize = []
            for vdim in list(var_dim_name):
                dim_id.append(file_dim_keys.index(vdim))
                shape_resize.append(dim_size[file_dim_keys.index(vdim)])
            shape_resize = tuple(shape_resize)
                
            # Valeurs
            if vtype != 'c':
                #print
                #print elname, var_val
                #print var_dim_name
                #print file_dim_keys, dim_size
                #print shape_resize
                var_val = np.reshape(var_val, shape_resize)
            else:
                var_val = vars[elname]['value']
                
            netCDFVar.assignValue(var_val)
                        
            del var_val

    file.close()
# END
# ------------------------------------------------------------------------------



# ===============================================================================
# Ecriture d'une variable dans un fichier FORTRAN UNFORMATTED
# Entrees :
#  - ID fichier
#  - variable (tableau)
# ===============================================================================
def var_writefor(f, swapendian, var):


    # -- Calcul du nombre d'octets de la variable
    
    ndims = len(var.shape)
    size = 1
    for i in range(ndims):
        size = size * var.shape[i]

    # -- Permette ordre des dimensions si tableau
    #   (ordre des dimensions inverse entre PYTHON vs FORTRAN)
    if ndims > 1: var = np.transpose(var)
    
    longueur_octets = np.array(size * var.itemsize, np.int32)

    # -- Ecriture

    #if swapendian == 0: 
    f.write(longueur_octets.tostring())
    f.write(var.tostring())
    f.write(longueur_octets.tostring())
    #else:                # swap bytes
    #    print 'on swap les bits'
    #    f.write(longueur_octets.byteswapped().tostring())
    #    f.write(var.byteswapped().tostring())
    #   f.write(longueur_octets.byteswapped().tostring())
# END
# ------------------------------------------------------------------------------



# ===============================================================================
# Lecture d'une variable dans un fichier FORTRAN UNFORMATTED
# Entrees :
#  - ID fichier
#  - variable (tableau)
# ===============================================================================
def var_readfor(f, swapendian, var):

    # -- Calcul du nombre d'octets de la variable
    
    lo = f.read(4)
    #if swapendian == 0:
    longueur_octets = np.fromstring(lo, np.int32)[0]
    #else:              # swap bytes
    #longueur_octets = np.fromstring(lo, np.int32).byteswapped()[0]
    #print 'longueur_octets',longueur_octets
    
    
    # -- Lecture de la variable
    shape = var.shape
    data = f.read(longueur_octets)

    #if swapendian == 0:
    data = np.fromstring(data, np.float64)
    
    #else:              # swap bytes
    #    data = np.fromstring(data, np.float64).byteswapped()
    data.shape = shape

    # -- Lecture du nombre d'octets
    lo = f.read(4)
    #if swapendian == 0:
    longueur_octets = np.fromstring(lo, np.int32)[0]
    #else:              # swap bytes
    #longueur_octets = np.fromstring(lo, np.int32).byteswapped()[0]
    #print     longueur_octets
    
    return data
    


# END
# ------------------------------------------------------------------------------
