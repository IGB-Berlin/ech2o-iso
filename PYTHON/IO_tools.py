# =============================================================================
# Creation recursive des repertoires
# =============================================================================
def mkdir(dirname):
    import os,glob

    path = os.getcwd()+'/'
    
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



#*******************************************************************************
# Retourne les positions de caracteres 'cible' dans une chaine de caractere
# Entree : chaine   STRING
#          cibles   LISTE de STRING
#
# Sortie : indices  LISTE
#*******************************************************************************
def retind(chaine, cibles):
    """indices = def retind(chaine, cibles):
    Cherche la position des caracteres 'cible' dans la chaine 'chaine'
    et retourne la position 'indices' des elements 'cible' dans 'chaine'"""

    import Numeric as N
    
    pos = -1
    indices = [-1]
    string = N.array(chaine)

    for icible in range(len(cibles)):
        pos = N.nonzero((string==cibles[icible]))
        if pos != -1:
            if icible == 0:
                indices = N.array(pos).tolist()
            else:
                indices = [indices, N.array(pos).tolist()]

    return indices
# END
# ------------------------------------------------------------------------------


# ==============================================================================
# Sort out a list according to a predifined order
# provided in the template
#
# ------------------------------------------------------------------------------
def triname(names, template):

    name = []
    for pname in template:
        try:
            test = names.index(pname)
        except:
            test = -1
        if test != -1: name.append(pname)
    return name
# ==============================================================================


# ===============================================================
# Genere l'ecriture d'un fichier texte
#
# Inputs : - chaine de caractere
#          - nom du fichier de sortie
#
# Outputs :
# ---------------------------------------------------------------
def text_write(config, fileout_name):
    fout = open(fileout_name, 'w')
    fout.write(config)
    fout.close()
# fin orchidee_write
# ---------------------------------------------------------------


#*******************************************************************************
# Lit les valeurs contenues dans un fichier ASCII les retourne sous forme
# de tableau (marche uniquement pour un fichier vecteur pour le moment)
#*******************************************************************************
## def get_array_from_ascii(file):
##     """array = def get_array_from_ascii(file):
##     Lit les valeurs contenues dans le fichier ascii 'file'
##     et les retourne sous forme de tableau

##     RQ : marche pour des vecteurs pour le moment"""

##     import os
    

##     fic = open(file,'r')
##     lignes = fic.readlines()
##     fic.close()

##     nlignes = len(lignes)

##     array = [1e20]
##     for lig in lignes:
##         array.append(float(lig.strip()))
##     array.remove(1e20)
    
##     return array
# END
# ------------------------------------------------------------------------------


#*******************************************************************************
# Lit les valeurs contenues dans un fichier ASCII les retourne sous forme
# de tableau 
#*******************************************************************************
def get_array_from_ascii(file):
    """array = def get_array_from_ascii(file):
    Lit les valeurs contenues dans le fichier ascii 'file'
    et les retourne sous forme de tableau
    """

    import os
    import Numeric as N
    
    fic = open(file,'r')

    # Initialisation du tableau
    lignes = fic.readlines()
    nlig = len(lignes)
    ncol = len(lignes[0].split())
    #print 'nlignes:',nlig
    #print 'ncolonnes:',ncol
    
    data = N.zeros((nlig,ncol),N.Float64)
    
    for i in range(nlig):
        ligne = lignes[i].split()
        for j in range(ncol):
            data[i,j] = float(ligne[j])
    fic.close()


    return data
    
# END
# ------------------------------------------------------------------------------


#*******************************************************************************
# Ecrit un tableau de valeurs sous forme de fichier ASCII
#*******************************************************************************
def write_array_to_ascii(array,file=None ,mode = 'w',sep= '  '):
    """def write_array_to_ascii(array,file,mode,sep):
    Ecrit les valeurs d'un tableau  dans le fichier ascii 'file'
    en fonction du mode d'ecriture 'w' (write) ou 'a' (append)
    Les valeurs sont separees par le separateur sep
    """

    import os

    if file == None : sys.exit('WRITE_ARRAY_TO_ASCII : you must provide a file name. STOP')

    # On passe le tableau en N.array si c'est une liste
    if type(array) == type([]):
        import Numeric as N
        array = N.array(array)

    Nlig = array.shape[0]
    
    # generation du format chaine de caracteres
    chaine = []
    for i in range(Nlig):
        buf = (str(array[i]))
        buf = buf.replace('[','')
        buf = buf.replace(']','')
        buf = buf.replace('\n','')
        buf = buf.strip(' ')
        chaine.append(buf + sep + '\n')


    # ecriture
    fic = open(file,mode)
    for i in range(Nlig):
        fic.write(chaine[i])
    fic.close
                       
    
# END
# ------------------------------------------------------------------------------



#*******************************************************************************
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
#*******************************************************************************
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

    from Scientific.IO import NetCDF
    import os, sys


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



#*******************************************************************************
# Ecriture fichiers NetCDF
# Entrees : filename, nom du fichier d'entree
#           vars : LISTE de dictionnaires (1 par variable)
#           attr : LISTE de dictionnaires (1 par attribut global)
#           dims : LISTE de dictionnaires (1 par dimension)
#           append : ajout ou non a 1 fichier existant
#*******************************************************************************
def writenc(filename, vars = None, gattr = None, dims = None, append = 0, dim_order = None, var_order = None):
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

    import Numeric as N
    from Scientific.IO import NetCDF
    import os, sys, copy


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
    if gattr != None:
                
        if type(gattr) == type({}):   # dictionnaire
            gattr = [gattr]
        
        for i in range(len(gattr)):
            elem = gattr[i]
            setattr(file,elem['name'] , elem['value'])

    # -- Dimensions --
    if dims != None:

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

                if vars != None:
                    var_name = vars.keys()
                    icnt = 0
                    while dim_unlimited == -999.: 
                        for vdim in list(vars[var_name[icnt]]['dim_name']):
                            if vdim == dim_name:
                                indice=list(vars[var_name[icnt]]['dim_name']).index(vdim)
                                dim_unlimited=int(vars[var_name[icnt]]['value'].shape[indice])
                                break
                        icnt = icnt+1

                else:
                    dim_unlimited = 1
                    
                elem['size'] = dim_unlimited

                
            else:
                elem['size'] = int(elem['size'])
           
            file.createDimension(elem['name'], elem['size'])


            
    # -- Variables ---
    if vars != None:

        
         # if a variable order is prescribed, change the order of the variables
        if var_order != None:
            if len(var_order) != len(vars.keys()): sys.exit('* WRITENC : you have provided VAR_ORDER with a length that does not match that of VAR')
            var_name = var_order
        else:
            var_name = vars.keys() # nom des elements

            
        # write the variables
        for elname in var_name:
          
            # Definition de la variable
            if vars[elname].has_key('datatype') == True:
                vtype = vars[elname]['datatype']
            else:
                if type(vars[elname]['value']) == type([]): #'list':
                    vars[elname]['value'] = N.array(vars[elname]['value'])
                vtype = vars[elname]['value'].typecode()

            if vtype == 'i' :
                var_val = N.array(vars[elname]['value'], N.Int32, savespace=1)
            if vtype == 'l' :
                var_val = N.array(vars[elname]['value'], N.Int32, savespace=1)
            if vtype == 'f':
                var_val = N.array(vars[elname]['value'], N.Float32, savespace=1)
            if vtype == 'd':
                var_val = N.array(vars[elname]['value'], N.Float64, savespace=1)

            ###ancienne originale
            ###var_val = N.array(vars[elname]['value'], N.Float64, savespace=1)
            ###
            ###if vars[elname].has_key('datatype') == True:
            ###    vtype = vars[elname]['datatype']
            ###else:
            ###    vtype = var_val.typecode()

        
            netCDFVar = file.createVariable(elname, vtype,
                                            (vars[elname]['dim_name'])) 

            # Attributs de la variable
            attnames = vars[elname]['attr_name']
            attvalues = vars[elname]['attr_value']
            for j in range(len(attnames)):
                setattr(netCDFVar, attnames[j], attvalues[j])
            
            # Redimensionnement pour etre en accord avec les dimensions
            file_dim_keys  = file.dimensions.keys()
            dim_size = []
            for elt in file_dim_keys:  dim_size.append(file.dimensions[elt])
                      
            var_dim_name = vars[elname]['dim_name']
            dim_id = []
            shape_resize = []
            for vdim in list(var_dim_name):
                dim_id.append(file_dim_keys.index(vdim))
                shape_resize.append(dim_size[file_dim_keys.index(vdim)])
            shape_resize = tuple(shape_resize)
                
            # Valeurs
            if vtype != 'c':
                #print elname
                #print var_val
                #print shape_resize
                var_val = N.reshape(var_val, shape_resize)
            else:
                var_val = vars[elname]['value']

            netCDFVar.assignValue(var_val)
                        
            del var_val

    file.close()
# END
# ------------------------------------------------------------------------------



#*******************************************************************************
# Ecriture d'une variable dans un fichier FORTRAN UNFORMATTED
# Entrees :
#  - ID fichier
#  - variable (tableau)
#*******************************************************************************
def write_var_to_fortran_file(f, var):

    import Numeric as N

    # Calcul du nombre d'octets de la variables
    ndims = len(var.shape)
    size = 0
    for i in range(ndims):
        size = size * var.shape[i]

        
    longueur_octets = N.array(size * var.itemsize(), N.Int32)
    
    # ecriture
    f.write(longueur_octets.tostring())
    f.write(var.tostring())
    f.write(longueur_octets.tostring())
    
# END
# ------------------------------------------------------------------------------
